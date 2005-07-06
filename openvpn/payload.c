/*
 ***** THIS IS WORK IN PROGRESS (and correctly ifdef'd out) --jjo Jun/2005 *****
 * 
 *  Payload conntrack optimizations for OpenVPN
 *    Intended for filtering TCP retransmissions over reliable links
 *
 * QUICK copy-n-paste for jjo:
 * make payload.o CFLAGS="-D USE_PAYLOAD_CONNTRACK=1 -Wall "
 *
 *  Author: JuanJo Ciarlante <jjo@mendoza.gov.ar>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (see the file COPYING included with this
 *  distribution); if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef WIN32
#include "config-win32.h"
#else
#include "config.h"
#endif

#ifdef USE_PAYLOAD_CONNTRACK

#include "syshead.h"
#include "forward.h"
#include "memdbg.h"

#include "list.h"
#include "forward-inline.h"
#include "payload.h"
#define PAYLOAD_HASH_SIZE 256 /* Should be configurable by option */

#define PAYLOAD_N_TCPSEGS 8
/* 
 * struct payload_tuple_id : uniq 5-upla (proto==TCP is implicit)
 */
struct payload_tuple_id {
  /* BEGIN uniq TCP 5-upla id */
  uint32_t ip_saddr, ip_daddr;
  uint16_t tcp_sport, tcp_dport;
  /* END   uniq TCP 5-upla id */
};

/* 
 * payload_tuple: 1 per TCP connection, currently only TX side hook
 */
struct payload_tuple {
  struct payload_tuple_id id;
  /* round robin array with PAYLOAD_N_TCPSEGS latest tcp segments: */
  struct {
    struct openvpn_tcphdr tcph;
    int tcp_len;
    int hits;
  } tcp_seg[PAYLOAD_N_TCPSEGS];
  int tcp_seg_idx; 	/* next slot to use */
  time_t last_used;
  time_t expires;
  int conn_hits;
  int deleted;
};

/* shortcuts macros for easier printing */
#define PAYLOAD_FMT_MIN "0x%08x:%d -> 0x%08x:%d"
#define PAYLOAD_FMT_MIN_ARGS(pip, ptcp) \
					ntohl(pip->saddr), ntohs(ptcp->source),\
					ntohl(pip->daddr), ntohs(ptcp->dest)

#define PAYLOAD_FMT_MED PAYLOAD_FMT_MIN " (seq=%08x, ack_seq=%08x)"
#define PAYLOAD_FMT_MED_ARGS(pip, ptcp) PAYLOAD_FMT_MIN_ARGS(pip, ptcp) ,\
					ntohl(ptcp->seq),ntohl(ptcp->ack_seq)

#define PAYLOAD_FMT_FULL PAYLOAD_FMT_MED " w=%d [%c%c%c%c%c%c]" /* doesn't show WSCALE */
#define PAYLOAD_FMT_FULL_ARGS(pip, ptcp) PAYLOAD_FMT_MED_ARGS(pip,ptcp) ,\
					ntohs(ptcp->window), \
					OPENVPN_TCPH_FIN_MASK & ptcp->flags ? 'F' : '.', \
					OPENVPN_TCPH_SYN_MASK & ptcp->flags ? 'S' : '.', \
					OPENVPN_TCPH_RST_MASK & ptcp->flags ? 'R' : '.', \
					OPENVPN_TCPH_PSH_MASK & ptcp->flags ? 'P' : '.', \
					OPENVPN_TCPH_ACK_MASK & ptcp->flags ? 'A' : '.', \
					OPENVPN_TCPH_URG_MASK & ptcp->flags ? 'P' : '.' 
/* hash related functions */
static inline const uint8_t *
payload_tuple_hash_ptr (const struct payload_tuple *pt)
{
  return (uint8_t *) &pt->id;
}

static inline uint32_t
payload_tuple_hash_len (const struct payload_tuple *pt)
{
  return (uint32_t) sizeof (pt->id);
}

uint32_t 
payload_tuple_hash_func (const void *key, uint32_t iv)
{
  return hash_func(payload_tuple_hash_ptr((const struct payload_tuple *)key),
		   payload_tuple_hash_len((const struct payload_tuple *)key),
		   iv);
}

bool
payload_tuple_compare_func(const void *key1, const void *key2)
{
  return memcmp(&((const struct payload_tuple *)key1)->id, 
		&((const struct payload_tuple *)key2)->id, 
		sizeof ((const struct payload_tuple *)key1)->id) == 0;
}


/* create a new conntrack entry */
static struct payload_tuple *
payload_tuple_new(void)
{
  struct payload_tuple *pt;
  ALLOC_OBJ_CLEAR(pt, struct payload_tuple);
  return pt;
}
void
payload_tuple_delete(struct payload_tuple *pt)
{
  free(pt);
}

/* initialize conntrack entry with {ip,tcp} hdr data */
static void
payload_tuple_id_init(struct payload_tuple_id *pt_id, const struct openvpn_iphdr *pip, const struct openvpn_tcphdr *ptcp)
{
  ASSERT(pt_id);
  CLEAR(*pt_id);
  pt_id->ip_saddr=pip->saddr;
  pt_id->ip_daddr=pip->daddr;
  pt_id->tcp_sport=ptcp->source;
  pt_id->tcp_dport=ptcp->dest;
}

/*
 * How many buckets in hash to gc per pass. (shamelessly stolen from multi.c)
 */
#define PAYLOAD_GC_DIVISOR     256  /* How many passes to cover whole hash table */
#define PAYLOAD_GC_MIN          16  /* Minimum number of buckets per pass */
#define PAYLOAD_GC_MAX        1024  /* Maximum number of buckets per pass */
static int
payload_gc_buckets_per_pass (int n_buckets)
{
  return constrain_int (n_buckets / PAYLOAD_GC_DIVISOR, PAYLOAD_GC_MIN, PAYLOAD_GC_MAX);
}

/*
 * Adjust next event wakeup based on how many conntrack entries present
 * currently it programs
 */
static void
payload_gc_adjust_timers(struct context *c)
{
  struct payload_context *pc=c->c2.payload_context;
  int n=hash_n_elements(pc->hash);
  if (n>0)
    {
      if (pc->gc.wakeup.n != 1)
	event_timeout_init (&pc->gc.wakeup, 1, now);
      context_reschedule_sec(c, 1);
      reset_coarse_timers(c);
      dmsg(D_PAYLOAD_CONNTRACK, "payload_gc_adjust_timers: n_elem=%d wakeup.n=%d", n, pc->gc.wakeup.n);
    }
  else
    {
      event_timeout_clear(&pc->gc.wakeup);
      dmsg(D_PAYLOAD_CONNTRACK, "payload_gc_adjust_timers: n_elem=%d timeout_cleared", n);
    }
}
/*
 * garbage collect old (pt->deleted) entries
 */
int
payload_gc_run(struct payload_context *pc, int start_bucket, int end_bucket)
{
  int n_deleted=0, n_total=0, n_total_hits=0;
  struct hash_element *he;
  struct hash_iterator hi;
  if (start_bucket < 0)
    {
      start_bucket = 0;
      end_bucket = hash_n_buckets (pc->hash);
    }
  hash_iterator_init_range (pc->hash, &hi, true, start_bucket, end_bucket);
  while ((he = hash_iterator_next (&hi)) != NULL)
    {
      struct payload_tuple *pt = (struct payload_tuple *) he->value;
      n_total++;
      n_total_hits+=pt->conn_hits;
      if (now > pt->expires)
	pt->deleted++;
      /*
       * could have been marked syncronously (from payload_tcp_dd ...) or
       * async (here) by above check
       */
      if (pt->deleted)
	{
	  n_deleted++;
	  dmsg(D_PAYLOAD_CONNTRACK, "payload_gc_run DELETED " PAYLOAD_FMT_MIN, 
	      ntohl(pt->id.ip_saddr), ntohs(pt->id.tcp_sport),
	      ntohl(pt->id.ip_daddr), ntohs(pt->id.tcp_dport));
	  payload_tuple_delete(pt);
	  hash_iterator_delete_element (&hi);
	}
    }
  hash_iterator_free (&hi);
  dmsg(D_PAYLOAD_CONNTRACK, "payload_gc_run(%d,%d) DELETED %d/%d entries, n_elem=%d, total_hits=%d" , 
      start_bucket, end_bucket,
      n_deleted, n_total, hash_n_elements(pc->hash), n_total_hits);
  return n_deleted;
}

static void
payload_gc_dowork(struct payload_context *pc)
{
  if (pc->gc.bucket_base >= hash_n_buckets(pc->hash))
    pc->gc.bucket_base = 0;
  payload_gc_run(pc, pc->gc.bucket_base, pc->gc.bucket_base+pc->gc.buckets_per_pass);
  pc->gc.bucket_base+= pc->gc.buckets_per_pass;
}

/*
 * conntrack cleaner from event ticks
 */
void
check_payload_gc_dowork (struct context *c)
{
  struct payload_context *pc=c->c2.payload_context;
  ASSERT(pc);
  payload_gc_dowork(pc);
  payload_gc_adjust_timers(c);
}

/*
 * Loop over TCP options, call callback() if matches passed TCP opt  (wildcard: OPENVPN_TCPOPT_ANY)
 *
 * Some examples:
 * - is SACK present?
 *    if(tcp_opt_process(buf, OPENVPN_TCPOPT_SACK, NULL, NULL, pip, ptcp)) ...
 * - if SACK is present, call myfunc(..., myarg, ...)
 *    if(tcp_opt_process(buf, OPENVPN_TCPOPT_SACK, myfunc, myarg, pip, ptcp)) ...
 * - is any option present? (except EOL, NOP)
 *    if(tcp_opt_process(buf, OPENVPN_TCPOPT_ANY, NULL, NULL, pip, ptcp)) ...
 * - for each option (except EOL, NOP) call myfunc(..., myarg, ...)
 *    if(tcp_opt_process(buf, OPENVPN_TCPOPT_ANY, myfunc, myarg, pip, ptcp)) ...
 *    
 */
static inline bool
tcp_opt_process (struct buffer *buf, int optnum,  bool (*callback)(uint8_t *opt, int optlen, void *callback_arg), void *callback_arg)
{
  int hlen, olen, optlen;
  uint8_t *opt;
  struct openvpn_tcphdr *ptcp;

  ASSERT (BLEN (buf) >= (int) sizeof (struct openvpn_tcphdr));

  verify_align_4 (buf);

  ptcp = (struct openvpn_tcphdr *) BPTR (buf);
  hlen = OPENVPN_TCPH_GET_DOFF (ptcp->doff_res);

  /* Invalid header length or header without options. */
  if (hlen <= (int) sizeof (struct openvpn_tcphdr)
      || hlen > BLEN (buf))
    return false;

  for (olen = hlen - sizeof (struct openvpn_tcphdr),
	 opt = (uint8_t *)(ptcp + 1);
       olen > 0;
       olen -= optlen, opt += optlen) {
    if (*opt == OPENVPN_TCPOPT_EOL)
      break;
    else if (*opt == OPENVPN_TCPOPT_NOP)
      optlen = 1;
    else {
      optlen = *(opt + 1);
      if (optlen <= 0 || optlen > olen)
        break;
      /* 
       * TCP opt found, callback() passed function and return true if callback succeded 
       */

      if (optnum == OPENVPN_TCPOPT_ANY || optnum == *opt) {
	if (callback) {
	  if (callback(opt, optlen, callback_arg))
	    return true;
	} else 
	  /* if no callback function, just return true (opt found) */
	  return true;
      }
    }
  }
  return false;
}

static inline bool
tcp_dd_opt_skip_segment(uint8_t *opt, int optlen, void *arg)
{
  switch(*opt)
    {
    case OPENVPN_TCPOPT_SACK:
    case OPENVPN_TCPOPT_WSCALE:
      return true;
    }
  return false;
}
/*
 * careful logic: TCP may update window, probe 0-window and alike
 * yet could be seen as a retransmission
 */

static inline int
payload_tcp_dd_drop_hit(struct context *c, const struct openvpn_iphdr *pip, const struct openvpn_tcphdr *ptcp, struct buffer *buf)
{
  bool result=false;
  struct hash_element *he;
  uint32_t hv;
  struct hash_bucket *bucket;
  struct payload_context *pc=c->c2.payload_context;
  struct payload_tuple_id pt_id;
  struct payload_tuple *pt;
  struct openvpn_tcphdr tcph = *ptcp; /* local copy: stored tcph is NOT equal to sniffed one */
  int ip_totlen=ntohs(pip->tot_len);
  int tcp_len=ip_totlen-OPENVPN_IPH_GET_LEN (pip->version_len); /* tcphdr+DATA */

  ASSERT(pc);
  tcph.check=0; /* TCP chksum will vary with eg. tstamp */
  payload_tuple_id_init(&pt_id, pip, ptcp);

  hv= hash_value (pc->hash, &pt_id);
  bucket= hash_bucket (pc->hash, hv);
  hash_bucket_lock (bucket);
  he = hash_lookup_fast (pc->hash, bucket, &pt_id, hv);
  if (he)
    {
      int i;
      pt = (struct payload_tuple *) he->value;
      /*
       * Avoid filtering out if:
       * - SYN or FIN
       * - zero window 
       * - zero window probe (data size=1)
       * - SACK or WSCALE option present
       */
      if ( (ptcp->flags & (OPENVPN_TCPH_SYN_MASK|OPENVPN_TCPH_FIN_MASK))
	   || tcph.window == 0
	   || (ip_totlen-OPENVPN_TCPH_GET_DOFF(tcph.doff_res))==1
	   || tcp_opt_process(buf, OPENVPN_TCPOPT_ANY, tcp_dd_opt_skip_segment, NULL))
	{
	  dmsg(D_PAYLOAD_CONNTRACK, "payload_tcp_dd_drop_hit: SKIP segment " PAYLOAD_FMT_FULL, 
	      PAYLOAD_FMT_FULL_ARGS(pip, ptcp));
	  goto done;
	}

      /* search for seen TCP header */
      for(i=0;i<PAYLOAD_N_TCPSEGS;i++)
	{
	  /* 
	   * if:
	   * 	same seq nums (optimization to avoid next checks if unequal)
	   * 	AND same tcp segment header (w/o options)
	   * 	AND same tcp header len (incl. options sizes)
	   *
	   */
	  if (pt->tcp_seg[i].tcph.seq==tcph.seq 
	      && memcmp(pt->tcp_seg+i, &tcph, sizeof tcph)==0
	      && pt->tcp_seg[i].tcp_len == tcp_len)
	    {
	      pt->tcp_seg[i].hits++;
	      pt->conn_hits++;
	      dmsg(D_PAYLOAD_CONNTRACK, "payload_tcp_dd_drop_hit HIT! conn_hits=%d timeleft=%ld " PAYLOAD_FMT_FULL, 
		  pt->conn_hits, pt->expires-now,
		  PAYLOAD_FMT_FULL_ARGS(pip, ptcp));
	      result=true;	/* _IS_ a dup */
	      goto done;
	    }
	}
      if (pt->expires < now)
	{
	  pt->deleted=true;
	  goto out;
	}
    }
  else
    {
      dmsg(D_PAYLOAD_CONNTRACK, "payload_tcp_dd_drop_hit CREATED " PAYLOAD_FMT_FULL, PAYLOAD_FMT_FULL_ARGS(pip, ptcp));
      pt = payload_tuple_new();
      pt->id = pt_id;
      hash_add_fast(pc->hash, bucket, &pt_id, hv, pt);
      payload_gc_adjust_timers(c);
    }

  /* Not found or new ... use next available slot */
  pt->tcp_seg[pt->tcp_seg_idx].tcph=tcph;
  pt->tcp_seg[pt->tcp_seg_idx].tcp_len=tcp_len;
  pt->tcp_seg[pt->tcp_seg_idx].hits=0;
  pt->tcp_seg_idx++;
  pt->tcp_seg_idx%=PAYLOAD_N_TCPSEGS;

done:
  /* refresh entry */
  pt->expires=now+pc->tcp_retrans;
  pt->last_used=now;

out:
  hash_bucket_unlock (bucket);
  return result;
}

/*
 * External interface: payload_tcp_retrans_drop() returns true
 * if this segment IS a dup (should be dropped by caller)
 */
bool
payload_tcp_retrans_drop (struct context *c, struct buffer *buf)
{
  const struct openvpn_iphdr *pip;
  int hlen;
  if (BLEN (buf) < (int) sizeof (struct openvpn_iphdr))
    return false;
  pip = (struct openvpn_iphdr *) BPTR (buf);

  hlen = OPENVPN_IPH_GET_LEN (pip->version_len);

  if (pip->protocol == OPENVPN_IPPROTO_TCP
      && ntohs (pip->tot_len) == BLEN (buf)
      && (ntohs (pip->frag_off) & OPENVPN_IP_OFFMASK) == 0
      && hlen <= BLEN (buf)
      && BLEN (buf) - hlen
      >= (int) sizeof (struct openvpn_tcphdr))
    {
      struct buffer newbuf = *buf;
      if (buf_advance (&newbuf, hlen))
	{
	  const struct openvpn_tcphdr *ptcp = (struct openvpn_tcphdr *) BPTR (&newbuf);
	  /*
	     dmsg(D_PAYLOAD_CONNTRACK, "payload_tcp_retrans_drop : " PAYLOAD_FMT_FULL, PAYLOAD_FMT_FULL_ARGS(pip, ptcp));
	     */
	  if (c->c2.payload_context && payload_tcp_dd_drop_hit(c, pip, ptcp, buf))
	    return true;
	}
    }
  return false;

}
/* Initialize payload (conntrack) hash table */
struct payload_context *
payload_new(int tcp_retrans)
{
  struct payload_context *pc;
  ASSERT(tcp_retrans);

  ALLOC_OBJ_CLEAR(pc, struct payload_context);
  pc->tcp_retrans = tcp_retrans;
  pc->hash= hash_init(PAYLOAD_HASH_SIZE, payload_tuple_hash_func, payload_tuple_compare_func);
  pc->gc.buckets_per_pass = payload_gc_buckets_per_pass(PAYLOAD_HASH_SIZE);
#if 0+WIP
  pc->payload_context>schedule = schedule_init();
  pc->payload_context>max_tcp_conns = 64; /* XXX: obviously must be configurable, test 4now */
#endif
  event_timeout_clear(&pc->gc.wakeup);
  return pc;
}
void payload_free(struct payload_context *pc)
{
  struct hash_iterator hi;
  struct hash_element *he;
  ASSERT(pc);

  event_timeout_clear(&pc->gc.wakeup);
  hash_iterator_init (pc->hash, &hi, true);
  while ((he = hash_iterator_next (&hi)))
    {
      struct payload_tuple *pt = (struct payload_tuple *) he->value;
      pt->deleted++;
    }
  hash_iterator_free (&hi);
  payload_gc_run(pc, -1, 0);
  free(pc);
}

#endif
/* 
vim: cino={.5s,\:.5s,+.5s,t0,g0,^-2,e-2,n-2,p2s,(0,=.5s:sw=4:ts=8:sts=4
 */

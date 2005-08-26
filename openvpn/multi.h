/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single TCP/UDP port, with support for SSL/TLS-based
 *             session authentication and key exchange,
 *             packet encryption, packet authentication, and
 *             packet compression.
 *
 *  Copyright (C) 2002-2005 OpenVPN Solutions LLC <info@openvpn.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2
 *  as published by the Free Software Foundation.
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

#ifndef MULTI_H
#define MULTI_H

#if P2MP_SERVER

#include "common.h"
#include "init.h"
#include "forward.h"
#include "mroute.h"
#include "fastlook.h"
#include "mbuf.h"
#include "list.h"
#include "schedule.h"
#include "pool.h"
#include "mudp.h"
#include "mtcp.h"
#include "perf.h"

/*
 * Walk (don't run) through the routing table,
 * deleting old entries, and possibly multi_instance
 * structs as well which have been marked for deletion.
 */
struct multi_reap
{
  int bucket_base;
  int buckets_per_pass;
  time_t last_call;
};

#ifdef FAST_IO
/*
 * Handle queuing of deferred MPP_PRE_SELECT actions
 */
struct multi_postprocess_defer_instance
{
  bool queued;
};
#endif

/*
 * One multi_instance object per client instance.
 */
struct multi_instance {
  struct schedule_entry se;    /* this must be the first element of the structure */
  struct gc_arena gc;
  //MUTEX_DEFINE (mutex);
  bool defined;
  bool halt;
  int refcount;
  time_t created;
  struct timeval wakeup;       /* absolute time */
  struct mroute_addr real;
  ifconfig_pool_handle vaddr_handle;
  const char *msg_prefix;

  /* queued outgoing data in Server/TCP mode */
  unsigned int tcp_rwflags;
  struct mbuf_set *tcp_link_out_deferred;
  bool socket_set_called;

  in_addr_t reporting_addr;       /* IP address shown in status listing */

#ifdef FAST_IO
  struct multi_postprocess_defer_instance mpdi;
#endif

  bool did_open_context;
  bool did_real_hash;
  bool did_iter;
  bool connection_established_flag;
  bool did_iroutes;

  struct context context;
};

#ifdef FAST_IO

struct multi_postprocess_defer
{
  int iter;

#ifdef FAST_IO_DEBUG
  int max;
#endif

  int n;
  struct multi_instance *list[MPD_MAX_QUEUED_INSTANCES];
};
#endif

/*
 * One multi_context object per server daemon thread.
 */
struct multi_context {
  struct hash *hash;   /* client instances indexed by real address */
  struct hash *vhash;  /* client instances indexed by virtual address */
  struct hash *iter;   /* like real address hash but optimized for iteration */
  struct schedule *schedule;
  struct mbuf_set *mbuf;
  struct multi_tcp *mtcp;
  struct ifconfig_pool *ifconfig_pool;
  struct frequency_limit *new_connection_limiter;
  struct mroute_helper *route_helper;
  struct multi_reap *reaper;
  struct mroute_addr local;
  bool enable_c2c;
  int max_clients;
  int tcp_queue_limit;
  int status_file_version;

#ifdef FAST_IO
  struct multi_postprocess_defer mpd;
#endif

#ifdef FAST_ADDR_LOOKUP
  struct fast_addr fast_addr;
  struct fast_addr fast_vaddr;
#endif

  time_t per_second_trigger;

  struct multi_instance *pending;
  struct multi_instance *earliest_wakeup;
  struct multi_instance **mpp_touched;

  bool io_order_toggle;

  struct context top;
};

/*
 * Host route
 */
struct multi_route
{
  struct mroute_addr addr;
  struct multi_instance *instance;

  /* must not collide with MGI_ (multi.c), or S_ (misc.h) flags */
# define MULTI_ROUTE_CACHE   (1<<8)
# define MULTI_ROUTE_AGEABLE (1<<9)
# define MULTI_ROUTE_MASK    (MULTI_ROUTE_CACHE|MULTI_ROUTE_AGEABLE)
# define MULTI_LOOKUP_CACHE  (1<<10)
  unsigned int flags;

  unsigned int cache_generation;
  time_t last_reference;
};

/*
 * top level function, called by openvpn.c
 */
void tunnel_server (struct context *top);

const char *multi_instance_string (const struct multi_instance *mi, bool null, struct gc_arena *gc);

void multi_bcast (struct multi_context *m,
		  const struct buffer *buf,
		  struct multi_instance *src,
		  const struct mroute_addr *srcaddr);

/*
 * Called by mtcp.c, mudp.c, or other (to be written) protocol drivers
 */

void multi_init (struct multi_context *m, struct context *t, bool tcp_mode);
void multi_uninit (struct multi_context *m);

void multi_top_init (struct multi_context *m, const struct context *top, const bool alloc_buffers);
void multi_top_free (struct multi_context *m);

struct multi_instance *multi_create_instance (struct multi_context *m, const struct mroute_addr *real);
void multi_close_instance (struct multi_context *m, struct multi_instance *mi, bool shutdown);

bool multi_process_timeout (struct multi_context *m, const unsigned int mpp_flags);

/* flags for multi_process_post */
#define MPP_PRE_SELECT                 (1<<0)
#define MPP_CLOSE_ON_SIGNAL            (1<<1)
#define MPP_RECORD_TOUCH               (1<<2)

#ifdef FAST_IO
#define MPP_POSTPROCESS_DEFER          (1<<3)
#endif

bool multi_process_post (struct multi_context *m, struct multi_instance *mi, const unsigned int flags);

bool multi_process_incoming_link (struct multi_context *m, struct multi_instance *instance, const unsigned int mpp_flags);
bool multi_process_incoming_tun (struct multi_context *m, const unsigned int mpp_flags);

void multi_process_drop_outgoing_tun (struct multi_context *m, const unsigned int mpp_flags);

void multi_print_status (struct multi_context *m, struct status_output *so, const int version);

struct multi_instance *multi_get_queue (struct mbuf_set *ms);

void multi_add_mbuf (struct multi_context *m,
		     struct multi_instance *mi,
		     struct mbuf_buffer *mb);

void multi_ifconfig_pool_persist (struct multi_context *m, bool force);

bool multi_process_signal (struct multi_context *m);

void multi_close_instance_on_signal (struct multi_context *m, struct multi_instance *mi);

void init_management_callback_multi (struct multi_context *m);
void uninit_management_callback_multi (struct multi_context *m);

/*
 * Is instance ready with respect to work thread locking?
 */
static inline bool
multi_instance_ready (const struct multi_instance *mi)
{
  return mi != NULL;
}

static inline struct multi_instance *
multi_instance_ref (struct multi_instance *mi)
{
  return mi;
}

/*
 * Return true if our output queue is not full
 */
static inline bool
multi_output_queue_ready (const struct multi_context *m,
			  const struct multi_instance *mi)
{
  if (mi->tcp_link_out_deferred)
    return mbuf_len (mi->tcp_link_out_deferred) <= m->tcp_queue_limit;
  else
    return true;
}

/*
 * Determine which instance has pending output
 * and prepare the output for sending in
 * the to_link buffer.
 */
static inline struct multi_instance *
multi_process_outgoing_link_pre (struct multi_context *m)
{
  struct multi_instance *mi = NULL;

  if (m->pending)
    mi = m->pending;
  else if (mbuf_defined (m->mbuf))
    mi = multi_get_queue (m->mbuf);
  return multi_instance_ref (mi);
}

/*
 * Instance reference counting
 */

static inline void
multi_instance_inc_refcount (struct multi_instance *mi)
{
  ++mi->refcount;
}

static inline void
multi_instance_dec_refcount (struct multi_instance *mi)
{
  if (--mi->refcount <= 0)
    {
      ASSERT (mi->halt);
      gc_free (&mi->gc);
      free (mi);
    }
}

static inline void
multi_route_del (struct multi_route *route)
{
  multi_instance_dec_refcount (route->instance);
  free (route);
}

static inline bool
multi_route_defined (const struct multi_context *m,
		     const struct multi_route *r)
{
  if (r->instance->halt || !multi_instance_ready (r->instance))
    return false;
  else if ((r->flags & MULTI_ROUTE_CACHE)
	   && r->cache_generation != m->route_helper->cache_generation)
    return false;
  else if ((r->flags & MULTI_ROUTE_AGEABLE)
	   && r->last_reference + m->route_helper->ageable_ttl_secs < now)
    return false;
  else
    return true;
}

/*
 * Set a msg() function prefix with our current client instance ID.
 */

static inline void
set_prefix (struct multi_instance *mi)
{
#ifdef MULTI_DEBUG_EVENT_LOOP
  if (mi->msg_prefix)
    printf ("[%s]\n", mi->msg_prefix);
#endif
  msg_set_prefix (mi->msg_prefix);
}

static inline void
clear_prefix (void)
{
#ifdef MULTI_DEBUG_EVENT_LOOP
  printf ("[NULL]\n");
#endif
  msg_set_prefix (NULL);
}

/*
 * Instance Reaper
 *
 * Reaper constants.  The reaper is the process where the virtual address
 * and virtual route hash table is scanned for dead entries which are
 * then removed.  The hash table could potentially be quite large, so we
 * don't want to reap in a single pass.
 */

#define REAP_MAX_WAKEUP   10  /* Do reap pass at least once per n seconds */
#define REAP_DIVISOR     256  /* How many passes to cover whole hash table */
#define REAP_MIN          16  /* Minimum number of buckets per pass */
#define REAP_MAX        1024  /* Maximum number of buckets per pass */

/*
 * Mark a cached host route for deletion after this
 * many seconds without any references.
 */
#define MULTI_CACHE_ROUTE_TTL 60

static inline void
multi_reap_process (struct multi_context *m)
{
  void multi_reap_process_dowork (struct multi_context *m);
  if (m->reaper->last_call != now)
    multi_reap_process_dowork (m);
}

static inline void
multi_process_per_second_timers (struct multi_context *m)
{
  if (m->per_second_trigger != now)
    {
      void multi_process_per_second_timers_dowork (struct multi_context *m);
      multi_process_per_second_timers_dowork (m);
      m->per_second_trigger = now;
    }
}

/*
 * Compute earliest timeout expiry from the set of
 * all instances.  Output:
 *
 * m->earliest_wakeup : instance needing the earliest service.
 * dest               : earliest timeout as a delta in relation
 *                      to current time.
 */
static inline void
multi_get_timeout (struct multi_context *m, struct timeval *dest)
{
  struct timeval tv, current;

  m->earliest_wakeup = (struct multi_instance *) schedule_get_earliest_wakeup (m->schedule, &tv);
  if (m->earliest_wakeup)
    {
      ASSERT (!gettimeofday (&current, NULL));
      tv_delta (dest, &current, &tv);
      if (dest->tv_sec >= REAP_MAX_WAKEUP)
	{
	  m->earliest_wakeup = NULL;
	  dest->tv_sec = REAP_MAX_WAKEUP;
	  dest->tv_usec = 0;
	}
    }
  else
    {
      dest->tv_sec = REAP_MAX_WAKEUP;
      dest->tv_usec = 0;
    }
}

/*
 * Send a packet to TUN/TAP interface.
 */
static inline bool
multi_process_outgoing_tun (struct multi_context *m, const unsigned int mpp_flags)
{
  struct multi_instance *mi = m->pending;
  bool ret = true;

  ASSERT (mi);
#ifdef MULTI_DEBUG_EVENT_LOOP
  printf ("%s -> TUN len=%d\n",
	  id(mi),
	  mi->context.c2.to_tun.len);
#endif
  set_prefix (mi);
  process_outgoing_tun (&mi->context);
  ret = multi_process_post (m, mi, mpp_flags);
  clear_prefix ();
  return ret;
}

static inline bool
multi_process_outgoing_link_dowork (struct multi_context *m, struct multi_instance *mi, const unsigned int mpp_flags)
{
  bool ret = true;
  set_prefix (mi);
  process_outgoing_link (&mi->context);
  ret = multi_process_post (m, mi, mpp_flags);
  clear_prefix ();
  return ret;
}

/*
 * Check for signals.
 */
#define MULTI_CHECK_SIG(m) EVENT_LOOP_CHECK_SIGNAL (&(m)->top, multi_process_signal, (m))

/*
 * Set currently pending instance
 */
static inline void
multi_set_pending (struct multi_context *m, struct multi_instance *mi)
{
  m->pending = multi_instance_ref (mi);
}

#ifdef FAST_IO

/*
 * Handle queuing of deferred MPP_PRE_SELECT actions
 * inline functions
 */

static inline void
multi_postprocess_defer_reset (struct multi_context *m)
{
  void multi_postprocess_defer_reset_dowork (struct multi_context *m);
  m->mpd.iter = 0;
  if (m->mpd.n > 0)
    multi_postprocess_defer_reset_dowork (m);
}

static inline void
multi_postprocess_defer_add (struct multi_context *m, struct multi_instance *mi)
{
#ifdef FAST_IO_DEBUG
  void multi_postprocess_defer_max_exceeded (struct multi_postprocess_defer *mpd);
#endif

  ++m->mpd.iter;

#ifdef FAST_IO_DEBUG
  if (m->mpd.iter > m->mpd.max)
    multi_postprocess_defer_max_exceeded (&m->mpd);
#endif

  if (!mi->mpdi.queued)
    {
      ASSERT (m->mpd.n < MPD_MAX_QUEUED_INSTANCES);
      m->mpd.list[m->mpd.n++] = mi;
      mi->mpdi.queued = true;
    }
}

static inline struct multi_instance *
multi_postprocess_defer_get (struct multi_context *m)
{
  struct multi_instance *mi = NULL;
  if (m->mpd.n > 0)
    {
      mi = m->mpd.list[--m->mpd.n];
      mi->mpdi.queued = false;
    }
  return mi;
}

static inline bool 
multi_postprocess_defer_must_flush (struct multi_context *m)
{
  return (!m->top.c2.event_set_status_hint)
    || (m->mpd.iter >= MPD_MAX_ITERATIONS)
    || (m->mpd.n == MPD_MAX_QUEUED_INSTANCES);
}
#endif /* FAST_IO */

#endif /* P2MP_SERVER */
#endif /* MULTI_H */

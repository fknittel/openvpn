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

#ifdef WIN32
#include "config-win32.h"
#else
#include "config.h"
#endif

#include "syshead.h"

#if P2MP_SERVER

#include "multi.h"

#include "memdbg.h"

#include "forward-inline.h"
#include "fastlook-inline.h"

/*
 * Get a client instance based on real address.  If
 * the instance doesn't exist, create it while
 * maintaining real address hash table atomicity.
 */

struct multi_instance *
multi_get_create_instance_udp (struct multi_context *m)
{
  struct multi_instance *mi = NULL;
  struct mroute_addr real;

  if (mroute_extract_openvpn_sockaddr (&real, &m->top.c2.from_addr, true))
    {
#ifdef FAST_ADDR_LOOKUP
      if ((mi = multi_fast_addr_lookup (&m->fast_addr, &real)) != NULL)
	return mi;
#endif
      {
	struct gc_arena gc = gc_new ();
	struct hash *hash = m->hash;
	const uint32_t hv = hash_value (hash, &real);
	struct hash_bucket *bucket = hash_bucket (hash, hv);
	struct hash_element *he;
  
	hash_bucket_lock (bucket);
	he = hash_lookup_fast (hash, bucket, &real, hv);

	if (he)
	  {
	    mi = (struct multi_instance *) he->value;
	  }
	else
	  {
	    if (!m->top.c2.tls_auth_standalone
		|| tls_pre_decrypt_lite (m->top.c2.tls_auth_standalone, &m->top.c2.from_addr, &m->top.c2.buf))
	      {
		if (frequency_limit_event_allowed (m->new_connection_limiter))
		  {
		    mi = multi_create_instance (m, &real);
		    if (mi)
		      {
			hash_add_fast (hash, bucket, &mi->real, hv, mi);
			mi->did_real_hash = true;
		      }
		  }
		else
		  {
		    msg (D_MULTI_ERRORS,
			 "MULTI: Connection from %s would exceed new connection frequency limit as controlled by --connect-freq",
			 mroute_addr_print (&real, &gc));
		  }
	      }
	  }

	hash_bucket_unlock (bucket);

#ifdef ENABLE_DEBUG
	if (check_debug_level (D_MULTI_DEBUG))
	  {
	    const char *status;

	    if (he && mi)
	      status = "[succeeded]";
	    else if (!he && mi)
	      status = "[created]";
	    else
	      status = "[failed]";
	
	    dmsg (D_MULTI_DEBUG, "GET INST BY REAL: %s %s",
		  mroute_addr_print (&real, &gc),
		  status);
	  }
#endif

#ifdef FAST_ADDR_LOOKUP
	multi_fast_addr_save (&m->fast_addr, &real, mi);
#endif

	gc_free (&gc);
      }
    }

  ASSERT (!(mi && mi->halt));
  return mi;
}

/*
 * Send a packet to TCP/UDP socket.
 */
static void
multi_process_outgoing_link (struct multi_context *m, const unsigned int mpp_flags)
{
  struct multi_instance *mi = multi_process_outgoing_link_pre (m);
  if (mi)
    multi_process_outgoing_link_dowork (m, mi, mpp_flags);
}

/*
 * Return the io_wait() flags appropriate for
 * a point-to-multipoint tunnel.
 */
static inline unsigned int
p2mp_iow_flags (const struct multi_context *m)
{
  unsigned int flags = IOW_WAIT_SIGNAL;
  if (m->pending)
    {
      if (TUN_OUT (&m->pending->context))
	flags |= IOW_TO_TUN;
      if (LINK_OUT (&m->pending->context))
	flags |= IOW_TO_LINK;
    }
  else if (mbuf_defined (m->mbuf))
    flags |= IOW_MBUF;
  else
    flags |= IOW_READ;

  return flags;
}

/*
 * Process an I/O event.
 */

static inline void
do_socket_read (struct multi_context *m, const unsigned int mpp_flags)
{
  if (!read_incoming_link (&m->top, &m->top.c2.from_addr))
    multi_process_incoming_link (m, NULL, mpp_flags);
}

static inline void
do_tun_read (struct multi_context *m, const unsigned int mpp_flags)
{
  if (!read_incoming_tun (&m->top))
    multi_process_incoming_tun (m, mpp_flags);
}

static void
multi_process_io_udp (struct multi_context *m)
{
  const unsigned int status = m->top.c2.event_set_status;
  unsigned int mpp_flags;

#ifdef MULTI_DEBUG_EVENT_LOOP
  char buf[16];
  buf[0] = 0;
  if (status & SOCKET_READ)
    strcat (buf, "SR/");
  else if (status & SOCKET_WRITE)
    strcat (buf, "SW/");
  else if (status & TUN_READ)
    strcat (buf, "TR/");
  else if (status & TUN_WRITE)
    strcat (buf, "TW/");
  printf ("IO %s\n", buf);
#endif

#ifdef ENABLE_MANAGEMENT
  if (status & (MANAGEMENT_READ|MANAGEMENT_WRITE))
    {
      ASSERT (management);
      management_io (management);
    }
#endif

  /*
   * When fast_io is set, we wait until input buffers are depleted before
   * performing an MPP_PRE_SELECT action on all touched instances.
   */
#ifdef FAST_IO
  if (m->top.c2.default_iow_flags & IOW_FAST_IO)
    mpp_flags = MPP_POSTPROCESS_DEFER | MPP_CLOSE_ON_SIGNAL;
  else
#endif
    mpp_flags = MPP_PRE_SELECT | MPP_CLOSE_ON_SIGNAL;

  /*
   * Process I/O
   */
  if (status & (SOCKET_WRITE|TUN_WRITE))
    {
      if (status & SOCKET_WRITE)
	{
	  multi_process_outgoing_link (m, mpp_flags);
	}
      else if (status & TUN_WRITE)
	{
	  multi_process_outgoing_tun (m, mpp_flags);
	}
    }
  else if (status & (SOCKET_READ|TUN_READ))
    {
      if (!m->io_order_toggle)
	{
	  if (status & SOCKET_READ)
	    do_socket_read (m, mpp_flags);
	  else if (status & TUN_READ)
	    do_tun_read (m, mpp_flags);
	  m->io_order_toggle = true;
	}
      else
	{
	  if (status & TUN_READ)
	    do_tun_read (m, mpp_flags);
	  else if (status & SOCKET_READ)
	    do_socket_read (m, mpp_flags);
	  m->io_order_toggle = false;
	}
    }

#ifdef FAST_IO
  /*
   * If must_flush condition, do MPP_PRE_SELECT action on deferred instances.
   */
  if ((m->top.c2.default_iow_flags & IOW_FAST_IO) && !m->pending && multi_postprocess_defer_must_flush (m))
    {
      struct multi_instance *mi;
      const unsigned int defer_mpp_flags = MPP_PRE_SELECT | MPP_CLOSE_ON_SIGNAL;
      while ((mi = multi_postprocess_defer_get (m)))
	{
	  while (true)
	    {
	      multi_process_post (m, mi, defer_mpp_flags);
	      if (!m->pending)
		break;
	      if (LINK_OUT (&m->pending->context))
		multi_process_outgoing_link (m, defer_mpp_flags);
	      else if (TUN_OUT (&m->pending->context))
		multi_process_outgoing_tun (m, defer_mpp_flags);
	    }
	}
      multi_postprocess_defer_reset (m);
    }
#endif
}

/*
 * Top level event loop for single-threaded operation.
 * UDP mode.
 */

static int
tunnel_server_udp_event_loop (void *arg)
{
  struct multi_context *m = (struct multi_context *) arg;
  int ret = 0;

  while (true)
    {
      unsigned int io_flags;

      perf_push (PERF_EVENT_LOOP);

      io_flags = p2mp_iow_flags (m);

      if (!is_io_wait_fast_path (&m->top, io_flags))
	{
	  /* set up and do the io_wait() */
	  if (!IS_SIG (&m->top))
	    {
	      multi_get_timeout (m, &m->top.c2.timeval);
	      io_wait_slow (&m->top, io_flags);
	    }
	  MULTI_CHECK_SIG (m);
	}

      /* timeout? */
      if (m->top.c2.event_set_status != ES_TIMEOUT)
	{
	  /* process I/O */
	  multi_process_io_udp (m);
	  MULTI_CHECK_SIG (m);
	}
      else
	{
	  multi_process_timeout (m, MPP_PRE_SELECT|MPP_CLOSE_ON_SIGNAL);
	}
      
      /* check on status of coarse timers */
      if (!m->pending && !ess_hint (&m->top))
	multi_process_per_second_timers (m);

      perf_pop ();
    }

  return ret;
}

void
tunnel_server_udp (struct context *top)
{
  struct multi_context multi;

  top->mode = CM_TOP;
  context_clear_2 (top);

  /* initialize top-tunnel instance */
  init_instance_handle_signals (top, top->es, CC_HARD_USR1_TO_HUP);
  if (IS_SIG (top))
    return;
  
  /* initialize global multi_context object */
  multi_init (&multi, top, false);

  /* initialize our cloned top object */
  multi_top_init (&multi, top, true);

  /* initialize management interface */
  init_management_callback_multi (&multi);

  /* finished with initialization */
  initialization_sequence_completed (top, ISC_SERVER); /* --mode server --proto udp */

  /* per-packet event loop */
  tunnel_server_udp_event_loop (&multi);

  /* shut down management interface */
  uninit_management_callback_multi (&multi);

  /* save ifconfig-pool */
  multi_ifconfig_pool_persist (&multi, true);

  /* tear down tunnel instance (unless --persist-tun) */
  multi_uninit (&multi);
  multi_top_free (&multi);
  close_instance (top);
}

#endif

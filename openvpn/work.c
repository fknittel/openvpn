/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single UDP port, with support for SSL/TLS-based
 *             session authentication and key exchange,
 *             packet encryption, packet authentication, and
 *             packet compression.
 *
 *  Copyright (C) 2002-2004 James Yonan <jim@yonan.net>
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

#include "syshead.h"

#ifdef USE_PTHREAD

#include "error.h"
#include "buffer.h"
#include "fdmisc.h"
#include "work.h"

#include "memdbg.h"

static void *
thread_func (void *arg)
{
  /*struct work_thread *wt = (struct work_thread *) arg;*/

  return NULL;
}

struct work_thread *
work_thread_init (const int n_threads, const int nice_work)
{
  struct work_thread *wt;

  ASSERT (n_threads >= 2);

  if (n_threads > 2)
    msg (M_INFO, "NOTE: currently a maximum of 2 threads is supported");

  ALLOC_OBJ_CLEAR (wt, struct work_thread);
  wt->n_threads = n_threads;
  wt->nice_work = nice_work;
  
  /* initialize pthread */
  openvpn_thread_init ();

  /*
   * Make a socket for foreground and background threads
   * to communicate.  The background thread will set its
   * end to blocking, while the foreground will set its
   * end to non-blocking.
   */
  if (socketpair (PF_UNIX, SOCK_DGRAM, 0, wt->sd) == -1)
    msg (M_ERR, "socketpair call failed");

  /* set socket properties */
  set_nonblock (wt->sd[MAIN_THREAD]);
  set_cloexec (wt->sd[MAIN_THREAD]);
  set_cloexec (wt->sd[WORKER_THREAD]);

  wt->thread_id = openvpn_thread_create (thread_func, (void*)wt);

  return wt;
}

void
work_thread_close (struct work_thread *wt)
{
  if (wt)
    {
      openvpn_thread_join (wt->thread_id);
      openvpn_thread_cleanup ();
      close (wt->sd[MAIN_THREAD]);
      close (wt->sd[WORKER_THREAD]);
      free (wt);
    }
}

void
work_thread_enable (struct work_thread *wt, void *arg, work_thread_event_loop_t event_loop)
{
  wt->event_loop_arg = arg;
  wt->event_loop = event_loop;
}

void
work_thread_disable (struct work_thread *wt)
{
  wt->event_loop_arg = NULL;
  wt->event_loop = NULL;
}

#else
static void dummy(void) {}
#endif

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

/*
 * Socket read/write functions.
 */

static int
recv_uchar (const int fd)
{
  unsigned char c;
  const ssize_t size = read (fd, &c, sizeof (c));
  if (size == sizeof (c))
    return c;
  else
    return -1;
}

static int
send_uchar (const int fd, const int value)
{
  ssize_t size;
  unsigned char c;

  ASSERT (value >= 0 && value <= 0xFF);
  c = (unsigned char) value;
  size = write (fd, &c, sizeof (c));
  if (size == sizeof (c))
    return (int) size;
  else
    return -1;
}

/*
 * This is the work thread function
 */
static void *
thread_func (void *arg)
{
  struct work_thread *wt = (struct work_thread *) arg;

  msg (D_WORK_THREAD, "WORK THREAD: starting");

  /* change thread priority if requested */
  set_nice (wt->nice_work);

  msg (D_WORK_THREAD, "WORK THREAD: exiting");

  return NULL;
}

struct work_thread *
work_thread_init (const int n_threads, const int nice_work)
{
  struct work_thread *wt;

  ASSERT (n_threads >= 2);

  ALLOC_OBJ_CLEAR (wt, struct work_thread);
  wt->n_threads = n_threads;
  
  /* initialize pthread */
  openvpn_thread_init ();

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

      free (wt);
    }
}

#else
static void dummy(void) {}
#endif

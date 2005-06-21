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

#ifndef WORK_H
#define WORK_H

#include "thread.h"
#include "ssl.h"
#include "misc.h"
#include "plugin.h"
#include "event.h"

#ifdef USE_PTHREAD

/*
 * Object that represents a work thread.
 */
struct work_thread
{
  int n_threads;
  int nice_work;
  openvpn_thread_t thread_id;
  int sd[2];                   /* PF_UNIX/SOCK_DGRAM communication handles */
  bool brk;
};

struct work_thread *work_thread_init (const int n_threads, const int nice_work);

void work_thread_close (struct work_thread *wt);

#endif /* USE_PTHREAD */
#endif /* WORK_H */

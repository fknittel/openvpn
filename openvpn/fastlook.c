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

#include "fastlook-inline.h"

struct multi_instance *
multi_fast_addr_lookup_foo (struct fast_addr *fa, const struct mroute_addr *addr)
{
  return multi_fast_addr_lookup (fa, addr);
}

void
multi_fast_addr_save_foo (struct fast_addr *fa, const struct mroute_addr *addr, struct multi_instance *mi)
{
  multi_fast_addr_save (fa, addr, mi);
}

bool
memeq6_foo (const void *m1, const void *m2)
{
  return memeq6 (m1, m2);

}

void
multi_fast_addr_reset (struct fast_addr *fa)
{
  CLEAR (*fa);
}

#endif

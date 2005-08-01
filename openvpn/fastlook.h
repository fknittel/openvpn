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

#ifndef FASTLOOK_H
#define FASTLOOK_H

#include "common.h"

#if P2MP_SERVER && defined(FAST_ADDR_LOOKUP)

struct multi_instance;

/*
 * Extremely fast small-sized cache for IPv4 -> struct multi_instance lookup
 */

#if !(FAST_ADDR_LOOKUP == 4 || FAST_ADDR_LOOKUP == 8)
#error FAST_ADDR_LOOKUP must be 4, 8, or undefined
#endif

struct fast_addr_4
{
  struct multi_instance *mi;
  in_addr_t addr4;
};

struct fast_addr_6 /* number of bytes in address, not IP version # */
{
  struct multi_instance *mi;
  uint8_t addr6[6];
  uint8_t dummy[6];
};

struct fast_addr
{
  int type;
  unsigned int index;
  union {
    struct fast_addr_4 list4[FAST_ADDR_LOOKUP];
    struct fast_addr_6 list6[FAST_ADDR_LOOKUP];
  } u;
};

void multi_fast_addr_reset (struct fast_addr *fa);

#endif

#endif

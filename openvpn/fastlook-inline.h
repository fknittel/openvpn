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

#ifndef FASTLOOK_INLINE_H
#define FASTLOOK_INLINE_H

#include "common.h"

#if P2MP_SERVER && defined(FAST_ADDR_LOOKUP)

#include "mroute.h"
#include "fastlook.h"

static inline bool
memeq6 (const void *m1, const void *m2)
{
#if USE_MEMCMP_6
  return memcmp (m1, m2, 6) == 0;
#else
  return *((uint32_t*)m1) == *((uint32_t*)m2)
    && *((uint16_t*)m1+2) == *((uint16_t*)m2+2);
#endif
}

#define MFAL_TEST4  { if (fa->u.list4[i].addr4 == a) return fa->u.list4[i].mi; }
#define MFAL_TEST6  { if (memeq6(fa->u.list6[i].addr6, addr->addr)) return fa->u.list6[i].mi; }
#define MFAL_INCR   { i = ((i + 1) & (FAST_ADDR_LOOKUP-1)); }

static inline struct multi_instance *
multi_fast_addr_lookup (struct fast_addr *fa, const struct mroute_addr *addr)
{
  if (fa->type == addr->type)
    {
      unsigned int i = fa->index;
      switch (fa->type)
	{
	case MR_ADDR_IPV4:
	  {
	    in_addr_t a;
	    memcpy (&a, addr->addr, 4);

#if FAST_ADDR_LOOKUP >= 4
	    MFAL_TEST4;
	    MFAL_INCR;
	    MFAL_TEST4;
	    MFAL_INCR;
	    MFAL_TEST4;
	    MFAL_INCR;
	    MFAL_TEST4;
#endif
#if FAST_ADDR_LOOKUP >= 8
	    MFAL_INCR;
	    MFAL_TEST4;
	    MFAL_INCR;
	    MFAL_TEST4;
	    MFAL_INCR;
	    MFAL_TEST4;
	    MFAL_INCR;
	    MFAL_TEST4;
#endif
	  }
	  break;
	case MR_ADDR_IPV4|MR_WITH_PORT:
	case MR_ADDR_ETHER:
	  {
#if FAST_ADDR_LOOKUP >= 4
	    MFAL_TEST6;
	    MFAL_INCR;
	    MFAL_TEST6;
	    MFAL_INCR;
	    MFAL_TEST6;
	    MFAL_INCR;
	    MFAL_TEST6;
#endif
#if FAST_ADDR_LOOKUP >= 8
	    MFAL_INCR;
	    MFAL_TEST6;
	    MFAL_INCR;
	    MFAL_TEST6;
	    MFAL_INCR;
	    MFAL_TEST6;
	    MFAL_INCR;
	    MFAL_TEST6;
#endif
	  }
	  break;
	}
    }
  return NULL;
}

#undef MFAL_TEST4
#undef MFAL_TEST6
#undef MFAL_INCR

static inline void
multi_fast_addr_save (struct fast_addr *fa, const struct mroute_addr *addr, struct multi_instance *mi)
{
  if (fa->type == addr->type)
    {
      switch (fa->type)
	{
	case MR_ADDR_IPV4:
	  {
	    struct fast_addr_4 *e;
	    --fa->index;
	    fa->index &= (FAST_ADDR_LOOKUP-1);
	    e = &fa->u.list4[fa->index];
	    memcpy (&e->addr4, addr->addr, 4);
	    e->mi = mi;
	  }
	  break;
	case MR_ADDR_IPV4|MR_WITH_PORT:
	case MR_ADDR_ETHER:
	  {
	    struct fast_addr_6 *e;
	    --fa->index;
	    fa->index &= (FAST_ADDR_LOOKUP-1);
	    e = &fa->u.list6[fa->index];
	    memcpy (&e->addr6, addr->addr, 6);
	    e->mi = mi;
	  }
	  break;
	}
    }
}

#endif
#endif

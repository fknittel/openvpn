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

#ifndef MROUTE_H
#define MROUTE_H

#if P2MP_SERVER

#include "buffer.h"
#include "list.h"
#include "route.h"

#define IP_MCAST_SUBNET_MASK  ((in_addr_t)240<<24)
#define IP_MCAST_NETWORK      ((in_addr_t)224<<24)

/* Return status values for mroute_extract_addr_from_packet */
#define MROUTE_EXTRACT_SUCCEEDED (1<<1)
#define MROUTE_EXTRACT_BCAST     (1<<2)
#define MROUTE_EXTRACT_MCAST     (1<<3)
#define MROUTE_EXTRACT_IGMP      (1<<4)

/*
 * Choose the largest address possible with
 * any of our supported types, which is IPv6
 * with port number.
 */
#define MR_MAX_ADDR_LEN 20

/*
 * Address Types
 */
#define MR_ADDR_NONE             0
#define MR_ADDR_ETHER            1
#define MR_ADDR_IPV4             2
#define MR_ADDR_IPV6             3
#define MR_ADDR_UNIX             4
#define MR_ADDR_MASK             7

/* Address type mask indicating that port # is part of address */
#define MR_WITH_PORT             8

/* Address type mask indicating that netbits is part of address */
#define MR_WITH_NETBITS          16

struct mroute_addr {
  uint8_t len;      /* length of address */
  uint8_t unused;
  uint8_t type;     /* MR_ADDR/MR_WITH flags */
  uint8_t netbits;  /* number of bits in network part of address,
		       valid if MR_WITH_NETBITS is set */
  uint8_t addr[MR_MAX_ADDR_LEN];  /* actual address */
};

/*
 * Number of bits in an address.  Should be raised for
 * IPv6.
 */
#define MR_HELPER_NET_LEN 32

/*
 * Used to help maintain CIDR routing table.
 */
struct mroute_helper {
  /*MUTEX_DEFINE (mutex);*/
  unsigned int cache_generation; /* incremented when route added */
  int ageable_ttl_secs;          /* host route cache entry time-to-live*/
  int n_net_len;                 /* length of net_len array */
  uint8_t net_len[MR_HELPER_NET_LEN];      /* CIDR netlengths in descending order */
  int net_len_refcount[MR_HELPER_NET_LEN]; /* refcount of each netlength */
};

unsigned int mroute_extract_addr_from_packet (struct mroute_addr *src,
					      struct mroute_addr *dest,
					      struct buffer *buf,
					      int tunnel_type,
					      struct mroute_addr *esrc,
					      struct mroute_addr *edest);

bool mroute_extract_sockaddr_in (struct mroute_addr *addr,
				 const struct sockaddr_in *saddr,
				 bool use_port);

struct openvpn_sockaddr;
int mroute_extract_openvpn_sockaddr(struct mroute_addr *addr, 
				 const struct openvpn_sockaddr *osaddr,
				 int use_port);

bool mroute_learnable_address (const struct mroute_addr *addr);

uint32_t mroute_addr_hash_function (const void *key, uint32_t iv);
bool mroute_addr_compare_function (const void *key1, const void *key2);

void mroute_addr_init (struct mroute_addr *addr);

const char *mroute_addr_print (const struct mroute_addr *ma,
			       struct gc_arena *gc);

void mroute_addr_mask_host_bits (struct mroute_addr *ma);

struct mroute_helper *mroute_helper_init (int ageable_ttl_secs);
void mroute_helper_free (struct mroute_helper *mh);
void mroute_helper_add_iroute (struct mroute_helper *mh, const struct iroute *ir);
void mroute_helper_del_iroute (struct mroute_helper *mh, const struct iroute *ir);

static inline void
mroute_helper_lock (struct mroute_helper *mh)
{
  /*mutex_lock (&mh->mutex);*/
}

static inline void
mroute_helper_unlock (struct mroute_helper *mh)
{
  /*mutex_unlock (&mh->mutex);*/
}

static inline bool
mroute_addr_equal (const struct mroute_addr *a1, const struct mroute_addr *a2)
{
  if (a1->type != a2->type)
    return false;
  if (a1->netbits != a2->netbits)
    return false;
  if (a1->len != a2->len)
    return false;
  return memcmp (a1->addr, a2->addr, a1->len) == 0;
}

static inline const uint8_t *
mroute_addr_hash_ptr (const struct mroute_addr *a)
{
  /* NOTE: depends on ordering of struct mroute_addr */
  return (uint8_t *) &a->type;
}

static inline uint32_t
mroute_addr_hash_len (const struct mroute_addr *a)
{
  return (uint32_t) a->len + 2;
}

static inline void
mroute_extract_in_addr_t (struct mroute_addr *dest, const in_addr_t src)
{
  dest->type = MR_ADDR_IPV4;
  dest->netbits = 0;
  dest->len = 4;
  *(in_addr_t*)dest->addr = htonl (src);
}

#endif /* P2MP_SERVER */
#endif /* MROUTE_H */

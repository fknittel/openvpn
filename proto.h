/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single TCP/UDP port, with support for SSL/TLS-based
 *             session authentication and key exchange,
 *             packet encryption, packet authentication, and
 *             packet compression.
 *
 *  Copyright (C) 2002-2009 OpenVPN Technologies, Inc. <sales@openvpn.net>
 *  Copyright (C) 2010      Fabian Knittel <fabian.knittel@avona.com>
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

#ifndef PROTO_H
#define PROTO_H

#include "common.h"
#include "buffer.h"

#pragma pack(1)

/*
 * Tunnel types
 */
#define DEV_TYPE_UNDEF 0
#define DEV_TYPE_NULL  1
#define DEV_TYPE_TUN   2    /* point-to-point IP tunnel */
#define DEV_TYPE_TAP   3    /* ethernet (802.3) tunnel */

/* TUN topologies */

#define TOP_UNDEF   0
#define TOP_NET30   1
#define TOP_P2P     2
#define TOP_SUBNET  3

/*
 * IP and Ethernet protocol structs.  For portability,
 * OpenVPN needs its own definitions of these structs, and
 * names have been adjusted to avoid collisions with
 * native structs.
 */

#define OPENVPN_ETH_ALEN 6            /* ethernet address length */
struct openvpn_ethhdr 
{
  uint8_t dest[OPENVPN_ETH_ALEN];     /* destination ethernet addr */
  uint8_t source[OPENVPN_ETH_ALEN];   /* source ethernet addr	*/

# define OPENVPN_ETH_P_IPV4   0x0800  /* IPv4 protocol */
# define OPENVPN_ETH_P_IPV6   0x86DD  /* IPv6 protocol */
# define OPENVPN_ETH_P_ARP    0x0806  /* ARP protocol */
# define OPENVPN_ETH_P_8021Q  0x8100  /* 802.1Q protocol */
  uint16_t proto;                     /* packet type ID field */
};

struct openvpn_8021qhdr
{
  uint8_t dest[OPENVPN_ETH_ALEN];     /* destination ethernet addr */
  uint8_t source[OPENVPN_ETH_ALEN];   /* source ethernet addr	*/

  uint16_t tpid;                      /* 802.1Q Tag Protocol Identifier */
# define OPENVPN_8021Q_MASK_VID htons (0x0FFF) /* mask VID out of pcp_cfi_vid */
# define OPENVPN_8021Q_MASK_PCP htons (0xE000) /* mask PCP out of pcp_cfi_vid */
# define OPENVPN_8021Q_MASK_CFI htons (0x1000) /* mask CFI out of pcp_cfi_vid */
  uint16_t pcp_cfi_vid;               /* bit fields, see IEEE 802.1Q */
  uint16_t proto;                     /* contained packet type ID field */
};

/*
 * Size difference between a regular Ethernet II header and an Ethernet II
 * header with additional IEEE 802.1Q tagging.
 */
#define SIZE_ETH_TO_8021Q_HDR (sizeof (struct openvpn_8021qhdr) - sizeof (struct openvpn_ethhdr))

struct openvpn_arp {
# define ARP_MAC_ADDR_TYPE 0x0001
  uint16_t mac_addr_type;       /* 0x0001 */

  uint16_t proto_addr_type;     /* 0x0800 */
  uint8_t  mac_addr_size;       /* 0x06 */
  uint8_t  proto_addr_size;     /* 0x04 */

# define ARP_REQUEST 0x0001
# define ARP_REPLY   0x0002
  uint16_t arp_command;         /* 0x0001 for ARP request, 0x0002 for ARP reply */

  uint8_t   mac_src[OPENVPN_ETH_ALEN];
  in_addr_t ip_src;
  uint8_t   mac_dest[OPENVPN_ETH_ALEN];
  in_addr_t ip_dest;
};

struct openvpn_iphdr {
# define OPENVPN_IPH_GET_VER(v) (((v) >> 4) & 0x0F)
# define OPENVPN_IPH_GET_LEN(v) (((v) & 0x0F) << 2)
  uint8_t    version_len;

  uint8_t    tos;
  uint16_t   tot_len;
  uint16_t   id;

# define OPENVPN_IP_OFFMASK 0x1fff
  uint16_t   frag_off;

  uint8_t    ttl;

# define OPENVPN_IPPROTO_IGMP 2 /* IGMP protocol */
# define OPENVPN_IPPROTO_TCP  6 /* TCP protocol */
# define OPENVPN_IPPROTO_UDP 17 /* UDP protocol */
  uint8_t    protocol;

  uint16_t   check;
  uint32_t   saddr;
  uint32_t   daddr;
  /*The options start here. */
};

/*
 * UDP header
 */
struct openvpn_udphdr {
  uint16_t   source;
  uint16_t   dest;
  uint16_t   len;
  uint16_t   check;
};

/*
 * TCP header, per RFC 793.
 */
struct openvpn_tcphdr {
  uint16_t      source;    /* source port */
  uint16_t      dest;      /* destination port */
  uint32_t      seq;       /* sequence number */
  uint32_t      ack_seq;   /* acknowledgement number */

# define OPENVPN_TCPH_GET_DOFF(d) (((d) & 0xF0) >> 2)
  uint8_t       doff_res;

# define OPENVPN_TCPH_FIN_MASK (1<<0)
# define OPENVPN_TCPH_SYN_MASK (1<<1)
# define OPENVPN_TCPH_RST_MASK (1<<2)
# define OPENVPN_TCPH_PSH_MASK (1<<3)
# define OPENVPN_TCPH_ACK_MASK (1<<4)
# define OPENVPN_TCPH_URG_MASK (1<<5)
# define OPENVPN_TCPH_ECE_MASK (1<<6)
# define OPENVPN_TCPH_CWR_MASK (1<<7)
  uint8_t       flags;

  uint16_t      window;
  uint16_t      check;
  uint16_t      urg_ptr;
};

#define	OPENVPN_TCPOPT_EOL     0
#define	OPENVPN_TCPOPT_NOP     1
#define	OPENVPN_TCPOPT_MAXSEG  2
#define OPENVPN_TCPOLEN_MAXSEG 4

#pragma pack()

/*
 * The following macro is used to update an
 * internet checksum.  "acc" is a 32-bit
 * accumulation of all the changes to the
 * checksum (adding in old 16-bit words and
 * subtracting out new words), and "cksum"
 * is the checksum value to be updated.
 */
#define ADJUST_CHECKSUM(acc, cksum) { \
  (acc) += (cksum); \
  if ((acc) < 0) { \
    (acc) = -(acc); \
    (acc) = ((acc) >> 16) + ((acc) & 0xffff); \
    (acc) += (acc) >> 16; \
    (cksum) = (uint16_t) ~(acc); \
  } else { \
    (acc) = ((acc) >> 16) + ((acc) & 0xffff); \
    (acc) += (acc) >> 16; \
    (cksum) = (uint16_t) (acc); \
  } \
}

/*
 * We are in a "liberal" position with respect to MSS,
 * i.e. we assume that MSS can be calculated from MTU
 * by subtracting out only the IP and TCP header sizes
 * without options.
 *
 * (RFC 879, section 7).
 */
#define MTU_TO_MSS(mtu) (mtu - sizeof(struct openvpn_iphdr) \
                             - sizeof(struct openvpn_tcphdr))

/*
 * If raw tunnel packet is IPv4, return true and increment
 * buffer offset to start of IP header.
 */
bool is_ipv4 (int tunnel_type, struct buffer *buf);

#ifdef PACKET_TRUNCATION_CHECK
void ipv4_packet_size_verify (const uint8_t *data,
			      const int size,
			      const int tunnel_type,
			      const char
			      *prefix,
			      counter_type *errors);
#endif

#ifdef ENABLE_VLAN_TAGGING
# define OPENVPN_8021Q_MIN_VID 1
# define OPENVPN_8021Q_MAX_VID 4094

/*
 * Retrieve the Priority Code Point (PCP) from the IEEE 802.1Q header.
 *
 * @param hdr Pointer to the Ethernet header with IEEE 802.1Q tagging.
 * @return    Returns the PCP in host byte order.
 */
static inline uint16_t
vlanhdr_get_pcp (const struct openvpn_8021qhdr *hdr)
{
  return ntohs (hdr->pcp_cfi_vid & OPENVPN_8021Q_MASK_PCP);
}
/*
 * Retrieve the Canonical Format Indicator (CFI) from the IEEE 802.1Q header.
 *
 * @param hdr Pointer to the Ethernet header with IEEE 802.1Q tagging.
 * @return    Returns the CFI in host byte order.
 */
static inline uint16_t
vlanhdr_get_cfi (const struct openvpn_8021qhdr *hdr)
{
  return ntohs (hdr->pcp_cfi_vid & OPENVPN_8021Q_MASK_CFI);
}
/*
 * Retrieve the VLAN Identifier (VID) from the IEEE 802.1Q header.
 *
 * @param hdr Pointer to the Ethernet header with IEEE 802.1Q tagging.
 * @return    Returns the VID in host byte order.
 */
static inline uint16_t
vlanhdr_get_vid (const struct openvpn_8021qhdr *hdr)
{
  return ntohs (hdr->pcp_cfi_vid & OPENVPN_8021Q_MASK_VID);
}

/*
 * Set the Priority Code Point (PCP) in an IEEE 802.1Q header.
 *
 * @param hdr Pointer to the Ethernet header with IEEE 802.1Q tagging.
 * @param pcp The PCP to set (in host byte order).
 */
static inline void
vlanhdr_set_pcp (struct openvpn_8021qhdr *hdr, const uint16_t pcp)
{
  hdr->pcp_cfi_vid = (hdr->pcp_cfi_vid & ~OPENVPN_8021Q_MASK_PCP) |
		     (htons (pcp) & OPENVPN_8021Q_MASK_PCP);
}
/*
 * Set the Canonical Format Indicator (CFI) in an IEEE 802.1Q header.
 *
 * @param hdr Pointer to the Ethernet header with IEEE 802.1Q tagging.
 * @param cfi The CFI to set (in host byte order).
 */
static inline void
vlanhdr_set_cfi (struct openvpn_8021qhdr *hdr, const uint16_t cfi)
{
  hdr->pcp_cfi_vid = (hdr->pcp_cfi_vid & ~OPENVPN_8021Q_MASK_CFI) |
		     (htons (cfi) & OPENVPN_8021Q_MASK_CFI);
}
/*
 * Set the VLAN Identifier (VID) in an IEEE 802.1Q header.
 *
 * @param hdr Pointer to the Ethernet header with IEEE 802.1Q tagging.
 * @param vid The VID to set (in host byte order).
 */
static inline void
vlanhdr_set_vid (struct openvpn_8021qhdr *hdr, const uint16_t vid)
{
  hdr->pcp_cfi_vid = (hdr->pcp_cfi_vid & ~OPENVPN_8021Q_MASK_VID) |
		     (htons (vid) & OPENVPN_8021Q_MASK_VID);
}
#endif

#endif

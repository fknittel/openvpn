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

#ifndef FORWARD_H
#define FORWARD_H

#include "openvpn.h"
#include "occ.h"
#include "ping.h"

#define TUN_OUT(c)      (BLEN(&(c)->c2.to_tun) > 0)
#define LINK_OUT(c)     (BLEN(&(c)->c2.to_link) > 0)
#define ANY_OUT(c)      (TUN_OUT(c) || LINK_OUT(c))

#ifdef ENABLE_FRAGMENT
#define TO_LINK_FRAG(c) ((c)->c2.fragment && fragment_outgoing_defined ((c)->c2.fragment))
#else
#define TO_LINK_FRAG(c) (false)
#endif

#define TO_LINK_DEF(c)  (LINK_OUT(c) || TO_LINK_FRAG(c))

/* the 1<<n bits below must be exclusive of ST_RW_MASK */

#define IOW_READ_LINK       ((1<<4) | SOCKET_READ)
#define IOW_TO_LINK         ((1<<5) | SOCKET_WRITE)

#define IOW_READ_TUN        ((1<<6) | TUN_READ)
#define IOW_TO_TUN          ((1<<7) | TUN_WRITE)

#define IOW_MBUF            ((1<<8) | SOCKET_WRITE)
#define IOW_READ_TUN_FORCE  ((1<<9) | TUN_READ)

#define IOW_CHECK_RESIDUAL  (1<<10)
#define IOW_FRAG            (1<<11)
#define IOW_SHAPER          (1<<12)
#define IOW_WAIT_SIGNAL     (1<<13)

#define IOW_FAST_IO         (1<<14)

#define IOW_READ            (IOW_READ_TUN|IOW_READ_LINK)

void p2p_iow_flags_init (struct context *c);

void io_wait_slow (struct context *c, const unsigned int flags);

void pre_select (struct context *c);

void process_io (struct context *c, bool *io_order_toggle);

void encrypt_sign (struct context *c, bool comp_frag);

const char *wait_status_string (struct context *c, struct gc_arena *gc);
void show_wait_status (struct context *c);

bool read_incoming_link (struct context *c, struct openvpn_sockaddr *from);
void process_incoming_link (struct context *c, struct openvpn_sockaddr *from);
bool read_incoming_tun (struct context *c);
void process_incoming_tun (struct context *c);
bool process_outgoing_link (struct context *c);
bool process_outgoing_tun (struct context *c);

bool send_control_channel_string (struct context *c, const char *str, int msglevel);

#define PIPV4_PASSTOS         (1<<0)
#define PIPV4_MSSFIX          (1<<1)

void process_ipv4_header (struct context *c, unsigned int flags, struct buffer *buf);

#if P2MP
void schedule_exit (struct context *c, const int n_seconds);
#endif

#endif /* FORWARD_H */

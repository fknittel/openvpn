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

#ifndef SOCKET_H
#define SOCKET_H

#include "buffer.h"
#include "common.h"
#include "error.h"
#include "proto.h"
#include "mtu.h"
#include "win32.h"
#include "event.h"
#include "proxy.h"
#include "socks.h"
#include "misc.h"

/*
 * OpenVPN's default port number as assigned by IANA.
 */
#define OPENVPN_PORT 1194

/*
 * Number of seconds that "resolv-retry infinite"
 * represents.
 */
#define RESOLV_RETRY_INFINITE 1000000000

#define REMOTE_LIST_SIZE 64

struct remote_entry
{
  const char *hostname;
  int port;
};

struct remote_list
{
  int len;
  int current;
  bool no_advance;
  struct remote_entry array[REMOTE_LIST_SIZE];
};

/* 
 * packet_size_type is used to communicate packet size
 * over the wire when stream oriented protocols are
 * being used
 */

typedef uint16_t packet_size_type;

/* convert a packet_size_type from host to network order */
#define htonps(x) htons(x)

/* convert a packet_size_type from network to host order */
#define ntohps(x) ntohs(x)

struct openvpn_sockaddr {
	union {
		struct sockaddr sa;
		struct sockaddr_in in;
#ifdef USE_PF_INET6
		struct sockaddr_in6 in6;
#endif
#ifdef USE_PF_UNIX
		struct sockaddr_un un;
#endif
	} addr;
#if ENABLE_IP_PKTINFO
	union {
		struct in_pktinfo in;
#ifdef USE_PF_INET6
		struct in6_pktinfo in6;
#endif
	} pi;	/* Multihome support for UDP */
#endif
};
/* IP addresses which are persistant across SIGUSR1s */
struct link_socket_addr
{
	struct openvpn_sockaddr local;
	struct openvpn_sockaddr remote;
	struct openvpn_sockaddr actual;
};

struct link_socket_info
{
  struct link_socket_addr *lsa;
  bool connection_established;
  const char *ipchange_command;
  const struct plugin_list *plugins;
  bool remote_float;  
  int proto;                    /* Protocol (PROTO_x defined below) */
  int mtu_changed;              /* Set to true when mtu value is changed */
};

/*
 * Used to extract packets encapsulated in streams into a buffer,
 * in this case IP packets embedded in a TCP stream.
 */
struct stream_buf
{
  struct buffer buf_init;
  struct buffer residual;
  int maxlen;
  bool residual_fully_formed;

  struct buffer buf;
  struct buffer next;
  int len;     /* -1 if not yet known */

  bool error;  /* if true, fatal TCP error has occurred,
		  requiring that connection be restarted */
};

/*
 * Used to set socket buffer sizes
 */
struct socket_buffer_size
{
  int rcvbuf;
  int sndbuf;
};

/*
 * This is the main socket structure used by OpenVPN.  The SOCKET_
 * defines try to abstract away our implementation differences between
 * using sockets on Posix vs. Win32.
 */
struct link_socket
{
  struct link_socket_info info;

  socket_descriptor_t sd;

#ifdef ENABLE_SOCKS
  socket_descriptor_t ctrl_sd;  /* only used for UDP over Socks */
#endif

#ifdef WIN32
  struct overlapped_io reads;
  struct overlapped_io writes;
  struct rw_handle rw_handle;
  struct rw_handle listen_handle; /* For listening on TCP socket in server mode */
#endif

  /* used for printing status info only */
  unsigned int rwflags_debug;

  /* used for long-term queueing of pre-accepted socket listen */
  bool listen_persistent_queued;

  /* set on initial call to init phase 1 */
  struct remote_list *remote_list;
  const char *remote_host;
  int remote_port;
  const char *local_host;
  int local_port;
  bool bind_local;

# define INETD_NONE   0
# define INETD_WAIT   1
# define INETD_NOWAIT 2
  int inetd;

# define LS_MODE_DEFAULT           0
# define LS_MODE_TCP_LISTEN        1
# define LS_MODE_TCP_ACCEPT_FROM   2
  int mode;

  int resolve_retry_seconds;
  int connect_retry_seconds;
  int mtu_discover_type;

  struct socket_buffer_size socket_buffer_sizes;

  int mtu;                      /* OS discovered MTU, or 0 if unknown */

  bool did_resolve_remote;

  /* for stream sockets */
  struct stream_buf stream_buf;
  struct buffer stream_buf_data;
  bool stream_reset;
# define SF_USE_IP_PKTINFO (1<<0)
  unsigned int socket_flags;

#ifdef ENABLE_HTTP_PROXY
  /* HTTP proxy */
  struct http_proxy_info *http_proxy;
#endif

#ifdef ENABLE_SOCKS
  /* Socks proxy */
  struct socks_proxy_info *socks_proxy;
  struct openvpn_sockaddr  socks_relay; /* Socks UDP relay address */
#endif

#if defined(ENABLE_HTTP_PROXY) || defined(ENABLE_SOCKS)
  /* The OpenVPN server we will use the proxy to connect to */
  const char *proxy_dest_host;
  int proxy_dest_port;
#endif

#if PASSTOS_CAPABILITY
  /* used to get/set TOS. */
  uint8_t ptos;
  bool ptos_defined;
#endif

#ifdef ENABLE_DEBUG
  int gremlin; /* --gremlin bits */
#endif
};

/*
 * Some Posix/Win32 differences.
 */

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

#ifdef WIN32

#define openvpn_close_socket(s) closesocket(s)

int socket_recv_queue (struct link_socket *sock, int maxsize);

int socket_send_queue (struct link_socket *sock,
		       struct buffer *buf,
		       const struct sockaddr_in *to);

int socket_finalize (
		     SOCKET s,
		     struct overlapped_io *io,
		     struct buffer *buf,
		     struct sockaddr_in *from);

#else

#define openvpn_close_socket(s) close(s)

#endif

struct link_socket *link_socket_new (void);

/*
 * Initialize link_socket object.
 */

void
link_socket_init_phase1 (struct link_socket *sock,
			 const char *local_host,
			 struct remote_list *remote_list,
			 int local_port,
			 int proto,
			 int mode,
			 const struct link_socket *accept_from,
#ifdef ENABLE_HTTP_PROXY
			 struct http_proxy_info *http_proxy,
#endif
#ifdef ENABLE_SOCKS
			 struct socks_proxy_info *socks_proxy,
#endif
#ifdef ENABLE_DEBUG
			 int gremlin,
#endif
			 bool bind_local,
			 bool remote_float,
			 int inetd,
			 struct link_socket_addr *lsa,
			 const char *ipchange_command,
			 const struct plugin_list *plugins,
			 int resolve_retry_seconds,
			 int connect_retry_seconds,
			 int mtu_discover_type,
			 int rcvbuf,
			 int sndbuf,
			 const unsigned int flags);

void link_socket_init_phase2 (struct link_socket *sock,
			      const struct frame *frame,
			      volatile int *signal_received);

void link_socket_post_fork (const struct link_socket *sock,
			    const struct sockaddr_in *remote);

void socket_adjust_frame_parameters (struct frame *frame, int proto);

void frame_adjust_path_mtu (struct frame *frame, int pmtu, int proto);

void link_socket_close (struct link_socket *sock);

#define PS_SHOW_PORT_IF_DEFINED (1<<0)
#define PS_SHOW_PORT            (1<<1)
#define PS_SHOW_PKTINFO         (1<<2)
const char *print_sockaddr_ex (const struct openvpn_sockaddr *addr,
			       const char* separator,
			       int flags,
			       struct gc_arena *gc);

const char *print_sockaddr (const struct openvpn_sockaddr *addr,
			    struct gc_arena *gc);
const char *print_link_sockaddr (const struct openvpn_sockaddr *act,
				      struct gc_arena *gc);


int addr_guess_type(int proto, const char *name);
#define IA_EMPTY_IF_UNDEF (1<<0)
#define IA_NET_ORDER      (1<<1)
const char *print_in_addr_t (in_addr_t addr, unsigned int flags, struct gc_arena *gc);

#define SA_IP_PORT        (1<<0)
#define SA_SET_IF_NONZERO (1<<1)
void setenv_sockaddr (struct env_set *es,
		      const char *name_prefix,
		      const struct openvpn_sockaddr *addr,
		      const bool flags);

void setenv_in_addr_t (struct env_set *es,
		       const char *name_prefix,
		       in_addr_t addr,
		       const bool flags);
void bad_address_length (int actual, int expected);

in_addr_t link_socket_current_remote (const struct link_socket_info *info);

void link_socket_connection_initiated (const struct buffer *buf,
				       struct link_socket_info *info,
				       const struct openvpn_sockaddr *addr,
				       const char *common_name,
				       struct env_set *es);

void link_socket_bad_incoming_addr (struct buffer *buf,
				    const struct link_socket_info *info,
				    const struct openvpn_sockaddr *from_addr);

void link_socket_bad_outgoing_addr (void);

void setenv_trusted (struct env_set *es, const struct link_socket_info *info);

void remote_list_randomize (struct remote_list *l);

/*
 * Low-level functions
 */

/* return values of openvpn_inet_aton */
#define OIA_HOSTNAME   0
#define OIA_IP         1
#define OIA_ERROR     -1
int openvpn_inet_aton (const char *dotted_quad, struct in_addr *addr);

socket_descriptor_t create_socket_tcp (void);

socket_descriptor_t socket_do_accept (socket_descriptor_t sd,
				      struct openvpn_sockaddr *act,
				      const bool nowait);

/*
 * proto related
 */

/* 
 * Use enum's instead of #define to allow for easier
 * optional proto support
 */
enum proto_num {
	PROTO_NONE, /* catch for uninitialized */
	PROTO_UDPv4,
	PROTO_TCPv4_SERVER,
	PROTO_TCPv4_CLIENT,
	PROTO_TCPv4,
	PROTO_UDPv6,
	PROTO_TCPv6_SERVER,
	PROTO_TCPv6_CLIENT,
	PROTO_TCPv6,
	PROTO_UNIX_DGRAM,
	PROTO_UNIX_STREAM,
	PROTO_N
};

struct proto_names {
  const char *short_form;
  const char *display_form;
  bool	is_dgram;
  bool	is_net;
  sa_family_t proto_af;
};

extern const struct proto_names proto_names[PROTO_N];

static inline bool
proto_is_net(int proto)
{
  ASSERT (proto >= 0 && proto < PROTO_N);
  return proto_names[proto].is_net;
}

static inline bool
proto_is_dgram(int proto)
{
  ASSERT (proto >= 0 && proto < PROTO_N);
  return proto_names[proto].is_dgram;
}

static inline bool
proto_is_udp(int proto)
{
  ASSERT (proto >= 0 && proto < PROTO_N);
  return proto_names[proto].is_dgram && proto_names[proto].is_net;
}

static inline bool
proto_is_tcp(int proto)
{
  ASSERT (proto >= 0 && proto < PROTO_N);
  return (!proto_names[proto].is_dgram) && proto_names[proto].is_net;
}

/*
 * DNS resolution
 */

#define GETADDR_RESOLVE               (1<<0)
#define GETADDR_FATAL                 (1<<1)
#define GETADDR_HOST_ORDER            (1<<2)
#define GETADDR_MENTION_RESOLVE_RETRY (1<<3)
#define GETADDR_FATAL_ON_SIGNAL       (1<<4)
#define GETADDR_WARN_ON_SIGNAL        (1<<5)
#define GETADDR_MSG_VIRT_OUT          (1<<6)
#define GETADDR_TRY_ONCE              (1<<7)

in_addr_t getaddr (unsigned int flags,
		   const char *hostname,
		   int resolve_retry_seconds,
		   bool *succeeded,
		   volatile int *signal_received);

/*
 * Transport protocol naming and other details.
 */

int ascii2proto (const char* proto_name);
const char *proto2ascii (int proto, bool display_form);
const char *proto2ascii_all (struct gc_arena *gc);
int proto_remote (int proto, bool remote);
const char *addr_family_name(int af);

/*
 * Overhead added to packets by various protocols.
 */
#define IPv4_UDP_HEADER_SIZE              28
#define IPv4_TCP_HEADER_SIZE              40
#define IPv6_UDP_HEADER_SIZE              48
#define IPv6_TCP_HEADER_SIZE              60

static const int proto_overhead[PROTO_N] = { /* indexed by PROTO_x */
  0,
  IPv4_UDP_HEADER_SIZE, /* IPv4 */
  IPv4_TCP_HEADER_SIZE,
  IPv4_TCP_HEADER_SIZE,
  IPv4_TCP_HEADER_SIZE,
#ifdef USE_PF_INET6
  IPv6_UDP_HEADER_SIZE, /* IPv6 */
  IPv6_TCP_HEADER_SIZE,
  IPv6_TCP_HEADER_SIZE,
  IPv6_TCP_HEADER_SIZE,
#endif
#ifdef USE_PF_UNIX
  0,			/* AF_UNIX proxies, assume no overhead */
  0,
#endif
};

static inline int
datagram_overhead (int proto)
{
  ASSERT (proto >= 0 && proto < PROTO_N);
  return proto_overhead [proto];
}

/*
 * Misc inline functions
 */

static inline int
remote_list_len (const struct remote_list *rl)
{
  if (rl)
    return rl->len;
  else
    return 0;
}

static inline bool
legal_ipv4_port (int port)
{
  return port > 0 && port < 65536;
}

static inline bool
link_socket_proto_connection_oriented (int proto)
{
  return !proto_is_dgram(proto);
}

static inline bool
link_socket_proto_stream_oriented (int proto)
{
  return proto_is_tcp(proto);
}

static inline bool
link_socket_connection_oriented (const struct link_socket *sock)
{
  if (sock)
    return link_socket_proto_connection_oriented (sock->info.proto);
  else
    return false;
}
static inline bool
addr_defined (const struct openvpn_sockaddr *addr)
{
  if (!addr) return 0;
  switch (addr->addr.sa.sa_family) {
    case AF_INET: return addr->addr.in.sin_addr.s_addr != 0;
#ifdef USE_PF_UNIX
    case AF_UNIX: return addr->addr.un.sun_path[0] != 0;
#endif
#ifdef USE_PF_INET6
    case AF_INET6: return !IN6_IS_ADDR_UNSPECIFIED(&addr->addr.in6.sin6_addr);
#endif
    default: return 0;
  }
}
static inline bool
addr_defined_ipi (const struct openvpn_sockaddr *addr)
{
#if ENABLE_IP_PKTINFO
  if (!addr) return 0;
  switch (addr->addr.sa.sa_family) {
    case AF_INET: return addr->pi.in.ipi_spec_dst.s_addr != 0;
#ifdef USE_PF_UNIX
    case AF_UNIX: ASSERT(0);
#endif
#ifdef USE_PF_INET6
    case AF_INET6: return !IN6_IS_ADDR_UNSPECIFIED(&addr->pi.in6.ipi6_addr);
#endif
    default: return 0;
  }
#else
  ASSERT(0);
  return 0; /* NOTREACHED */
#endif
}
static inline bool
addr_defined_sa (const struct sockaddr *addr)
{
	struct openvpn_sockaddr osa;
	memcpy(&osa.addr.sa, addr, sizeof (osa.addr));
	return addr_defined(&osa);
}

static inline bool
link_addr_defined (const struct openvpn_sockaddr *act)
{
  return addr_defined (act);
}
static inline bool
addr_match (const struct openvpn_sockaddr *a1, const struct openvpn_sockaddr *a2)
{
  switch(a1->addr.sa.sa_family) {
    case AF_INET:
      return a1->addr.in.sin_addr.s_addr == a2->addr.in.sin_addr.s_addr;
#ifdef USE_PF_UNIX
    case AF_UNIX:
      return strncmp(a1->addr.un.sun_path, a2->addr.un.sun_path, sizeof (a1->addr.un.sun_path)) == 0;
#endif
#ifdef USE_PF_INET6
    case AF_INET6:
      return IN6_ARE_ADDR_EQUAL(&a1->addr.in6.sin6_addr, &a2->addr.in6.sin6_addr);
#endif
  }
  ASSERT(0);
  return false;
}

static inline in_addr_t
addr_host (const struct openvpn_sockaddr *addr)
{
  /* 
   * "public" addr returned is checked against ifconfig for
   * possible clash: non sense for now given
   * that we do ifconfig only IPv4
   */
#if defined(USE_PF_INET6) || defined(USE_PF_UNIX)
  if(addr->addr.sa.sa_family != AF_INET)
    return 0;
#else 
  ASSERT(addr->addr.sa.sa_family == AF_INET);
#endif
  return ntohl (addr->addr.in.sin_addr.s_addr);
}

static inline bool
addr_port_match (const struct openvpn_sockaddr *a1, const struct openvpn_sockaddr *a2)
{
  switch(a1->addr.sa.sa_family) {
    case AF_INET:
      return a1->addr.in.sin_addr.s_addr == a2->addr.in.sin_addr.s_addr
	&& a1->addr.in.sin_port == a2->addr.in.sin_port;
#ifdef USE_PF_UNIX
    case AF_UNIX:
      return strncmp(a1->addr.un.sun_path, a2->addr.un.sun_path, sizeof (a1->addr.un.sun_path)) == 0;
#endif
#ifdef USE_PF_INET6
    case AF_INET6:
      return IN6_ARE_ADDR_EQUAL(&a1->addr.in6.sin6_addr, &a2->addr.in6.sin6_addr) 
	&& a1->addr.in6.sin6_port == a2->addr.in6.sin6_port;
#endif
  }
  ASSERT(0);
  return false;
}

static inline bool
addr_match_proto (const struct openvpn_sockaddr *a1,
		  const struct openvpn_sockaddr *a2,
		  const int proto)
{
  return link_socket_proto_connection_oriented (proto)
    ? addr_match (a1, a2)
    : addr_port_match (a1, a2);
}

static inline void
addr_copy_sa(struct openvpn_sockaddr *dst, const struct openvpn_sockaddr *src)
{
  dst->addr = src->addr;
}

static inline void
addr_copy_host(struct openvpn_sockaddr *dst, const struct openvpn_sockaddr *src)
{
   switch(src->addr.sa.sa_family) {
     case AF_INET:
       dst->addr.in.sin_addr.s_addr = src->addr.in.sin_addr.s_addr;
       break;
#ifdef USE_PF_UNIX
     case AF_UNIX:
       strncpynt(dst->addr.un.sun_path, src->addr.un.sun_path, sizeof dst->addr.un.sun_path);
       break;
#endif
#ifdef USE_PF_INET6
     case AF_INET6: 
       dst->addr.in6.sin6_addr = src->addr.in6.sin6_addr;
       break;
#endif
   }
}


static inline void
addr_zero_host(struct openvpn_sockaddr *addr)
{
   switch(addr->addr.sa.sa_family) {
     case AF_INET:
       addr->addr.in.sin_addr.s_addr = 0;
       break;
#ifdef USE_PF_UNIX
     case AF_UNIX:
       *addr->addr.un.sun_path=0;
       break;
#endif
#ifdef USE_PF_INET6
     case AF_INET6: 
       memset(&addr->addr.in6.sin6_addr, 0, sizeof (struct in6_addr));
       break;
#endif
   }
}

static inline bool
addr_inet4or6(struct sockaddr *addr)
{
	return addr->sa_family == AF_INET || addr->sa_family == AF_INET6;
}
int
addr_guess_family(int proto, const char *name);

static inline int
af_addr_size(sa_family_t af)
{
#if defined(USE_PF_INET6) || defined (USE_PF_UNIX)
   switch(af) {
     case AF_INET: return sizeof (struct sockaddr_in);
#ifdef USE_PF_UNIX
     case AF_UNIX: return sizeof (struct sockaddr_un);
#endif
#ifdef USE_PF_INET6
     case AF_INET6: return sizeof (struct sockaddr_in6);
#endif
     default: 
#if 0
      /* could be called from socket_do_accept() with empty addr */
      msg (M_ERR, "Bad address family: %d\n", addr->sa_family);
      ASSERT(0);
#endif
     	return 0;
   }
#else /* only AF_INET */
   return sizeof(struct sockaddr_in);
#endif
}


static inline bool
link_addr_port_match (const struct openvpn_sockaddr *a1, const struct openvpn_sockaddr *a2)
{
  return addr_port_match (a1, a2);
}

static inline bool
socket_connection_reset (const struct link_socket *sock, int status)
{
  if (link_socket_connection_oriented (sock))
    {
      if (sock->stream_reset || sock->stream_buf.error)
	return true;
      else if (status < 0)
	{
	  const int err = openvpn_errno_socket ();
#ifdef WIN32
	  return err == WSAECONNRESET || err == WSAECONNABORTED;
#else
	  return err == ECONNRESET;
#endif
	}
    }
  return false;
}

static inline bool
link_socket_verify_incoming_addr (struct buffer *buf,
				  const struct link_socket_info *info,
				  const struct openvpn_sockaddr *from_addr)
{
  if (buf->len > 0)
    {
      switch (from_addr->addr.sa.sa_family) {
#ifdef USE_PF_UNIX
	case AF_UNIX:
#endif
#ifdef USE_PF_INET6
	case AF_INET6:
#endif
	case AF_INET:
	  if (!addr_defined (from_addr))
	    return false;
	  if (info->remote_float || !addr_defined (&info->lsa->remote))
	    return true;
	  if (addr_match_proto (from_addr, &info->lsa->remote, info->proto))
	    return true;
      }
    }
  return false;
}

static inline void
link_socket_get_outgoing_addr (struct buffer *buf,
			      const struct link_socket_info *info,
			      struct openvpn_sockaddr **act)
{
  if (buf->len > 0)
    {
      struct link_socket_addr *lsa = info->lsa;
      if (link_addr_defined (&lsa->actual))
	{
	  //addr_copy(addr, &lsa->actual.addr);
	  *act = &lsa->actual;
	}
      else
	{
	  link_socket_bad_outgoing_addr ();
	  buf->len = 0;
	}
    }
}

static inline void
link_socket_set_outgoing_addr (const struct buffer *buf,
			       struct link_socket_info *info,
			       const struct openvpn_sockaddr *addr,
			       const char *common_name,
			       struct env_set *es)
{
  if (!buf || buf->len > 0)
    {
      struct link_socket_addr *lsa = info->lsa;
      if (
	  /* new or changed address? */
	  (!info->connection_established
	   || !addr_match_proto (addr, &lsa->actual, info->proto))
	  /* address undef or address == remote or --float */
	  && (info->remote_float
	      || !addr_defined (&lsa->remote)
	      || addr_match_proto (addr, &lsa->remote, info->proto))
	  )
	{
	  link_socket_connection_initiated (buf, info, addr, common_name, es);
	}
    }
}

/*
 * Stream buffer handling -- stream_buf is a helper class
 * to assist in the packetization of stream transport protocols
 * such as TCP.
 */

void stream_buf_init (struct stream_buf *sb, struct buffer *buf);
void stream_buf_close (struct stream_buf* sb);
bool stream_buf_added (struct stream_buf *sb, int length_added);

static inline bool
stream_buf_read_setup (struct link_socket* sock)
{
  bool stream_buf_read_setup_dowork (struct link_socket* sock);
  if (link_socket_connection_oriented (sock))
    return stream_buf_read_setup_dowork (sock);
  else
    return true;
}

/*
 * Socket Read Routines
 */

int link_socket_read_tcp (struct link_socket *sock,
			  struct buffer *buf);

#ifdef WIN32

static inline int
link_socket_read_udp_win32 (struct link_socket *sock,
			    struct buffer *buf,
			    struct sockaddr_in *from)
{
  return socket_finalize (sock->sd, &sock->reads, buf, from);
}

#else

int link_socket_read_udp_posix (struct link_socket *sock,
				struct buffer *buf,
				int maxsize,
				struct openvpn_sockaddr *from);

#endif
#ifdef USE_PF_UNIX
int link_socket_read_unix_dgram (struct link_socket *sock,
				struct buffer *buf,
				int maxsize,
				struct sockaddr_un *from);
#endif

/* read a TCP or UDP packet from link */
static inline int
link_socket_read (struct link_socket *sock,
		  struct buffer *buf,
		  int maxsize,
		  struct openvpn_sockaddr *from)
{
  if (proto_is_udp(sock->info.proto)) /* unified UDPv4 and UDPv6 */
    {
      int res;

#ifdef WIN32
      res = link_socket_read_udp_win32 (sock, buf, &from->addr.in);
#else
      res = link_socket_read_udp_posix (sock, buf, maxsize, from);
#endif
      return res;
    }
  else if (proto_is_tcp(sock->info.proto)) /* unified TCPv4 and TCPv6 */
    {
      /* from address was returned by accept */
      addr_copy_sa(from, &sock->info.lsa->actual);
      return link_socket_read_tcp (sock, buf);
    }
#ifdef USE_PF_UNIX
  else if (sock->info.proto == PROTO_UNIX_DGRAM)
    {
      int res;
      res = link_socket_read_unix_dgram (sock, buf, maxsize, &from->addr.un);
      return res;
    }
#endif
  else
    {
      ASSERT (0);
      return -1; /* NOTREACHED */
    }
}

/*
 * Socket Write routines
 */

int link_socket_write_tcp (struct link_socket *sock,
			   struct buffer *buf,
			   struct openvpn_sockaddr *to);

#ifdef WIN32

static inline int
link_socket_write_win32 (struct link_socket *sock,
			 struct buffer *buf,
			 struct openvpn_sockaddr *to)
{
  int err = 0;
  int status = 0;
  if (overlapped_io_active (&sock->writes))
    {
      status = socket_finalize (sock->sd, &sock->writes, NULL, NULL);
      if (status < 0)
	err = WSAGetLastError ();
    }
  socket_send_queue (sock, buf, &to->addr.in);
  if (status < 0)
    {
      WSASetLastError (err);
      return status;
    }
  else
    return BLEN (buf);
}

#else

static inline int
link_socket_write_udp_posix (struct link_socket *sock,
			     struct buffer *buf,
			     struct openvpn_sockaddr *to)
{
#if ENABLE_IP_PKTINFO
  int link_socket_write_udp_posix_sendmsg (struct link_socket *sock,
					   struct buffer *buf,
					   struct openvpn_sockaddr *to);

  /*
  if (sock->info.proto == PROTO_UDPv4 && (sock->socket_flags & SF_USE_IP_PKTINFO)
	  && to->pi.in.ipi_spec_dst.s_addr)

  */
  if (proto_is_udp(sock->info.proto) && (sock->socket_flags & SF_USE_IP_PKTINFO)
	  && addr_defined_ipi(to))
    return link_socket_write_udp_posix_sendmsg (sock, buf, to);
  else
#endif
    return sendto (sock->sd, BPTR (buf), BLEN (buf), 0,
		   &to->addr.sa,
		   (socklen_t) af_addr_size(to->addr.sa.sa_family));
}

static inline int
link_socket_write_tcp_posix (struct link_socket *sock,
			     struct buffer *buf,
			     struct openvpn_sockaddr *to)
{
  return send (sock->sd, BPTR (buf), BLEN (buf), MSG_NOSIGNAL);
}

#endif
#ifdef USE_PF_UNIX
static inline int
link_socket_write_unix_dgram (struct link_socket *sock,
			     struct buffer *buf,
			     struct sockaddr_un *to)
{
  return sendto (sock->sd, BPTR (buf), BLEN (buf), 0,
		 (struct sockaddr *) to,
		 (socklen_t) sizeof (*to));
}
#endif


static inline int
link_socket_write_udp (struct link_socket *sock,
		       struct buffer *buf,
		       struct openvpn_sockaddr *to)
{
#ifdef WIN32
  return link_socket_write_win32 (sock, buf, to);
#else
  return link_socket_write_udp_posix (sock, buf, to);
#endif
}

/* write a TCP or UDP packet to link */
static inline int
link_socket_write (struct link_socket *sock,
		   struct buffer *buf,
		   struct openvpn_sockaddr *to)
{
  if (proto_is_udp(sock->info.proto)) /* unified UDPv4 and UDPv6 */
    {
      return link_socket_write_udp (sock, buf, to);
    }
  else if (proto_is_tcp(sock->info.proto)) /* unified TCPv4 and TCPv6 */
    {
      return link_socket_write_tcp (sock, buf, to);
    }
#ifdef USE_PF_UNIX
  else if (sock->info.proto == PROTO_UNIX_DGRAM)
    {
      return link_socket_write_unix_dgram (sock, buf, &to->addr.un);
    }
#endif
  else
    {
      ASSERT (0);
      return -1; /* NOTREACHED */
    }
}

#if PASSTOS_CAPABILITY

/*
 * Extract TOS bits.  Assumes that ipbuf is a valid IPv4 packet.
 */
static inline void
link_socket_extract_tos (struct link_socket *ls, const struct buffer *ipbuf)
{
  if (ls && ipbuf)
    {
      struct openvpn_iphdr *iph = (struct openvpn_iphdr *) BPTR (ipbuf);
      ls->ptos = iph->tos;
      ls->ptos_defined = true;
    }
}

/*
 * Set socket properties to reflect TOS bits which were extracted
 * from tunnel packet.
 */
static inline void
link_socket_set_tos (struct link_socket *ls)
{
  if (ls && ls->ptos_defined)
    setsockopt (ls->sd, IPPROTO_IP, IP_TOS, &ls->ptos, sizeof (ls->ptos));
}

#endif

/*
 * Socket I/O wait functions
 */

static inline bool
socket_read_residual (const struct link_socket *s)
{
  return s && s->stream_buf.residual_fully_formed;
}

static inline event_t
socket_event_handle (const struct link_socket *s)
{
#ifdef WIN32
  return &s->rw_handle;
#else
  return s->sd;
#endif
}

event_t socket_listen_event_handle (struct link_socket *s);

unsigned int
socket_set (struct link_socket *s,
	    struct event_set *es,
	    unsigned int rwflags,
	    void *arg,
	    unsigned int *persistent);

static inline void
socket_set_listen_persistent (struct link_socket *s,
			      struct event_set *es,
			      void *arg)
{
  if (s && !s->listen_persistent_queued)
    {
      event_ctl (es, socket_listen_event_handle (s), EVENT_READ, arg);
      s->listen_persistent_queued = true;
    }
}

static inline void
socket_reset_listen_persistent (struct link_socket *s)
{
#ifdef WIN32
  reset_net_event_win32 (&s->listen_handle, s->sd);
#endif
}

const char *socket_stat (const struct link_socket *s, unsigned int rwflags, struct gc_arena *gc);

#endif /* SOCKET_H */

/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single TCP/UDP port, with support for SSL/TLS-based
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

#include "forward.h"
#include "helper.h"
#include "pool.h"
#include "push.h"

#include "memdbg.h"

static const char *
print_netmask (int netbits, struct gc_arena *gc)
{
  struct buffer out = alloc_buf_gc (128, gc);
  const in_addr_t netmask = netbits_to_netmask (netbits);

  buf_printf (&out, "%s (/%d)", print_in_addr_t (netmask, 0, gc), netbits);

  return BSTR (&out);
}

static const char *
print_opt_route_gateway (const in_addr_t route_gateway, struct gc_arena *gc)
{
  struct buffer out = alloc_buf_gc (128, gc);
  ASSERT (route_gateway);
  buf_printf (&out, "route-gateway %s", print_in_addr_t (route_gateway, 0, gc));
  return BSTR (&out);
}

static const char *
print_opt_route (const in_addr_t network, const in_addr_t netmask, struct gc_arena *gc)
{
  struct buffer out = alloc_buf_gc (128, gc);
  ASSERT (network);

  if (netmask)
    buf_printf (&out, "route %s %s",
		print_in_addr_t (network, 0, gc),
		print_in_addr_t (netmask, 0, gc));
  else
    buf_printf (&out, "route %s",
		print_in_addr_t (network, 0, gc));

  return BSTR (&out);
}

static const char *
print_str_int (const char *str, const int i, struct gc_arena *gc)
{
  struct buffer out = alloc_buf_gc (128, gc);
  buf_printf (&out, "%s %d", str, i);
  return BSTR (&out);
}

static void
helper_add_route (const in_addr_t network, const in_addr_t netmask, struct options *o)
{
  rol_check_alloc (o);
  add_route_to_option_list (o->routes,
			    print_in_addr_t (network, 0, &o->gc),
			    print_in_addr_t (netmask, 0, &o->gc),
			    NULL,
			    NULL);
}

static void
verify_common_subnet (const char *opt, const in_addr_t a, const in_addr_t b, const in_addr_t subnet)
{
  struct gc_arena gc = gc_new ();
  if ((a & subnet) != (b & subnet))
    msg (M_USAGE, "Options Error: %s IP addresses %s and %s are not in the same %s subnet",
	 opt,
	 print_in_addr_t (a, 0, &gc),
	 print_in_addr_t (b, 0, &gc),
	 print_in_addr_t (subnet, 0, &gc));
  gc_free (&gc);
}

/*
 * Process server, server-bridge, and client helper
 * directives after the parameters themselves have been
 * parsed and placed in struct options.
 */
void
helper_client_server (struct options *o)
{
  struct gc_arena gc = gc_new ();
#if P2MP
  /*
   * Get tun/tap/null device type
   */
  const int dev = dev_type_enum (o->dev, o->dev_type);

  /*
   *
   * HELPER DIRECTIVE:
   *
   * server 10.8.0.0 255.255.255.0
   *
   * EXPANDS TO:
   *
   * mode server
   * tls-server
   *
   * if tun:
   *   ifconfig 10.8.0.1 10.8.0.2 
   *   ifconfig-pool 10.8.0.4 10.8.0.251
   *   route 10.8.0.0 255.255.255.0
   *   if client-to-client:
   *     push "route 10.8.0.0 255.255.255.0"
   *   else if !linear-addr:
   *     push "route 10.8.0.1"
   *
   * if tap:
   *   ifconfig 10.8.0.1 255.255.255.0
   *   ifconfig-pool 10.8.0.2 10.8.0.254 255.255.255.0
   *   push "route-gateway 10.8.0.1"
   */
  if (o->server_defined)
    {
      int netbits = -2;
      bool status = false;

      if (o->client)
	msg (M_USAGE, "Options Error: --server and --client cannot be used together");

      if (o->server_bridge_defined)
	msg (M_USAGE, "Options Error: --server and --server-bridge cannot be used together");

      if (o->shared_secret_file)
	msg (M_USAGE, "Options Error: --server and --secret cannot be used together (you must use SSL/TLS keys)");

      if (o->ifconfig_pool_defined)
	msg (M_USAGE, "Options Error: --server already defines an ifconfig-pool, so you can't also specify --ifconfig-pool explicitly");

      if (!(dev == DEV_TYPE_TAP || dev == DEV_TYPE_TUN))
	msg (M_USAGE, "Options Error: --server directive only makes sense with --dev tun or --dev tap");

      status = netmask_to_netbits (o->server_network, o->server_netmask, &netbits);
      if (!status)
	msg (M_USAGE, "Options Error: --server directive network/netmask combination is invalid");

      if (netbits < 0)
	msg (M_USAGE, "Options Error: --server directive netmask is invalid");

      if (netbits < IFCONFIG_POOL_MIN_NETBITS)
	msg (M_USAGE, "Options Error: --server directive netmask allows for too many host addresses (subnet must be %s or higher)",
	     print_netmask (IFCONFIG_POOL_MIN_NETBITS, &gc));

      if (dev == DEV_TYPE_TUN)
	{
	  int pool_end_reserve = 4;

	  if (netbits > 29)
	    msg (M_USAGE, "Options Error: --server directive when used with --dev tun must define a subnet of %s or lower",
		 print_netmask (29, &gc));

	  if (netbits == 29)
	    pool_end_reserve = 0;

	  o->mode = MODE_SERVER;
	  o->tls_server = true;
	  o->ifconfig_local = print_in_addr_t (o->server_network + 1, 0, &o->gc);
	  o->ifconfig_remote_netmask = print_in_addr_t (o->server_network + 2, 0, &o->gc);
	  o->ifconfig_pool_defined = true;
	  o->ifconfig_pool_start = o->server_network + 4;
	  o->ifconfig_pool_end = (o->server_network | ~o->server_netmask) - pool_end_reserve;
	  helper_add_route (o->server_network, o->server_netmask, o);
	  if (o->enable_c2c)
	    push_option (o, print_opt_route (o->server_network, o->server_netmask, &o->gc), M_USAGE);
	  else if (!o->ifconfig_pool_linear)
	    push_option (o, print_opt_route (o->server_network + 1, 0, &o->gc), M_USAGE);
	}
      else if (dev == DEV_TYPE_TAP)
	{
	  if (netbits >= 30)
	    msg (M_USAGE, "Options Error: --server directive when used with --dev tap must define a subnet of %s or lower",
		 print_netmask (30, &gc));

	  o->mode = MODE_SERVER;
	  o->tls_server = true;
	  o->ifconfig_local = print_in_addr_t (o->server_network + 1, 0, &o->gc);
	  o->ifconfig_remote_netmask = print_in_addr_t (o->server_netmask, 0, &o->gc);
	  o->ifconfig_pool_defined = true;
	  o->ifconfig_pool_start = o->server_network + 2;
	  o->ifconfig_pool_end = (o->server_network | ~o->server_netmask) - 1;
	  o->ifconfig_pool_netmask = o->server_netmask;
	  push_option (o, print_opt_route_gateway (o->server_network + 1, &o->gc), M_USAGE);
	}
      else
	{
	  ASSERT (0);
	}

      if (o->proto == PROTO_TCPv4)
	o->proto = PROTO_TCPv4_SERVER;
    }

  /*
   * HELPER DIRECTIVE:
   *
   * server-bridge 10.8.0.4 255.255.255.0 10.8.0.128 10.8.0.254
   *
   * EXPANDS TO:
   *
   * mode server
   * tls-server
   *
   * ifconfig-pool 10.8.0.128 10.8.0.254 255.255.255.0
   * push "route-gateway 10.8.0.4"
   */
  else if (o->server_bridge_defined)
    {
      if (o->client)
	msg (M_USAGE, "Options Error: --server-bridge and --client cannot be used together");

      if (o->ifconfig_pool_defined)
	msg (M_USAGE, "Options Error: --server-bridge already defines an ifconfig-pool, so you can't also specify --ifconfig-pool explicitly");

      if (o->shared_secret_file)
	msg (M_USAGE, "Options Error: --server-bridge and --secret cannot be used together (you must use SSL/TLS keys)");

      if (dev != DEV_TYPE_TAP)
	msg (M_USAGE, "Options Error: --server-bridge directive only makes sense with --dev tap");

      verify_common_subnet ("--server-bridge", o->server_bridge_ip, o->server_bridge_pool_start, o->server_bridge_netmask); 
      verify_common_subnet ("--server-bridge", o->server_bridge_pool_start, o->server_bridge_pool_end, o->server_bridge_netmask); 
      verify_common_subnet ("--server-bridge", o->server_bridge_ip, o->server_bridge_pool_end, o->server_bridge_netmask); 

      o->mode = MODE_SERVER;
      o->tls_server = true;
      o->ifconfig_pool_defined = true;
      o->ifconfig_pool_start = o->server_bridge_pool_start;
      o->ifconfig_pool_end = o->server_bridge_pool_end;
      o->ifconfig_pool_netmask = o->server_bridge_netmask;
      push_option (o, print_opt_route_gateway (o->server_bridge_ip, &o->gc), M_USAGE);

      if (o->proto == PROTO_TCPv4)
	o->proto = PROTO_TCPv4_SERVER;
    }

  /*
   * HELPER DIRECTIVE:
   *
   * client
   *
   * EXPANDS TO:
   *
   * pull
   * tls-client
   */
  else if (o->client)
    {
      o->pull = true;
      o->tls_client = true;

      if (o->proto == PROTO_TCPv4)
	o->proto = PROTO_TCPv4_CLIENT;
    }
#endif

  if (o->proto == PROTO_TCPv4)
    msg (M_USAGE, "Options Error: --proto tcp is ambiguous in this context.  Please specify --proto tcp-server or --proto tcp-client");

  gc_free (&gc);
}

/*
 *
 * HELPER DIRECTIVE:
 *
 * keepalive 10 60
 *
 * EXPANDS TO:
 *
 * if mode server:
 *   ping 10
 *   ping-restart 120
 *   push "ping 10"
 *   push "ping-restart 60"
 * else
 *   ping 10
 *   ping-restart 60
 */
void
helper_keepalive (struct options *o)
{
  if (o->keepalive_ping || o->keepalive_timeout)
    {
      /*
       * Sanity checks.
       */
      if (o->keepalive_ping <= 0 || o->keepalive_timeout <= 0)
	msg (M_USAGE, "Options Error: --keepalive parameters must be > 0");
      if (o->keepalive_ping * 2 > o->keepalive_timeout)
	msg (M_USAGE, "Options Error: the second parameter to --keepalive (restart timeout=%d) must be at least twice the value of the first parameter (ping interval=%d).  A ratio of 1:5 or 1:6 would be even better.  Recommended setting is --keepalive 10 60.",
	     o->keepalive_timeout,
	     o->keepalive_ping);
      if (o->ping_send_timeout || o->ping_rec_timeout)
	msg (M_USAGE, "Options Error: --keepalive conflicts with --ping, --ping-exit, or --ping-restart.  If you use --keepalive, you don't need any of the other --ping directives.");

      /*
       * Expand.
       */
      if (o->mode == MODE_POINT_TO_POINT)
	{
	  o->ping_rec_timeout_action = PING_RESTART;
	  o->ping_send_timeout = o->keepalive_ping;
	  o->ping_rec_timeout = o->keepalive_timeout;
	}
#if P2MP
      else if (o->mode == MODE_SERVER)
	{
	  o->ping_rec_timeout_action = PING_RESTART;
	  o->ping_send_timeout = o->keepalive_ping;
	  o->ping_rec_timeout = o->keepalive_timeout * 2;
	  push_option (o, print_str_int ("ping", o->keepalive_ping, &o->gc), M_USAGE);
	  push_option (o, print_str_int ("ping-restart", o->keepalive_timeout, &o->gc), M_USAGE);
	}
#endif
      else
	{
	  ASSERT (0);
	}
    }
}

/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single UDP port, with support for TLS-based
 *             session authentication and key exchange,
 *             packet encryption, packet authentication, and
 *             packet compression.
 *
 *  Copyright (C) 2002 James Yonan <jim@yonan.net>
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

/*
 * This file is based on the TUN/TAP driver interface routines
 * from VTun by Maxim Krasnyansky <max_mk@yahoo.com>.
 */

#include "config.h"

#include "syshead.h"

#include "tun.h"
#include "fdmisc.h"
#include "error.h"
#include "buffer.h"
#include "common.h"
#include "misc.h"

#include "memdbg.h"

static bool
is_dev_type (const char *dev, const char *dev_type, const char *match_type)
{
  ASSERT (dev);
  ASSERT (match_type);
  if (dev_type)
    return !strcmp (dev_type, match_type);
  else
    return !strncmp (dev, match_type, strlen (match_type));
}

const char *
dev_type_string(const char *dev, const char *dev_type)
{
  if (is_dev_type (dev, dev_type, "tun"))
    return "tun";
  else if (is_dev_type (dev, dev_type, "tap"))
    return "tap";
  else if (is_dev_type (dev, dev_type, "null"))
    return "null";
  else
    return "[unknown-dev-type]";
}

const char *
dev_component_in_dev_node (const char *dev_node)
{
  const char *ret;
  const int dirsep = '/';

  if (dev_node)
    {
      ret = rindex (dev_node, dirsep);
      if (ret && *ret)
	++ret;
      else
	ret = dev_node;
      if (*ret)
	return ret;
    }
  return NULL;
}

/*
 * Called by the open_tun function of OSes to check if we
 * explicitly support IPv6.
 *
 * In this context, explicit means that the OS expects us to
 * do something special to the tun socket in order to support
 * IPv6, i.e. it is not transparent.
 *
 * ipv6_explicitly_supported should be set to false if we don't
 * have any explicit IPv6 code in the tun device handler.
 *
 * If ipv6_explicitly_supported is true, then we have explicit
 * OS-specific tun dev code for handling IPv6.  If so, tt->ipv6
 * is set according to the --tun-ipv6 command line option.
 */
static void
ipv6_support (bool ipv6, bool ipv6_explicitly_supported, struct tuntap* tt)
{
  tt->ipv6 = false;
  if (ipv6_explicitly_supported)
    tt->ipv6 = ipv6;
  else if (ipv6)
    msg (M_WARN, "NOTE: explicit support for IPv6 tun devices is not provided for this OS");
}

/* do ifconfig */
void
do_ifconfig (const char *dev, const char *dev_type,
	     const char *ifconfig_local, const char *ifconfig_remote,
	     int tun_mtu)
{
  if (ifconfig_local && ifconfig_remote)
    {
      char command_line[256];

      if (!is_dev_type (dev, dev_type, "tun"))
	msg (M_FATAL, "%s is not a tun device.  The --ifconfig option works only for tun devices.  You should use an --up script to ifconfig a tap device.", dev);

#if defined(TARGET_LINUX)

      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s %s pointopoint %s mtu %d",
		dev,
		ifconfig_local,
		ifconfig_remote,
		tun_mtu
		);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, "Linux ifconfig failed", true);

#elif defined(TARGET_SOLARIS)

      /* example: ifconfig tun2 10.2.0.2 10.2.0.1 mtu 1450 netmask 255.255.255.255 up */
      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s %s %s mtu %d netmask 255.255.255.255 up",
		dev,
		ifconfig_local,
		ifconfig_remote,
		tun_mtu
		);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, "Solaris ifconfig failed", true);

#elif defined(TARGET_OPENBSD)

      /*
       * OpenBSD tun devices appear to be persistent by default.  It seems in order
       * to make this work correctly, we need to delete the previous instance
       * (if it exists), and re-ifconfig.  Let me know if you know a better way.
       */

      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s delete",
		dev);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, NULL, false);
      msg (M_INFO, "NOTE: Tried to delete pre-existing tun instance -- No Problem if failure");


      /* example: ifconfig tun2 10.2.0.2 10.2.0.1 mtu 1450 netmask 255.255.255.255 up */
      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s %s %s mtu %d netmask 255.255.255.255 up",
		dev,
		ifconfig_local,
		ifconfig_remote,
		tun_mtu
		);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, "OpenBSD ifconfig failed", true);

#elif defined(TARGET_NETBSD)

      snprintf (command_line, sizeof (command_line),
	        IFCONFIG_PATH " %s %s %s mtu %d netmask 255.255.255.255 up",
	        dev,
		ifconfig_local,
                ifconfig_remote,
                tun_mtu
                );
      msg (M_INFO, "%s", command_line);
      system_check (command_line, "NetBSD ifconfig failed", true);

#elif defined(TARGET_DARWIN)

      /*
       * Darwin seems to exibit similar behaviour to OpenBSD...
       */

      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s delete",
		dev);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, NULL, false);
      msg (M_INFO, "NOTE: Tried to delete pre-existing tun instance -- No Problem if failure");


      /* example: ifconfig tun2 10.2.0.2 10.2.0.1 mtu 1450 netmask 255.255.255.255 up */
      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s %s %s mtu %d netmask 255.255.255.255 up",
		dev,
		ifconfig_local,
		ifconfig_remote,
		tun_mtu
		);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, "Darwin ifconfig failed", true);

#elif defined(TARGET_FREEBSD)

      /* example: ifconfig tun2 10.2.0.2 10.2.0.1 mtu 1450 netmask 255.255.255.255 up */
      snprintf (command_line, sizeof (command_line),
		IFCONFIG_PATH " %s %s %s mtu %d netmask 255.255.255.255 up",
		dev,
		ifconfig_local,
		ifconfig_remote,
		tun_mtu
		);
      msg (M_INFO, "%s", command_line);
      system_check (command_line, "FreeBSD ifconfig failed", true);
#else
      msg (M_FATAL, "Sorry, but I don't know how to do 'ifconfig' commands on this operating system.  You should ifconfig your tun/tap device manually or use an --up script.");
#endif
    }
}

void
clear_tuntap (struct tuntap *tuntap)
{
  tuntap->fd = -1;
#ifdef TARGET_SOLARIS
  tuntap->ip_fd = -1;
#endif
  tuntap->ipv6 = false;
  CLEAR (tuntap->actual);
}

static void
open_null (struct tuntap *tt)
{
  clear_tuntap (tt);
  strncpynt (tt->actual, "null", sizeof (tt->actual));
}

static void
open_tun_generic (const char *dev, const char *dev_node,
		  bool ipv6, bool ipv6_explicitly_supported,
		  struct tuntap *tt)
{
  char tunname[64];

  clear_tuntap (tt);

  ipv6_support (ipv6, ipv6_explicitly_supported, tt);

  if (!strcmp(dev, "null"))
    {
      open_null (tt);
    }
  else
    {
      if (dev_node)
	snprintf (tunname, sizeof (tunname), "%s", dev_node);
      else
	snprintf (tunname, sizeof (tunname), "/dev/%s", dev);
      if ((tt->fd = open (tunname, O_RDWR)) < 0)
	msg (M_ERR, "Cannot open tun/tap dev %s", tunname);
      set_nonblock (tt->fd);
      msg (M_INFO, "tun/tap device %s opened", tunname);
      strncpynt (tt->actual, dev, sizeof (tt->actual));
    }
}

static void
close_tun_generic (struct tuntap *tt)
{
  if (tt->fd >= 0)
    close (tt->fd);
  clear_tuntap (tt);
}

#if defined(TARGET_LINUX)

#ifdef HAVE_LINUX_IF_TUN_H	/* New driver support */

#if defined(HAVE_TUN_PI) && defined(HAVE_IPHDR) && defined(HAVE_IOVEC) && defined(ETH_P_IPV6) && defined(ETH_P_IP) && defined(HAVE_READV) && defined(HAVE_WRITEV)
#define LINUX_IPV6 1
/* #warning IPv6 ON */
#else
#define LINUX_IPV6 0
/* #warning IPv6 OFF */
#endif

void
open_tun (const char *dev, const char *dev_type, const char *dev_node, bool ipv6,
	  struct tuntap *tt)
{
  struct ifreq ifr;

  clear_tuntap (tt);

  ipv6_support (ipv6, LINUX_IPV6, tt);

  if (!strcmp(dev, "null"))
    {
      open_null (tt);
    }
  else
    {
      if (!dev_node)
	dev_node = "/dev/net/tun";
      if ((tt->fd = open (dev_node, O_RDWR)) < 0)
	msg (M_ERR, "Cannot open tun/tap dev %s", dev_node);

      CLEAR (ifr);
      if (!tt->ipv6)
	ifr.ifr_flags = IFF_NO_PI;

      if (is_dev_type (dev, dev_type, "tun"))
	{
	  ifr.ifr_flags |= IFF_TUN;
	}
      else if (is_dev_type (dev, dev_type, "tap"))
	{
	  ifr.ifr_flags |= IFF_TAP;
	}
      else
	{
	  msg (M_FATAL, "I don't recognize device %s as a tun or tap device",
	       dev);
	}
      if (strlen (dev) > 3)		/* unit number specified? */
	strncpynt (ifr.ifr_name, dev, IFNAMSIZ);

      if (ioctl (tt->fd, TUNSETIFF, (void *) &ifr) < 0)
	msg (M_ERR, "Cannot ioctl TUNSETIFF %s", dev);

      set_nonblock (tt->fd);
      msg (M_INFO, "tun/tap device %s opened", ifr.ifr_name);
      strncpynt (tt->actual, ifr.ifr_name, sizeof (tt->actual));
    }
}

#ifdef TUNSETPERSIST

void
tuncfg (const char *dev, const char *dev_type, const char *dev_node, bool ipv6, int persist_mode)
{
  struct tuntap tt;

  open_tun (dev, dev_type, dev_node, ipv6, &tt);
  if (ioctl (tt.fd, TUNSETPERSIST, persist_mode) < 0)
    msg (M_ERR, "Cannot ioctl TUNSETPERSIST(%d) %s", persist_mode, dev);
  close_tun (&tt);
  msg (M_INFO, "Persist state set to: %s", (persist_mode ? "ON" : "OFF"));
}

#endif /* TUNSETPERSIST */

#else

void
open_tun (const char *dev, const char *dev_type, const char *dev_node, bool ipv6,
	  struct tuntap *tt)
{
  open_tun_generic (dev, dev_node, ipv6, false, tt);
}

#endif /* HAVE_LINUX_IF_TUN_H */

void
close_tun (struct tuntap *tt)
{
  close_tun_generic (tt);
}

int
write_tun (struct tuntap* tt, uint8_t *buf, int len)
{
#if LINUX_IPV6
  if (tt->ipv6)
    {
      struct tun_pi pi;
      struct iphdr *iph;
      struct iovec vect[2];
      int ret;

      iph = (struct iphdr *)buf;

      pi.flags = 0;

      if(iph->version == 6)
	pi.proto = htons(ETH_P_IPV6);
      else
	pi.proto = htons(ETH_P_IP);

      vect[0].iov_len = sizeof(pi);
      vect[0].iov_base = &pi;
      vect[1].iov_len = len;
      vect[1].iov_base = buf;

      ret = writev(tt->fd, vect, 2);
      return(ret - sizeof(pi));
    }
  else
#endif
    return write (tt->fd, buf, len);
}

int
read_tun (struct tuntap* tt, uint8_t *buf, int len)
{
#if LINUX_IPV6
  if (tt->ipv6)
    {
      struct iovec vect[2];
      struct tun_pi pi;
      int ret;

      vect[0].iov_len = sizeof(pi);
      vect[0].iov_base = &pi;
      vect[1].iov_len = len;
      vect[1].iov_base = buf;

      ret = readv(tt->fd, vect, 2);
      return(ret - sizeof(pi));
    }
  else
#endif
    return read (tt->fd, buf, len);
}

#elif defined(TARGET_SOLARIS)

void
open_tun (const char *dev, const char *dev_type, const char *dev_node, bool ipv6,
	  struct tuntap *tt)
{
  int if_fd, muxid, ppa = -1;
  struct ifreq ifr;
  const char *ptr;
  const char *ip_node;
  const char *dev_tuntap_type;
  int link_type;
  bool is_tun;

  clear_tuntap (tt);

  ipv6_support (ipv6, false, tt);

  if (!strcmp(dev, "null"))
    {
      open_null (tt);
      return;
    }

  if (is_dev_type (dev, dev_type, "tun"))
    {
      ip_node = "/dev/udp";
      if (!dev_node)
	dev_node = "/dev/tun";
      dev_tuntap_type = "tun";
      link_type = I_PLINK;
      is_tun = true;
    }
  else if (is_dev_type (dev, dev_type, "tap"))
    {
      ip_node = "/dev/ip";
      if (!dev_node)
	dev_node = "/dev/tap";
      dev_tuntap_type = "tap";
      link_type = I_PLINK; /* was: I_LINK */
      is_tun = false;
    }
  else
    {
      msg (M_FATAL, "I don't recognize device %s as a tun or tap device",
	   dev);
    }
  
  /* get unit number */
  if (*dev)
    {
      ptr = dev;
      while (*ptr && !isdigit ((int) *ptr))
	ptr++;
      ppa = atoi (ptr);
    }

  if ((tt->ip_fd = open (ip_node, O_RDWR, 0)) < 0)
    msg (M_ERR, "Can't open %s", ip_node);

  if ((tt->fd = open (dev_node, O_RDWR, 0)) < 0)
    msg (M_ERR, "Can't open %s", dev_node);

  /* Assign a new PPA and get its unit number. */
  if ((ppa = ioctl (tt->fd, TUNNEWPPA, ppa)) < 0)
    msg (M_ERR, "Can't assign new interface");

  if ((if_fd = open (dev_node, O_RDWR, 0)) < 0)
    msg (M_ERR, "Can't open %s (2)", dev_node);

  if (ioctl (if_fd, I_PUSH, "ip") < 0)
    msg (M_ERR, "Can't push IP module");

  /* Assign ppa according to the unit number returned by tun device */
  if (ioctl (if_fd, IF_UNITSEL, (char *) &ppa) < 0)
    msg (M_ERR, "Can't set PPA %d", ppa);

  if ((muxid = ioctl (tt->ip_fd, link_type, if_fd)) < 0)
    msg (M_ERR, "Can't link %s device to IP", dev_tuntap_type);

  close (if_fd);

  snprintf (tt->actual, sizeof (tt->actual), "%s%d", dev_tuntap_type, ppa);

  CLEAR (ifr);
  strncpynt (ifr.ifr_name, tt->actual, sizeof (ifr.ifr_name));
  ifr.ifr_ip_muxid = muxid;

  if (ioctl (tt->ip_fd, SIOCSIFMUXID, &ifr) < 0)
    {
      ioctl (tt->ip_fd, I_PUNLINK, muxid);
      msg (M_ERR, "Can't set multiplexor id");
    }

  set_nonblock (tt->fd);
}

/*
 * Close TUN device. 
 */
void
close_tun (struct tuntap* tt)
{
  if (tt->fd >= 0)
    {
      struct ifreq ifr;

      CLEAR (ifr);
      strncpynt (ifr.ifr_name, tt->actual, sizeof (ifr.ifr_name));

     if (ioctl (tt->ip_fd, SIOCGIFFLAGS, &ifr) < 0)
	msg (M_ERR, "Can't get iface flags");

      if (ioctl (tt->ip_fd, SIOCGIFMUXID, &ifr) < 0)
	msg (M_ERR, "Can't get multiplexor id");

      if (ioctl (tt->ip_fd, I_PUNLINK, ifr.ifr_ip_muxid) < 0)
	msg (M_ERR, "Can't unlink interface");

      close (tt->ip_fd);
      close (tt->fd);
    }
  clear_tuntap (tt);
}

int
write_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  struct strbuf sbuf;
  sbuf.len = len;
  sbuf.buf = buf;
  return putmsg (tt->fd, NULL, &sbuf, 0) >= 0 ? sbuf.len : -1;
}

int
read_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  struct strbuf sbuf;
  int f = 0;

  sbuf.maxlen = len;
  sbuf.buf = buf;
  return getmsg (tt->fd, NULL, &sbuf, &f) >= 0 ? sbuf.len : -1;
}

#elif defined(TARGET_OPENBSD)

void
open_tun (const char *dev, const char *dev_type, const char *dev_node, bool ipv6,
	  struct tuntap *tt)
{
  open_tun_generic (dev, dev_node, ipv6, false, tt);
}

void
close_tun (struct tuntap* tt)
{
  close_tun_generic (tt);
}

static inline int
openbsd_modify_read_write_return (int len)
{
 if (len > 0)
    return len > sizeof (u_int32_t) ? len - sizeof (u_int32_t) : 0;
  else
    return len;
}

int
write_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  u_int32_t type = htonl (AF_INET);
  struct iovec iv[2];

  iv[0].iov_base = &type;
  iv[0].iov_len = sizeof (type);
  iv[1].iov_base = buf;
  iv[1].iov_len = len;

  return openbsd_modify_read_write_return (writev (tt->fd, iv, 2));
}

int
read_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  u_int32_t type;
  struct iovec iv[2];

  iv[0].iov_base = &type;
  iv[0].iov_len = sizeof (type);
  iv[1].iov_base = buf;
  iv[1].iov_len = len;

  return openbsd_modify_read_write_return (readv (tt->fd, iv, 2));
}

#elif defined(TARGET_FREEBSD)

void
open_tun (const char *dev, const char *dev_type, const char *dev_node, bool ipv6,
	  struct tuntap *tt)
{
  open_tun_generic (dev, dev_node, ipv6, false, tt);

  if (tt->fd >= 0)
    {
      int i = 0;

      /* Disable extended modes */
      ioctl (tt->fd, TUNSLMODE, &i);
      ioctl (tt->fd, TUNSIFHEAD, &i);
    }
}

void
close_tun (struct tuntap* tt)
{
  close_tun_generic (tt);
}

int
write_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  return write (tt->fd, buf, len);
}

int
read_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  return read (tt->fd, buf, len);
}

#else /* generic */

void
open_tun (const char *dev, const char *dev_type, const char *dev_node, bool ipv6,
	  struct tuntap *tt)
{
  open_tun_generic (dev, dev_node, ipv6, false, tt);
}

void
close_tun (struct tuntap* tt)
{
  close_tun_generic (tt);
}

int
write_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  return write (tt->fd, buf, len);
}

int
read_tun (struct tuntap* tt, uint8_t *buf, int len)
{
  return read (tt->fd, buf, len);
}

#endif

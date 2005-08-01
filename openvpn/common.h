/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single UDP port, with support for SSL/TLS-based
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

#ifndef COMMON_H
#define COMMON_H

/*
 * This parameter controls the TLS channel buffer size.  Among
 * other things, this buffer must be large enough to contain
 * the full --push/--pull list.  If you increase it, do so
 * on both server and client.
 * Recommended default: 1024
 */
#define TLS_CHANNEL_BUF_SIZE 1024

/*
 * Maximum number of --route directives.
 * Recommended default: 100
 */
#define MAX_ROUTES 100

/*
 * Maximum number of groups.
 * Recommended default: 64
 */
#define MAX_GROUPS 64

/*
 * Extremely fast small-sized cache for IPv4 -> struct multi_instance lookup.
 * Currently supported values are 4 or 8.  To disable, leave undefined.
 * Recommended default: 4
 */
#define FAST_ADDR_LOOKUP 4

/*
 * OpenVPN needs to be able to do extremely fast equality comparisons of
 * 6-byte strings.  If set to 1, use memcmp.  Otherwise, do a separate
 * uint32 and uint16 comparison.  gcc 3.x on x86 seems to generate better
 * code with a 0 setting.
 */
#define USE_MEMCMP_6 0

/*
 * (Experimental) Define to 1 to change OpenVPN protocol to improve buffer
 * alignment efficiency.  Define to 0 otherwise.  Both server and all
 * client executables must use the same ALIGN_OPTIMIZE setting.
 * Note that setting this to 1 changes the OpenVPN protocol.
 * Recommended default: 0
 */
#define ALIGN_OPTIMIZE 0

/*
 * Verify buffer alignment at run-time.  Generally should only be
 * used by developers or testers.
 * Recommended default: 0
 */
#define VERIFY_ALIGNMENT 0

/*
 * Preferred alignment of buffer data in bytes.
 * Must be a power of 2.
 * Recommended default: 4
 */
#define BUFFER_ALIGN 4

/*
 * Tuning parameters for --fast-io option.
 */

/*
 * Maximum number of client instances which can be queued at one time.
 * Recommended default: 16
 */
#define MPD_MAX_QUEUED_INSTANCES 16

/*
 * Maximum number of I/O iterations before queue flush is forced.
 * Recommended default: 128
 */
#define MPD_MAX_ITERATIONS       128

/*
 * the --client-config-dir default file.
 */
#define CCD_DEFAULT "DEFAULT"

/*
 * MTU parameters
 */

/*
 * Standard ethernet MTU
 * Recommended default: 1500
 */
#define ETHERNET_MTU       1500

/*
 * It is a fatal error if mtu is less than
 * this value for tun device.
 * Recommended default: 100
 */
#define TUN_MTU_MIN        100

/*
 * Default MTU of network over which tunnel data will pass by TCP/UDP.
 * Recommended default: 1500
 */
#define LINK_MTU_DEFAULT   1500

/*
 * Default MTU of tunnel device.
 * Recommended default: 1500
 */
#define TUN_MTU_DEFAULT    1500

/*
 * MTU Defaults for TAP devices
 * Recommended default: 32
 */
#define TAP_MTU_EXTRA_DEFAULT  32

/*
 * Default MSSFIX value, used for reducing TCP MTU size
 * Recommended default: 1450
 */
#define MSSFIX_DEFAULT     1450

/*
 * Statistics counters.
 */
typedef unsigned long counter_type;

/*
 * Time intervals
 */
typedef int interval_t;

/*
 * Used as an upper bound for timeouts.
 */
#define BIG_TIMEOUT  (60*60*24*7)  /* one week (in seconds) */

/*
 * Printf formats for special types
 */
#define counter_format          "%lu"
#define ptr_format              "0x%08lx"
#define time_format             "%lu"
#define fragment_header_format  "0x%08x"

/* these are used to cast the arguments
 * and MUST match the formats above */
typedef unsigned long time_type;
typedef unsigned long ptr_type;

#endif

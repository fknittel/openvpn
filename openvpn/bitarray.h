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

/*
 * Bit arrays
 */

#ifndef BITARRAY_H
#define BITARRAY_H

#define BA_TYPE       uint32_t
#define BA_TYPE_BITS  (sizeof(BA_TYPE)*8)

#define BA_DIM(x,y)             (((x)+((y)-1))/(y))
#define BA_ARRAY_INDEX(array,n) ((array)[(n) / BA_TYPE_BITS])
#define BA_BIT_INDEX(n)         (1 << ((n) % BA_TYPE_BITS))

#define BA_DEFINE(array,size) BA_TYPE array[BA_DIM(size, BA_TYPE_BITS)]

#define BA_ZERO(array,size)   memset((array), 0x00, sizeof(array))
#define BA_ONE(array,size)    memset((array), 0xFF, sizeof(array))
#define BA_SET(array, n)      (BA_ARRAY_INDEX(array, n) |= BA_BIT_INDEX(n))
#define BA_CLR(array, n)      (BA_ARRAY_INDEX(array, n) &= ~BA_BIT_INDEX(n))
#define BA_ISSET(array, n)    (BA_ARRAY_INDEX(array, n) & BA_BIT_INDEX(n))

#endif

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
 * Each session is identified by a random 8-byte session identifier.
 *
 * For efficiency, the session id is only transmitted over the control
 * channel (which only sees traffic occasionally when keys are being
 * negotiated).
 */

#if defined(USE_CRYPTO) && defined(USE_SSL)

#ifndef SESSION_ID_H
#define SESSION_ID_H

#include "basic.h"
#include "buffer.h"

struct session_id
{
  unsigned char id[8];
};

extern const struct session_id _session_id_zero;

static inline bool
session_id_equal (const struct session_id *sid1,
		  const struct session_id *sid2)
{
  return !memcmp (sid1, sid2, sizeof (struct session_id));
}

static inline bool
session_id_defined (const struct session_id *sid1)
{
  return memcmp (sid1, &_session_id_zero, sizeof (struct session_id) != 0);
}

static inline bool
session_id_read (struct session_id *sid, struct buffer *buf)
{
  return buf_read (buf, sid, sizeof (struct session_id));
}

static inline bool
session_id_write_prepend (const struct session_id *sid, struct buffer *buf)
{
  return buf_write_prepend (buf, sid, sizeof (struct session_id));
}

static inline bool
session_id_write (const struct session_id *sid, struct buffer *buf)
{
  return buf_write (buf, sid, sizeof (struct session_id));
}

void session_id_random (struct session_id *sid);

const char *session_id_print (const struct session_id *sid);

#endif /* SESSION_ID_H */
#endif /* USE_CRYPTO && USE_SSL */

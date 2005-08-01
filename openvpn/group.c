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

#ifdef WIN32
#include "config-win32.h"
#else
#include "config.h"
#endif

#include "syshead.h"

#if GROUPS

#include "ssl-x509.h"
#include "group.h"

static inline bool
group_ip_range_equal (const struct group_ip_range *r1, const struct group_ip_range *r2)
{
  return r1->start == r2->start && r1->end == r2->end;
}

static bool
group_pf_add_range (struct group_pf_list *list, const struct group_ip_range *range, struct gc_arena *gc)
{
  struct group_pf_entry **store = &list->list;
  while (*store)
    {
      struct group_pf_entry *e = *store;
      if (group_ip_range_equal (&e->range, range))
	return false;
      *store = e->next;
    }
  
  ALLOC_OBJ_CLEAR_GC (*store, struct group_pf_entry, gc);
  (*store)->range = *range;
  return true;
}

static bool
group_inherit_acl_group (struct group_entry *dest, const struct group_entry *src, const int n_groups)
{
  int i;
  bool modified = false;
  for (i = 0; i < n_groups; ++i)
    {
      if (BA_ISSET (src->acl_group, i) && !BA_ISSET (dest->acl_group, i))
	{
	  BA_SET (dest->acl_group, i);
	  modified = true;
	}
    }
  return modified;
}

static bool
group_inherit_acl_addr (struct group_entry *dest, const struct group_entry *src, struct gc_arena *gc)
{
  bool modified = false;
  const struct group_pf_entry *s = src->acl_addr.list;
  while (s)
    {
      modified |= group_pf_add_range (&dest->acl_addr, &s->range, gc);
      s = s->next;
    }
}

void
group_inherit_compile (struct group_info *gi, struct gc_arena *gc)
{
  int si, di;
  bool modified = true;
  int count = 0;

  while (modified)
    {
      modified = false;
      if (count >= MAX_INHERIT_PASSES)
	msg (M_FATAL, "Group Error: group inheritance graph failed to converge after %d passes", MAX_INHERIT_PASSES);
      for (di = 0; di < gi->n_groups; ++di)
	{
	  struct group_entry *dest = &gi->groups[di];
	  for (si = 0; si < gi->n_groups; ++si)
	    {
	      struct group_entry *src = &gi->groups[si];
	      if (di != si && BA_ISSET (dest->acl_inherit, si)) {
		modified |= group_inherit_acl_group (dest, src, gi->n_groups);
		modified |= group_inherit_acl_addr (dest, src, gc);
	      }
	    }
	}
      ++count;
    }
}

#ifdef ENABLE_DEBUG

/* dump group_info struct in human readable format */
void
group_info_print (const struct group_info *gi, const char *prefix, struct status_output *so)
{
}

#endif /* ENABLE_DEBUG */

#endif

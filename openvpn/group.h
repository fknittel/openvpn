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

#ifndef OPENVPN_GROUP_H
#define OPENVPN_GROUP_H

#if GROUPS

#include "bitarray.h"
#include "buffer.h"
#include "list.h"
#include "status.h"
#include "mroute.h"

#define MAX_INHERIT_PASSES 32

/*
 * Start and end IP addresses of pool range (endpoints included)
 */
struct group_ip_range
{
  in_addr_t start;
  in_addr_t end;
};

/*
 * Packet filter entry.  Allow packet if dest address is a member
 * of range.
 */
struct group_pf_entry
{
  struct group_ip_range range;
  struct group_pf_entry *next;
};

/*
 * List of packet filter entries.
 */
struct group_pf_list
{
  struct group_pf_entry *list;
};

/*
 * A single group definition
 */
struct group_entry
{
  /* group name */
  const char *name;

  /* set if at least one --group-ip-pool directive is specified */
  bool pool_enabled;

  /* set if at least one --group-acl directive is specified */
  bool acl_enabled;

  /* how to allocate IP addresses for this group */
  struct group_ip_range pool;

  /* which IP addresses can members of this group interact with */
  struct group_pf_list acl_addr;

  /* which groups can this group interact with if --client-to-client is enabled */
  BA_DEFINE(acl_group, MAX_GROUPS);

  /* which groups does this group inherit access permissions from */
  BA_DEFINE(acl_inherit, MAX_GROUPS);
};

/*
 * A single term in an AND equality test against
 * an x509 subject string.
 */
struct group_x509_match_term
{
  const char *name;  /* e.g. "C" */
  const char *value; /* e.g. "US" */
  struct group_x509_match_term *next;
};

/*
 * Represents a --group-by-x509 directive
 */
struct group_x509_match_list
{
  int group_index;                             /* assign client to this group if ... */
  struct group_x509_match_term *match_terms;   /* ... x509 subject matches against match_terms */
  struct group_x509_match *next;
};

struct group_info
{
# define GI_DEFINED          (1<<0)
# define GI_ENABLE_POOL      (1<<1)
# define GI_ACL_CC_BY_GROUP  (1<<2)
# define GI_ACL_CC_BY_ADDR   (1<<3)
# define GI_ACL_CS_BY_ADDR   (1<<4)
  unsigned int flags;

  /* these two parameters are mutually exclusive -- only
     one should be defined */
  struct group_x509_match_list *group_by_x509; /* --group-by-x509 directive */

  /* set by --group */
  int default_group;

  /* group list */
  int n_groups;
  struct group_entry groups[MAX_GROUPS];
};

/*
 * Represents server state for one group
 */
struct group_server_context_entry
{
  struct hash *pf_cache; /* packet filter cache */
};

/*
 * Represents server-wide state
 */
struct group_server_context
{
  struct group_server_context_entry groups[MAX_GROUPS];
};

/*
 * Valid source addresses for a given client.  Used in
 * TAP mode only.
 */
struct client_valid_src
{
  struct group_pf_list *list;
  struct hash *cache;
};

/*
 * Group state per-client, placed in client's context_2 state
 */
struct group_context
{
  /* inherited from group_info::flags */
# define GC_OWNED             (1<<16)
# define GC_VALID_SRC_DEFINED (1<<17)
  unsigned int flags;

  const struct group_info *info;
  struct group_server_context *server_context;
  int group_index;

  struct client_valid_src valid_src; /* per-client */
};

/*
 * Public functions
 */

struct group_info *group_info_new (struct gc_arena *gc);

/* process config file directives */

void group_ip_pool    (struct group_info *gi, struct gc_arena *gc, const int msglevel,
		      const char *group_name, const char *start_ip, const char *end_ip);

void group_by_x509    (struct group_info *gi, struct gc_arena *gc, const int msglevel,
		      const char *group_name, char *parms[]);

int  group_default    (struct group_info *gi, struct gc_arena *gc, const int msglevel,
		      const char *group_name);

void group_acl_enable (struct group_info *gi, struct gc_arena *gc, const int msglevel,
		      char *parms[]);

void group_acl        (struct group_info *gi, struct gc_arena *gc, const int msglevel,
		      const char *group_name, char *parms[]);

/* compile group inheritance graph */
void group_inherit_compile (struct group_info *gi, struct gc_arena *gc);

/* struct group_context */

void group_context_init (struct group_context *g, const struct group_info *gi, const bool server);
void group_context_inherit (struct group_context *dest, const struct group_context *src);
void group_context_detach (struct group_context *g);
void group_context_close (struct group_context *g);

bool group_context_get_ip_range (const struct group_context *g, struct group_ip_range *range);
void group_context_set_x509_subject (struct group_context *g, const char *subject);
bool group_context_allow_acl_by_group (const struct group_context *g, const struct group_context *group);
bool group_context_allow_acl_by_addr (struct group_context *g, const struct mroute_addr *addr);

bool group_context_allow_acl (struct group_context *src, const struct mroute_addr *srcaddr,
			      struct group_context *dest, const struct mroute_addr *destaddr);

/* client source address packet filter, used for --iroute in TAP mode */

void client_src_subnet_add (struct group_context *g, const in_addr_t network, const int netbits);
bool client_src_addr_validate (struct group_context *g, const struct mroute_addr *srcaddr);

#ifdef ENABLE_DEBUG

/* dump group_info struct in human readable format */
void group_info_print (const struct group_info *gi, const char *prefix, struct status_output *so);

#endif /* ENABLE_DEBUG */

/*
 * Inline functions
 */

static inline bool
group_context_defined (const struct group_context *gc)
{
  return gc && (gc->flags & GI_DEFINED) != 0;
}

static inline bool
group_context_id_defined (const struct group_context *gc)
{
  return group_context_defined (gc) && gc->group_index >= 0;
}

#endif
#endif

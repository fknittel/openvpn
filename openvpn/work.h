/*
 *  OpenVPN -- An application to securely tunnel IP networks
 *             over a single UDP port, with support for SSL/TLS-based
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

#ifndef WORK_H
#define WORK_H

#include "thread.h"
#include "ssl.h"
#include "misc.h"
#include "plugin.h"
#include "event.h"

/* Maximum recursive nesting of work_thread_event_loop_t virtual function */
#define WORK_THREAD_MAX_RECURSION 16 /* if you increase to 255 or more, code changes required */

/* status return from work_thread_event_loop_t virtual function */
#define WT_EVENT_LOOP_NORMAL 0
#define WT_EVENT_LOOP_BREAK  1

/* the event loop */
typedef int (*work_thread_event_loop_t) (void *arg);

/*
 * Thread lock levels
 */
#define TL_INACTIVE 0 /* work thread is inactive on this instance */
#define TL_LIGHT    1 /* instance is partially locked -- tunnel packet forwarding is still allowed */
#define TL_FULL     2 /* instance is fully locked by work thread */

#ifdef USE_PTHREAD

#define WORK_THREAD_EXIT_MESSAGE 0xFF

#if WORK_THREAD_MAX_RECURSION >= WORK_THREAD_EXIT_MESSAGE
#error WORK_THREAD_MAX_RECURSION is too large
#endif

/*
 * Thread Context flags.  Must not collide with 
 * S_ flags in misc.h, MULTI_ flags in multi.h, or
 * SP_/MGI_ flags in multi.c.
 */
#define WTF_LIGHT (1<<30) /* put a "light" lock on instance, still maintaining tunnel forwarding */

struct thread_context
{
  int thread_level;     /* TL_ levels */
  unsigned int flags;   /* WTF_, S_, MULTI_, SP_, or MGI_ flags */
  void *arg1;           /* mode specific -- in server mode will be a struct multi_context */
  void *arg2;           /* mode specific -- in server mode will be a struct multi_instance */

  int event_loop_retval;

  /* save state prior to recursive call into event loop, return state * */
  void *(*save)(struct thread_context *context);

  /* recursive return state restore, passed saved state */
  void (*restore)(struct thread_context *context, void *save_data);
};

/* Index into sd */
#define MAIN_THREAD   0
#define WORKER_THREAD 1

/* Message types */
#define TMS_BIO_READ  0
#define TMS_SCRIPT    1
#define TMS_PLUGIN    2

struct thread_message_bio_read {
  BIO *b;
  void *buf;
  int len;
  int ret;
};

struct thread_message_execute_script {
  const char *command;
  const struct env_set *es;
  unsigned int flags;
  int ret;
};

struct thread_message_execute_plugin {
  const struct plugin_list *pl;
  int type;
  const char *args;
  struct env_set *es;
  int ret;
};

union thread_message_union {
  struct thread_message_bio_read bio_read;
  struct thread_message_execute_script script;
  struct thread_message_execute_plugin plugin;
};


struct thread_message {
# define TMS_UNDEF    0
# define TMS_PENDING  1
# define TMS_FINISHED 2
  int state;

  int type;
  union thread_message_union u;
};

/*
 * Object that represents a work thread.
 */
struct work_thread
{
  int n_threads;
  int nice_work;
  openvpn_thread_t thread_id;
  int sd[2];                   /* PF_UNIX/SOCK_DGRAM communication handles */
  bool brk;
  void *event_loop_arg;
  work_thread_event_loop_t event_loop;
  int stack_index;
  struct thread_message stack[WORK_THREAD_MAX_RECURSION];
};

struct work_thread *work_thread_init (const int n_threads, const int nice_work);

void work_thread_close (struct work_thread *wt);

void work_thread_enable (struct work_thread *wt, void *arg, work_thread_event_loop_t event_loop);

void work_thread_disable (struct work_thread *wt);

bool work_thread_break (struct work_thread *wt);

/*
 * Functions implemented by work thread
 */

int work_thread_bio_read (struct work_thread *wt,
			  struct thread_context *tc,
			  BIO *b,
			  void *buf,
			  const int len);

int work_thread_system (struct work_thread *wt,
			struct thread_context *tc,
			const char *command,
			const struct env_set *es,
			const unsigned int flags);

bool work_thread_system_check (struct work_thread *wt,
			       struct thread_context *tc,
			       const char *command,
			       const struct env_set *es,
			       const unsigned int flags,
			       const char *error_message);

int work_thread_plugin_call (struct work_thread *wt,
			     struct thread_context *tc,
			     const struct plugin_list *pl,
			     const int type,
			     const char *args,
			     struct env_set *es);

/*
 * Inline functions
 */

static inline bool
work_thread_ready_level (const struct thread_context *tc, const int thread_level)
{
  return tc->thread_level <= thread_level;
}

static inline bool
work_thread_ready (const struct thread_context *tc)
{
  return work_thread_ready_level (tc, TL_INACTIVE);
}

static inline void
work_thread_socket_set (struct work_thread *wt,
			struct event_set *es,
			void *arg,
			unsigned int *persistent)
{
  void work_thread_socket_set_dowork (struct work_thread *wt,
				      struct event_set *es,
				      void *arg,
				      unsigned int *persistent);

  if (wt->stack_index > 0 || wt->brk || persistent)
    work_thread_socket_set_dowork (wt, es, arg, persistent);
}

#endif /* USE_PTHREAD */
#endif /* WORK_H */

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

#define WT_EVENT_LOOP_NORMAL 0
#define WT_EVENT_LOOP_BREAK  1

typedef int (*work_thread_event_loop_t) (void *arg);

#ifdef USE_PTHREAD

#define WORK_THREAD_MAX_RECURSION 16

/* Index into sd */
#define MAIN_THREAD   0
#define WORKER_THREAD 1

/* Message types */
#define TMS_SHUTDOWN  0
#define TMS_BIO_READ  1
#define TMS_SCRIPT    2
#define TMS_PLUGIN    3

struct thread_message_bio_read {
  BIO *b;
  void *buf;
  int len;
  int ret;
};

struct thread_message_execute_script {
  char *command;
  const char **envp;
  unsigned int flags;
  int ret;
};

struct thread_message_execute_plugin {
  openvpn_plugin_func_v1 func;
  openvpn_plugin_handle_t handle;
  int type;
  const char **argv;
  const char **envp;
  int ret;
};

union thread_message_union {
  struct thread_message_bio_read bio_read;
  struct thread_message_execute_script script;
  struct thread_message_execute_plugin plugin;
};

struct thread_message {
  int type;
  int sd[2];                   /* PF_UNIX/SOCK_DGRAM communication handles */
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
  void *event_loop_arg;
  work_thread_event_loop_t event_loop;
  int stack_index;
  struct thread_message stack[WORK_THREAD_MAX_RECURSION];
};

struct work_thread *work_thread_init (const int n_threads, const int nice_work);

void work_thread_close (struct work_thread *wt);

void work_thread_socket_set (struct work_thread *wt,
			     struct event_set *es,
			     void *arg,
			     unsigned int *persistent);

void work_thread_enable (struct work_thread *wt, void *arg, work_thread_event_loop_t event_loop);

void work_thread_disable (struct work_thread *wt);

/*
 * Functions implemented by work thread
 */

int work_thread_bio_read (struct work_thread *wt,
			  BIO *b,
			  void *buf,
			  const int len);

int work_thread_system (struct work_thread *wt,
			const char *command,
			const struct env_set *es,
			const unsigned int flags);

int work_thread_system_check (struct work_thread *wt,
			      const char *command,
			      const struct env_set *es,
			      const unsigned int flags,
			      const char *error_message);

int work_thread_plugin_call (struct work_thread *wt,
			     const struct plugin_list *pl,
			     const int type,
			     const char *args,
			     struct env_set *es);

#endif /* USE_PTHREAD */
#endif /* WORK_H */

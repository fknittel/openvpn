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

#ifdef WIN32
#include "config-win32.h"
#else
#include "config.h"
#endif

#include "syshead.h"

#ifdef USE_PTHREAD

#include "error.h"
#include "buffer.h"
#include "fdmisc.h"
#include "work.h"

#include "memdbg.h"

/*
 * Socket read/write functions.
 */

static int
recv_uchar (const int fd)
{
  unsigned char c;
  const ssize_t size = read (fd, &c, sizeof (c));
  if (size == sizeof (c))
    return c;
  else
    return -1;
}

static int
send_uchar (const int fd, const int value)
{
  ssize_t size;
  unsigned char c;

  ASSERT (value >= 0 && value <= 0xFF);
  c = (unsigned char) value;
  size = write (fd, &c, sizeof (c));
  if (size == sizeof (c))
    return (int) size;
  else
    return -1;
}

static inline struct thread_message *
thread_message_index (struct work_thread *wt, const int i)
{
  if (i < 0 || i >= WORK_THREAD_MAX_RECURSION)
    return NULL;
  return &wt->stack[i];
}

static struct thread_message *
alloc_thread_message (struct work_thread *wt)
{
  struct thread_message *ret = thread_message_index (wt, wt->stack_index);
  if (ret)
    CLEAR (*ret);
  return ret;
}

static inline struct thread_message *
current_thread_message (struct work_thread *wt, const int si_offset)
{
  return thread_message_index (wt, wt->stack_index + si_offset);
}

static void
post_message_to_work_thread (struct work_thread *wt)
{
  struct thread_message *tm = thread_message_index (wt, wt->stack_index);

  ASSERT (tm);
  tm->state = TMS_PENDING;
  dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: sending task code %d, si=%d", tm->type, wt->stack_index);
  send_uchar (wt->sd[MAIN_THREAD], wt->stack_index);
  ++wt->stack_index;
}

static void
post_exit_message_to_work_thread (struct work_thread *wt)
{
  if (wt->sd[MAIN_THREAD] >= 0)
    send_uchar (wt->sd[MAIN_THREAD], WORK_THREAD_EXIT_MESSAGE);
}

static void
pop_thread_message (struct work_thread *wt)
{
  ASSERT (--wt->stack_index >= 0);
  dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: pop thread message, si=%d", wt->stack_index);
}

/*
 * This is the work thread function
 */
static void *
thread_func (void *arg)
{
  struct work_thread *wt = (struct work_thread *) arg;

  msg (D_WORK_THREAD, "WORK THREAD: starting");

  /* change thread priority if requested */
  set_nice (wt->nice_work);

  while (true)
    {
      struct thread_message *tm;

      /*
       * Receive message from parent.
       */
      const int cmd = recv_uchar (wt->sd[WORKER_THREAD]);

      dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: received command in stack index %d", cmd);

      if (cmd == WORK_THREAD_EXIT_MESSAGE)
	break;
      if (cmd == -1)
	{
	  msg (D_WORK_THREAD_ERRORS, "WORK THREAD: read error on parent socket");
	  continue;
	}

      /*
       * Dispatch action.
       */
      tm = thread_message_index (wt, cmd);
      ASSERT (tm);
      ASSERT (tm->state == TMS_PENDING);
      switch (tm->type)
	{
	case TMS_BIO_READ:
	  sleep_milliseconds (250); // JYFIXME
	  tm->u.bio_read.ret = BIO_read (tm->u.bio_read.b,
					 tm->u.bio_read.buf,
					 tm->u.bio_read.len);
	  break;
	case TMS_SCRIPT:
	  tm->u.script.ret = openvpn_system (tm->u.script.command,
					     tm->u.script.es,
					     tm->u.script.flags);
	  break;
	case TMS_PLUGIN:
	  tm->u.plugin.ret = plugin_call (tm->u.plugin.pl,
					  tm->u.plugin.type,
					  tm->u.plugin.args,
					  tm->u.plugin.es);
	  break;
	default:
	  ASSERT (0);
	}

      /*
       * Tell parent we're finished.
       */
      if (send_uchar (wt->sd[WORKER_THREAD], cmd) == -1)
	msg (D_WORK_THREAD_ERRORS, "WORK THREAD: write error on parent socket");

      dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: task completed");
    }

  msg (D_WORK_THREAD, "WORK THREAD: exiting");

  return NULL;
}

static void
open_socketpair (int *sd)
{
  if (socketpair (PF_UNIX, SOCK_DGRAM, 0, sd) == -1)
    msg (M_ERR, "socketpair call failed");

  /* set socket properties */
  set_nonblock (sd[MAIN_THREAD]);
  set_cloexec (sd[MAIN_THREAD]);
  set_cloexec (sd[WORKER_THREAD]);
}

static void
close_socketpair (int *sd)
{
  if (sd[0] >= 0)
    {
      close (sd[0]);
      sd[0] = -1;
    }
  if (sd[1] >= 0)
    {
      close (sd[1]);
      sd[1] = -1;
    }
}

struct work_thread *
work_thread_init (const int n_threads, const int nice_work)
{
  struct work_thread *wt;

  ASSERT (n_threads >= 2);

  if (n_threads > 2)
    msg (M_INFO, "NOTE: currently a maximum of 2 threads is supported");

  ALLOC_OBJ_CLEAR (wt, struct work_thread);
  wt->n_threads = n_threads;
  wt->nice_work = nice_work;
  
  /* initialize pthread */
  openvpn_thread_init ();

  /*
   * Make a socket for foreground and background threads
   * to communicate.  The background thread will set its
   * end to blocking, while the foreground will set its
   * end to non-blocking.
   */
  open_socketpair (wt->sd);

  wt->thread_id = openvpn_thread_create (thread_func, (void*)wt);

  return wt;
}

void
work_thread_close (struct work_thread *wt)
{
  if (wt)
    {
      post_exit_message_to_work_thread (wt);
      openvpn_thread_join (wt->thread_id);
      openvpn_thread_cleanup ();
      close_socketpair (wt->sd);

      free (wt);
    }
}

void
work_thread_enable (struct work_thread *wt, void *arg, work_thread_event_loop_t event_loop)
{
  wt->event_loop_arg = arg;
  wt->event_loop = event_loop;
}

void
work_thread_disable (struct work_thread *wt)
{
  wt->event_loop_arg = NULL;
  wt->event_loop = NULL;
}

/*
 * Functions to synchronize main thread event loop
 * with worker thread.
 */

static void
set_es_override (struct event_set *es, void *arg)
{
  struct event_set_return esr;
  esr.rwflags = EVENT_READ;
  esr.arg = arg;
  event_ctl_override (es, &esr); /* force event_wait to return with EVENT_READ status */
}

void
work_thread_socket_set_dowork (struct work_thread *wt,
			       struct event_set *es,
			       void *arg,
			       unsigned int *persistent)
{
  if (!wt->brk)
    {
      unsigned int enable = 0;
      struct thread_message *tm = current_thread_message (wt, -1);

      if (tm)
	{
	  switch (tm->state)
	    {
	    case TMS_PENDING:
	      /*
	       * If work thread task is still pending,
	       * tell event loop to wait on work thread
	       * communication socket.
     	       */
	      enable = 1;
	      break;

	    case TMS_FINISHED:
	      /*
	       * If work thread task is finished, override
	       * event_wait so that it will return an
	       * immediate EVENT_READ on the work
	       * thread socket.  This in turn will
	       * cause work_thread_break to be called
	       * which can then return a true status
	       * to the calling event loop, causing it
	       * to return to the recursive event
	       * loop level which queued the work thread
	       * task.
	       */
	      set_es_override (es, arg);
	      break;
	    default:
	      ASSERT (0);
	    }
	}

      if (persistent)
	{
	  if (enable != *persistent)
	    {
	      if (enable)
		event_ctl (es, wt->sd[MAIN_THREAD], EVENT_READ, arg);
	      else
		event_del (es, wt->sd[MAIN_THREAD]);
	      *persistent = enable;
	    }
	}
      else if (enable)
	{
	  event_ctl (es, wt->sd[MAIN_THREAD], EVENT_READ, arg);
	}
    }
  else
    set_es_override (es, arg);
}

static bool
work_thread_break_dowork (struct work_thread *wt, const int si_offset)
{
  struct thread_message *tm;
  int msg;

  /*
   * Read data from work thread socket.
   * A received uchar indicates that the
   * corresponding stack index has transitioned
   * to TMS_FINISHED.
   */
  while (true)
    {
      msg = recv_uchar (wt->sd[MAIN_THREAD]);
      if (msg == -1)
	break;
      tm = thread_message_index (wt, msg);
      if (tm)
	{
	  ASSERT (tm->state == TMS_PENDING);
	  tm->state = TMS_FINISHED;
	}
    }

  /*
   * Determine if the event loop which called us
   * should break.
   */
  tm = current_thread_message (wt, si_offset);
  if (tm)
    {
      switch (tm->state)
	{
	case TMS_PENDING:
	  break;
	case TMS_FINISHED:
	  return true;
	default:
	  ASSERT (0);
	}
    }
  return false;
}

bool
work_thread_break (struct work_thread *wt)
{
  if (wt->brk)
    return true;
  else
    return work_thread_break_dowork (wt, -1);
}

static void
work_thread_flush (struct work_thread *wt)
{
  while (!work_thread_break_dowork (wt, 0))
    {
      dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: Break event, waiting for work thread to finish");
      openvpn_sleep (1);
    }
}

static int
work_thread_event_loop (struct work_thread *wt)
{
  int ret;
  dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: event loop ENTER, si=%d", wt->stack_index);
  ret = (*wt->event_loop)(wt->event_loop_arg);
  dmsg (D_WORK_THREAD_DEBUG, "WORK THREAD: event loop EXIT, si=%d", wt->stack_index);
  return ret;
}

/*
 * Functions implemented by work thread
 */

int
work_thread_bio_read (struct work_thread *wt,
		      struct thread_context *tc,
		      BIO *b,
		      void *buf,
		      const int len)
{
  if (!wt->brk)
    {
      if (wt->event_loop)
	{
	  struct thread_message *tm = alloc_thread_message (wt);
	  if (tm)
	    {
	      void *save_data = (*tc->save)(tc);
	      
	      tm->type = TMS_BIO_READ;
	      tm->u.bio_read.b = b;
	      tm->u.bio_read.buf = buf;
	      tm->u.bio_read.len = len;
	      
	      post_message_to_work_thread (wt);
	      tc->event_loop_retval = work_thread_event_loop (wt);
	      pop_thread_message (wt);
	      (*tc->restore)(tc, save_data);

	      switch (tc->event_loop_retval)
		{
		case WT_EVENT_LOOP_BREAK:
		  return tm->u.bio_read.ret;
		case WT_EVENT_LOOP_NORMAL:
		  wt->brk = true;
		  work_thread_flush (wt);
		  return -1;
		default:
		  ASSERT (0);
		}
	    }
	}
      return BIO_read (b, buf, len);
    }
  return -1;
}

int
work_thread_system (struct work_thread *wt,
		    struct thread_context *tc,
		    const char *command,
		    const struct env_set *es,
		    const unsigned int flags)
{
  if (!wt->brk)
    {
      if (wt->event_loop)
	{
	  struct thread_message *tm = alloc_thread_message (wt);
	  if (tm)
	    {
	      void *save_data = (*tc->save)(tc);

	      tm->type = TMS_SCRIPT;
	      tm->u.script.command = command;
	      tm->u.script.es = es;
	      tm->u.script.flags = flags;

	      post_message_to_work_thread (wt);
	      tc->event_loop_retval = work_thread_event_loop (wt);
	      pop_thread_message (wt);
	      (*tc->restore)(tc, save_data);

	      switch (tc->event_loop_retval)
		{
		case WT_EVENT_LOOP_BREAK:
		  return tm->u.script.ret;
		case WT_EVENT_LOOP_NORMAL:
		  wt->brk = true;
		  work_thread_flush (wt);
		  return -1;
		default:
		  ASSERT (0);
		}
	    }
	}
      return openvpn_system (command, es, flags);
    }
  return -1;
}

bool
work_thread_system_check (struct work_thread *wt,
			  struct thread_context *tc,
			  const char *command,
			  const struct env_set *es,
			  const unsigned int flags,
			  const char *error_message)
{
  const int wts_ret = work_thread_system (wt, tc, command, es, flags);
  return system_check_error (wts_ret, flags, error_message);
}

int
work_thread_plugin_call (struct work_thread *wt,
			 struct thread_context *tc,
			 const struct plugin_list *pl,
			 const int type,
			 const char *args,
			 struct env_set *es)
{
  if (!wt->brk)
    {
      if (wt->event_loop)
	{
	  struct thread_message *tm = alloc_thread_message (wt);
	  if (tm)
	    {
	      void *save_data = (*tc->save)(tc);

	      tm->type = TMS_PLUGIN;
	      tm->u.plugin.pl = pl;
	      tm->u.plugin.type = type;
	      tm->u.plugin.args = args;
	      tm->u.plugin.es = es;
	  
	      post_message_to_work_thread (wt);
	      tc->event_loop_retval = work_thread_event_loop (wt);
	      pop_thread_message (wt);
	      (*tc->restore)(tc, save_data);

	      switch (tc->event_loop_retval)
		{
		case WT_EVENT_LOOP_BREAK:
		  return tm->u.plugin.ret;
		case WT_EVENT_LOOP_NORMAL:
		  wt->brk = true;
		  work_thread_flush (wt);
		  return 1;
		default:
		  ASSERT (0);
		}
	    }
	}
      return plugin_call (pl, type, args, es);
    }
  return 1;
}

#else
static void dummy(void) {}
#endif

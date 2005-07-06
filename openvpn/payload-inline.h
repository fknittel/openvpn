#ifndef _PAYLOAD_INLINE_H
#define _PAYLOAD_INLINE_H

#ifdef USE_PAYLOAD_CONNTRACK
/*
 * Should run conntrack GC?
 */
static inline void
check_payload_gc (struct context *c)
{
  void check_payload_gc_dowork (struct context *c);
  if (!c->c2.payload_context)
    return;
  if (!event_timeout_defined (&c->c2.payload_context->gc.wakeup))
    return;
  if (event_timeout_trigger (&c->c2.payload_context->gc.wakeup, &c->c2.timeval, ETT_DEFAULT))
    check_payload_gc_dowork (c);
}

#endif
#endif /* _PAYLOAD_INLINE_H */
/* 
vim: cino={.5s,\:.5s,+.5s,t0,g0,^-2,e-2,n-2,p2s,(0,=.5s:sw=4:ts=8:sts=4
 */

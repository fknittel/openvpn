#ifndef _PAYLOAD_H
#define _PAYLOAD_H
struct payload_context {
  struct hash *hash;
  struct schedule *schedule; /* unused by now */
  int max_tcp_conns; /* unused by now */
  int tcp_retrans;
  struct {
    time_t per_second_trigger;
    int bucket_base;
    int buckets_per_pass;
    struct event_timeout wakeup;
  } gc;
};

struct payload_context * payload_new(int tcp_retrans);
void payload_free(struct payload_context *);
int payload_tcp_retrans_drop(struct context *c, struct buffer *buf);
#endif /* _PAYLOAD_H */
/* 
vim: cino={.5s,\:.5s,+.5s,t0,g0,^-2,e-2,n-2,p2s,(0,=.5s:sw=4:ts=8:sts=4
 */

/* arith_bingcd.c
 * by Joppe W. Bos and Michael Naehrig (mnaehrig@microsoft.com), 
 * Cryptography Research Group, MSR Redmond, 2014
 * 
 * This file is part of the ARITH library version 0.01
 *
 * DISCLAIMER:
 * This is RESEARCH code.
 * Parts of the code might be incomplete.
 * Please use at your own risk.
 * DO NOT USE IN PRODUCTION CODE!!!
 * This code is still under active development.
 */

#include "arith.h"
#include "sint.h"
#include "uint.h"

/* Binary extended GCD from Knuth, section 4.5.2 exercise 39.
u1 = 1;
u2 = 0;
u3 = u;
v1 = v;
v2 = 1-u;
v3 = v;
if (u is odd) {
  t1 = 0;
  t2 = -1;
  t3 = -v;
  goto Y4;
} else {
  t1 = 1;
  t2 = 0;
  t3 = u;
}

Y3:
if (t1 is even and t2 is even) {
  t1 = t1/2;
  t2 = t2/2;
  t3 = t3/2;
} else {
  t1 = (t1+v) / 2;
  t2 = (t2-u)/2;
  t3 = t3 / 2;
}
Y4:
if (t3 is even) goto Y3;

if (t3 > 0) {
  u1 = t1;
  u2 = t2;
  u3 = t3;
} else {
  v1 = v - t1;
  v2 = -u - t2;
  v3 = -t3;
}

t1 = u1 - v1;
t2 = u2 - v2;
t3 = u3 - v3;

if (t1 <= 0) {
  t1 = t1 + v;
  t2 = t2 - u;
}
if (t3 != 0) goto Y3;

return (u1, u2, u3*2^k);
*/

/* S = u^-1 mod v, This only works for odd v */
void binary_gcd (uint64_t *S, uint64_t *U, uint64_t *V, uint64_t *mem, int n) {
  int i;
  sint_t u1, u3, v1, v3, t1, t3, zero, u, v;

  /* Increase the limb size with one since addition might produce a carry. */
  u1->n   = mem + 0*(n+1);
  u3->n   = mem + 1*(n+1);
  v1->n   = mem + 2*(n+1);
  v3->n   = mem + 3*(n+1);
  t1->n   = mem + 4*(n+1);
  t3->n   = mem + 5*(n+1);
  u->n    = mem + 6*(n+1);
  v->n    = mem + 7*(n+1);
  zero->n = mem + 8*(n+1);
  
  for (i=0; i < n; i++) {
    u->n[i] = U[i];
    v->n[i] = V[i];
  }
  u->n[n] = (uint64_t) 0;
  v->n[n] = (uint64_t) 0;
  v->sign = 1;
  u->sign = 1;

  for (i=0; i < n+1; i++) zero->n[i] = (uint64_t) 0;
  u1->sign = 1; u3->sign = 1;
  v1->sign = 1; v3->sign = 1;

  u1->n[0] = (uint64_t) 1; for (i=1; i < n+1; i++) u1->n[i] = (uint64_t) 0;  

  for (i=0; i < n+1; i++) u3->n[i] = u->n[i];

  for (i=0; i < n+1; i++) v1->n[i] = v->n[i];
  for (i=0; i < n+1; i++) v3->n[i] = v->n[i];

  if (u->n[0]%2 == 1) {
    for (i=0; i < n+1; i++) t1->n[i] = (uint64_t) 0;
    t1->sign = 1;

    for (i=0; i < n+1; i++) t3->n[i] = v->n[i];
    t3->sign = -1;
  } else {
    t1->n[0] = (uint64_t) 1; for (i=1; i < n+1; i++) t1->n[i] = (uint64_t) 0;
    t1->sign = 1;

    for (i=0; i < n+1; i++) t3->n[i] = u->n[i];
    t3->sign = 1;

    if ((t1->n[0]%2==0)) {
      sint_div2 (t1,t1,n+1);
      sint_div2 (t3,t3,n+1);
    } else {
      sint_add (t1, t1, v, n+1);
      sint_div2 (t1,t1,n+1);
      sint_div2 (t3,t3,n+1);
    }
  }

  do {
    while (t3->n[0]%2 == 0) {    
      if ((t1->n[0]%2 == 0)) {
        sint_div2 (t1,t1,n+1);
        sint_div2 (t3,t3,n+1);
      } else {
        sint_add (t1, t1, v, n+1);
        sint_div2 (t1,t1,n+1);
        sint_div2 (t3,t3,n+1);
      }
    }

    if (uint_cmp (t3->n,zero->n,n+1) != 0 && t3->sign > 0) {
      for (i=0; i < n+1; i++) u1->n[i] = t1->n[i];
      for (i=0; i < n+1; i++) u3->n[i] = t3->n[i];
    } else {
      sint_sub (v1,v,t1,n+1);
      for (i=0; i < n+1; i++) v3->n[i] = t3->n[i];
      v3->sign = -t3->sign;
    }

    sint_sub (t1, u1, v1, n+1);
    sint_sub (t3, u3, v3, n+1);

    if (t1->sign < 0 || uint_cmp (t1->n, zero->n, n+1) == 0) {
      sint_add (t1, t1, v, n+1);
    } 
  } while (uint_cmp (t3->n,zero->n,n+1) != 0);

  for (i=0; i < n; i++) S[i] = u1->n[i];
}

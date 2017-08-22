/* mont_arith.c
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

#include "mont_arith.h"

void mont2normal (uint64_t *c, uint64_t *a, montmul_t *m, void *mem) {
  memset (c, 0, m->n * sizeof(uint64_t));
  c[0] = 1;
  montmul_ct (c,a,c,m,mem);
}

void normal2mont (uint64_t *c, uint64_t *a, montmul_t *m, void *mem) {
  montmul_ct (c, a, m->R2, m, mem);
}


/* Perform the final Montgomery subtraction.
 * We do this in constant-time.
 * 'a' can be a limb longer compared to the modulus used. 
 * Return 1 if a >= Fp_config->p
 * return 0 if a <  Fp_config->p
 */
static int mont_cmp (uint64_t *a, montmul_t *m) {
  int ret = 0;

  /* 'a' can be a limb longer than the prime */ 
  ret = (a[m->n] != 0);
  ret |= uint_cmpge_ct (a, m->p, m->n);
  return ret;
}

void montmul_ct (uint64_t *c, uint64_t *a, uint64_t *b, montmul_t *m, void *tmp) { 
  int i;
  uint64_t q, m1, m2, *mem;
  mem = (uint64_t *) tmp;
  
  mul_1_n (mem, a[0], b, m->n);
  q = m->mu * (*mem);
  mul_1_n_add_m_div1 (mem, q, m->p, m->n, m->n+1);
  
  for (i=1; i < m->n; i++) {
    mem[m->n] += mul_1_n_add_n (mem, a[i], b, m->n);
    q = m->mu * (*mem);
    mul_1_n_add_m_div1 (mem, q, m->p, m->n, m->n+1);
  }
  
  /* We want to compute the conditional subtraction:
   *   if (c > 0) {
   *     c = c - p
   *   }
   * But we remove the branch in order to run in constant-time.
   */
  m1 = (mont_cmp (mem, m) > 0);
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_sub (c, mem, m->p, m->n);
  for (i=0; i < m->n; i++) {
    c[i] = (mem[i] & m2) | (c[i] & m1);
  }
}

void montred_ct (uint64_t *c, uint64_t *a, montmul_t *m, void *tmp) { 
  int i, size;
  uint64_t q, m1, m2, *mem;
  mem = (uint64_t *) tmp;
  
  size = 2*m->n;
  q = m->mu * a[0];
  for (i=0; i < size; i++) mem[i] = a[i];
  mem[size] = 0;
  size++;

  mul_1_n_add_m_div1 (mem, q, m->p, m->n, size);
  size--;
  
  for (i=1; i < m->n; i++) {
    q = m->mu * (*mem);
    mul_1_n_add_m_div1 (mem, q, m->p, m->n, size);
    size--;
  }
  
  /* We want to compute the conditional subtraction:
   *   if (c > 0) {
   *     c = c - p
   *   }
   * But we remove the branch to protect against side-channel resistance.
   */
  m1 = (mont_cmp (mem, m) > 0);
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_sub (c, mem, m->p, m->n);
  for (i=0; i < m->n; i++) {
    c[i] = (mem[i] & m2) | (c[i] & m1);
  }
}

/* tmp should contain enough space for n 64-bit limbs */
void mod_add (uint64_t *c, uint64_t *a, uint64_t *b, montmul_t *m, void *tmp) {
  uint64_t carry, m1, m2, *mem;
  int i;
 
  mem = (uint64_t *) tmp;
  carry = uint_add (c, a, b, m->n);

  /* Compute the conditional-subtraction:
   * if (carry || uint_cmpgt (c, Fp_config->p) > 0) {
   *   uint_sub (c, c, Fp_config->p, Fp_config->n);
   * }
   * using straight-line code
   */  
  m1 = carry;
  m1 |= uint_cmpge_ct (c, m->p, m->n);  
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_sub (mem, c, m->p, m->n);
  for (i=0; i < m->n; i++) {
    c[i] = (mem[i] & m1) | (c[i] & m2);  
  }
}

/* tmp should contain enough space for n 64-bit limbs */
void mod_sub (uint64_t *c, uint64_t *a, uint64_t *b, montmul_t *m, void *tmp) {
  uint64_t borrow, m1, m2, *mem;
  int i;
  
  mem = (uint64_t *) tmp;
  borrow = uint_sub (c, a, b, m->n);
  
  /* Compute the conditional-addition:
   * if (borrow) {
   *   c = c + p;
   * }
   * using straight-line code
   */  
  m1 = borrow;
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_add (mem, c, m->p, m->n);
  for (i=0; i < m->n; i++) {
    c[i] = (mem[i] & m1) | (c[i] & m2);  
  }
}

void mont_sqr (uint64_t *c, uint64_t *a, montmul_t *m, void *tmp) {
  uint64_t u, v, *p, t, t0, t1, t2, q, c1, c2, m1, m2;
  int i, j, n = m->n;

  p = (uint64_t *) tmp;
  for (i=0; i < n+1; i++) p[i] = (uint64_t) 0;

  for (i=0; i < n; i++) {
    v = p[i];
    u=0;
    muladd (v,u,a[i],a[i]);
    p[i] = v;
    t = 0;
    for (j=i+1; j < n; j++) {
      mul (t0,t1,a[i],a[j]);
      t2 = (t1 >> (uint64_t) 63);
      t1 = ((t1 << (uint64_t) 1) | (t0 >> (uint64_t) 63));
      t0 = (t0 << (uint64_t) 1);
      t0 = t0 + u;
      c1 = (t0 < u);
      t1 = t1 + c1;
      c2 = (t1 < c1);
      t2 = t2 + c2;

      t0 = t0 + p[j];
      c1 = (t0 < p[j]);
      t1 = t1 + c1;
      c2 = (t1 < c1);
      t2 = t2 + c2;

      t1 = t1 + t;
      c1 = (t1 < t);
      t2 = t2 + c1;

      p[j] = t0;
      u = t1;
      t = t2;
    }
    v = u + p[m->n];
    u = t + (v<u);
    p[n] = v;
    p[n+1] = u;
    q = p[0] * m->mu;
    v=p[0];
    u=0;
    muladd (v,u,m->p[0],q);
    for (j=1; j < n; j++) {
      v=p[j];
      muladdadd (v,u,m->p[j],q);
      p[j-1] = v;
    }
    v = p[n] + u;
    u = (v < u);
    p[n-1] = v;
    p[n] = p[m->n+1] + u;
  }

  /* We want to compute the conditional subtraction:
   *   if (c > 0) {
   *     c = c - p
   *   }
   * But we remove the branch to protect against side-channel resistance.
   */
  m1 = (mont_cmp (p, m) > 0);
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_sub (c, p, m->p, m->n);
  for (i=0; i < m->n; i++) {
    c[i] = (p[i] & m2) | (c[i] & m1);
  }

}

#if 0
void montred_sl_ct_n (uint64_t *c, uint64_t *a, montmul_t *m, void *tmp) { 
  int i, size;
  uint64_t q, m1, m2, *mem;
  mem = (uint_64_t *) tmp;
  
  q = m->mu * a[0];                                   /* q has one limb  */
  size = 2*m->n;
  mul_1_n_add_m_div1_s (mem, a, q, m->p, m->n, size); 
  // (a + q*p) / 2^64
  size--;
  
  for (i=1; i < m->n-1; i++) {    
    q = mu * (*mem);                               /* q has one limb  */
    mul_1_n_add_m_div1 (mem, q, m->p, m->n, size);
    size--;
  }

  /* Final loop which stores the result in c */
  q = mu * (*mem); /* q has one limb  */
  mul_1_n_add_m_div1_s (c, mem, q, m->p, m->n, size); 
}
#endif

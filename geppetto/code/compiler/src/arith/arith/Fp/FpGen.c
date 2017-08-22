/* Fp/Fp.c
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
 
#include "FpGen.h"
#include "../uint.h"
#include "../mont_arith.h"
#include <math.h>

/* Allocate mod_size 64-bit limbs of space in dst 
 * Return > 0 on success and
 * return < 0 on error. 
 */
static int FpGen_alloc (FpGen_t dst, int n) {
  dst->limbs = NULL;
  dst->limbs = (uint64_t *) malloc (n * sizeof (uint64_t));
  memset(dst->limbs, 0, n * sizeof (uint64_t));
  if (dst->limbs == NULL) return ERR_OUT_OF_MEMORY;
  return ERR_SUCCESS;
}

int FpGen_init (FpGen_config_t* cfg, FpGen_t a) {
   return FpGen_alloc (a, cfg->m->n);
}

int FplGen_init (FpGen_config_t* cfg, FpGen_t a) {
   return FpGen_alloc (a, 2*cfg->m->n+1);
}

void FpGen_free (FpGen_config_t* cfg, FpGen_t src) {
  free ((src)->limbs);
  src = NULL;
}

/* Compute a^-1 mod 2^64  */
static uint64_t modinv64 (uint64_t a) {
  uint64_t c= (uint64_t) 1, b;
  int i;

  b=0;
  for (i=0; i < 64; i++) {
    uint64_t m = ((uint64_t) 1 << i);
    if ((c&m) != (uint64_t) 0) {
      c = c - a;
      b |= m;
    }
    a <<= (uint64_t) 1;
  }
  return b;
}


void FpGen_free_config (FpGen_config_t* cfg) {
  free (cfg->m->p);
  free (cfg->m);
  free (cfg->mem);
  free (cfg->p);
}

void FpGen_modinv (FpGen_config_t* cfg, FpGen_t dst, FpGen_t src) {
  binary_gcd (dst->limbs, src->limbs, cfg->m->p, (uint64_t*) cfg->mem, cfg->m->n);
  cfg->montmul (dst->limbs, dst->limbs, cfg->m->R3, cfg->m, cfg->mem);
}

// FpGen_config->mem need +1 mem
void FpGen_mul (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b) {
  /* Redirect to the current version of Montgomery multiplication. */
  cfg->montmul (c->limbs, a->limbs, b->limbs, cfg->m, cfg->mem);
}

#if 0
void FpGen_sqr (FpGen_config_t* cfg, FpGen_t c, FpGen_t a) {
  mont_sqr (c->limbs, a->limbs, cfg->m, cfg->mem);
}
#endif

void FpGen_add (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b) {
  mod_add (c->limbs, a->limbs, b->limbs, cfg->m, cfg->mem); 
}

/* Addition without reduction. No check for overflow! */
void FpGen_add_no_red (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b) {
  uint_add (c->limbs, a->limbs, b->limbs, cfg->m->n);
}

/* Addition without reduction on long elements. */
void FplGen_add_no_red (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b) {
   uint_add (c->limbs, a->limbs, b->limbs, 2*cfg->m->n);
}

void FpGen_sub (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b) {
  mod_sub (c->limbs, a->limbs, b->limbs, cfg->m, cfg->mem); 
}

/* Addition without reduction. No check for overflow! */
void FpGen_sub_no_red (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b) {
  uint_sub (c->limbs, a->limbs, b->limbs, cfg->m->n);
}

/* TODO FIX We assume here that we can add without overflow. */
void FpGen_div2 (FpGen_config_t* cfg, FpGen_t c, FpGen_t a) {
  uint64_t *t, m1, m2;
  int lsb, i;
  t = (uint64_t *) cfg->mem;
  uint_add (t, a->limbs, cfg->m->p, cfg->m->n);
  lsb = (a->limbs[0]&1);
  m1 = (uint64_t) 0 - lsb;
  m2 = ~m1;
  for (i=0; i < cfg->m->n; i++) {
    c->limbs[i] = ((t[i] & m1) | (a->limbs[i] & m2));
  }
  uint_div2 (c->limbs, c->limbs, cfg->m->n);
}

/* Addition without reduction. No check for overflow! */
void FplGen_sub_no_red (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b) {
  uint_sub (c->limbs, a->limbs, b->limbs, 2*cfg->m->n);
}

/* Subtraction of Fpl elements, option 1 (used with h=1 or h=2) in Aranha et al. */
void FplGen_sub_o1_ct (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b, int h, void *mem) {
  uint64_t *tmp1, *tmp2;
  int i;

  tmp1 = (uint64_t *) mem;
  tmp2 = (uint64_t *) mem + 2*cfg->m->n;

  for (i=0; i < 2*cfg->m->n; i++) tmp1[i] = a->limbs[i];

  // TODO fix this division
  uint_div2 (tmp2, cfg->p2N, 2*cfg->m->n);
  for (i=1; i < h; i++) uint_div2 (tmp2, tmp2, 2*cfg->m->n);
  
  // TODO fix this addition with fewer limbs
  uint_add (tmp1, tmp1, tmp2, 2*cfg->m->n);
  uint_sub(c->limbs, tmp1, b->limbs, 2*cfg->m->n);
}


/* Subtraction of Fpl elements, constant-time version of option 2 in Aranha et al. */
void FplGen_sub_o2_ct (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b, void *mem) {
  uint64_t borrow, m1, m2, *tmp;
  int i;

  tmp = (uint64_t *) mem;
  borrow = uint_sub(c->limbs, a->limbs, b->limbs, 2*cfg->m->n);
  /* Compute the conditional-addition:
   * if (borrow) {
   *   c = c + p;
   * }
   * using straight-line code as in mod_sub.
   */  
  m1 = borrow;
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_add (tmp, c->limbs + cfg->m->n, cfg->m->p, cfg->m->n);
  for (i=cfg->m->n; i < 2*cfg->m->n; i++) {
    c->limbs[i] = (tmp[i - cfg->m->n] & m1) | (c->limbs[i] & m2);  
  }
}

/* Addition of Fpl elements, constant-time version of "option 2" addition in Aranha et al. */
void FplGen_add_o2_ct (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b, void *mem) {
  uint64_t carry, m1, m2, *tmp;
  int i;
 
  tmp = (uint64_t *) mem;
  carry = uint_add (c->limbs, a->limbs, b->limbs, 2*cfg->m->n);

  m1 = carry;
  m1 |= uint_cmpge_ct (c->limbs, cfg->p2N, 2*cfg->m->n);  
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_sub (tmp, c->limbs + cfg->m->n, cfg->m->p, cfg->m->n);
  for (i=cfg->m->n; i < 2*cfg->m->n; i++) {
    c->limbs[i] = (tmp[i - cfg->m->n] & m1) | (c->limbs[i] & m2);  
  }
}

int FpGen_cmp (FpGen_config_t* cfg, FpGen_t a, FpGen_t b) {
  return uint_cmp (a->limbs, b->limbs, cfg->m->n);
}

/* Just the multiplication part, no reduction.
 * Make sure C holds enough space for 2*FpGen_config->n limbs. 
 */
void FpGen_mul_no_red (FpGen_config_t* cfg, FplGen_t c, FpGen_t a, FpGen_t b) {
  uint_mul (c->limbs, a->limbs, b->limbs, (uint64_t *) cfg->mem, cfg->m->n); 
}

void FplGen_red (FpGen_config_t* cfg, FpGen_t c, FplGen_t a) {
  montred_ct (c->limbs, a->limbs, cfg->m, cfg->mem);
}

void FpGen_mul3 (FpGen_config_t* cfg, FpGen_t c, FpGen_t a) {
  uint64_t *t;
  t = (uint64_t *) cfg->mem;  
  mod_add (t, a->limbs, a->limbs, cfg->m, (uint64_t *) cfg->mem+cfg->m->n); 
  mod_add (c->limbs, t, a->limbs, cfg->m, (uint64_t *) cfg->mem+cfg->m->n); 
}

/* 1 <= n < 64 */
/* Does not work: ToDo Fix later */
/*void FpGen_mul_2exp (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, int n) {
  int i;
  uint64_t *m = (uint64_t *) cfg->mem;
  m[cfg->m->n] = uint_mul_2exp (m, a->limbs, n, cfg->m->n);
  for (i=cfg->m->n+1; i < 2*cfg->m->n; i++) {
    m[i] = (uint64_t) 0;
  }
  montred_ct (c->limbs, m, cfg->m, (uint64_t *) cfg->mem + 2*cfg->m->n);
}*/

int FpGen_initialize_config (FpGen_config_t* cfg, uint64_t *m, int n) {
  int i, s1, s2;
  uint64_t *tmp;

  if (n <= 0) return ERR_INVALID_LIMB_SIZE;
  if ((m[0]&1) == 0) return ERR_EVEN_MODULUS;
  
  /* Allocate max (3*n,2(n-4+3floor( log_2(n-3) ))) limbs of additional memory.
   * Where the latter is upper-bound for Karatsuba 
   * NOTE: currently the inversion requires a lot of memory: 
   * 9*(n+1) limbs, this is what we currently allocate. 
   * TODO: Fix this later.
   */
  s1 = (int) (2*(n-4+3*((uint64_t)(log((double)n-3.0)/log(2.0)))));
  s2 = 3*n;
  if (s1 < s2) s1 = s2;
  s1 = 9*(n+1);

  cfg->mem = NULL;
  cfg->mem = (void *) malloc (s1 * sizeof (uint64_t));
  if (cfg->mem == NULL) return ERR_OUT_OF_MEMORY;
  
  cfg->m = NULL;
  cfg->m = (montmul_t *) malloc (sizeof (montmul_t));
  if (cfg->m == NULL) return ERR_OUT_OF_MEMORY;

  cfg->m->n = n;
  
  cfg->m->p = NULL;
  cfg->m->p = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (cfg->m->p == NULL) return ERR_OUT_OF_MEMORY;
  for (i=0; i < n; i++) cfg->m->p[i] = m[i];
  cfg->m->mu =0 -modinv64 (m[0]);


  /* Compute R^2 mod p by repeatedly doubling the number. */
  cfg->m->R2 = NULL;
  cfg->m->R2 = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (cfg->m->R2 == NULL) return ERR_OUT_OF_MEMORY;

  tmp = (uint64_t *) cfg->mem;

#if 0
  for (i=0; i < n-1; i++) {
    tmp[i] = (uint64_t) 0;
  }
  tmp[n-1] = (uint64_t) 1;
  tmp[n] = (uint64_t) 0;

  i = 0;
  while (tmp[n] == 0 && uint_cmp (tmp, cfg->m->p, n) <= 0) {
    uint_add (tmp, tmp, tmp, n+1);
    i++;
  }
  uint_sub_mn (tmp, tmp, cfg->m->p, n+1, n);
  for (j=0; j < n; j++) cfg->m->R2[j] = tmp[j];
  for (i=(n-1)*64+i; i < 2*n*64; i++) {  
    mod_add (cfg->m->R2, cfg->m->R2, cfg->m->R2, cfg->m, cfg->mem);
  }
#else  
  // R^2
  for (i=0; i < 2*n; i++) tmp[i] = (uint64_t) 0;
  tmp[2*n] = 1;
  uint_sub_mn (tmp, tmp, cfg->m->p, 2*n+1, n);
  divrem (cfg->m->R2, tmp, cfg->m->p, n);
#endif  


  cfg->m->R3 = NULL;
  cfg->m->R3 = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (cfg->m->R3 == NULL) return ERR_OUT_OF_MEMORY;
  montmul_ct (cfg->m->R3, cfg->m->R2, cfg->m->R2, cfg->m, cfg->mem);

  /* Regular Montgomery arithmetic. */
  cfg->p = NULL;
  cfg->montmul = montmul_ct;

  cfg->p2N = NULL;
  cfg->p2N = (uint64_t *) malloc (2 * n * sizeof (uint64_t));
  if (cfg->p2N == NULL) return ERR_OUT_OF_MEMORY;
  for (i=0; i < n; i++) cfg->p2N[i] = (uint64_t) 0;
  for (i=n; i < 2*n; i++) cfg->p2N[i] = m[i-n];

  return ERR_SUCCESS;
}

void FpGen_get (FpGen_config_t* cfg, uint64_t *t, FpGen_t a) {
  mont2normal (t, a->limbs, cfg->m, (void*)((uint64_t *)cfg->mem+cfg->m->n));
}

void FpGen_print (FpGen_config_t* cfg, FpGen_t a) {
  uint64_t *t = (uint64_t *) cfg->mem;
  FpGen_get (cfg, t, a);
  Print (t, cfg->m->n);
}

void FpGen_set_ui (FpGen_config_t* cfg, FpGen_t c, uint64_t x) {
  int i;
  uint64_t *mem = (uint64_t *) cfg->mem;
  mem[0] = x;
  for (i=1; i < cfg->m->n; i++) mem[i] = (uint64_t) 0;
  normal2mont (c->limbs, mem, cfg->m, (uint64_t *) cfg->mem + cfg->m->n);  
}

void FpGen_set (FpGen_config_t* cfg, FpGen_t c, uint64_t *a) {
  normal2mont (c->limbs, a, cfg->m, cfg->mem);
}

void FpGen_copy (FpGen_config_t* cfg, FpGen_t c, FpGen_t a) {
  int i;
  for (i=0; i < cfg->m->n; i++) c->limbs[i] = a->limbs[i];
}

void FplGen_copy (FpGen_config_t* cfg, FplGen_t c, FplGen_t a) {
  int i;
  for (i=0; i < 2*cfg->m->n; i++) c->limbs[i] = a->limbs[i];
}

void FpGen_neg (FpGen_config_t* cfg, FpGen_t c, FpGen_t a) {
	mod_sub (c->limbs, cfg->m->p, a->limbs, cfg->m, cfg->mem);
}

void FpGen_rand (FpGen_config_t* cfg, FpGen_t c) {
  int i;
  for (i=0; i < cfg->m->n; i++) {
    c->limbs[i] =  random_uint64_t();
  }
  if ((cfg->m->p[i-1]-1) != 0) c->limbs[i-1] %= (cfg->m->p[i-1]-1);
  else c->limbs[i-1] = 0;
}

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
 
#include "Fp.h"
#include "../uint.h"
#include "../mont_arith.h"
#include <math.h>

#define BIT(a,i) (((a)[i/64] >> (uint64_t) (i%64))&1)

Fp_config_t Fp_config;

/* Allocate mod_size 64-bit limbs of space in dst 
 * Return > 0 on success and
 * return < 0 on error. 
 */
static int Fp_alloc (Fp_t dst, int n) {
  dst->limbs = NULL;
  dst->limbs = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (dst->limbs == NULL) return ERR_OUT_OF_MEMORY;
  return ERR_SUCCESS;
}

int Fp_init (Fp_t a) {
   return Fp_alloc (a, Fp_config.m->n);
}

int Fpl_init (Fp_t a) {
   return Fp_alloc (a, 2*Fp_config.m->n+1);
}

void Fp_free (Fp_t src) {
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


void Fp_free_config (void) {
  Fp_free(Fp_config.tmp);
  free (Fp_config.m->p);
  free(Fp_config.m->R2);
  free(Fp_config.m->R3);
  free (Fp_config.m);
  free (Fp_config.mem);
  free (Fp_config.p);
  free (Fp_config.sqrt_exp_p3mod4);  
  free (Fp_config.p2N);
}

void Fp_modinv (Fp_t dst, Fp_t src) {
  binary_gcd (dst->limbs, src->limbs, Fp_config.m->p, (uint64_t*) Fp_config.mem, Fp_config.m->n);
  Fp_config.montmul (dst->limbs, dst->limbs, Fp_config.m->R3, Fp_config.m, Fp_config.mem);
}

// Fp_config->mem need +1 mem
void Fp_mul (Fp_t c, Fp_t a, Fp_t b) {
  /* Redirect to the current version of Montgomery multiplication. */
  Fp_config.montmul (c->limbs, a->limbs, b->limbs, Fp_config.m, Fp_config.mem);
}

#if 0
void Fp_sqr (Fp_t c, Fp_t a) {
  mont_sqr (c->limbs, a->limbs, Fp_config.m, Fp_config.mem);
}
#endif

void Fp_add (Fp_t c, Fp_t a, Fp_t b) {
  mod_add (c->limbs, a->limbs, b->limbs, Fp_config.m, Fp_config.mem); 
}

/* Addition without reduction. No check for overflow! */
void Fp_add_no_red (Fp_t c, Fp_t a, Fp_t b) {
  uint_add (c->limbs, a->limbs, b->limbs, Fp_config.m->n);
}

/* Addition without reduction on long elements. */
void Fpl_add_no_red (Fpl_t c, Fpl_t a, Fpl_t b) {
   uint_add (c->limbs, a->limbs, b->limbs, 2*Fp_config.m->n);
}

void Fp_sub (Fp_t c, Fp_t a, Fp_t b) {
  mod_sub (c->limbs, a->limbs, b->limbs, Fp_config.m, Fp_config.mem); 
}

/* Addition without reduction. No check for overflow! */
void Fp_sub_no_red (Fp_t c, Fp_t a, Fp_t b) {
  uint_sub (c->limbs, a->limbs, b->limbs, Fp_config.m->n);
}

/* TODO FIX We assume here that we can add without overflow. */
void Fp_div2 (Fp_t c, Fp_t a) {
  uint64_t *t, m1, m2;
  int lsb, i;
  t = (uint64_t *) Fp_config.mem;
  uint_add (t, a->limbs, Fp_config.m->p, Fp_config.m->n);
  lsb = (a->limbs[0]&1);
  m1 = (uint64_t) 0 - lsb;
  m2 = ~m1;
  for (i=0; i < Fp_config.m->n; i++) {
    c->limbs[i] = ((t[i] & m1) | (a->limbs[i] & m2));
  }
  uint_div2 (c->limbs, c->limbs, Fp_config.m->n);
}

/* Addition without reduction. No check for overflow! */
void Fpl_sub_no_red (Fpl_t c, Fpl_t a, Fpl_t b) {
  uint_sub (c->limbs, a->limbs, b->limbs, 2*Fp_config.m->n);
}

/* Subtraction of Fpl elements, option 1 (used with h=1 or h=2) in Aranha et al. */
void Fpl_sub_o1_ct (Fpl_t c, Fpl_t a, Fpl_t b, int h, void *mem) {
  uint64_t *tmp1, *tmp2;
  int i;

  tmp1 = (uint64_t *) mem;
  tmp2 = (uint64_t *) mem + 2*Fp_config.m->n;

  for (i=0; i < 2*Fp_config.m->n; i++) tmp1[i] = a->limbs[i];

  // TODO fix this division
  uint_div2 (tmp2, Fp_config.p2N, 2*Fp_config.m->n);
  for (i=1; i < h; i++) uint_div2 (tmp2, tmp2, 2*Fp_config.m->n);
  
  // TODO fix this addition with fewer limbs
  uint_add (tmp1, tmp1, tmp2, 2*Fp_config.m->n);
  uint_sub(c->limbs, tmp1, b->limbs, 2*Fp_config.m->n);
}


/* Subtraction of Fpl elements, constant-time version of option 2 in Aranha et al. */
void Fpl_sub_o2_ct (Fpl_t c, Fpl_t a, Fpl_t b, void *mem) {
  uint64_t borrow, m1, m2, *tmp;
  int i;

  tmp = (uint64_t *) mem;
  borrow = uint_sub(c->limbs, a->limbs, b->limbs, 2*Fp_config.m->n);
  /* Compute the conditional-addition:
   * if (borrow) {
   *   c = c + p;
   * }
   * using straight-line code as in mod_sub.
   */  
  m1 = borrow;
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_add (tmp, c->limbs + Fp_config.m->n, Fp_config.m->p, Fp_config.m->n);
  for (i=Fp_config.m->n; i < 2*Fp_config.m->n; i++) {
    c->limbs[i] = (tmp[i - Fp_config.m->n] & m1) | (c->limbs[i] & m2);  
  }
}

/* Addition of Fpl elements, constant-time version of "option 2" addition in Aranha et al. */
void Fpl_add_o2_ct (Fpl_t c, Fpl_t a, Fpl_t b, void *mem) {
  uint64_t carry, m1, m2, *tmp;
  int i;
 
  tmp = (uint64_t *) mem;
  carry = uint_add (c->limbs, a->limbs, b->limbs, 2*Fp_config.m->n);

  m1 = carry;
  m1 |= uint_cmpge_ct (c->limbs, Fp_config.p2N, 2*Fp_config.m->n);  
  m1 = (uint64_t) 0 - m1;
  m2 = ~m1;
  uint_sub (tmp, c->limbs + Fp_config.m->n, Fp_config.m->p, Fp_config.m->n);
  for (i=Fp_config.m->n; i < 2*Fp_config.m->n; i++) {
    c->limbs[i] = (tmp[i - Fp_config.m->n] & m1) | (c->limbs[i] & m2);  
  }
}

int Fp_cmp (Fp_t a, Fp_t b) {
  return uint_cmp (a->limbs, b->limbs, Fp_config.m->n);
}

/* Just the multiplication part, no reduction.
 * Make sure C holds enough space for 2*Fp_config->n limbs. 
 */
void Fp_mul_no_red (Fpl_t c, Fp_t a, Fp_t b) {
  uint_mul (c->limbs, a->limbs, b->limbs, (uint64_t *) Fp_config.mem, Fp_config.m->n); 
}

void Fpl_red (Fp_t c, Fpl_t a) {
  montred_ct (c->limbs, a->limbs, Fp_config.m, Fp_config.mem);
}

void Fp_mul3 (Fp_t c, Fp_t a) {
  uint64_t *t;
  t = (uint64_t *) Fp_config.mem;  
  mod_add (t, a->limbs, a->limbs, Fp_config.m, (uint64_t *) Fp_config.mem+Fp_config.m->n); 
  mod_add (c->limbs, t, a->limbs, Fp_config.m, (uint64_t *) Fp_config.mem+Fp_config.m->n); 
}

/* 1 <= n < 64 */
/* Does not work: ToDo Fix later */
/*void Fp_mul_2exp (Fp_t c, Fp_t a, int n) {
  int i;
  uint64_t *m = (uint64_t *) Fp_config.mem;
  m[Fp_config.m->n] = uint_mul_2exp (m, a->limbs, n, Fp_config.m->n);
  for (i=Fp_config.m->n+1; i < 2*Fp_config.m->n; i++) {
    m[i] = (uint64_t) 0;
  }
  montred_ct (c->limbs, m, Fp_config.m, (uint64_t *) Fp_config.mem + 2*Fp_config.m->n);
}*/

int Fp_initialize_config (uint64_t *m, int n) {
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

  Fp_config.mem = NULL;
  Fp_config.mem = (void *) malloc (s1 * sizeof (uint64_t));
  if (Fp_config.mem == NULL) return ERR_OUT_OF_MEMORY;
  
  Fp_config.m = NULL;
  Fp_config.m = (montmul_t *) malloc (sizeof (montmul_t));
  if (Fp_config.m == NULL) return ERR_OUT_OF_MEMORY;

  Fp_config.m->n = n;
  
  Fp_config.m->p = NULL;
  Fp_config.m->p = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (Fp_config.m->p == NULL) return ERR_OUT_OF_MEMORY;
  for (i=0; i < n; i++) Fp_config.m->p[i] = m[i];
  Fp_config.m->mu =0 -modinv64 (m[0]);

  /* Compute R^2 mod p by repeatedly doubling the number. */
  Fp_config.m->R2 = NULL;
  Fp_config.m->R2 = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (Fp_config.m->R2 == NULL) return ERR_OUT_OF_MEMORY;

  tmp = (uint64_t *) Fp_config.mem;

#if 0
  for (i=0; i < n-1; i++) {
    tmp[i] = (uint64_t) 0;
  }
  tmp[n-1] = (uint64_t) 1;
  tmp[n] = (uint64_t) 0;

  i = 0;
  while (tmp[n] == 0 && uint_cmp (tmp, Fp_config.m->p, n) <= 0) {
    uint_add (tmp, tmp, tmp, n+1);
    i++;
  }
  uint_sub_mn (tmp, tmp, Fp_config.m->p, n+1, n);
  for (j=0; j < n; j++) Fp_config.m->R2[j] = tmp[j];
  for (i=(n-1)*64+i; i < 2*n*64; i++) {  
    mod_add (Fp_config.m->R2, Fp_config.m->R2, Fp_config.m->R2, Fp_config.m, Fp_config.mem);
  }
#else  
  // R^2
  for (i=0; i < 2*n; i++) tmp[i] = (uint64_t) 0;
  tmp[2*n] = 1;
  uint_sub_mn (tmp, tmp, Fp_config.m->p, 2*n+1, n);
  divrem (Fp_config.m->R2, tmp, Fp_config.m->p, n);
#endif  


  Fp_config.m->R3 = NULL;
  Fp_config.m->R3 = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (Fp_config.m->R3 == NULL) return ERR_OUT_OF_MEMORY;
  montmul_ct (Fp_config.m->R3, Fp_config.m->R2, Fp_config.m->R2, Fp_config.m, Fp_config.mem);

  /* Regular Montgomery arithmetic. */
  Fp_config.p = NULL;
  Fp_config.montmul = montmul_ct;

  Fp_config.p2N = NULL;
  Fp_config.p2N = (uint64_t *) malloc (2 * n * sizeof (uint64_t));
  if (Fp_config.p2N == NULL) return ERR_OUT_OF_MEMORY;
  for (i=0; i < n; i++) Fp_config.p2N[i] = (uint64_t) 0;
  for (i=n; i < 2*n; i++) Fp_config.p2N[i] = m[i-n];

  tmp[0] = (uint64_t) 1;

  Fp_config.sqrt_exp_p3mod4 = NULL;
  Fp_config.sqrt_exp_p3mod4 = (uint64_t *)malloc((n+1) * sizeof (uint64_t));
  Fp_config.sqrt_exp_p3mod4[n] = 0;
  if (Fp_config.sqrt_exp_p3mod4 == NULL) return ERR_OUT_OF_MEMORY;
  uint_add_mn(Fp_config.sqrt_exp_p3mod4, Fp_config.m->p, tmp, n, 1);
  uint_div2(Fp_config.sqrt_exp_p3mod4, Fp_config.sqrt_exp_p3mod4, n + 1);
  uint_div2(Fp_config.sqrt_exp_p3mod4, Fp_config.sqrt_exp_p3mod4, n + 1);
  Fp_config.sqrt_exp_len_p3mod4 = uint_nb_bits(Fp_config.sqrt_exp_p3mod4, n + 1);

  Fp_config.Fp_sqrt = p3mod4_Fp_sqrt;

  if (Fp_init(Fp_config.tmp) < 0) printf("Fp_config memory error.\n");

  return ERR_SUCCESS;
}

void Fp_get (uint64_t *t, Fp_t a) {
  mont2normal (t, a->limbs, Fp_config.m, (void*)((uint64_t *)Fp_config.mem+Fp_config.m->n));
}

void Fp_print (Fp_t a) {
  uint64_t *t = (uint64_t *) Fp_config.mem;
  Fp_get (t, a);
  Print (t, Fp_config.m->n);
}

void Fp_set_ui (Fp_t c, uint64_t x) {
  int i;
  uint64_t *mem = (uint64_t *) Fp_config.mem;
  mem[0] = x;
  for (i=1; i < Fp_config.m->n; i++) mem[i] = (uint64_t) 0;
  normal2mont (c->limbs, mem, Fp_config.m, (uint64_t *) Fp_config.mem + Fp_config.m->n);  
}

void Fp_set (Fp_t c, uint64_t *a) {
  normal2mont (c->limbs, a, Fp_config.m, Fp_config.mem);
}

void Fp_copy (Fp_t c, Fp_t a) {
  int i;
  for (i=0; i < Fp_config.m->n; i++) c->limbs[i] = a->limbs[i];
}

void Fpl_copy (Fpl_t c, Fpl_t a) {
  int i;
  for (i=0; i < 2*Fp_config.m->n; i++) c->limbs[i] = a->limbs[i];
}

void Fp_neg (Fp_t c, Fp_t a) {
	mod_sub (c->limbs, Fp_config.m->p, a->limbs, Fp_config.m, Fp_config.mem);
}

void Fp_rand (Fp_t c) {
  int i;
  for (i=0; i < Fp_config.m->n; i++) {
    c->limbs[i] =  random_uint64_t();
  }
  if ((Fp_config.m->p[i-1]-1) != 0) c->limbs[i-1] %= (Fp_config.m->p[i-1]-1);
  else c->limbs[i-1] = 0;
}

void p3mod4_Fp_sqrt(Fp_t c, Fp_t a) {
	int i;
	Fp_t tmp; 
	tmp->limbs  = Fp_config.tmp->limbs;
	Fp_copy(c, a);
	Fp_copy(tmp, a);
	for (i = Fp_config.sqrt_exp_len_p3mod4 - 2; i >= 0; i--) {
		Fp_mul (c, c, c);
		if (BIT(Fp_config.sqrt_exp_p3mod4,i) == 1) Fp_mul(c, c, tmp);
	}
}

void Fp_sqrt(Fp_t c, Fp_t a) {
	Fp_config.Fp_sqrt(c, a);
}
/* uint.h
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

#ifndef __UINT_H

#define __UINT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#pragma intrinsic(_umul128)

/* Subtract n-limb integers
 * c = a - b
 * No size check on n
 * Return borrow-out
 */
int uint_sub (uint64_t *c, uint64_t *a, uint64_t *b, int n);

/* Addition of n-limb integers
 * c = a - b
 * No size check on n
 * Return carry-out
 */
int uint_add (uint64_t *c, uint64_t *a, uint64_t *b, int n);

/* Multiplication of the two n-limb integers a and b 
 * to the result dst which should have enough space for 2n limbs.
 * In- and output may not overlap.
 */
void uint_mul (uint64_t *dst, uint64_t *a, uint64_t *b, uint64_t *mem, int n);

/* Compare the two n-limb integers a and b
 * Return > 0, if a > b
 *          0, if a==b
 *        < 0, if a < b
 */
int uint_cmp (uint64_t *a, uint64_t *b, int n);

/* Returns 1 if a > b
* Returns 0 if a <= b
* Where a and b are both 'len'-limb 64-bit integers.
* This function runs in constant time.
*/
int uint_cmpgt_ct (uint64_t *a, uint64_t *b, int n);

/* Returns 1 if a >= b
 * Returns 0 if a < b
 * Where a and b are both 'len'-limb 64-bit integers.
 * This function runs in constant time.
 */
int uint_cmpge_ct (uint64_t *a, uint64_t *b, int n);

/* Compute c = floor(a/2) = (a>>1) */
void uint_div2 (uint64_t *c, uint64_t *a, int n);

int uint_add_mn (uint64_t *c, uint64_t *a, uint64_t *b, int m, int n);

int uint_sub_mn (uint64_t *c, uint64_t *a, uint64_t *b, int m, int n);

/* Shift a by s (1 <= s < 64) positions and put the result in c,
 * return the bits shifted out.
 */
uint64_t uint_mul_2exp (uint64_t *c, uint64_t *a, uint64_t s, int n);

/************************* Internal functions ************************/

static void Print (uint64_t *a, int n) {
  int i;

  for (i=n-1; i>= 0; i--) {
	  printf ("%08X%08X", (uint32_t) (a[i]>>(uint64_t)32), (uint32_t) a[i]);
	  //printf ("%u %u ", (uint32_t) (a[i]>>(uint64_t)32), (uint32_t) a[i]);
  }
  //printf ("\n");
 }


/* 'Good-old' schoolbook (aka textbook) multiplication */
void mul_schoolbook (uint64_t *c, uint64_t *a, uint64_t *b, unsigned int n);

void mul_schoolbook_mn (uint64_t *c, uint64_t *a, uint64_t *b, unsigned int m, unsigned int n);

/* Macro which determines when to switch from schoolbook to Karatsuba.
 * Initial test show this is reasonable value for modern 64-bit machines.
 */
extern uint32_t SCHOOLBOOK_THRESHOLD; 

/* Karatsuba algorithm which computes on n-limb integers a, b
 * dst need to have space for 2n-limbs and should not overlap with a or b
 */
void mul_karatsuba (uint64_t *dst, uint64_t *a, uint64_t *b, int n);

/* Platform specific macro which implement some basic functionality */

/* Currently these macros work on 64-bit platforms
 * when compiling with Visual Studio.
 */

/* 64 x 64 --> 128-bit multiplication
 * (c1,c0) = a * b
 */
#define mul(c0,c1,a,b) c0 = _umul128 (a,b,&c1)

/* Multiply-and-accumulate
 * (c1,c0) = a*b+c0 
 */
#define muladd(c0,c1,a,b) \
do { \
  uint64_t _c = c0; \
  mul (c0,c1,a,b); \
  c0 = c0 + _c; \
  c1 = c1 + (c0 < _c); \
} while (0)

/* Multiply-and-accumulate-accumulate
 * (c1,c0) = a*b+c0+c1
 */
#define muladdadd(c0,c1,a,b) \
do { \
  uint64_t _C0 = c0, _C1 = c1; \
  mul (c0,c1,a,b); \
  c0 = c0 + _C0; \
  c1 = c1 + (c0 < _C0); \
  c0 = c0 + _C1; \
  c1 = c1 + (c0 < _C1); \
} while (0)

/* mul_1_n: Compute c = a * b where a is a single limb, 
 * b is n-limbs and c has room for n+1 limbs.
 * Note that c cannot overlap with either a or b.
 * We assume that n > 0.
 */
static void mul_1_n (uint64_t *c, const uint64_t a, const uint64_t * const b, int n) {
  int i;
  mul (c[0],c[1],a,b[0]);
  for (i=1; i < n; i++) {
    muladd (c[i], c[i+1], a, b[i]);
  }
}

/* mul_1_n_add_n: Compute
 * c = c + a * b
 * where the input a has a single limb and b,c have n limbs 
 * return the most-significant (n+1)th limb
 */
static uint64_t mul_1_n_add_n (uint64_t *c, const uint64_t a, const uint64_t * const b, int n) {
  int i;
  uint64_t t;

  muladd (c[0],t,a,b[0]);
  for (i=1; i < n; i++) {
    muladdadd (c[i], t, a, b[i]);
  }
  return t;
}

void uint_karatsuba (uint64_t *a, uint64_t *b, uint64_t *c, uint64_t *w, int la, int lb);

/*
 * Count number of bits of uint64_t* and uint64_t
 */
unsigned int uint_nb_bits(uint64_t *c, unsigned int n);
unsigned int uint64_nb_bits(uint64_t c);

/*
 * Generate a random uint64_t element
 */
uint64_t random_uint64_t();

#if 0
/* mul_1_n_add_n: Compute
 * c = (c + a * b)/2^64
 * where the input a has a single limb and b,c has n limbs 
 * and the result c has n limbs.
 */
static void mul_1_n_add_n_div1 (uint64_t *c, const uint64_t a, const uint64_t * const b, int n) {
  int i;
  uint64_t t;

  muladd (c[0],t,a,b[0]);
  for (i=1; i < n; i++) {
    muladdadd (c[i-1], t, a, b[i]);
  }
  c[i-1] = t;
}
#endif

/* mul_1_n_add_m_div1: Compute
 * c = (c + a * b)/2^64
 * where the input a has a single limb, b has n limbs and c has m limbs (m>n)
 * storing the m-limb result in c.
 */
static void mul_1_n_add_m_div1 (uint64_t *c, const uint64_t a, const uint64_t * const b, int n, int m) {
  int i;
  uint64_t t;

  muladd (c[0],t,a,b[0]);
  for (i=1; i < n; i++) {
    muladdadd (c[i], t, a, b[i]);
    c[i-1] = c[i];
  }
  c[i-1] = c[i] + t;
  t = (c[i-1] < t);

  for (i++; i < m; i++) {
    c[i-1] = c[i] + t;
    t = (c[i-1] < t);
  }
  c[i-1] = t;
}

#endif /* __UINT_H */

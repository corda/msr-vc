/* uint.c
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

#include "uint.h"
#include "arith_error.h"
#include "Fp/Fp.h"

uint32_t SCHOOLBOOK_THRESHOLD = 30; /* Initial test show this is reasonable value for modern 64-bit machines. */

int uint_cmp (uint64_t *a, uint64_t *b, int n) {
  int i;
  for (i=n-1; i >= 0; i--) {
    if (a[i] > b[i]) return 1;
    if (a[i] < b[i]) return -1;
  }
  return 0;
}

/* Returns 1 if a >= b
 * Returns 0 if a < b
 * Where a and b are both 'len'-limb 64-bit integers.
 * This function runs in constant time.
 */
int uint_cmpge_ct (uint64_t *a, uint64_t *b, int len) {
  int i, m;

  /* We compute this using recursive definition:
   * F(a, b, n) := a_0 >= b_0                                        , if n == 0
   *               (a_i >= b_i) /\ (not(a_i == b_i) \/ F(a, b, n-1)  , else
   */
  m = (a[0] >= b[0]);
  for (i=1; i < len; i++) {
    m = ((a[i] >= b[i]) && (!(a[i] == b[i]) || m));
  }
  return m;
}

int uint_sub (uint64_t *c, uint64_t *a, uint64_t *b, int n) {
  int i, b1, b2;

  b1 = (b[0] > a[0]);
  c[0] = a[0] - b[0];
  
  for (i=1; i < n; i++) {
    b2 = (b[i] > a[i]);
    c[i] = a[i] - b[i];

    b2 += (b1 > c[i]);
    c[i] -= b1;
    b1 = b2;
  }
  return b1;
}

/* c = a + b
* No size check on n
* Return carry-out
*/
int uint_add (uint64_t *c, uint64_t *a, uint64_t *b, int n) {
  int i, c1;
  uint64_t C;
  
  C = a[0] + b[0];
  c1 = (C < a[0]);
  c[0] = C;
  
  for (i=1; i < n; i++) {
    C = a[i] + c1;
    c1 = (C < c1);
    C = C + b[i];
    c1 += (C < b[i]);
    c[i] = C;
  }
  return c1;
}

void mul_schoolbook (uint64_t *c, uint64_t *a, uint64_t *b, unsigned int n) {
  uint64_t A;
  unsigned int j;

  A = a[0];
  mul_1_n (c, A, b, n);

  for (j=1; j < n; j++) {
    A = a[j];
    c[j+n] = mul_1_n_add_n (&c[j], A, b, n);
  }
}

void mul_schoolbook_mn (uint64_t *c, uint64_t *a, uint64_t *b, unsigned int m, unsigned int n) {
  uint64_t A;
  unsigned int j;

  A = a[0];
  mul_1_n (c, A, b, n);

  for (j=1; j < m; j++) {
    A = a[j];
    c[j+n] = mul_1_n_add_n (&c[j], A, b, n);
  }
}

void uint_mul (uint64_t *dst, uint64_t *a, uint64_t *b, uint64_t *mem, int n) {
  if ((uint32_t) n <= SCHOOLBOOK_THRESHOLD) {
    mul_schoolbook (dst, a, b, n);
  } else {
    //memset (dst,0,2*n*sizeof(uint64_t)+1);
    int i;
    for (i=0; i < 2*n; i++) dst[i] = 0;
    uint_karatsuba (dst, a, b, mem, n, n);
  }
}

/* c = a + b where
 * a has m-limbs and 
 * b has n-limbs (m > n)
 */
int uint_add_mn (uint64_t *c, uint64_t *a, uint64_t *b, int m, int n) {
  int i, c1;
  uint64_t C;
  
  C = a[0] + b[0];
  c1 = (C < a[0]);
  c[0] = C;
  
  for (i=1; i < n; i++) {
    C = a[i] + c1;
    c1 = (C < c1);
    C = C + b[i];
    c[i] = C;
    c1 += (C < b[i]);
  }
  for (; i < m; i++) {
    C = a[i] + c1;
    c1 = (C < c1);
    c[i] = C;
  }
  return c1;
}

int uint_sub_mn (uint64_t *c, uint64_t *a, uint64_t *b, int m, int n) {
  int i, b1, b2;

  b1 = (b[0] > a[0]);
  c[0] = a[0] - b[0];
  
  for (i=1; i < n; i++) {
    b2 = (b[i] > a[i]);
    c[i] = a[i] - b[i];
    b2 += (b1 > c[i]);
    c[i] -= b1;
    b1 = b2;
  }
  for (; i < m; i++) {
    b2 = (b1 > a[i]);
    c[i] = a[i] - b1;
    b1 = b2;
  }
  return b1;
}

/* Karatsuba multiplication following the efficient storage ideas from
 * Roman E. Maeder 
 * Design and Implementation of Symbolic Computation Systems 
 * Lecture Notes in Computer Science Volume 722,  1993,   pp 59-65 
 * This requires a maximum of 2(n-4+3floor( log_2(n-3) )) storage
 * for n-limb input when recursively going to 4-limbs (which we don't do.
 */
void uint_karatsuba (uint64_t *c, uint64_t *a, uint64_t *b, uint64_t *w, int la, int lb) {
  int m, lt;

  if ((uint32_t) la <= SCHOOLBOOK_THRESHOLD) {
    mul_schoolbook (c, a, b, la);
    return;
  }
  m = (la+1)/2;
  memcpy (w, a, m*sizeof(uint64_t));
  
  //w[m] = 0;
  //w[la-m] += uint_add (w, w, &a[m], la-m);
  w[m] = uint_add_mn (w, w, &a[m], m, la-m);

  memcpy (&w[m+1], b, m*sizeof(uint64_t));

  //w[m+1+m] = 0;
  //w[m+1+(lb-m)] += uint_add (&w[m+1], &w[m+1], &b[m], lb-m);
  w[m+1+m] = uint_add_mn (&w[m+1], &w[m+1], &b[m], m, lb-m);

  uint_karatsuba (&c[m], w, &w[m+1], &w[2*(m+1)], m+1, m+1);
  lt = (la-m) + (lb-m) + 1;
  memset (w, 0, lt*sizeof(uint64_t));
  uint_karatsuba (w, &a[m], &b[m], &w[lt], la-m, lb-m);
  
  //c[2*m+(la-m)+(lb-m)] += uint_add (&c[2*m], &c[2*m], w, (la-m)+(lb-m));
  uint_add_mn (&c[2*m], &c[2*m], w, 2*m, (la-m)+(lb-m));
  
  //c[m+(la-m)+(lb-m)] -= uint_sub (&c[m], &c[m], w, (la-m)+(lb-m));
  uint_sub_mn (&c[m], &c[m], w, la+lb-m, (la-m)+(lb-m));

  lt = m+m+1;
  memset (w, 0, lt * sizeof(uint64_t));
  uint_karatsuba (w, a, b, &w[lt], m, m);
  
  //c[m+m] += uint_add (c, c, w, m+m);
  uint_add_mn (c, c, w, la+lb, m+m);

  //c[m+m+m] -= uint_sub (&c[m], &c[m], w, m+m);
  uint_sub_mn (&c[m], &c[m], w, la+lb-m, m+m);
}

/* Compute c = floor(a/2) = (a>>1) */
void uint_div2 (uint64_t *c, uint64_t *a, int n) {
  int i;  
  for (i=0; i < n-1; i++) {
    c[i] = (((uint64_t) a[i] >> (uint64_t) 1) | ((uint64_t)a[i+1] << (uint64_t) 63)); 
  }
  c[n-1] = ((uint64_t) a[n-1] >> (uint64_t) 1);
}

uint64_t uint_mul_2exp (uint64_t *c, uint64_t *a, uint64_t s, int n) {
  int i;
  uint64_t ret;
  ret = (a[n-1] >> (uint64_t) (64-s));
  for (i=n-1; i >= 1; i--) {
    c[i] = ((a[i] << s) | (a[i-1] >> (uint64_t) (64-s)));
  }
  c[0] = (a[0] << s);
  return ret;
}

/*
 * Return the bit-size of c
 */
unsigned int uint_nb_bits(uint64_t *c, unsigned int n)
{
	unsigned int nb_bits = n*64;
	int i;
	for (i=n-1; i>=0; i--)
	{
		if ( c[i] != (uint64_t) 0 )
			break;
		else
			nb_bits -= 64;
	}
	nb_bits -= (64-uint64_nb_bits(c[i]));
	return(nb_bits);
}

/*
 * Return the bit-size of a uint64_t
 */
unsigned int uint64_nb_bits(uint64_t c)
{
	unsigned int nb_bits=64, i;
	for (i=63; i>=0; i--)
	{
		if ( (c>>i)&1 )
			break;
		else
			nb_bits--;
	}
	return(nb_bits);
}

/*
 * Random uint64_t element
 */
uint64_t random_uint64_t()
{
	return ((uint64_t) rand () ^ ((uint64_t) rand () << (uint64_t) 15) ^ ((uint64_t) rand () << (uint64_t) 30) ^ ((uint64_t) rand () << (uint64_t) 45) ^ ((uint64_t) rand () << (uint64_t) 60));
}
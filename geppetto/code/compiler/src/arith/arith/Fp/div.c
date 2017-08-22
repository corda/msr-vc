/* Fp/div.c
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

#include "../arith.h"
#include "../uint.h"
#include <math.h>

/* From Zimmermann's book, Algorithm 1.6 with beta=2^64 
 * in the case m=n. */
/* This code was quickly made to compute modular reduction, it does 
 * *not* compute the quotient.
 * ToDo: Many things can be optimized.
 */
void divrem (uint64_t *r, uint64_t *a, uint64_t *b, int n) {
  int m, i, j, k;
  uint64_t *A, q, *T1, *T2, *B;

  m = n+1;
  A = (uint64_t *) malloc ((m+n)*sizeof (uint64_t));
  B = (uint64_t *) malloc (n*sizeof (uint64_t));
  T1 = (uint64_t *) malloc ((m+n)*sizeof (uint64_t));
  T2 = (uint64_t *) malloc ((m+n)*sizeof (uint64_t));
 
  for (i=0; i < m+n-1; i++) A[i] = a[i];
  A[i] = 0;
  for (i=0; i < n; i++) B[i] = b[i];

  //printf("divrem step 1\n");

  /* Normalize B */
  k=0;
  while (B[n-1] <= ((uint64_t) 1<< (uint64_t) 63)) {
    uint_mul_2exp (B, B, 1, n);
	k++;
  }
  uint_mul_2exp (A, A, k, m+n);

  //printf("divrem step 2\n");
  
  for (i=0; i < m; i++) T2[i] = 0;  
  for (i=m; i < m+n; i++) T2[i] = B[i-m];
  if (uint_cmp (A, T2, m+n) >= 0) {  
    uint_sub (A, A, T2, m+n);
  }

  //printf("divrem step 3\n");
    
  for (j=m-1; j >= 0; j--) {
    uint64_t b2[2], Q[2], a2[2], tmp[2];

	Q[0] = 0;
	Q[1] = 0;
    b2[0] = B[n-1];
	b2[1] = 0;
    a2[0] = A[n+j-1];
	a2[1] = A[n+j];

	//if (j%100==0) printf("divrem step 4, j=%d\n", j);

    while (uint_cmp (a2, b2, 2) >= 0) {
	  uint64_t out;
      i=0;
	  while (uint_cmp (b2, a2, 2) <= 0) {
        out = uint_mul_2exp (b2, b2, 1, 2);
	    i++;
	    if (out == 1) break;
      }
      i--;
	  uint_div2 (b2, b2, 2);
	  b2[1] |= ((uint64_t) out << (uint64_t) 63);
	  if (i < 64) {
        tmp[0] = ((uint64_t) 1 << (uint64_t) i);
	    tmp[1] = 0;
      } else {
	    tmp[0] = 0;
	    tmp[1] = ((uint64_t) 1 << (uint64_t) (i-64));
      }
      uint_add (Q, Q, tmp, 2);
	  uint_sub (a2, a2, b2, 2);

	  b2[0] = B[n-1];
	  b2[1] = 0;
    }

	//if (j%100==0) printf("divrem step 4, j=%d (while ok)\n", j);
	
    if (Q[1] != 0) {
	  q = (uint64_t) 0xFFFFFFFFFFFFFFFF;
    } else {
      q = Q[0];
    }

    //for (i=0; i < j; i++) T2[i] = 0;
    //for (i=j; i < j+n; i++) T2[i] = B[i-j];
	//mul_schoolbook_mn (T1, T2, &q, j+n, 1);
    
	mul_schoolbook_mn (T2, B, &q, n, 1);
    for (i=0; i < j; i++)     T1[i] = 0;
    for (i=j; i < j+n+1; i++) T1[i] = T2[i-j];
	for (i=0; i < j; i++) T2[i] = 0;
    for (i=j; i < j+n; i++) T2[i] = B[i-j];

	//if (j%100==0) printf("divrem step 4 mult_schoolbook, j=%d, m=%d, n=%d\n", j, m, n);

	for (i=j+n+1; i < m+n; i++) T1[i] = 0;
    while (uint_cmp (T1, A, m+n) > 0) {
	  uint_sub_mn (T1, T1, T2, m+n, j+n);
    }
    uint_sub (A, A, T1, m+n);
  }
  for (i=0; i < k; i++) uint_div2 (A, A, m+n);
  for (i=0; i < n; i++) r[i] = A[i];
  
  free (A);
  free (B);
  free (T1);
  free (T2);
}

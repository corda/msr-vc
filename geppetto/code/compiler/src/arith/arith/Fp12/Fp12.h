/* Fp12/Fp12.h
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


#ifndef __FP12_H
#define __FP12_H

#include "../arith_error.h"
#include "../Fp6/Fp6.h"


/* This contructions realizes F_{p^12} as a quadratic extension over F_{p^6}
 * by an irreducible polynomial X^2 - v, where v is in F_{p^6}.
 */

/* An element in F_{p^12} consists of two coefficients from F_{p^6}, 
 * i.e. a = a0 + a1*w in F_{p^12} for w in F_{p^6} with w^2 = v. 
 */
typedef struct {
  Fp6_t a0;
  Fp6_t a1;
} fp12_t;

/* Our basic type is a pointer to the struct. */
typedef fp12_t Fp12_t[1];


int Fp12_initialize_config (void);

void Fp12_free_config (void);


/************** Regular Fp12 functions ***************/

/* Initialize an element in Fp12 by initializing two Fp6 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp12_init(Fp12_t dst);

/* Free the memory allocated in src. */
void Fp12_free (Fp12_t src);

/* Get a random element in Fp12. */
void Fp12_rand (Fp12_t c);

/* Print the (regular non-Montgomery form) value of a */
void Fp12_print (Fp12_t a);

/* Set c to 0. */
void Fp12_set_zero (Fp12_t c);

/* Set c to 1. */
void Fp12_set_one (Fp12_t c);

/* Copy an Fp12 element a to an Fp12 element c. */
void Fp12_copy(Fp12_t c, Fp12_t a);

/* Compare two Fp12 elements a and b for equality. */
int Fp12_cmpeq(Fp12_t a, Fp12_t b);

/* Negate an Fp12 element. */
void Fp12_neg(Fp12_t c, Fp12_t a);

/* Add two Fp12 elements coefficient wise. */
void Fp12_add (Fp12_t c, Fp12_t a, Fp12_t b);

/* Subtract an Fp12 element b from an Fp12 element a coefficient wise. */
void Fp12_sub (Fp12_t c, Fp12_t a, Fp12_t b);

/* Multiply two Fp12 elements.*/
void Fp12_mul(Fp12_t c, Fp12_t a, Fp12_t b);

/* Square a general Fp12 element. */
void Fp12_squ(Fp12_t c, Fp12_t a);

/* Multiply an Fp12 element by an Fp6 element. */
void Fp12_mulFp6(Fp12_t c, Fp12_t a, Fp6_t b);

/* Multiply an Fp12 element by an Fp2 element. */
void Fp12_mulFp2(Fp12_t c, Fp12_t a, Fp2_t b);

/* Multiply an Fp12 element by an Fp1 element. */
void Fp12_mulFp(Fp12_t c, Fp12_t a, Fp_t b);

/* Compute the p-power Frobenius of an Fp12 element. */
void Fp12_ppow(Fp12_t c, Fp12_t a);

/* Compute the p^2-power Frobenius of an Fp12 element. */
void Fp12_p2pow(Fp12_t c, Fp12_t a);

/* Compute the p^3-power Frobenius of an Fp12 element. */
void Fp12_p3pow(Fp12_t c, Fp12_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void Fp12_p6pow(Fp12_t c, Fp12_t a);

/* Invert an Fp12 element. */
void Fp12_inv(Fp12_t c, Fp12_t a);

/* Multiply two Fp12 elements using lazy reduction.*/
void Fp12_mul_lazy(Fp12_t c, Fp12_t a, Fp12_t b);

/* Square an Fp12 element using lazy reduction.*/
void Fp12_squ_lazy(Fp12_t c, Fp12_t a);

/* Multiply two Fp12 elements (the second one sparse) using lazy reduction.*/
void Fp12_mul_sparse_twist_lazy(Fp12_t c, Fp12_t a, Fp2_t b00, Fp2_t b01, Fp2_t b11);

void Fp12_mul_sparse_untwist_lazy(Fp12_t c, Fp12_t a, Fp2_t b01, Fp2_t b11, Fp2_t b12);

/* Multiply two Fp12 elements.*/
void w2minusv_Fp12_mul(Fp12_t c, Fp12_t a, Fp12_t b, void *mem);

/* Square a general Fp12 element. */
void w2minusv_Fp12_squ(Fp12_t c, Fp12_t a, void *mem);

/* Compute the p-power Frobenius of an Fp12 element. */
void w2minusv_v3minusxi_p1mod6_Fp12_ppow(Fp12_t c, Fp12_t a);

/* Compute the p^2-power Frobenius of an Fp12 element. */
void w2minusv_v3minusxi_p1mod6_Fp12_p2pow(Fp12_t c, Fp12_t a);

/* Compute the p^3-power Frobenius of an Fp12 element. */
void w2minusv_v3minusxi_p1mod6_Fp12_p3pow(Fp12_t c, Fp12_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void w2minusv_Fp12_p6pow(Fp12_t c, Fp12_t a);

/* Invert an Fp12 element. */
void w2minusv_Fp12_inv(Fp12_t c, Fp12_t a, void *mem);

/* Multiply two Fp12 elements using lazy reduction. */
void w2minusv_Fp12_mul_lazy(Fp12_t c, Fp12_t a, Fp12_t b, void *mem);

/* Square an Fp12 element using lazy reduction. */
void w2minusv_Fp12_squ_lazy(Fp12_t c, Fp12_t a, void *mem);

/* Multiply an Fp12 element a and a sparse Fp12 element b using lazy reduction. */
void w2minusv_Fp12_mul_sparse_twist_lazy(Fp12_t c, Fp12_t a, Fp2_t b00, Fp2_t b01, Fp2_t b11, void *mem);

/* Multiply an Fp12 element a and a sparse Fp12 element b using lazy reduction. */
void w2minusv_Fp12_mul_sparse_untwist_lazy(Fp12_t c, Fp12_t a, Fp2_t b01, Fp2_t b11, Fp2_t b12, void *mem);

/***************** Internal stuff *******************/

typedef struct {
	void *mem;

  void (*Fp12_mul)(Fp12_t, Fp12_t, Fp12_t, void *);
  void (*Fp12_squ)(Fp12_t, Fp12_t, void *);
  void (*Fp12_inv)(Fp12_t, Fp12_t, void *);
  void (*Fp12_ppow)(Fp12_t, Fp12_t);
  void (*Fp12_p2pow)(Fp12_t, Fp12_t);
  void (*Fp12_p3pow)(Fp12_t, Fp12_t);
  void (*Fp12_p6pow)(Fp12_t, Fp12_t);
  void (*Fp12_mul_lazy)(Fp12_t, Fp12_t, Fp12_t, void *);
  void (*Fp12_squ_lazy)(Fp12_t, Fp12_t, void *);
  void (*Fp12_mul_sparse_twist_lazy)(Fp12_t, Fp12_t, Fp2_t, Fp2_t, Fp2_t, void *);
  void (*Fp12_mul_sparse_untwist_lazy)(Fp12_t, Fp12_t, Fp2_t, Fp2_t, Fp2_t, void *);

  Fp2_t wppow; // w^(p-1)
  Fp_t wp2pow; // w^(p^2-1)
  Fp_t wp2pow3; // (w^3)^(p^2-1)
  Fp2_t wp3pow;// w^(p^3-1)
  Fp2_t wppow3;// (w^3)^(p-1)
  Fp2_t wppow3inv; // wpppow3^(-1)
} Fp12_config_t;

extern Fp12_config_t Fp12_config; /* Defined in Fp12.c */

#endif /* __FP12_H */
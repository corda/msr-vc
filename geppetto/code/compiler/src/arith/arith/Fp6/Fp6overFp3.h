/* Fp6/Fp6overFp3.h
 * by Joppe W. Bos and Michael Naehrig (mnaehrig@microsoft.com), 
 * Cryptography Research Group, MSR Redmond
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

#ifndef __FP6VIAFP3_H
#define __FP6VIAFP3_H

#include "../arith_error.h"
#include "../Fp3/Fp3.h"


/* This contructions realizes F_{p^6} as a quadratic extension over F_{p^3}
 * by an irreducible polynomial X^2 - v, where v is in F_{p^3}.
 */

/* An element in F_{p^6} consists of two coefficients from F_{p^3}, 
 * i.e. a = a0 + a1*w in F_{p^6} for w in F_{p^6} with w^2 = v.
 * Note that sometimes v is denoted u here because the implementation of 
 * Fp3 uses u as its primitive element.
 */
typedef struct {
  Fp3_t a0;
  Fp3_t a1;
} fp6q_t;

/* Our basic type is a pointer to the struct. */
typedef fp6q_t Fp6q_t[1];


int Fp6q_initialize_config (void);

void Fp6q_free_config (void);


/************** Regular Fp6q functions ***************/

/* Initialize an element in Fp6overFp3 by initializing two Fp3 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp6q_init(Fp6q_t dst);

/* Free the memory allocated in src. */
void Fp6q_free (Fp6q_t src);

/* Get a random element in Fp6overFp3. */
void Fp6q_rand (Fp6q_t c);

/* Print the (regular non-Montgomery form) value of a */
void Fp6q_print (Fp6q_t a);

/* Set c to 0. */
void Fp6q_set_zero (Fp6q_t c);

/* Set c to 1. */
void Fp6q_set_one (Fp6q_t c);

/* Copy an Fp6q element a to an Fp6q element c. */
void Fp6q_copy(Fp6q_t c, Fp6q_t a);

/* Compare two Fp6q elements a and b for equality. */
int Fp6q_cmpeq(Fp6q_t a, Fp6q_t b);

/* Negate an Fp6q element. */
void Fp6q_neg(Fp6q_t c, Fp6q_t a);

/* Add two Fp6q elements coefficient wise. */
void Fp6q_add (Fp6q_t c, Fp6q_t a, Fp6q_t b);

/* Subtract an Fp6q element b from an Fp6q element a coefficient wise. */
void Fp6q_sub (Fp6q_t c, Fp6q_t a, Fp6q_t b);

/* Multiply two Fp6q elements.*/
void Fp6q_mul(Fp6q_t c, Fp6q_t a, Fp6q_t b);

/* Square a general Fp6q element. */
void Fp6q_squ(Fp6q_t c, Fp6q_t a);

/* Multiply an Fp6q element by an Fp3 element. */
void Fp6q_mulFp6(Fp6q_t c, Fp6q_t a, Fp3_t b);

/* Multiply an Fp6q element by an Fp1 element. */
void Fp6q_mulFp(Fp6q_t c, Fp6q_t a, Fp_t b);

/* Compute the p-power Frobenius of an Fp6q element. */
void Fp6q_ppow(Fp6q_t c, Fp6q_t a);

/* Compute the p^2-power Frobenius of an Fp6q element. */
void Fp6q_p2pow(Fp6q_t c, Fp6q_t a);

/* Compute the p^3-power Frobenius of an Fp6q element. */
void Fp6q_p3pow(Fp6q_t c, Fp6q_t a);

/* Invert an Fp6q element. */
void Fp6q_inv(Fp6q_t c, Fp6q_t a);

/* Multiply two Fp6q elements using lazy reduction.*/
void Fp6q_mul_lazy(Fp6q_t c, Fp6q_t a, Fp6q_t b);

/* Square an Fp6q element using lazy reduction.*/
void Fp6q_squ_lazy(Fp6q_t c, Fp6q_t a);

/* Multiply two Fp6q elements (the second one sparse) using lazy reduction.*/
void Fp6q_mul_sparse_twist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b00, Fp_t b01, Fp_t b11);

void Fp6q_mul_sparse_untwist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b01, Fp_t b11, Fp_t b12);

/* Multiply two Fp6q elements.*/
void w2minusv_Fp6q_mul(Fp6q_t c, Fp6q_t a, Fp6q_t b, void *mem);

/* Square a general Fp6q element. */
void w2minusv_Fp6q_squ(Fp6q_t c, Fp6q_t a, void *mem);

/* Compute the p-power Frobenius of an Fp6q element. */
void w2minusv_v3minusxi_p1mod6_Fp6q_ppow(Fp6q_t c, Fp6q_t a);

/* Compute the p^2-power Frobenius of an Fp6q element. */
void w2minusv_v3minusxi_p1mod6_Fp6q_p2pow(Fp6q_t c, Fp6q_t a);

/* Compute the p^3-power Frobenius of an Fp6q element. */
void w2minusv_Fp6q_p3pow(Fp6q_t c, Fp6q_t a);

/* Invert an Fp6q element. */
void w2minusv_Fp6q_inv(Fp6q_t c, Fp6q_t a, void *mem);

/* Multiply two Fp6q elements using lazy reduction. */
void w2minusv_Fp6q_mul_lazy(Fp6q_t c, Fp6q_t a, Fp6q_t b, void *mem);

/* Square an Fp6q element using lazy reduction. */
void w2minusv_Fp6q_squ_lazy(Fp6q_t c, Fp6q_t a, void *mem);

/* Multiply an Fp6q element a and a sparse Fp6q element b using lazy reduction. */
void w2minusv_Fp6q_mul_sparse_twist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b00, Fp_t b01, Fp_t b11, void *mem);

/* Multiply an Fp6q element a and a sparse Fp6q element b using lazy reduction. */
void w2minusv_Fp6q_mul_sparse_untwist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b01, Fp_t b11, Fp_t b12, void *mem);

/***************** Internal stuff *******************/

typedef struct {
	void *mem;

  void (*Fp6q_mul)(Fp6q_t, Fp6q_t, Fp6q_t, void *);
  void (*Fp6q_squ)(Fp6q_t, Fp6q_t, void *);
  void (*Fp6q_inv)(Fp6q_t, Fp6q_t, void *);
  void (*Fp6q_ppow)(Fp6q_t, Fp6q_t);
  void (*Fp6q_p2pow)(Fp6q_t, Fp6q_t);
  void (*Fp6q_p3pow)(Fp6q_t, Fp6q_t);
  void (*Fp6q_mul_lazy)(Fp6q_t, Fp6q_t, Fp6q_t, void *);
  void (*Fp6q_squ_lazy)(Fp6q_t, Fp6q_t, void *);
  void (*Fp6q_mul_sparse_twist_lazy)(Fp6q_t, Fp6q_t, Fp_t, Fp_t, Fp_t, void *);
  void (*Fp6q_mul_sparse_untwist_lazy)(Fp6q_t, Fp6q_t, Fp_t, Fp_t, Fp_t, void *);

  Fp_t wppow; // w^(p-1)
  Fp_t wp2pow; // w^(p^2-1)
  Fp_t wp2pow3; // (w^3)^(p^2-1)
  Fp_t wp3pow;// w^(p^3-1)
  Fp_t wppow3;// (w^3)^(p-1)
  Fp_t wppow3inv; // wpppow3^(-1)
} Fp6q_config_t;

extern Fp6q_config_t Fp6q_config; /* Defined in Fp6overFp3.c */

#endif /* __FP6VIAFP3_H */
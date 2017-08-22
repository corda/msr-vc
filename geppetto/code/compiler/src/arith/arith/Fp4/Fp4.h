/* Fp4/Fp4.h
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

#ifndef __FP4_H
#define __FP4_H

#include "../arith_error.h"
#include "../Fp2/Fp2.h"

/* This contructions realizes F_{p^4} as a quadratic extension over F_{p^2}
 * by an irreducible polynomial X^2 - xi, where xi is in F_{p^2}.
 */

/* An element in F_{p^4} consists of two coefficients from F_{p^2}, 
 * i.e. a = a0 + a1*s in F_{p^4} for s in F_{p^4} with s^2 = xi. 
 */
typedef struct {
  Fp2_t a0;
  Fp2_t a1;
} fp4_t;

/* Our basic type is a pointer to the struct. */
typedef fp4_t Fp4_t[1];

int Fp4_initialize_config (void);

void Fp4_free_config (void);


/* Initialize an element in Fp4 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp4_init(Fp4_t dst);

/* Free the memory allocated in src. */
void Fp4_free (Fp4_t src);

/* Get a random element in Fp4. */
void Fp4_rand (Fp4_t c);

/* Set c to 0. */
void Fp4_set_zero (Fp4_t c);

/* Set c to 1. */
void Fp4_set_one (Fp4_t c);

/* Copy an Fp4 element a to an Fp4 element c. */
void Fp4_copy(Fp4_t c, Fp4_t a);

/* Print the (regular non-Montgomery form) value of a */
void Fp4_print (Fp4_t a);

/* Compare two Fp4 elements a and b for equality. */
int Fp4_cmpeq(Fp4_t a, Fp4_t b);

/* Negate an Fp4 element. */
void Fp4_neg(Fp4_t c, Fp4_t a);

/* Add two Fp4 elements coefficient wise. */
void Fp4_add (Fp4_t c, Fp4_t a, Fp4_t b);

/* Subtract an Fp4 element b from an Fp4 element a coefficient wise. */
void Fp4_sub (Fp4_t c, Fp4_t a, Fp4_t b);

/* Divide an Fp4 element a by 2 coefficient wise. */
void Fp4_div2 (Fp4_t c, Fp4_t a);

/* Multiply an Fp4 element by an Fp2 element. */
void Fp4_mulFp2(Fp4_t c, Fp4_t a, Fp2_t b);

/* Multiply an Fp4 element by an Fp element. */
void Fp4_mulFp(Fp4_t c, Fp4_t a, Fp_t b);

/* Multiply two Fp4 elements. */
void Fp4_mul (Fp4_t c, Fp4_t a, Fp4_t b);

/* Multiply two Fp4 elements using lazy reduction. */
void Fp4_mul_lazy (Fp4_t c, Fp4_t a, Fp4_t b);

/* Multiply compute 3*a. */
void Fp4_mul3 (Fp4_t c, Fp4_t a);

/* Square an Fp4 element. */
void Fp4_squ (Fp4_t c, Fp4_t a);

/* Square an Fp4 element using lazy reduction. */
void Fp4_squ_lazy (Fp4_t c, Fp4_t a);

/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void Fp4_squ_sep (Fp2_t c0, Fp2_t c1, Fp2_t a0, Fp2_t a1);

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2. */
void Fp4_muls (Fp4_t c, Fp4_t a);

/* Multiply an Fp4 element by the special element i in Fp2.*/
void Fp4_muli(Fp4_t c, Fp4_t a);

/* Multiply an Fp4 element by the special element xi in Fp2.*/
void Fp4_mulxi(Fp4_t c, Fp4_t a);

/* Compute the p^2-power Frobenius of an Fp4 element. */
void Fp4_p2pow(Fp4_t c, Fp4_t a);

/* Compute the p-power Frobenius of an Fp4 element. */
void Fp4_ppow(Fp4_t c, Fp4_t a);

/* Invert an Fp4 element. */
void Fp4_inv (Fp4_t c, Fp4_t a);

/* Multiply two Fp4 elements. */
void s2minusxi_Fp4_mul(Fp4_t c, Fp4_t a, Fp4_t b, void *mem);

/* Multiply two Fp4 elements using lazy reduction. */
void s2minusxi_Fp4_mul_lazy(Fp4_t c, Fp4_t a, Fp4_t b, void *mem);

/* Square an Fp4 element. */
void s2minusxi_Fp4_squ(Fp4_t c, Fp4_t a, void *mem);

/* Square an Fp4 element using lazy reduction. */
void s2minusxi_Fp4_squ_lazy(Fp4_t c, Fp4_t a, void *mem);

/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void s2minusxi_Fp4_squ_sep(Fp2_t c0, Fp2_t c1, Fp2_t a0, Fp2_t a1, void *mem);

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2. */
void s2minusxi_Fp4_muls(Fp4_t c, Fp4_t a, void *mem);

/* Inversion of an Fp4 element. This uses that Fp4 is constructed over Fp2 via s^2 - xi. */
void s2minusxi_Fp4_inv(Fp4_t c, Fp4_t a, void *mem);

/* Compute the p^2-power Frobenius of an Fp4 element. */
void binomial_Fp4_p2pow(Fp4_t c, Fp4_t a);

/* Compute the p-power Frobenius of an Fp4 element. */
void binomial_Fp4_ppow(Fp4_t c, Fp4_t a);

/* Multiply an Fp4 element by the special element i in Fp2. */
void p3mod4_Fp4_muli(Fp4_t c, Fp4_t a);

/* Multiply an Fp4 element by the special element xi in Fp2. */
void p3mod4_Fp4_mulxi(Fp4_t c, Fp4_t a);

/***************** Internal stuff *******************/

typedef struct {
	void *mem;

  void (*Fp4_mul)(Fp4_t, Fp4_t, Fp4_t, void *);
  void (*Fp4_squ)(Fp4_t, Fp4_t, void *);
  void (*Fp4_squ_sep)(Fp2_t, Fp2_t, Fp2_t, Fp2_t, void *);
  void (*Fp4_muls)(Fp4_t, Fp4_t, void *);
  void (*Fp4_muli)(Fp4_t, Fp4_t);
  void (*Fp4_mulxi)(Fp4_t, Fp4_t);
  void (*Fp4_p2pow)(Fp4_t, Fp4_t);
  void (*Fp4_ppow)(Fp4_t, Fp4_t);
  void (*Fp4_inv)(Fp4_t, Fp4_t, void *);
  void (*Fp4_mul_lazy)(Fp4_t, Fp4_t, Fp4_t, void *);
  void (*Fp4_squ_lazy)(Fp4_t, Fp4_t, void *);

  Fp2_t sppow;
} Fp4_config_t;

extern Fp4_config_t Fp4_config; /* Defined in Fp4.c */

#endif /* __FP4_H */
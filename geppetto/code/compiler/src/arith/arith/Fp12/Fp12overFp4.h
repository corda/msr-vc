/* Fp12/Fp12overFp4.h
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

#ifndef __FP12VIAFP4_H
#define __FP12VIAFP4_H

#include "../arith_error.h"
#include "../Fp4/Fp4.h"

/* This contructions realizes F_{p^12} as a cubic extension over F_{p^4}
 * by an irreducible polynomial X^3 - s, where s is in F_{p^4}.
 */

/* An element in F_{p^12} consists of three coefficients from F_{p^4}, 
 * i.e. a = a0 + a1*v + a2*v^2 in F_{p^12} for v in F_{p^12} with v^3 = s. 
 */
typedef struct {
  Fp4_t a0;
  Fp4_t a1;
  Fp4_t a2;
} fp12c_t;

/* Our basic type is a pointer to the struct. */
typedef fp12c_t Fp12c_t[1];

int Fp12c_initialize_config (void);

void Fp12c_free_config (void);


/************** Regular Fp12 functions ***************/

/* Initialize an element in Fp12 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp12c_init(Fp12c_t dst);

/* Free the memory allocated in src. */
void Fp12c_free (Fp12c_t src);

/* Get a random element in Fp12. */
void Fp12c_rand (Fp12c_t c);

/* Set c to 0. */
void Fp12c_set_zero (Fp12c_t c);

/* Set c to 1. */
void Fp12c_set_one (Fp12c_t c);

/* Copy an Fp12 element a to an Fp12 element c. */
void Fp12c_copy(Fp12c_t c, Fp12c_t a);

/* Print the (regular non-Montgomery form) value of a */
void Fp12c_print (Fp12c_t a);

/* Compare two Fp12 elements a and b for equality. */
int Fp12c_cmpeq(Fp12c_t a, Fp12c_t b);

/* Negate an Fp12 element. */
void Fp12c_neg(Fp12c_t c, Fp12c_t a);

/* Add two Fp12 elements coefficient wise. */
void Fp12c_add (Fp12c_t c, Fp12c_t a, Fp12c_t b);

/* Subtract an Fp12 element b from an Fp12 element a coefficient wise. */
void Fp12c_sub (Fp12c_t c, Fp12c_t a, Fp12c_t b);

/* Multiply two Fp12 elements. */
void Fp12c_mul (Fp12c_t c, Fp12c_t a, Fp12c_t b);

/* Square an Fp12 element. */
void Fp12c_squ (Fp12c_t c, Fp12c_t a);

/* Multiply an Fp12 element by an Fp element. */
void Fp12c_mulFp (Fp12c_t c, Fp12c_t a, Fp_t b);

/* Multiply an Fp12 element by an Fp2 element. */
void Fp12c_mulFp2 (Fp12c_t c, Fp12c_t a, Fp2_t b);

/* Multiply an Fp12 element by an Fp2 element. */
void Fp12c_mulFp4 (Fp12c_t c, Fp12c_t a, Fp4_t b);

/* Invert an Fp12 element. */
void Fp12c_inv (Fp12c_t c, Fp12c_t a);

/* Multiply an Fp12 element by the special element v. */
void Fp12c_mulv(Fp12c_t c, Fp12c_t a);

/* Multiply an Fp12 element by the special element i. */
void Fp12c_muli(Fp12c_t c, Fp12c_t a);

/* Multiply an Fp12 element by the special element s. */
void Fp12c_muls(Fp12c_t c, Fp12c_t a);

/* Compute the p-power Frobenius of an Fp12 element. */
void Fp12c_ppow(Fp12c_t c, Fp12c_t a);

/* Compute the p^2-power Frobenius of an Fp12 element. */
void Fp12c_p2pow(Fp12c_t c, Fp12c_t a);

/* Compute the p^4-power Frobenius of an Fp12 element. */
void Fp12c_p4pow(Fp12c_t c, Fp12c_t a);

/* Multiply an Fp12c element with a sparse Fp12c element. */
void Fp12c_mul_sparse01 (Fp12c_t c, Fp12c_t a, Fp4_t b0, Fp4_t b1);
void Fp12c_mul_sparse12 (Fp12c_t c, Fp12c_t a, Fp4_t b1, Fp4_t b2);

/* Multiply two Fp12 elements. */
void v3minusxi_Fp12c_mul(Fp12c_t c, Fp12c_t a, Fp12c_t b, void *mem);

/* Square an Fp12 element. */
void v3minuss_Fp12c_squ(Fp12c_t c, Fp12c_t a, void *mem);

/* Multiply an Fp12 element by the special element v. */
void v3minuss_Fp12c_mulv(Fp12c_t c, Fp12c_t a, void *mem);

/* Multiply an Fp12 element by the special element i. */
void p3mod4_Fp12c_muli(Fp12c_t c, Fp12c_t a);

/* Multiply an Fp12 element by the special element s. */
void s2minusxi_Fp12c_muls(Fp12c_t c, Fp12c_t a);

/* Compute the p-power Frobenius of an Fp12 element. */
void p1mod3_v3minuss_Fp12c_ppow(Fp12c_t c, Fp12c_t a);

/* Compute the p^2-power Frobenius of an Fp12 element. */
void p1mod3_v3minuss_Fp12c_p2pow(Fp12c_t c, Fp12c_t a);

/* Compute the p^4-power Frobenius of an Fp12 element. */
void p1mod3_v3minuss_Fp12c_p4pow(Fp12c_t c, Fp12c_t a);

/* Inversion of an Fp12 element. */
void v3minuss_Fp12c_inv(Fp12c_t c, Fp12c_t a, void *mem);

/* Multiply an Fp12c element with a 12-sparse Fp12c element. */
void v3minuss_Fp12c_mul_sparse01 (Fp12c_t c, Fp12c_t a, Fp4_t b1, Fp4_t b2, void *mem);
void v3minuss_Fp12c_mul_sparse12 (Fp12c_t c, Fp12c_t a, Fp4_t b1, Fp4_t b2, void *mem);


/***************** Internal stuff *******************/

typedef struct {
	void *mem;
  void (*Fp12c_mul)(Fp12c_t, Fp12c_t, Fp12c_t, void *);
  void (*Fp12c_squ)(Fp12c_t, Fp12c_t, void *);
  void (*Fp12c_inv)(Fp12c_t, Fp12c_t, void *);
  void (*Fp12c_mulv)(Fp12c_t, Fp12c_t, void *);
  void (*Fp12c_ppow)(Fp12c_t, Fp12c_t);
  void (*Fp12c_p2pow)(Fp12c_t, Fp12c_t);
  void (*Fp12c_p4pow)(Fp12c_t, Fp12c_t);
  void (*Fp12c_muli)(Fp12c_t, Fp12c_t);
  void (*Fp12c_muls)(Fp12c_t, Fp12c_t);
  void (*Fp12c_mul_sparse01)(Fp12c_t, Fp12c_t, Fp4_t, Fp4_t, void *);
  void (*Fp12c_mul_sparse12)(Fp12c_t, Fp12c_t, Fp4_t, Fp4_t, void *);

  Fp2_t vppow;  // v^(p-1)  
  Fp2_t vppow2; // vppow^2
  Fp_t vp2pow;  // v^(p^2-1)
  Fp_t vp2pow2; // vp2pow^2
  Fp_t vp4pow;  // v^(p^4-1)
  Fp_t vp4pow2; // vp4pow^2
} Fp12c_config_t;

extern Fp12c_config_t Fp12c_config; /* Defined in Fp12viaFp4.c */

#endif /* __Fp12VIAFP4_H */
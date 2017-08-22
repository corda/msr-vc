/* Fp6/Fp6.h
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

#ifndef __FP6_H
#define __FP6_H

#include "../arith_error.h"
#include "../Fp2/Fp2.h"

/* This contructions realizes F_{p^6} as a cubic extension over F_{p^2}
 * by an irreducible polynomial X^3 - xi, where xi is in F_{p^2}.
 */

/* An element in F_{p^6} consists of three coefficients from F_{p^2}, 
 * i.e. a = a0 + a1*v + a2*v^2 in F_{p^6} for v in F_{p^6} with v^3 = xi. 
 */
typedef struct {
  Fp2_t a0;
  Fp2_t a1;
  Fp2_t a2;
} fp6_t;

/* Our basic type is a pointer to the struct. */
typedef fp6_t Fp6_t[1];

/* A long element in F_{p^6} of double the size of a regular F_{p^6} element.
 * It consists of three long (double size) coefficients from F_{p^2} (type fp2l_t), 
 * i.e. A = A0 + A1*v + A2*v^2 in F_{p^6} for v in F_{p^6} with v^3 = xi.
 * Such elements are used for lazy reduction techniques where reduction is 
 * separated from multiplication and is postponed until after possibly several
 * other operations.
 */
typedef struct {
  Fp2l_t A0;
  Fp2l_t A1;
  Fp2l_t A2;
} fp6l_t;

/* Our basic type is a pointer to the struct. */
typedef fp6l_t Fp6l_t[1];

int Fp6_initialize_config (void);

void Fp6_free_config (void);


/************** Regular Fp6 functions ***************/

/* Initialize an element in Fp6 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp6_init(Fp6_t dst);

/* Free the memory allocated in src. */
void Fp6_free (Fp6_t src);

/* Get a random element in Fp6. */
void Fp6_rand (Fp6_t c);

/* Set c to 0. */
void Fp6_set_zero (Fp6_t c);

/* Set c to 1. */
void Fp6_set_one (Fp6_t c);

/* Copy an Fp6 element a to an Fp6 element c. */
void Fp6_copy(Fp6_t c, Fp6_t a);

/* Print the (regular non-Montgomery form) value of a */
void Fp6_print (Fp6_t a);

/* Compare two Fp6 elements a and b for equality. */
int Fp6_cmpeq(Fp6_t a, Fp6_t b);

/* Negate an Fp6 element. */
void Fp6_neg(Fp6_t c, Fp6_t a);

/* Add two Fp6 elements coefficient wise. */
void Fp6_add (Fp6_t c, Fp6_t a, Fp6_t b);

/* Subtract an Fp6 element b from an Fp6 element a coefficient wise. */
void Fp6_sub (Fp6_t c, Fp6_t a, Fp6_t b);

/* Multiply two Fp6 elements. */
void Fp6_mul (Fp6_t c, Fp6_t a, Fp6_t b);

/* Square an Fp6 element. */
void Fp6_squ (Fp6_t c, Fp6_t a);

/* Multiply an Fp6 element by an Fp element. */
void Fp6_mulFp (Fp6_t c, Fp6_t a, Fp_t b);

/* Multiply an Fp6 element by an Fp2 element. */
void Fp6_mulFp2 (Fp6_t c, Fp6_t a, Fp2_t b);

/* Invert an Fp6 element. */
void Fp6_inv (Fp6_t c, Fp6_t a);

/* Multiply an fp6 element by the special element v. */
void Fp6_mulv(Fp6_t c, Fp6_t a);

/* Compute the p-power Frobenius of an Fp6 element. */
void Fp6_ppow(Fp6_t c, Fp6_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void Fp6_p2pow(Fp6_t c, Fp6_t a);

/* Compute the p^3-power Frobenius of an Fp6 element. */
void Fp6_p3pow(Fp6_t c, Fp6_t a);

/* Multiply two Fp6 elements. */
void v3minusxi_Fp6_mul(Fp6_t c, Fp6_t a, Fp6_t b, void *mem);

/* Square an Fp6 element. */
void v3minusxi_Fp6_squ(Fp6_t c, Fp6_t a, void *mem);

/* Multiply an fp6 element by the special element v. */
void v3minusxi_Fp6_mulv(Fp6_t c, Fp6_t a, void *mem);

/* Compute the p-power Frobenius of an Fp6 element. */
void p1mod3_v3minusxi_Fp6_ppow(Fp6_t c, Fp6_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void p1mod3_v3minusxi_Fp6_p2pow(Fp6_t c, Fp6_t a);

/* Compute the p^3-power Frobenius of an Fp6 element. */
void p1mod3_v3minusxi_Fp6_p3pow(Fp6_t c, Fp6_t a);

/* Inversion of an Fp6 element. */
void v3minusxi_Fp6_inv(Fp6_t c, Fp6_t a, void *mem);


/************** Fp6 functions for lazy reduction ***************/

/* Initialize a long Fp6 element. */
int Fp6l_init(Fp6l_t dst);

/* Free the memory allocated in src. */
void Fp6l_free (Fp6l_t src);

/* Copy an Fp6l element a to an Fp6l element c. */
void Fp6l_copy(Fp6l_t C, Fp6l_t A);

/* Separate reduction function calling the Fp2 reduction for each coefficient. */
void Fp6l_red(Fp6_t c, Fp6l_t A);

/* Subtraction of Fp6l elements, constant-time version of option 2 in Aranha et al. */
void Fp6l_sub_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b);

/* Addition of Fp6l elements, constant-time version of "option 2" addition in Aranha et al. */
void Fp6l_add_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b);

/* Multiply an Fp6 element by an Fp element. */
void Fp6_mulFp2_no_red (Fp6l_t c, Fp6_t a, Fp2_t b);

/* Multiply two Fp6 elements without reduction. */
void Fp6_mul_no_red (Fp6l_t c, Fp6_t a, Fp6_t b);

/* Multiply an Fp6 element and a sparse Fp6 element without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void Fp6_mul_sparse01_no_red (Fp6l_t c, Fp6_t a, Fp2_t b0, Fp2_t b1);
void Fp6_mul_sparse12_no_red (Fp6l_t c, Fp6_t a, Fp2_t b0, Fp2_t b1);

/* Multiply two Fp6 elements without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void v3minusxi_Fp6_mul_no_red (Fp6l_t c, Fp6_t a, Fp6_t b, void *mem);

/* Multiply an Fp6 element and a sparse Fp6 element without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void v3minusxi_Fp6_mul_sparse01_no_red (Fp6l_t c, Fp6_t a, Fp2_t b0, Fp2_t b1, void *mem);
void v3minusxi_Fp6_mul_sparse12_no_red (Fp6l_t c, Fp6_t a, Fp2_t b1, Fp2_t b2, void *mem);


/***************** Internal stuff *******************/

typedef struct {
	void *mem;
  void (*Fp6_mul)(Fp6_t, Fp6_t, Fp6_t, void *);
  void (*Fp6_squ)(Fp6_t, Fp6_t, void *);
  void (*Fp6_inv)(Fp6_t, Fp6_t, void *);
  void (*Fp6_mulv)(Fp6_t, Fp6_t, void *);
  void (*Fp6_ppow)(Fp6_t, Fp6_t);
  void (*Fp6_p2pow)(Fp6_t, Fp6_t);
  void (*Fp6_p3pow)(Fp6_t, Fp6_t);
  void (*Fp6_mul_no_red)(Fp6l_t, Fp6_t, Fp6_t, void *);
  void (*Fp6_mul_sparse01_no_red)(Fp6l_t, Fp6_t, Fp2_t, Fp2_t, void *);
  void (*Fp6_mul_sparse12_no_red)(Fp6l_t, Fp6_t, Fp2_t, Fp2_t, void *);

  Fp_t vppow;   // v^(p-1)/i  
  Fp_t vppow2;  // (i*vppow)^2 = -vppow^2
  Fp_t vp2pow;  // v^(p^2-1)
  Fp_t vp2pow2; // vp2pow^2
  Fp_t vppowinv; // -i/v^(p-1) = 1/v^(p-1)/i
} Fp6_config_t;

extern Fp6_config_t Fp6_config; /* Defined in Fp6.c */

#endif /* __FP6_H */
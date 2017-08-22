/* Fp3/Fp3.h
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

#ifndef __FP3_H
#define __FP3_H

#include "../arith_error.h"
#include "../Fp/Fp.h"

/* This contruction realizes F_{p^3} over F_{p}
 * by an irreducible polynomial X^3 - xi, where xi is in F_{p}.
 */

/* An element in F_{p^3} consists of three coefficients from F_{p}, 
 * i.e. a = a0 + a1*u + a2*u^2 in F_{p^3} for u in F_{p^3} with u^3 = xi. 
 */
typedef struct {
  Fp_t a0;
  Fp_t a1;
  Fp_t a2;
} fp3_t;

/* Our basic type is a pointer to the struct. */
typedef fp3_t Fp3_t[1];

/* A long element in F_{p^3} of double the size of a regular F_{p^3} element.
 * It consists of three long (double size) coefficients from F_p (type fpl_t), 
 * i.e. A = A0 + A1*u + A2*u^2 in F_{p^3} for u in F_{p^3} with u^3 = xi.
 * Such elements are used for lazy reduction techniques where reduction is 
 * separated from multiplication and is postponed until after possibly several
 * other operations.
 */
typedef struct {
  Fpl_t A0;
  Fpl_t A1;
  Fpl_t A2;
} fp3l_t;

/* Our basic type is a pointer to the struct. */
typedef fp3l_t Fp3l_t[1];



int Fp3_initialize_config (void);

void Fp3_free_config (void);


/************** Regular Fp3 functions ***************/

/* Initialize an element in Fp3 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp3_init(Fp3_t dst);

/* Free the memory allocated in src. */
void Fp3_free (Fp3_t src);

/* Get a random element in Fp3. */
void Fp3_rand (Fp3_t c);

/* Set c to 0. */
void Fp3_set_zero (Fp3_t c);

/* Set c to 1. */
void Fp3_set_one (Fp3_t c);

/* Copy an Fp3 element a to an Fp3 element c. */
void Fp3_copy(Fp3_t c, Fp3_t a);

/* Print the (regular non-Montgomery form) value of a */
void Fp3_print (Fp3_t a);

/* Compare two Fp3 elements a and b for equality. */
int Fp3_cmpeq(Fp3_t a, Fp3_t b);

/* Negate an Fp3 element. */
void Fp3_neg(Fp3_t c, Fp3_t a);

/* Add two Fp3 elements coefficient wise. */
void Fp3_add (Fp3_t c, Fp3_t a, Fp3_t b);

/* Subtract an Fp3 element b from an Fp3 element a coefficient wise. */
void Fp3_sub (Fp3_t c, Fp3_t a, Fp3_t b);

/* Multiply two Fp3 elements. */
void Fp3_mul (Fp3_t c, Fp3_t a, Fp3_t b);

/* Compute 3*a. */
void Fp3_mul3 (Fp3_t c, Fp3_t a);

/* Square an Fp3 element. */
void Fp3_squ (Fp3_t c, Fp3_t a);

/* Multiply an Fp3 element by an Fp element. */
void Fp3_mulFp (Fp3_t c, Fp3_t a, Fp_t b);

/* Invert an Fp3 element. */
void Fp3_inv (Fp3_t c, Fp3_t a);

/* Multiply an Fp3 element by the special element u. */
void Fp3_mulu(Fp3_t c, Fp3_t a);

/* Compute the p-power Frobenius of an Fp3 element. */
void Fp3_ppow(Fp3_t c, Fp3_t a);

/* Compute the p-power Frobenius of an Fp3 element. */
void p1mod3_Fp3_ppow(Fp3_t c, Fp3_t a);

/* Multiply two Fp3 elements. */
void u3plus2_Fp3_mul(Fp3_t c, Fp3_t a, Fp3_t b, void *mem);

/* Square an Fp3 element. */
void u3plus2_Fp3_squ(Fp3_t c, Fp3_t a, void *mem);

/* Multiply an Fp3 element by the special element u. */
void u3plus2_Fp3_mulu(Fp3_t c, Fp3_t a, void *mem);

/* Inversion of an Fp3 element. */
void u3plus2_Fp3_inv(Fp3_t c, Fp3_t a, void *mem);

/* Multiply two Fp3 elements. */
void u3minus3_Fp3_mul(Fp3_t c, Fp3_t a, Fp3_t b, void *mem);

/* Square an Fp3 element. */
void u3minus3_Fp3_squ(Fp3_t c, Fp3_t a, void *mem);

/* Multiply an Fp3 element by the special element u. */
void u3minus3_Fp3_mulu(Fp3_t c, Fp3_t a, void *mem);

/* Inversion of an Fp3 element. */
void u3minus3_Fp3_inv(Fp3_t c, Fp3_t a, void *mem);


/************** Fp3 functions for lazy reduction ***************/

/* Initialize a long Fp3 element. */
int Fp3l_init(Fp3l_t dst);

/* Free the memory allocated in src. */
void Fp3l_free (Fp3l_t src);

/* Copy an Fp3l element a to an Fp3l element c. */
void Fp3l_copy(Fp3l_t C, Fp3l_t A);

/* Separate reduction function calling the Fp reduction for each coefficient. */
void Fp3l_red(Fp3_t c, Fp3l_t A);

/* Add two Fp3 elements coefficient wise without reducing mod p for lazy reduction. */ 
void Fp3_add_no_red (Fp3_t c, Fp3_t a, Fp3_t b);

/* Add two long Fp3 elements coefficient wise without reducing mod p for lazy reduction. */
void Fp3l_add_no_red (Fp3l_t C, Fp3l_t A, Fp3l_t B);

/* Subtraction of Fp3l elements, option 1 in Aranha et al. */
void Fp3l_sub_o1 (Fp3l_t c, Fp3l_t a, Fp3l_t b, int h);

/* Subtraction of Fp3l elements, constant-time version of option 2 in Aranha et al. */
void Fp3l_sub_o2 (Fp3l_t c, Fp3l_t a, Fp3l_t b);

/* Addition of Fp3l elements, constant-time version of "option 2" addition in Aranha et al. */
void Fp3l_add_o2 (Fp3l_t c, Fp3l_t a, Fp3l_t b);

/* Multiplication by the special element v in Fp3, with v^3 = 3 in Fp. */
void Fp3l_mulu (Fp3l_t c, Fp3l_t a);

/* Multiply an Fp3 element by an Fp element. */
void Fp3_mulFp_no_red (Fp3l_t c, Fp3_t a, Fp_t b);

/* Multiply two Fp3 elements without reduction. */
void Fp3_mul_no_red(Fp3l_t c, Fp3_t a, Fp3_t b);
void Fp3_mul_no_red_o1(Fp3l_t c, Fp3_t a, Fp3_t b, int h);
void Fp3_mul_no_red_o2(Fp3l_t c, Fp3_t a, Fp3_t b);

/* Multiply in Fp3 using the lazy reduction multiplications and 3 reductions only. */
void Fp3_mul_lazy (Fp3_t c, Fp3_t a, Fp3_t b);

/* Square an Fp3 element without reduction. */
void Fp3_squ_no_red(Fp3l_t c, Fp3_t a);
void Fp3_squ_no_red_o1(Fp3l_t c, Fp3_t a, int h);
void Fp3_squ_no_red_o2(Fp3l_t c, Fp3_t a);

/* Square in Fp3 using the lazy reduction multiplications and 3 reductions only. */
void Fp3_squ_lazy (Fp3_t c, Fp3_t a);

/* Multiplication by the special element v in Fp3, with v^3 = 3 in Fp.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = 3*a2 + a0*v + a1*v^2.
 */
void u3minus3_Fp3l_mulu(Fp3l_t c, Fp3l_t a, void *mem);

/* Multiply two Fp3 elements without reduction resulting in a long Fp3l element.
 * This function assumes that u^3 = -2.
 */
void u3plus2_Fp3_mul_no_red_o1(Fp3l_t C, Fp3_t a, Fp3_t b, int h, void *mem);
void u3plus2_Fp3_mul_no_red_o2(Fp3l_t C, Fp3_t a, Fp3_t b, void *mem);

/* Squaring an Fp3 element without reduction resulting in a long Fp3l element. */
void u3plus2_Fp3_squ_no_red_o1(Fp3l_t C, Fp3_t a, int h, void *mem);
void u3plus2_Fp3_squ_no_red_o2(Fp3l_t C, Fp3_t a, void *mem);

/* Multiply two Fp3 elements without reduction resulting in a long Fp3l element.
 * This function assumes that u^3 = 3.
 */
void u3minus3_Fp3_mul_no_red(Fp3l_t C, Fp3_t a, Fp3_t b, void *mem);

/* Squaring an Fp3 element without reduction resulting in a long Fp3l element. */
void u3minus3_Fp3_squ_no_red (Fp3l_t C, Fp3_t a, void *mem);

/* Multiply an Fp3 element and a sparse Fp3 element without reduction.
   This function assumes that Fp3 is constructed via v^3 - xi.
 */
void Fp3_mul_sparse01_no_red (Fp3l_t c, Fp3_t a, Fp_t b0, Fp_t b1);
void Fp3_mul_sparse12_no_red (Fp3l_t c, Fp3_t a, Fp_t b0, Fp_t b1);

void v3minusxi_Fp3_mul_sparse01_no_red (Fp3l_t c, Fp3_t a, Fp_t b0, Fp_t b1, void *mem);
void v3minusxi_Fp3_mul_sparse12_no_red (Fp3l_t c, Fp3_t a, Fp_t b1, Fp_t b2, void *mem);

/***************** Internal stuff *******************/

typedef struct {
	void *mem;
  void (*Fp3_mul)(Fp3_t, Fp3_t, Fp3_t, void *);
  void (*Fp3_squ)(Fp3_t, Fp3_t, void *);
  void (*Fp3_inv)(Fp3_t, Fp3_t, void *);
  void (*Fp3_mulu)(Fp3_t, Fp3_t, void *);
  void (*Fp3l_mulu)(Fp3l_t, Fp3l_t, void *);
  void (*Fp3_ppow)(Fp3_t, Fp3_t);
  void (*Fp3_mul_no_red)(Fp3l_t, Fp3_t, Fp3_t, void *);
  void (*Fp3_mul_no_red_o1)(Fp3l_t, Fp3_t, Fp3_t, int, void *);
  void (*Fp3_mul_no_red_o2)(Fp3l_t, Fp3_t, Fp3_t, void *);
  void (*Fp3_squ_no_red)(Fp3l_t, Fp3_t, void *);
  void (*Fp3_squ_no_red_o1)(Fp3l_t, Fp3_t, int, void *);
  void (*Fp3_squ_no_red_o2)(Fp3l_t, Fp3_t, void *);
  void (*Fp3_mul_sparse01_no_red)(Fp3l_t, Fp3_t, Fp_t, Fp_t, void *);
  void (*Fp3_mul_sparse12_no_red)(Fp3l_t, Fp3_t, Fp_t, Fp_t, void *);

  Fp_t uppow;   // v^(p-1)  
  Fp_t uppow2;  // vppow^2

  Fp3l_t T0;    // tmp element
} Fp3_config_t;

extern Fp3_config_t Fp3_config; /* Defined in Fp3.c */

#endif /* __FP3_H */
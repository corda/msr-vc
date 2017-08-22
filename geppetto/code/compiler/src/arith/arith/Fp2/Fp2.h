/* Fp2/Fp2.h
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

#ifndef __FP2_H
#define __FP2_H

#include "../arith_error.h"
#include "../Fp/Fp.h"

/* An element in F_{p^2}, consists of two coefficients from F_p, 
 * i.e. a = a0 + a1*i in F_{p^2} for i in F_{p^2} with i^2 = -1. 
 */
typedef struct {
  Fp_t a0;
  Fp_t a1;
} fp2_t;

/* Our basic type is a pointer to the struct. */
typedef fp2_t Fp2_t[1];

/* A long element in F_{p^2} of double the size of a regular F_{p^2} element.
 * It consists of two long (double size) coefficients from F_p (type fpl_t), 
 * i.e. A = A0 + A1*i in F_{p^2} for i in F_{p^2} with i^2 = -1.
 * Such elements are used for lazy reduction techniques where reduction is 
 * separated from multiplication and is postponed until after possibly several
 * other operations.
 */
typedef struct {
  Fpl_t A0;
  Fpl_t A1;
} fp2l_t;

/* Our basic type is a pointer to the struct. */
typedef fp2l_t Fp2l_t[1];


/* Initialize the global parameters for the quadratic extension
 * given the base field Fp.
 * Return > 0 on success and
 * Return < 0 on error.
 */ 
int Fp2_initialize_config (void);

void Fp2_free_config (void);


/* Initialize an element in Fp2 by initializing two Fp elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp2_init(Fp2_t dst);

/* Free the memory allocated in src. */
void Fp2_free (Fp2_t src);

/* Get a random element in Fp2. */
void Fp2_rand (Fp2_t c);

/* Set c to the value of a0 + a1*i. */
void Fp2_set (Fp2_t c, uint64_t *a0, uint64_t *a1);

/* Set c to the value x. */
void Fp2_set_ui (Fp2_t c, uint64_t x);

/* Set c to 0. */
void Fp2_set_zero (Fp2_t c);

/* Set c to 1. */
void Fp2_set_one (Fp2_t c);

/* Print the (regular non-Montgomery form) value of a */
void Fp2_print (Fp2_t a);

/* Copy an Fp2 element a to an Fp2 element c. */
void Fp2_copy(Fp2_t c, Fp2_t a);

/* Compare two Fp2 elements a and b for equality. */
int Fp2_cmpeq(Fp2_t a, Fp2_t b);

/* Negate an Fp2 element. */
void Fp2_neg(Fp2_t c, Fp2_t a);

/* Add two Fp2 elements coefficient wise. */
void Fp2_add (Fp2_t c, Fp2_t a, Fp2_t b);

/* Subtract an Fp2 element b from an Fp2 element a coefficient wise. */
void Fp2_sub (Fp2_t c, Fp2_t a, Fp2_t b);

/* Divide an Fp2 element a by 2 coefficient wise. */
void Fp2_div2 (Fp2_t c, Fp2_t a);

/* Multiply an Fp2 element by 3. */
void Fp2_mul3 (Fp2_t c, Fp2_t a);

/* Multiply an Fp2 element by an Fp element. */
void Fp2_mulFp(Fp2_t c, Fp2_t a, Fp_t b);

/* Multiply an Fp2 element by an element i*b, where b is in Fp. */
void Fp2_muliFp(Fp2_t c, Fp2_t a, Fp_t b);

/* Multiply two Fp2 elements. */
void Fp2_mul(Fp2_t c, Fp2_t a, Fp2_t b);

/* Square a Fp2 element. */
void Fp2_squ(Fp2_t c, Fp2_t a);

/* Multiply an Fp2 element by the special element xi. */
void Fp2_mulxi(Fp2_t c, Fp2_t a);

/* Multiply an Fp2 element by the special element i. */
void Fp2_muli(Fp2_t c, Fp2_t a);

/* Invert an Fp2 element. */
void Fp2_inv(Fp2_t c, Fp2_t a);

/* Compute the p-power Frobenius of an Fp2 element. */
void Fp2_ppow(Fp2_t c, Fp2_t a);

/* Multiply two Fp2 elements. This function assumes that p = 3 (mod4).*/
void p3mod4_Fp2_mul(Fp2_t c, Fp2_t a, Fp2_t b, void *mem);

/* Square an Fp2 element. This function assumes that p = 3 (mod4).*/
void p3mod4_Fp2_squ(Fp2_t c, Fp2_t a, void *mem);

/* Multiply an Fp2 element by the special element xi = 1+i. This function assumes that p = 3 (mod4). */
void p3mod4_Fp2_mulxi(Fp2_t c, Fp2_t a, void *mem);

/* Multiply an Fp2 element by the special element i. This function assumes that p = 3 (mod4). */
void p3mod4_Fp2_muli(Fp2_t c, Fp2_t a, void *mem);

/* Multiply an Fp2 element by an element i*b, where b is in Fp. This function assumes that p = 3 (mod4). */
void p3mod4_Fp2_muliFp(Fp2_t c, Fp2_t a, Fp_t b, void *mem);

/* Compute the norm of an Fp2 element. This function assumes that p = 3 (mod4).*/
void p3mod4_Fp2_norm(Fp_t c, Fp2_t a, void *mem);

/* Invert an Fp2 element. This function assumes that p = 3 (mod4).*/
void p3mod4_Fp2_inv(Fp2_t c, Fp2_t a, void *mem);

/* Compute the p-power Frobenius of an Fp2 element. This assumes binomial to construct the extension. */
void binomial_Fp2_ppow(Fp2_t c, Fp2_t a);

/************** Fp2 functions for lazy reduction ***************/

/* Initialize a long Fp2 element. */
int Fp2l_init(Fp2l_t dst);

/* Free the memory allocated in src. */
void Fp2l_free (Fp2l_t src);

/* Copy an Fp2l element a to an Fp2l element c. */
void Fp2l_copy(Fp2l_t C, Fp2l_t A);

/* Separate reduction function calling the Fp reduction for each coefficient. */
void Fp2l_red(Fp2_t c, Fp2l_t A);

/* Add two Fp2 elements coefficient wise without reducing mod p for lazy reduction. */ 
void Fp2_add_no_red (Fp2_t c, Fp2_t a, Fp2_t b);
	
/* Add two long Fp2 elements coefficient wise without reducing mod p for lazy reduction. */
void Fp2l_add_no_red (Fp2l_t C, Fp2l_t A, Fp2l_t B);

/* Subtraction of Fp2l elements, option 1 in Aranha et al. */
void Fp2l_sub_o1 (Fp2l_t c, Fp2l_t a, Fp2l_t b, int h);

/* Subtraction of Fp2l elements, constant-time version of option 2 in Aranha et al. */
void Fp2l_sub_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b);

/* Addition of Fp2l elements, constant-time version of "option 2" addition in Aranha et al. */
void Fp2l_add_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b);

/* Multiply an Fp2 element by an Fp element. */
void Fp2_mulFp_no_red (Fp2l_t c, Fp2_t a, Fp_t b);

/* Multiply two Fp2 elements without reduction. */
void Fp2_mul_no_red_o1(Fp2l_t c, Fp2_t a, Fp2_t b, int h);

/* Multiply two Fp2 elements without reduction. */
void Fp2_mul_no_red_o2(Fp2l_t c, Fp2_t a, Fp2_t b);

/* Multiply in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void Fp2_mul_lazy (Fp2_t c, Fp2_t a, Fp2_t b);

/* Square a Fp2 element without reduction. */
void Fp2_squ_no_red(Fp2l_t c, Fp2_t a);

/* Square in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void Fp2_squ_lazy (Fp2_t c, Fp2_t a);

/* Multiply an Fp2l element by the special element xi. */
void p3mod4_Fp2l_mulxi(Fp2l_t c, Fp2l_t a, void *mem);

/* Multiply an Fp2l element by the special element xi. */
void Fp2l_mulxi(Fp2l_t c, Fp2l_t a);

/* Multiply two Fp2 elements without reduction resulting in a long Fp2l element.
 * This function assumes that p = 3 (mod4).
 */
void p3mod4_Fp2_mul_no_red_o1(Fp2l_t C, Fp2_t a, Fp2_t b, int h, void *mem);
void p3mod4_Fp2_mul_no_red_o2(Fp2l_t C, Fp2_t a, Fp2_t b, void *mem);

/* Squaring an Fp2 element without reduction resulting in a long Fp2l element. */
void p3mod4_Fp2_squ_no_red(Fp2l_t C, Fp2_t a, void *mem);

/***************** Internal stuff *******************/

typedef struct {
	void *mem;
  void (*Fp2_mul)(Fp2_t, Fp2_t, Fp2_t, void *);
  void (*Fp2_squ)(Fp2_t, Fp2_t, void *);
  void (*Fp2_inv)(Fp2_t, Fp2_t, void *);
  void (*Fp2_mulxi)(Fp2_t, Fp2_t, void *);
  void (*Fp2_muli)(Fp2_t, Fp2_t, void *);
  void (*Fp2_muliFp)(Fp2_t, Fp2_t, Fp_t, void *);
  void (*Fp2_norm)(Fp_t, Fp2_t, void *);
  void (*Fp2_ppow)(Fp2_t, Fp2_t);
  void (*Fp2_mul_no_red_o1)(Fp2l_t, Fp2_t, Fp2_t, int, void *);
  void (*Fp2_mul_no_red_o2)(Fp2l_t, Fp2_t, Fp2_t, void *);
  void (*Fp2_squ_no_red)(Fp2l_t, Fp2_t, void *);
  void (*Fp2l_mulxi)(Fp2l_t, Fp2l_t, void *);

   Fp2l_t T0;    // tmp element
} Fp2_config_t;

extern Fp2_config_t Fp2_config; /* Defined in Fp2.c */

#endif /* __FP2_H */
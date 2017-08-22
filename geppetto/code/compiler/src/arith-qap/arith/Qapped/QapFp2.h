#ifndef __FP2_H
#define __FP2_H

#include "Qaparith_error.h"
#include "QapFp.h"

/* An element in F_{p^2}, consists of two coefficients from F_p, 
 * i.e. a = a0 + a1*i in F_{p^2} for i in F_{p^2} with i^2 = -1. 
 */
typedef struct {
  QapFp_t a0;
  QapFp_t a1;
} Qapfp2_t;

/* Our basic type is a pointer to the struct. */
typedef Qapfp2_t QapFp2_t[1];

/* A long element in F_{p^2} of double the size of a regular F_{p^2} element.
 * It consists of two long (double size) coefficients from F_p (type fpl_t), 
 * i.e. A = A0 + A1*i in F_{p^2} for i in F_{p^2} with i^2 = -1.
 * Such elements are used for lazy reduction techniques where reduction is 
 * separated from multiplication and is postponed until after possibly several
 * other operations.
 */
typedef struct {
  QapFpl_t A0;
  QapFpl_t A1;
} fp2l_t;

/* Our basic type is a pointer to the struct. */
typedef fp2l_t Fp2l_t[1];


/* Initialize the global parameters for the quadratic extension
 * given the base field Fp.
 * Return > 0 on success and
 * Return < 0 on error.
 */ 
int QapFp2_initialize_config (void);

void QapFp2_free_config (void);


/* Initialize an element in Fp2 by initializing two Fp elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp2_init(QapFp2_t dst);

/* Free the memory allocated in src. */
void QapFp2_free (QapFp2_t src);

/* Get a random element in Fp2. */
void QapFp2_rand (QapFp2_t c);

/* Set c to the value of a0 + a1*i. */
//void QapFp2_set (QapFp2_t c, uint64_t *a0, uint64_t *a1);
void QapFp2_set (QapFp2_t c, QapFp_t a0, QapFp_t a1);

/* Set c to the value x. */
void QapFp2_set_ui (QapFp2_t c, uint64_t x);

/* Set c to 0. */
void QapFp2_set_zero (QapFp2_t c);

/* Set c to 1. */
void QapFp2_set_one (QapFp2_t c);

/* Print the (regular non-Montgomery form) value of a */
void QapFp2_print (QapFp2_t a);

/* Copy an Fp2 element a to an Fp2 element c. */
void QapFp2_copy(QapFp2_t c, QapFp2_t a);

/* Copy an Fp2 element (bit ? a : b) to an Fp2 element c */
void QapFp2_select(QapFp2_t c, int bit, QapFp2_t a, QapFp2_t b);

/* Compare two Fp2 elements a and b for equality. */
int QapFp2_cmpeq(QapFp2_t a, QapFp2_t b);

/* Negate an Fp2 element. */
void QapFp2_neg(QapFp2_t c, QapFp2_t a);

/* Add two Fp2 elements coefficient wise. */
void QapFp2_add (QapFp2_t c, QapFp2_t a, QapFp2_t b);

/* Subtract an Fp2 element b from an Fp2 element a coefficient wise. */
void QapFp2_sub (QapFp2_t c, QapFp2_t a, QapFp2_t b);

/* Divide an Fp2 element a by 2 coefficient wise. */
void QapFp2_div2 (QapFp2_t c, QapFp2_t a);

/* Multiply an Fp2 element by 3. */
void QapFp2_mul3 (QapFp2_t c, QapFp2_t a);

/* Multiply an Fp2 element by an Fp element. */
void QapFp2_mulFp(QapFp2_t c, QapFp2_t a, QapFp_t b);

/* Multiply an Fp2 element by an element i*b, where b is in Fp. */
void QapFp2_muliFp(QapFp2_t c, QapFp2_t a, QapFp_t b);

/* Multiply two Fp2 elements. */
void QapFp2_mul(QapFp2_t c, QapFp2_t a, QapFp2_t b);

/* Square a Fp2 element. */
void QapFp2_squ(QapFp2_t c, QapFp2_t a);

/* Multiply an Fp2 element by the special element xi. */
void QapFp2_mulxi(QapFp2_t c, QapFp2_t a);

/* Multiply an Fp2 element by the special element i. */
void QapFp2_muli(QapFp2_t c, QapFp2_t a);

/* Invert an Fp2 element. */
void QapFp2_inv(QapFp2_t c, QapFp2_t a);

/* Compute the p-power Frobenius of an Fp2 element. */
void QapFp2_ppow(QapFp2_t c, QapFp2_t a);

/* Multiply two Fp2 elements. This function assumes that p = 3 (mod4).*/
void Qapp3mod4_Fp2_mul(QapFp2_t c, QapFp2_t a, QapFp2_t b, TEMP_FORMAL);

/* Square an Fp2 element. This function assumes that p = 3 (mod4).*/
void Qapp3mod4_Fp2_squ(QapFp2_t c, QapFp2_t a, TEMP_FORMAL);

/* Multiply an Fp2 element by the special element xi = 1+i. This function assumes that p = 3 (mod4). */
void Qapp3mod4_Fp2_mulxi(QapFp2_t c, QapFp2_t a, TEMP_FORMAL);

/* Multiply an Fp2 element by the special element i. This function assumes that p = 3 (mod4). */
void Qapp3mod4_Fp2_muli(QapFp2_t c, QapFp2_t a, TEMP_FORMAL);

/* Multiply an Fp2 element by an element i*b, where b is in Fp. This function assumes that p = 3 (mod4). */
void Qapp3mod4_Fp2_muliFp(QapFp2_t c, QapFp2_t a, QapFp_t b, TEMP_FORMAL);

/* Compute the norm of an Fp2 element. This function assumes that p = 3 (mod4).*/
void Qapp3mod4_Fp2_norm(QapFp_t c, QapFp2_t a, TEMP_FORMAL);

/* Invert an Fp2 element. This function assumes that p = 3 (mod4).*/
void Qapp3mod4_Fp2_inv(QapFp2_t c, QapFp2_t a, TEMP_FORMAL);

/* Compute the p-power Frobenius of an Fp2 element. This assumes binomial to construct the extension. */
void Qapbinomial_Fp2_ppow(QapFp2_t c, QapFp2_t a);

/************** Fp2 functions for lazy reduction ***************/

/* Initialize a long Fp2 element. */
int QapFp2l_init(Fp2l_t dst);

/* Free the memory allocated in src. */
void QapFp2l_free (Fp2l_t src);

/* Copy an Fp2l element a to an Fp2l element c. */
void QapFp2l_copy(Fp2l_t C, Fp2l_t A);

/* Separate reduction function calling the Fp reduction for each coefficient. */
void QapFp2l_red(QapFp2_t c, Fp2l_t A);

/* Add two Fp2 elements coefficient wise without reducing mod p for lazy reduction. */ 
void QapFp2_add_no_red (QapFp2_t c, QapFp2_t a, QapFp2_t b);
	
/* Add two long Fp2 elements coefficient wise without reducing mod p for lazy reduction. */
void QapFp2l_add_no_red (Fp2l_t C, Fp2l_t A, Fp2l_t B);

/* Subtraction of Fp2l elements, option 1 in Aranha et al. */
void QapFp2l_sub_o1 (Fp2l_t c, Fp2l_t a, Fp2l_t b, int h);

/* Subtraction of Fp2l elements, constant-time version of option 2 in Aranha et al. */
void QapFp2l_sub_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b);

/* Addition of Fp2l elements, constant-time version of "option 2" addition in Aranha et al. */
void QapFp2l_add_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b);

/* Multiply an Fp2 element by an Fp element. */
void QapFp2_mulFp_no_red (Fp2l_t c, QapFp2_t a, QapFp_t b);

/* Multiply two Fp2 elements without reduction. */
void QapFp2_mul_no_red_o1(Fp2l_t c, QapFp2_t a, QapFp2_t b, int h);

/* Multiply two Fp2 elements without reduction. */
void QapFp2_mul_no_red_o2(Fp2l_t c, QapFp2_t a, QapFp2_t b);

/* Multiply in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void QapFp2_mul_lazy (QapFp2_t c, QapFp2_t a, QapFp2_t b);

/* Square a Fp2 element without reduction. */
void QapFp2_squ_no_red(Fp2l_t c, QapFp2_t a);

/* Square in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void QapFp2_squ_lazy (QapFp2_t c, QapFp2_t a);

/* Multiply an Fp2l element by the special element xi. */
void Qapp3mod4_Fp2l_mulxi(Fp2l_t c, Fp2l_t a, TEMP_FORMAL);

/* Multiply an Fp2l element by the special element xi. */
void QapFp2l_mulxi(Fp2l_t c, Fp2l_t a);

/* Multiply two Fp2 elements without reduction resulting in a long Fp2l element.
 * This function assumes that p = 3 (mod4).
 */
void Qapp3mod4_Fp2_mul_no_red_o1(Fp2l_t C, QapFp2_t a, QapFp2_t b, int h, TEMP_FORMAL);
void Qapp3mod4_Fp2_mul_no_red_o2(Fp2l_t C, QapFp2_t a, QapFp2_t b, TEMP_FORMAL);

/* Squaring an Fp2 element without reduction resulting in a long Fp2l element. */
void Qapp3mod4_Fp2_squ_no_red(Fp2l_t C, QapFp2_t a, TEMP_FORMAL);

/***************** Internal stuff *******************/

typedef struct {
	Dummy *mem;
   Fp2l_t T0;    // tmp element
} Fp2_config_t;

extern Fp2_config_t Fp2_config; /* Defined in Fp2.c */

#endif /* __FP2_H */

#ifndef __FP4_H
#define __FP4_H

#include "Qaparith_error.h"
#include "QapFp2.h"

/* This contructions realizes F_{p^4} as a quadratic extension over F_{p^2}
 * by an irreducible polynomial X^2 - xi, where xi is in F_{p^2}.
 */

/* An element in F_{p^4} consists of two coefficients from F_{p^2}, 
 * i.e. a = a0 + a1*s in F_{p^4} for s in F_{p^4} with s^2 = xi. 
 */
typedef struct {
  QapFp2_t a0;
  QapFp2_t a1;
} fp4_t;

/* Our basic type is a pointer to the struct. */
typedef fp4_t QapFp4_t[1];

int QapFp4_initialize_config (void);

void QapFp4_free_config (void);


/* Initialize an element in Fp4 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp4_init(QapFp4_t dst);

/* Free the memory allocated in src. */
void QapFp4_free (QapFp4_t src);

/* Get a random element in Fp4. */
void QapFp4_rand (QapFp4_t c);

/* Set c to 0. */
void QapFp4_set_zero (QapFp4_t c);

/* Set c to 1. */
void QapFp4_set_one (QapFp4_t c);

/* Copy an Fp4 element a to an Fp4 element c. */
void QapFp4_copy(QapFp4_t c, QapFp4_t a);

/* Print the (regular non-Montgomery form) value of a */
void QapFp4_print (QapFp4_t a);

/* Compare two Fp4 elements a and b for equality. */
int QapFp4_cmpeq(QapFp4_t a, QapFp4_t b);

/* Negate an Fp4 element. */
void QapFp4_neg(QapFp4_t c, QapFp4_t a);

/* Add two Fp4 elements coefficient wise. */
void QapFp4_add (QapFp4_t c, QapFp4_t a, QapFp4_t b);

/* Subtract an Fp4 element b from an Fp4 element a coefficient wise. */
void QapFp4_sub (QapFp4_t c, QapFp4_t a, QapFp4_t b);

/* Divide an Fp4 element a by 2 coefficient wise. */
void QapFp4_div2 (QapFp4_t c, QapFp4_t a);

/* Multiply an Fp4 element by an Fp2 element. */
void QapFp4_mulFp2(QapFp4_t c, QapFp4_t a, QapFp2_t b);

/* Multiply an Fp4 element by an Fp element. */
void QapFp4_mulFp(QapFp4_t c, QapFp4_t a, QapFp_t b);

/* Multiply two Fp4 elements. */
void QapFp4_mul (QapFp4_t c, QapFp4_t a, QapFp4_t b);

/* Multiply two Fp4 elements using lazy reduction. */
void QapFp4_mul_lazy (QapFp4_t c, QapFp4_t a, QapFp4_t b);

/* Multiply compute 3*a. */
void QapFp4_mul3 (QapFp4_t c, QapFp4_t a);

/* Square an Fp4 element. */
void QapFp4_squ (QapFp4_t c, QapFp4_t a);

/* Square an Fp4 element using lazy reduction. */
void QapFp4_squ_lazy (QapFp4_t c, QapFp4_t a);

/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void QapFp4_squ_sep (QapFp2_t c0, QapFp2_t c1, QapFp2_t a0, QapFp2_t a1);

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2. */
void QapFp4_muls (QapFp4_t c, QapFp4_t a);

/* Multiply an Fp4 element by the special element i in Fp2.*/
void QapFp4_muli(QapFp4_t c, QapFp4_t a);

/* Multiply an Fp4 element by the special element xi in Fp2.*/
void QapFp4_mulxi(QapFp4_t c, QapFp4_t a);

/* Compute the p^2-power Frobenius of an Fp4 element. */
void QapFp4_p2pow(QapFp4_t c, QapFp4_t a);

/* Compute the p-power Frobenius of an Fp4 element. */
void QapFp4_ppow(QapFp4_t c, QapFp4_t a);

/* Invert an Fp4 element. */
void QapFp4_inv (QapFp4_t c, QapFp4_t a);

/* Multiply two Fp4 elements. */
void Qaps2minusxi_Fp4_mul(QapFp4_t c, QapFp4_t a, QapFp4_t b, TEMP_FORMAL);

/* Multiply two Fp4 elements using lazy reduction. */
void Qaps2minusxi_Fp4_mul_lazy(QapFp4_t c, QapFp4_t a, QapFp4_t b, TEMP_FORMAL);

/* Square an Fp4 element. */
void Qaps2minusxi_Fp4_squ(QapFp4_t c, QapFp4_t a, TEMP_FORMAL);

/* Square an Fp4 element using lazy reduction. */
void Qaps2minusxi_Fp4_squ_lazy(QapFp4_t c, QapFp4_t a, TEMP_FORMAL);

/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void Qaps2minusxi_Fp4_squ_sep(QapFp2_t c0, QapFp2_t c1, QapFp2_t a0, QapFp2_t a1, TEMP_FORMAL);

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2. */
void Qaps2minusxi_Fp4_muls(QapFp4_t c, QapFp4_t a, TEMP_FORMAL);

/* Inversion of an Fp4 element. This uses that Fp4 is constructed over Fp2 via s^2 - xi. */
void Qaps2minusxi_Fp4_inv(QapFp4_t c, QapFp4_t a, TEMP_FORMAL);

/* Compute the p^2-power Frobenius of an Fp4 element. */
void Qapbinomial_Fp4_p2pow(QapFp4_t c, QapFp4_t a);

/* Compute the p-power Frobenius of an Fp4 element. */
void Qapbinomial_Fp4_ppow(QapFp4_t c, QapFp4_t a);

/* Multiply an Fp4 element by the special element i in Fp2. */
void Qapp3mod4_Fp4_muli(QapFp4_t c, QapFp4_t a);

/* Multiply an Fp4 element by the special element xi in Fp2. */
void Qapp3mod4_Fp4_mulxi(QapFp4_t c, QapFp4_t a);

/***************** Internal stuff *******************/

typedef struct {
  QapFp2_t sppow;
} Fp4_config_t;

extern Fp4_config_t Fp4_config; /* Defined in Fp4.c */

#endif /* __FP4_H */

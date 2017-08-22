#ifndef __FP12_H
#define __FP12_H

#include "Qaparith_error.h"
#include "QapFp6.h"


/* This contructions realizes F_{p^12} as a quadratic extension over F_{p^6}
 * by an irreducible polynomial X^2 - v, where v is in F_{p^6}.
 */

/* An element in F_{p^12} consists of two coefficients from F_{p^6}, 
 * i.e. a = a0 + a1*w in F_{p^12} for w in F_{p^6} with w^2 = v. 
 */
typedef struct {
  QapFp6_t a0;
  QapFp6_t a1;
} fp12_t;

/* Our basic type is a pointer to the struct. */
typedef fp12_t QapFp12_t[1];


int QapFp12_initialize_config (void);

void QapFp12_free_config (void);


/************** Regular Fp12 functions ***************/

/* Initialize an element in Fp12 by initializing two Fp6 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp12_init(QapFp12_t dst);

/* Free the memory allocated in src. */
void QapFp12_free (QapFp12_t src);

/* Get a random element in Fp12. */
void QapFp12_rand (QapFp12_t c);

/* Print the (regular non-Montgomery form) value of a */
void QapFp12_print (QapFp12_t a);

/* Set c to 0. */
void QapFp12_set_zero (QapFp12_t c);

/* Set c to 1. */
void QapFp12_set_one (QapFp12_t c);

/* Copy an Fp12 element a to an Fp12 element c. */
void QapFp12_copy(QapFp12_t c, QapFp12_t a);

/* Compare two Fp12 elements a and b for equality. */
int QapFp12_cmpeq(QapFp12_t a, QapFp12_t b);

/* Negate an Fp12 element. */
void QapFp12_neg(QapFp12_t c, QapFp12_t a);

/* Add two Fp12 elements coefficient wise. */
void QapFp12_add (QapFp12_t c, QapFp12_t a, QapFp12_t b);

/* Subtract an Fp12 element b from an Fp12 element a coefficient wise. */
void QapFp12_sub (QapFp12_t c, QapFp12_t a, QapFp12_t b);

/* Multiply two Fp12 elements.*/
void QapFp12_mul(QapFp12_t c, QapFp12_t a, QapFp12_t b);

/* Square a general Fp12 element. */
void QapFp12_squ(QapFp12_t c, QapFp12_t a);

/* Multiply an Fp12 element by an Fp6 element. */
void QapFp12_mulFp6(QapFp12_t c, QapFp12_t a, QapFp6_t b);

/* Multiply an Fp12 element by an Fp2 element. */
void QapFp12_mulFp2(QapFp12_t c, QapFp12_t a, QapFp2_t b);

/* Multiply an Fp12 element by an Fp1 element. */
void QapFp12_mulFp(QapFp12_t c, QapFp12_t a, QapFp_t b);

/* Compute the p-power Frobenius of an Fp12 element. */
void QapFp12_ppow(QapFp12_t c, QapFp12_t a);

/* Compute the p^2-power Frobenius of an Fp12 element. */
void QapFp12_p2pow(QapFp12_t c, QapFp12_t a);

/* Compute the p^3-power Frobenius of an Fp12 element. */
void QapFp12_p3pow(QapFp12_t c, QapFp12_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void QapFp12_p6pow(QapFp12_t c, QapFp12_t a);

/* Invert an Fp12 element. */
void QapFp12_inv(QapFp12_t c, QapFp12_t a);

/* Multiply two Fp12 elements using lazy reduction.*/
void QapFp12_mul_lazy(QapFp12_t c, QapFp12_t a, QapFp12_t b);

/* Square an Fp12 element using lazy reduction.*/
void QapFp12_squ_lazy(QapFp12_t c, QapFp12_t a);

/* Multiply two Fp12 elements (the second one sparse) using lazy reduction.*/
void QapFp12_mul_sparse_twist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b00, QapFp2_t b01, QapFp2_t b11);

void QapFp12_mul_sparse_untwist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b01, QapFp2_t b11, QapFp2_t b12);

/* Multiply two Fp12 elements.*/
void Qapw2minusv_Fp12_mul(QapFp12_t c, QapFp12_t a, QapFp12_t b, TEMP_FORMAL);

/* Square a general Fp12 element. */
void Qapw2minusv_Fp12_squ(QapFp12_t c, QapFp12_t a, TEMP_FORMAL);

/* Compute the p-power Frobenius of an Fp12 element. */
void Qapw2minusv_v3minusxi_p1mod6_Fp12_ppow(QapFp12_t c, QapFp12_t a);

/* Compute the p^2-power Frobenius of an Fp12 element. */
void Qapw2minusv_v3minusxi_p1mod6_Fp12_p2pow(QapFp12_t c, QapFp12_t a);

/* Compute the p^3-power Frobenius of an Fp12 element. */
void Qapw2minusv_v3minusxi_p1mod6_Fp12_p3pow(QapFp12_t c, QapFp12_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void Qapw2minusv_Fp12_p6pow(QapFp12_t c, QapFp12_t a);

/* Invert an Fp12 element. */
void Qapw2minusv_Fp12_inv(QapFp12_t c, QapFp12_t a, TEMP_FORMAL);

/* Multiply two Fp12 elements using lazy reduction. */
void Qapw2minusv_Fp12_mul_lazy(QapFp12_t c, QapFp12_t a, QapFp12_t b, TEMP_FORMAL);

/* Square an Fp12 element using lazy reduction. */
void Qapw2minusv_Fp12_squ_lazy(QapFp12_t c, QapFp12_t a, TEMP_FORMAL);

/* Multiply an Fp12 element a and a sparse Fp12 element b using lazy reduction. */
void Qapw2minusv_Fp12_mul_sparse_twist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b00, QapFp2_t b01, QapFp2_t b11, TEMP_FORMAL);

/* Multiply an Fp12 element a and a sparse Fp12 element b using lazy reduction. */
void Qapw2minusv_Fp12_mul_sparse_untwist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b01, QapFp2_t b11, QapFp2_t b12, TEMP_FORMAL);

/***************** Internal stuff *******************/

typedef struct {
  QapFp2_t wppow; // w^(p-1)
  QapFp_t wp2pow; // w^(p^2-1)
  QapFp_t wp2pow3; // (w^3)^(p^2-1)
  QapFp2_t wp3pow;// w^(p^3-1)
  QapFp2_t wppow3;// (w^3)^(p-1)
  QapFp2_t wppow3inv; // wpppow3^(-1)
} Fp12_config_t;

extern Fp12_config_t Fp12_config; /* Defined in Fp12.c */

#endif /* __FP12_H */

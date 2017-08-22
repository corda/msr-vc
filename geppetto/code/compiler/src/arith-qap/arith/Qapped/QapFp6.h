#ifndef __FP6_H
#define __FP6_H

#include "Qaparith_error.h"
#include "QapFp2.h"

/* This contructions realizes F_{p^6} as a cubic extension over F_{p^2}
 * by an irreducible polynomial X^3 - xi, where xi is in F_{p^2}.
 */

/* An element in F_{p^6} consists of three coefficients from F_{p^2}, 
 * i.e. a = a0 + a1*v + a2*v^2 in F_{p^6} for v in F_{p^6} with v^3 = xi. 
 */
typedef struct {
  QapFp2_t a0;
  QapFp2_t a1;
  QapFp2_t a2;
} fp6_t;

/* Our basic type is a pointer to the struct. */
typedef fp6_t QapFp6_t[1];

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

int QapFp6_initialize_config (void);

void QapFp6_free_config (void);


/************** Regular Fp6 functions ***************/

/* Initialize an element in Fp6 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp6_init(QapFp6_t dst);

/* Free the memory allocated in src. */
void QapFp6_free (QapFp6_t src);

/* Get a random element in Fp6. */
void QapFp6_rand (QapFp6_t c);

/* Set c to 0. */
void QapFp6_set_zero (QapFp6_t c);

/* Set c to 1. */
void QapFp6_set_one (QapFp6_t c);

/* Copy an Fp6 element a to an Fp6 element c. */
void QapFp6_copy(QapFp6_t c, QapFp6_t a);

/* Copy an Fp6 element (bit ? a : b) to an Fp6 element c */
void QapFp6_select(QapFp6_t c, int bit, QapFp6_t a, QapFp6_t b);

/* Print the (regular non-Montgomery form) value of a */
void QapFp6_print (QapFp6_t a);

/* Compare two Fp6 elements a and b for equality. */
int QapFp6_cmpeq(QapFp6_t a, QapFp6_t b);

/* Negate an Fp6 element. */
void QapFp6_neg(QapFp6_t c, QapFp6_t a);

/* Add two Fp6 elements coefficient wise. */
void QapFp6_add (QapFp6_t c, QapFp6_t a, QapFp6_t b);

/* Subtract an Fp6 element b from an Fp6 element a coefficient wise. */
void QapFp6_sub (QapFp6_t c, QapFp6_t a, QapFp6_t b);

/* Multiply two Fp6 elements. */
void QapFp6_mul (QapFp6_t c, QapFp6_t a, QapFp6_t b);

/* Square an Fp6 element. */
void QapFp6_squ (QapFp6_t c, QapFp6_t a);

/* Multiply an Fp6 element by an Fp element. */
void QapFp6_mulFp (QapFp6_t c, QapFp6_t a, QapFp_t b);

/* Multiply an Fp6 element by an Fp2 element. */
void QapFp6_mulFp2 (QapFp6_t c, QapFp6_t a, QapFp2_t b);

/* Invert an Fp6 element. */
void QapFp6_inv (QapFp6_t c, QapFp6_t a);

/* Multiply an fp6 element by the special element v. */
void QapFp6_mulv(QapFp6_t c, QapFp6_t a);

/* Compute the p-power Frobenius of an Fp6 element. */
void QapFp6_ppow(QapFp6_t c, QapFp6_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void QapFp6_p2pow(QapFp6_t c, QapFp6_t a);

/* Compute the p^3-power Frobenius of an Fp6 element. */
void QapFp6_p3pow(QapFp6_t c, QapFp6_t a);

/* Multiply two Fp6 elements. */
void Qapv3minusxi_Fp6_mul(QapFp6_t c, QapFp6_t a, QapFp6_t b, TEMP_FORMAL);

/* Square an Fp6 element. */
void Qapv3minusxi_Fp6_squ(QapFp6_t c, QapFp6_t a, TEMP_FORMAL);

/* Multiply an fp6 element by the special element v. */
void Qapv3minusxi_Fp6_mulv(QapFp6_t c, QapFp6_t a, TEMP_FORMAL);

/* Compute the p-power Frobenius of an Fp6 element. */
void Qapp1mod3_v3minusxi_Fp6_ppow(QapFp6_t c, QapFp6_t a);

/* Compute the p^2-power Frobenius of an Fp6 element. */
void Qapp1mod3_v3minusxi_Fp6_p2pow(QapFp6_t c, QapFp6_t a);

/* Compute the p^3-power Frobenius of an Fp6 element. */
void Qapp1mod3_v3minusxi_Fp6_p3pow(QapFp6_t c, QapFp6_t a);

/* Inversion of an Fp6 element. */
void Qapv3minusxi_Fp6_inv(QapFp6_t c, QapFp6_t a, TEMP_FORMAL);


/************** Fp6 functions for lazy reduction ***************/

/* Initialize a long Fp6 element. */
int QapFp6l_init(Fp6l_t dst);

/* Free the memory allocated in src. */
void QapFp6l_free (Fp6l_t src);

/* Copy an Fp6l element a to an Fp6l element c. */
void QapFp6l_copy(Fp6l_t C, Fp6l_t A);

/* Separate reduction function calling the Fp2 reduction for each coefficient. */
void QapFp6l_red(QapFp6_t c, Fp6l_t A);

/* Subtraction of Fp6l elements, constant-time version of option 2 in Aranha et al. */
void QapFp6l_sub_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b);

/* Addition of Fp6l elements, constant-time version of "option 2" addition in Aranha et al. */
void QapFp6l_add_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b);

/* Multiply an Fp6 element by an Fp element. */
void QapFp6_mulFp2_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b);

/* Multiply two Fp6 elements without reduction. */
void QapFp6_mul_no_red (Fp6l_t c, QapFp6_t a, QapFp6_t b);

/* Multiply an Fp6 element and a sparse Fp6 element without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void QapFp6_mul_sparse01_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b0, QapFp2_t b1);
void QapFp6_mul_sparse12_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b0, QapFp2_t b1);

/* Multiply two Fp6 elements without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void Qapv3minusxi_Fp6_mul_no_red (Fp6l_t c, QapFp6_t a, QapFp6_t b, TEMP_FORMAL);

/* Multiply an Fp6 element and a sparse Fp6 element without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void Qapv3minusxi_Fp6_mul_sparse01_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b0, QapFp2_t b1, TEMP_FORMAL);
void Qapv3minusxi_Fp6_mul_sparse12_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b1, QapFp2_t b2, TEMP_FORMAL);


/***************** Internal stuff *******************/

typedef struct {
  QapFp_t vppow;   // v^(p-1)/i  
  QapFp_t vppow2;  // (i*vppow)^2 = -vppow^2
  QapFp_t vp2pow;  // v^(p^2-1)
  QapFp_t vp2pow2; // vp2pow^2
  QapFp_t vppowinv; // -i/v^(p-1) = 1/v^(p-1)/i
} Fp6_config_t;

extern Fp6_config_t Fp6_config; /* Defined in Fp6.c */

#endif /* __FP6_H */

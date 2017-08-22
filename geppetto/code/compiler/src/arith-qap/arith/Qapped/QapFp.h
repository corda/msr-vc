// Pinocchio-specialized version, to avoid dynamic allocation

#ifndef __QAPFP_H
#define __QAPFP_H

//#include "Qaparith.h"
//#include "mont_arith.h"
#include "Qapuint.h"	// Print
#include "QapTempMem.h"

#pragma intrinsic(__rdtsc)

#include "PrimitiveIfc.h"

/* An element in F_p consists simply of an element of
  the aggregator's Pinocchio-Q */
typedef struct {
	Elem e;
} Qapfp_t;

/* Our basic type is a pointer to the structure. */
typedef Qapfp_t QapFp_t[1];
typedef QapFp_t QapFpl_t;

/* Initialize the global parameters for the modulus
 * given the modulus mod consisting of n 64-bit words.
 */ 
void QapFp_initialize_config ();

void QapFp_free_config (void);

/* Initialize an element in Fp given the parameters passed
 * to Fp_modulus_init.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp_init(QapFp_t dst);
int QapFpl_init(QapFpl_t dst);

/* Free the memory allocated in src,
 * set the number of limbs to zero and point to NULL. 
 */
void QapFp_free (QapFp_t src);
void QapFpl_free (QapFpl_t src);

/* Compute the modular inverse
 * dst = src^-1 mod p
 */
void QapFp_modinv (QapFp_t dst, QapFp_t src);

/* Compute c = a * b mod p */
void QapFp_mul (QapFp_t c, QapFp_t a, QapFp_t b);

/* Print the (regular non-Montgomery form) value of a */
void QapFp_print (QapFp_t a);

/* Set c to the value of a */
void QapFp_set(QapFp_t c, QapFp_t a);

void QapFp_set_ui (QapFp_t c, uint64_t x);

typedef uint32_t pin_int;

/* Set c to a literal value */
void QapFp_set_u32 (QapFp_t c,
	pin_int a7, pin_int a6, pin_int a5, pin_int a4,
	pin_int a3, pin_int a2, pin_int a1, pin_int a0);

void QapFp_set_literal32 (QapFp_t c, pin_int lit);

/* c = a - b mod p */
void QapFp_sub (QapFp_t c, QapFp_t a, QapFp_t b);

/* TODO FIX We assume here that we can add without overflow. */
void QapFp_div2 (QapFp_t c, QapFp_t a);

/* Addition without reduction. */
void QapFp_sub_no_red (QapFp_t c, QapFp_t a, QapFp_t b);

/* Addition without reduction. No check for overflow! */
void QapFpl_sub_no_red (QapFpl_t c, QapFpl_t a, QapFpl_t b);

/* Subtraction of Fpl elements, option 1 with h=1 or h=2 in Aranha et al. */
void QapFpl_sub_o1_ct (QapFpl_t c, QapFpl_t a, QapFpl_t b, int h, void *mem);

/* Subtraction of Fpl elements, option 2 in Aranha et al. */
void QapFpl_sub_o2_ct (QapFpl_t c, QapFpl_t a, QapFpl_t b, void *mem);

/* Addition of Fpl elements, constant-time version of "option 2" addition in Aranha et al. */
void QapFpl_add_o2_ct (QapFpl_t c, QapFpl_t a, QapFpl_t b, void *mem);

/* c = a + b mod p */
void QapFp_add (QapFp_t c, QapFp_t a, QapFp_t b);

void QapFp_select(QapFp_t c, int bit, QapFp_t a, QapFp_t b);

/* Addition without reduction. */
void QapFp_add_no_red (QapFp_t c, QapFp_t a, QapFp_t b);

/* Addition without reduction on long elements. */
void QapFpl_add_no_red (QapFpl_t c, QapFpl_t a, QapFpl_t b);

/* return true if a==b */
int QapFp_cmpeq (QapFp_t a, QapFp_t b);

/* Return +1 if a > b
 *         0 if a = b
 *        -1 if a < b
 */
int QapFp_cmp (QapFp_t a, QapFp_t b);

/* Just the multiplication part, no reduction.
 * c= a * b
 */
void QapFp_mul_no_red (QapFpl_t c, QapFp_t a, QapFp_t b);

/* Reduce a double-length number a
 * c = a mod p
 */
void QapFpl_red (QapFp_t c, QapFpl_t a);

void QapFp_copy (QapFp_t c, QapFp_t a);

void QapFpl_copy (QapFpl_t c, QapFpl_t a);

void QapFp_neg (QapFp_t c, QapFp_t a);

void QapFp_rand (QapFp_t c);

int QapFp_get_bit(QapFp_t a, int bit);

void QapFp_mul3 (QapFp_t c, QapFp_t a);

//void Fp_mul_2exp (QapFp_t c, QapFp_t a, int n);

#define QapFp_sqr(c,a) QapFp_mul (c, a, a)
//void QapFp_sqr (QapFp_t c, QapFp_t a);

/***************** Internal stuff *******************/

/* Our configuration structure. */
typedef struct {
	QapFp_t m;	// modulus
//  montmul_t m;
  int n; // limbs
} QapFp_config_t;

extern QapFp_config_t QapFp_config; /* Defined in Fp.c */

#define FP_CONFIG_N()	(QapFp_config.n)

Dummy* g_dummy;

/* Macros redirecting the memory allocating functions. */
#define QapFpl_free(x) QapFp_free (x)

/* Internal function.
 * Returns a^-1 mod 2^64
 */
uint64_t modinv64 (uint64_t a);

//void copy_limbs(m64_t* out, void* in, int fp_offset);

void QapDbgAssert(int condition);

#endif /* __QAPFP_H */

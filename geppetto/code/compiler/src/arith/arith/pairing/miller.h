/* miller.h
 * by Joppe W. Bos, Michael Naehrig, Cryptography Research Group, XCG, Microsoft Research 2013
 *
 * Header file for computing the Miller loop part of the pairing algorithm.
 * This file is part of the Pairing library.
 * We assume the target platform has access to 64-bit data-types.
 */

#ifndef __MILLER_H
#define __MILLER_H

#include "../arith_error.h"
#include "../Fp12/Fp12.h"
#include "../Curve/Fp_weierstrass.h"
#include "../Curve/Fp2_weierstrass.h"

int bn_miller_initialize_config (Fp2_t bt);

void bn_miller_free_config (void);

/* Line function for a doubling step on BN curves. */
void bn_line_dbl (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Fp2_t b3);

/* Line function for an addition step on BN curves. */
void bn_line_add (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Point_wp2_t P);

/* Compute the twisted point of untwist(Q)^p. */
void bn_Fp2_wp_ppow_untwist (Point_wp2_t Qp, Point_wp2_t Q);

/* Compute the twisted point of untwist(Q)^(p^2). */
void bn_Fp2_wp_p2pow_untwist (Point_wp2_t Qp, Point_wp2_t Q);

/* The Miller loop part of the optimal ate pairing on BN curves excluding the final exponentiation. */
void bn_optimal_ate_miller (Fp12_t f, Point_wp2_t Q, Point_wp_t P);

/***************** Internal stuff *******************/

typedef struct {
  int len6u2;
  int *naf6u2;
  int signu;
  Fp2_t bt3;

} bn_miller_config_t;

extern bn_miller_config_t bn_miller_config; /* Defined in miller.c */

typedef struct {
  int lenu;
  int *nafu;
  int signu;
  Fp2_t bt3;

} bls12_miller_config_t;

extern bls12_miller_config_t bls12_miller_config; /* Defined in miller.c */

#endif /* __MILLER_H */
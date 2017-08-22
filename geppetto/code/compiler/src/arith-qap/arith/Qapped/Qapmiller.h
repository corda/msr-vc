#ifndef __MILLER_H
#define __MILLER_H

#include "Qaparith_error.h"
#include "QapFp12.h"
#include "QapFp_weierstrass.h"
#include "QapFp2_weierstrass.h"

int bn_miller_initialize_config (QapFp2_t bt);

void bn_miller_free_config (void);

/* Line function for a doubling step on BN curves. */
void bn_line_dbl (QapFp2_t l00, QapFp2_t l01, QapFp2_t l10, QapPoint_wp2_t R, QapFp2_t b3);

/* Line function for an addition step on BN curves. */
void bn_line_add (QapFp2_t l00, QapFp2_t l01, QapFp2_t l10, QapPoint_wp2_t R, QapPoint_wp2_t P);

/* Compute the twisted point of untwist(Q)^p. */
void bn_Fp2_wp_ppow_untwist (QapPoint_wp2_t Qp, QapPoint_wp2_t Q);

/* Compute the twisted point of untwist(Q)^(p^2). */
void bn_Fp2_wp_p2pow_untwist (QapPoint_wp2_t Qp, QapPoint_wp2_t Q);

/* The Miller loop part of the optimal ate pairing on BN curves excluding the final exponentiation. */
void Qapbn_optimal_ate_miller (QapFp12_t f, QapPoint_wp2_t Q, QapPoint_wp_t P);


#endif /* __MILLER_H */

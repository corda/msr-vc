#ifndef __pairing_H
#define __pairing_H

#include "Qaparith_error.h"
#include "QapFp12.h"
#include "QapFp_weierstrass.h"
#include "QapFp2_weierstrass.h"

#define BN12        1     /* 254-bit BN curve, k=12  */
#define BLS12       2     /* 635-bit BLS curve, k=12 */
#define BN12CP      5     /* 254-bit BN curve, k=12, base field characteristic p divides group order of CP6 below. */
#define CP6         6     /* 507-bit Cocks-Pinch curve, k=6, group order over base field is divisible by base field characteristic of BN12CP. */
#define CP6b        8     /* 509-bit Cocks-Pinch curve, k=6, group order over base field is divisible by base field characteristic of BN12. */
#define BN12tiny    9     /* 53-bit toy example BN curve, k=12 */
#define CP6tiny    10     /* 107-bit toy example CP curve, k=6, group order divisible by BN12tiny base field prime*/

#ifndef TINY
#define PAIR_CURVE  BN12   /* Sets the curve to be used. */
#else 
#define PAIR_CURVE  BN12tiny
#endif // TINY

int Qapbn_bls12_miller_initialize_config (QapFp2_t bt);
void Qapbn_bls12_miller_free_config (void);

int Qapbn_bls12_finalexpo_initialize_config (void);
void Qapbn_bls12_finalexpo_free_config (void);

/* Line function for a doubling step on BN and BLS12 curves. */
void Qapbn_bls12_line_affdbl(QapFp2_t lambda, QapFp2_t nu, QapPoint_wp2_t R);
/* Line function for an addition step on BN and BLS12 curves. */
void Qapbn_bls12_line_affadd(QapFp2_t lambda, QapFp2_t nu, QapPoint_wp2_t R, QapPoint_wp2_t P);

/* The Miller loop part of the optimal ate pairing on BN curves excluding the final exponentiation. */
void Qapbn_optimal_ate_miller (QapFp12_t f, QapPoint_wp2_t Q, QapPoint_wp_t P);
void Qapbn_optimal_ate_miller_aff(QapFp12_t f, QapPoint_wp2_t Q, QapPoint_wp_t P);

/* The final exponentiation for BN curves with a negative parameter. */
void Qapbn_finalexpo_neg (QapFp12_t c, QapFp12_t a);

/* The final exponentiation for BN curves with a positive parameter. */
void Qapbn_finalexpo_pos (QapFp12_t c, QapFp12_t a);

/***************** Internal stuff *******************/

typedef struct {
  int lenu;
  int nafu[63];
  int len6u2;
  int naf6u2[66];
  int signu;
  QapFp2_t bt3;

} bn_config_t;

extern bn_config_t bn_config; /* Defined in miller.c */

#endif /* __PAIRING_H */

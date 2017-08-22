/* pairing/pairing.h
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

#ifndef __pairing_H
#define __pairing_H

#include "../arith_error.h"
#include "../Fp6/Fp6overFp3.h"
#include "../Fp12/Fp12.h"
#include "../Fp12/Fp12overFp4.h"
#include "../Curve/Fp_weierstrass.h"
#include "../Curve/Fp2_weierstrass.h"


#define BN12        1     /* 254-bit BN curve, k=12  */
#define BN12CP      5     /* 254-bit BN curve, k=12, base field characteristic p divides group order of CP6 below. */
#define CP6         6     /* 507-bit Cocks-Pinch curve, k=6, group order over base field is divisible by base field characteristic of BN12CP. */
#define CP6b        8     /* 509-bit Cocks-Pinch curve, k=6, group order over base field is divisible by base field characteristic of BN12. */
#define BN12tiny    9     /* 53-bit toy example BN curve, k=12 */
#define CP6tiny    10     /* 107-bit toy example CP curve, k=6, group order divisible by BN12tiny base field prime*/
#define CP3        11     /* 1023-bit Cocks-Pinch curve, k=3, group order over base field is divisible by base field characteristic of CP6b. */

extern int PAIR_CURVE;

int bn_bls12_miller_initialize_config (Fp2_t bt);
void bn_bls12_miller_free_config (void);

int bn_bls12_finalexpo_initialize_config (void);
void bn_bls12_finalexpo_free_config (void);

int cp6_miller_initialize_config (Fp_t bt);
void cp6_miller_free_config (void);

int cp6_finalexpo_initialize_config (void);
void cp6_finalexpo_free_config (void);

int cp3_miller_initialize_config(Fp_t bt);
void cp3_miller_free_config(void);

int cp3_finalexpo_initialize_config(void);
void cp3_finalexpo_free_config(void);

/* Line function for a doubling step on BN and BLS12 curves. */
void bn_bls12_line_dbl (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Fp2_t b3);
void bn_bls12_line_affdbl(Fp2_t lambda, Fp2_t nu, Point_wp2_t R);
/* Line function for an addition step on BN and BLS12 curves. */
void bn_bls12_line_add (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Point_wp2_t P);
void bn_bls12_line_affadd(Fp2_t lambda, Fp2_t nu, Point_wp2_t R, Point_wp2_t P);

/* Line function for a doubling step on CP6 curve. */
void cp6_line_dbl (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Fp_t b3);
void cp6_line_affdbl(Fp_t lambda, Fp_t nu, Point_wp_t R);
/* Line function for an addition step on CP6 curve. */
void cp6_line_add (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Point_wp_t P);
void cp6_line_add_proj (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Point_wp_t P);
void cp6_line_affadd(Fp_t lambda, Fp_t nu, Point_wp_t R, Point_wp_t P);

/* Line function for a doubling step on CP3 curve. */
void cp3_line_dbl(Fp_t l00, Fp_t l01, Fp_t l10, Fp_t l20, Point_wp_t R, Fp_t b);
/* Line function for an addition step on CP3 curve. */
void cp3_line_add(Fp_t l00, Fp_t l01, Fp_t l10, Fp_t l20, Point_wp_t R, Point_wp_t Q);

/* The Miller loop part of the optimal ate pairing on BN curves excluding the final exponentiation. */
void bn_optimal_ate_miller (Fp12_t f, Point_wp2_t Q, Point_wp_t P);
void bn_optimal_ate_miller_aff(Fp12_t f, Point_wp2_t Q, Point_wp_t P);
void bn_optimal_ate_miller_precompute_aff(Point_wp2_t Q, uint64_t *premem);
void bn_optimal_ate_miller_useprecomputed_aff(Fp12_t f, Point_wp_t P, uint64_t *premem);

/* The Miller loop part of the optimal ate pairing on the CP6 curve excluding the final exponentiation. */
void cp6_optimal_ate_miller (Fp6q_t f, Point_wp_t Q, Point_wp_t P);
/* The Miller loop part of the optimal ate pairing on the CP6b curve excluding the final exponentiation. */
void cp6b_optimal_ate_miller (Fp6q_t f, Point_wp_t Q, Point_wp_t P);
void cp6b_optimal_ate_miller_aff(Fp6q_t f, Point_wp_t Q, Point_wp_t P);
void cp6b_optimal_ate_miller_precompute(Point_wp_t Q, uint64_t *premem);
void cp6b_optimal_ate_miller_useprecomputed(Fp6q_t f, Point_wp_t P, uint64_t *premem);

/* The Miller loop part of the ate pairing on the CP3 curve excluding the final exponentiation. */
void cp3_ate_miller(Fp3_t f, Point_wp_t Q, Point_wp_t P);

/* The final exponentiation for BN curves with a negative parameter. */
void bn_finalexpo_neg (Fp12_t c, Fp12_t a);

/* The final exponentiation for BN curves with a positive parameter. */
void bn_finalexpo_pos (Fp12_t c, Fp12_t a);

/* The final exponentiation for the CP6 curve with a positive parameter. */
void cp6_finalexpo (Fp6q_t c, Fp6q_t a);
/* The final exponentiation for the CP6b curve with a positive parameter. */
void cp6b_finalexpo (Fp6q_t c, Fp6q_t a);

/* The final exponentiation for the CP3 curve with a positive parameter. */
void cp3_finalexpo(Fp3_t c, Fp3_t a);

void Fp12_squ_cyclotomic (Fp12_t c, Fp12_t a);
void Fp6q_squ_cyclotomic (Fp6q_t c, Fp6q_t a);


/***************** Internal stuff *******************/

typedef struct {
  int lenu;
  int *nafu;  
  int len6u2;
  int *naf6u2;
  int nafwt6u2;
  int signu;
  Fp2_t bt3;

} bn_config_t;

extern bn_config_t bn_config; /* Defined in miller.c */

typedef struct {
	int lenm;
	int len3m;
	int *nafm;
	int *naf3m;
	int nafwtm;
	int nafwt3m;
	Fp_t bt3;

} cp6_config_t;

extern cp6_config_t cp6_config; /* Defined in miller.c */

typedef struct {
	int lenT;
	int lena0;
	int lena1;
	int lena2;
	int naflena2;
	int *nafT;
	int *nafa2;
	int nafwtT;
	uint64_t *T;
	uint64_t *a0;
	uint64_t *a1;
	uint64_t *a2;
	Fp_t c1;
	Fp_t s11;
	Fp_t bt;

} cp3_config_t;

extern cp3_config_t cp3_config; /* Defined in miller.c */

extern cp6_config_t cp6_config; /* Defined in miller.c */

#endif /* __PAIRING_H */

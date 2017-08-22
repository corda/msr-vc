/* finalexpo.h
 * by Joppe W. Bos, Michael Naehrig, Cryptography Research Group, XCG, Microsoft Research 2013
 *
 * Header file for computing the final exponentiation part of the pairing algorithm.
 * This file is part of the Pairing library.
 * We assume the target platform has access to 64-bit data-types.
 */

#ifndef __FINALEXPO_H
#define __FINALEXPO_H

#include "../arith_error.h"
#include "../Fp4/Fp4.h"
#include "../Fp12/Fp12.h"


/***************** Internal stuff *******************/

typedef struct {
  int lenu;
  int *nafu;
  int signu;
} bn_finalexpo_config_t;

extern bn_finalexpo_config_t bn_finalexpo_config; /* Defined in finalexpo.c */

void bn_finalexpo_neg (Fp12_t c, Fp12_t a);

int bn_finalexpo_initialize_config (void);

void bn_finalexpo_free_config (void);

#endif /* __FINALEXPO_H */
// FpGen is a copy-pasta of Fp, "genericized" to take a config as a parameter
// rather than a global, so that it's a class that can be multiply
// instantiated with different constants. This is probably a reasonable way
// to maintain this code, if we want to push it upstream.
/* Fp/Fp.h
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
 
#ifndef __FPGEN_H
#define __FPGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../arith.h"
#include "../mont_arith.h"

#pragma intrinsic(__rdtsc)

/* Our configuration structure. */
typedef struct {
  montmul_t *m;
  void *mem;
  uint64_t *p; /* This pointer either contain NULL, if regular Montgomery 
                * multiplication is used, or the original modulus (and 2*p
                * is stored in the Montgomery structure). This is useful when
                * the result needs to be completely reduced modulo p.
                */
  uint64_t *p2N; /* contains 2^N*p for reduction after FplGen_t addition in the lazy
                  * lazy reduction techniques.
                  */
  /* The function pointer to the appropriate Montgomery multiplication routine */                
  void (*montmul)(uint64_t*,uint64_t*,uint64_t*,montmul_t*,void*); 
} FpGen_config_t;

/* An element in F_p consists simply of an array of 64-bit limbs. */
typedef struct {
  uint64_t *limbs;
} fpGen_t;

/* Our basic type is a pointer to the structure. */
typedef fpGen_t FpGen_t[1];
typedef FpGen_t FplGen_t;

void FpGen_free_config (FpGen_config_t* cfg);

/* Initialize an element in Fp given the parameters passed
 * to FpGen_modulus_init.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int FpGen_init(FpGen_config_t* cfg, FpGen_t dst);
int FplGen_init(FpGen_config_t* cfg, FplGen_t dst);

/* Free the memory allocated in src,
 * set the number of limbs to zero and point to NULL. 
 */
void FpGen_free (FpGen_config_t* cfg, FpGen_t src);
void FplGen_free (FpGen_config_t* cfg, FplGen_t src);

/* Compute the modular inverse
 * dst = src^-1 mod p
 */
void FpGen_modinv (FpGen_config_t* cfg, FpGen_t dst, FpGen_t src);

/* Compute c = a * b mod p */
void FpGen_mul (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b);

/* Print the (regular non-Montgomery form) value of a */
void FpGen_print (FpGen_config_t* cfg, FpGen_t a);

/* Set c to the value of a */
void FpGen_set (FpGen_config_t* cfg, FpGen_t c, uint64_t *a);

/* c = a - b mod p */
void FpGen_sub (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b);

/* TODO FIX We assume here that we can add without overflow. */
void FpGen_div2 (FpGen_config_t* cfg, FpGen_t c, FpGen_t a);

/* Addition without reduction. */
void FpGen_sub_no_red (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b);

/* Addition without reduction. No check for overflow! */
void FplGen_sub_no_red (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b);

/* Subtraction of Fpl elements, option 1 with h=1 or h=2 in Aranha et al. */
void FplGen_sub_o1_ct (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b, int h, void *mem);

/* Subtraction of Fpl elements, option 2 in Aranha et al. */
void FplGen_sub_o2_ct (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b, void *mem);

/* Addition of Fpl elements, constant-time version of "option 2" addition in Aranha et al. */
void FplGen_add_o2_ct (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b, void *mem);

/* c = a + b mod p */
void FpGen_add (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b);

/* Addition without reduction. */
void FpGen_add_no_red (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, FpGen_t b);

/* Addition without reduction on long elements. */
void FplGen_add_no_red (FpGen_config_t* cfg, FplGen_t c, FplGen_t a, FplGen_t b);

/* Return +1 if a > b
 *         0 if a = b
 *        -1 if a < b
 */
int FpGen_cmp (FpGen_config_t* cfg, FpGen_t a, FpGen_t b);

/* Just the multiplication part, no reduction.
 * c= a * b
 */
void FpGen_mul_no_red (FpGen_config_t* cfg, FplGen_t c, FpGen_t a, FpGen_t b);

/* Reduce a double-length number a
 * c = a mod p
 */
void FplGen_red (FpGen_config_t* cfg, FpGen_t c, FplGen_t a);

void FpGen_copy (FpGen_config_t* cfg, FpGen_t c, FpGen_t a);

void FplGen_copy (FpGen_config_t* cfg, FplGen_t c, FplGen_t a);

void FpGen_neg (FpGen_config_t* cfg, FpGen_t c, FpGen_t a);

void FpGen_rand (FpGen_config_t* cfg, FpGen_t c);

void FpGen_get (FpGen_config_t* cfg, uint64_t *t, FpGen_t a);

void FpGen_mul3 (FpGen_config_t* cfg, FpGen_t c, FpGen_t a);

//void FpGen_mul_2exp (FpGen_config_t* cfg, FpGen_t c, FpGen_t a, int n);

#define FpGen_sqr(cfg, c,a) FpGen_mul (cfg, c, a, a)
//void FpGen_sqr (FpGen_t c, FpGen_t a);

void FpGen_set_ui (FpGen_config_t* cfg, FpGen_t c, uint64_t x);

/***************** Internal stuff *******************/

/* Initialize the global parameters for the modulus
 * given the modulus mod consisting of n 64-bit words.
 * Return > 0 on success and
 * Return < 0 on error.
 */ 
int FpGen_initialize_config (FpGen_config_t* cfg, uint64_t *m, int n);

/* Macros redirecting the memory allocating functions. */
#define FplGen_free(x) FpGen_free (x)

/* Internal function.
 * Returns a^-1 mod 2^64
 */
uint64_t modinv64 (uint64_t a);

#endif /* __FPGEN_H */

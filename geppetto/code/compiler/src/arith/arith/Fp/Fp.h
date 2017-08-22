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
 
#ifndef __FP_H
#define __FP_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "../arith.h"
#include "../mont_arith.h"

#pragma intrinsic(__rdtsc)

/* An element in F_p consists simply of an array of 64-bit limbs. */
typedef struct {
  uint64_t *limbs;
} fp_t;

/* Our basic type is a pointer to the structure. */
typedef fp_t Fp_t[1];
typedef Fp_t Fpl_t;

/* Initialize the global parameters for the modulus
 * given the modulus mod consisting of n 64-bit words.
 * Return > 0 on success and
 * Return < 0 on error.
 */ 
int Fp_initialize_config (uint64_t *m, int n);

void Fp_free_config (void);

/* Initialize an element in Fp given the parameters passed
 * to Fp_modulus_init.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp_init(Fp_t dst);
int Fpl_init(Fpl_t dst);

/* Free the memory allocated in src,
 * set the number of limbs to zero and point to NULL. 
 */
void Fp_free (Fp_t src);
void Fpl_free (Fpl_t src);

/* Compute the modular inverse
 * dst = src^-1 mod p
 */
void Fp_modinv (Fp_t dst, Fp_t src);

/* Compute c = a * b mod p */
void Fp_mul (Fp_t c, Fp_t a, Fp_t b);

/* Print the (regular non-Montgomery form) value of a */
void Fp_print (Fp_t a);

/* Set c to the value of a */
void Fp_set (Fp_t c, uint64_t *a);

/* c = a - b mod p */
void Fp_sub (Fp_t c, Fp_t a, Fp_t b);

/* TODO FIX We assume here that we can add without overflow. */
void Fp_div2 (Fp_t c, Fp_t a);

/* Addition without reduction. */
void Fp_sub_no_red (Fp_t c, Fp_t a, Fp_t b);

/* Addition without reduction. No check for overflow! */
void Fpl_sub_no_red (Fpl_t c, Fpl_t a, Fpl_t b);

/* Subtraction of Fpl elements, option 1 with h=1 or h=2 in Aranha et al. */
void Fpl_sub_o1_ct (Fpl_t c, Fpl_t a, Fpl_t b, int h, void *mem);

/* Subtraction of Fpl elements, option 2 in Aranha et al. */
void Fpl_sub_o2_ct (Fpl_t c, Fpl_t a, Fpl_t b, void *mem);

/* Addition of Fpl elements, constant-time version of "option 2" addition in Aranha et al. */
void Fpl_add_o2_ct (Fpl_t c, Fpl_t a, Fpl_t b, void *mem);

/* c = a + b mod p */
void Fp_add (Fp_t c, Fp_t a, Fp_t b);

/* Addition without reduction. */
void Fp_add_no_red (Fp_t c, Fp_t a, Fp_t b);

/* Addition without reduction on long elements. */
void Fpl_add_no_red (Fpl_t c, Fpl_t a, Fpl_t b);

/* Return +1 if a > b
 *         0 if a = b
 *        -1 if a < b
 */
int Fp_cmp (Fp_t a, Fp_t b);

/* Just the multiplication part, no reduction.
 * c= a * b
 */
void Fp_mul_no_red (Fpl_t c, Fp_t a, Fp_t b);

/* Reduce a double-length number a
 * c = a mod p
 */
void Fpl_red (Fp_t c, Fpl_t a);

void Fp_copy (Fp_t c, Fp_t a);

void Fpl_copy (Fpl_t c, Fpl_t a);

void Fp_neg (Fp_t c, Fp_t a);

void Fp_rand (Fp_t c);

void Fp_get (uint64_t *t, Fp_t a);

void Fp_mul3 (Fp_t c, Fp_t a);

//void Fp_mul_2exp (Fp_t c, Fp_t a, int n);

#define Fp_sqr(c,a) Fp_mul (c, a, a)
//void Fp_sqr (Fp_t c, Fp_t a);

void Fp_set_ui (Fp_t c, uint64_t x);

void p3mod4_Fp_sqrt(Fp_t c, Fp_t a);

void Fp_sqrt(Fp_t c, Fp_t a);

/***************** Internal stuff *******************/

/* Our configuration structure. */
typedef struct {
  montmul_t *m;
  void *mem;
  uint64_t *p; /* This pointer either contain NULL, if regular Montgomery 
                * multiplication is used, or the original modulus (and 2*p
                * is stored in the Montgomery structure). This is useful when
                * the result needs to be completely reduced modulo p.
                */
  uint64_t *p2N; /* contains 2^N*p for reduction after Fpl_t addition in the
                  * lazy reduction techniques.
                  */

  /* The exponent (p+1)/4 used to compute square roots in Fp as a^((p+1)/4) if p = 3 (mod 4).*/
  uint64_t *sqrt_exp_p3mod4;
  unsigned int sqrt_exp_len_p3mod4;

  /* The function pointer to the appropriate Montgomery multiplication routine */                
  void (*montmul)(uint64_t*,uint64_t*,uint64_t*,montmul_t*,void*);
  
  void(*Fp_sqrt)(Fp_t, Fp_t); /*Square root function for Fp.*/

  Fp_t tmp; // Temp element
} Fp_config_t;

extern Fp_config_t Fp_config; /* Defined in Fp.c */

/* Macros redirecting the memory allocating functions. */
#define Fpl_free(x) Fp_free (x)

/* Internal function.
 * Returns a^-1 mod 2^64
 */
uint64_t modinv64 (uint64_t a);

#endif /* __FP_H */

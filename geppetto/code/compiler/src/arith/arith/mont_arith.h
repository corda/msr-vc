/* mont_arith.h
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

#ifndef __MONT_ARITH_H
#define __MONT_ARITH_H

#include "uint.h"

/* Data used for Montgomery multiplication. */
typedef struct {
  uint64_t *p;    // the prime modulus used
  int n;          // the number of limbs
  uint64_t mu;    // -p^-1 mod 2^64  
  uint64_t *R2;   // R^2 mod p = (2^(64*n))^2 mod p
  uint64_t *R3;   // R^3 mod p = (2^(64*n))^3 mod p
} montmul_t;

/* Convert a to its Montgomery residue. */
void normal2mont (uint64_t *c, uint64_t *a, montmul_t *m, void *mem);

/* Convert the Montgomery residue a back to its regular form. */
void mont2normal (uint64_t *c, uint64_t *a, montmul_t *m, void *mem);

/* Montgomery multiplication in constant-time. */
void montmul_ct (uint64_t *c, uint64_t *a, uint64_t *b, montmul_t *m, void *tmp);

/* Montgomery reduction of a double-length integer c in constant-time */
void montred_ct (uint64_t *c, uint64_t *a, montmul_t *m, void *tmp);

/* Montgomery addition (which is the same as modular addition) */
void mod_add (uint64_t *c, uint64_t *a, uint64_t *b, montmul_t *m, void *tmp);

/* Montgomery subtraction (which is the same as modular subtraction) */
void mod_sub (uint64_t *c, uint64_t *a, uint64_t *b, montmul_t *m, void *tmp);

#endif /* __MONT_ARITH_H */

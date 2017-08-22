/* sint.c
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

#ifndef __SINT_H
#define __SINT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct {
  uint64_t *n;
  int sign;
} sint;

typedef sint sint_t[1];

int sint_init (sint_t a, int n);
void sint_free (sint_t b);

void sint_sub (sint_t c, sint_t a, sint_t b, int n);
void sint_add (sint_t c, sint_t a, sint_t b, int n);
void sint_neg (sint_t c, sint_t a, int n);
void sint_mul (sint_t c, sint_t a, sint_t b, uint64_t *mem, int n);
void sint_set (sint_t c, sint_t a, int n);
void sint_set_mn (sint_t c, sint_t a, int m, int n);
void sint_div2 (sint_t c, sint_t a, int n);

#define sint_neg_single(c) (c)->sign = -(c)->sign

#endif /* __SINT_H */


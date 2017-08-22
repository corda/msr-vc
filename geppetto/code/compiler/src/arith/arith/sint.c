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

#include "sint.h"
#include "uint.h"
#include "arith_error.h"

int sint_init (sint_t a, int n) {
  a->n = NULL;
  a->n = (uint64_t *) malloc (n * sizeof (uint64_t));
  if (a->n == NULL) return ERR_OUT_OF_MEMORY;
  a->sign = 0;
  return ERR_SUCCESS;
}

void sint_free (sint_t b) {
  free (b->n);
}

void sint_neg (sint_t c, sint_t a, int n) {
  int i;

  for (i=0; i < n; i++) {
    c->n[i] = a->n[i];
  }
  c->sign = -a->sign;
}

/*void sint_set_ui (snum_t *c, uint64_t *a) {
  int i;

  for (i=0; i < 2; i++) {
    c->n[i] = a[i];
  }
  for (; i < c->alen; i++) {
    c->n[i] = 0;
  }
  c->sign = 1;
}*/

/* a and b have the same length. */
void sint_add (sint_t c, sint_t a, sint_t b, int n) {
  if (a->sign > 0 && b->sign > 0) { 
    uint_add (c->n, a->n, b->n, n); 
    c->sign = 1; 
  } else if (a->sign > 0 && b->sign < 0) {
    if (uint_cmp (a->n, b->n, n) > 0) {
      uint_sub (c->n, a->n, b->n, n);
      c->sign = 1;
    } else {
      uint_sub (c->n, b->n, a->n, n);
      c->sign = -1;
    }
  } else if (a->sign < 0 && b->sign > 0) {
    if (uint_cmp (b->n, a->n, n) > 0) {
      uint_sub (c->n, b->n, a->n, n); 
      c->sign = 1;
    } else {
      uint_sub (c->n, a->n, b->n, n); 
      c->sign = -1;
    }
  } else if (a->sign < 0 && b->sign < 0) { 
    uint_add (c->n, a->n, b->n, n); 
    c->sign = -1;
  }
}

void sint_sub (sint_t c, sint_t a, sint_t b, int n) {
  if (a->sign > 0 && b->sign > 0) {
    if (uint_cmp (a->n, b->n, n) > 0) {
      uint_sub (c->n, a->n, b->n, n); 
      c->sign = 1;
    } else {
      uint_sub (c->n, b->n, a->n, n); 
      c->sign = -1;
    }
  } else if (a->sign > 0 && b->sign < 0) { 
    uint_add (c->n, a->n, b->n, n); 
    c->sign = a->sign; 
  } else if (a->sign < 0 && b->sign > 0) { 
    uint_add (c->n, a->n, b->n, n); 
    c->sign = a->sign;
  } else if (a->sign < 0 && b->sign < 0) {
    if (uint_cmp (b->n, a->n, n) > 0) {
      uint_sub (c->n, b->n, a->n, n); 
      c->sign = 1;
    } else {
      uint_sub (c->n, a->n, b->n, n); 
      c->sign = -1;
    }
  }
}

void sint_mul (sint_t c, sint_t a, sint_t b, uint64_t *mem, int n) {
  uint_mul (c->n, a->n, b->n, mem, n);
  c->sign = a->sign * b->sign;
}

void sint_set (sint_t c, sint_t a, int n) {
  int i;
  for (i=0; i < n; i++) c->n[i] = a->n[i];
  c->sign = a->sign;
}

void sint_set_mn (sint_t c, sint_t a, int m, int n) {
  int i;
  for (i=0; i < n; i++) c->n[i] = a->n[i];
  for (; i < m; i++) c->n[i] = 0;
  c->sign = a->sign;
}

void sint_div2 (sint_t c, sint_t a, int n) {
  uint_div2 (c->n, a->n, n);
  c->sign = a->sign;
}

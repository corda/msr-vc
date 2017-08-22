/* Fp12/Fp12.c
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

#include "Fp12.h"
#include "../pairing/pairing.h"

Fp12_config_t Fp12_config;

/* Initialize an element in Fp12 by initializing two Fp6 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp12_init(Fp12_t dst) {
	int err = 0;
	err = Fp6_init(dst->a0);
	if (err < 0) return err;
	return Fp6_init(dst->a1);
}

/* Free the memory allocated in src. */
void Fp12_free (Fp12_t src) {
	Fp6_free(src->a0);
	Fp6_free(src->a1);
}

/* Get a random element in Fp12. */
void Fp12_rand (Fp12_t c) {
  Fp6_rand(c->a0);
  Fp6_rand(c->a1);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp12_print (Fp12_t a) {
  Fp6_print(a->a0);
  //printf("\n");
  Fp6_print(a->a1);
  //printf("\n");
}

/* Set c to 0. */
void Fp12_set_zero (Fp12_t c) {
  Fp6_set_zero (c->a0);
  Fp6_set_zero (c->a1);
}

/* Set c to 1. */
void Fp12_set_one (Fp12_t c) {
  Fp6_set_one (c->a0);
  Fp6_set_zero (c->a1);
}

/* Copy an Fp12 element a to an Fp12 element c. */
void Fp12_copy(Fp12_t c, Fp12_t a) {
	Fp6_copy(c->a0, a->a0);
	Fp6_copy(c->a1, a->a1);
}

/* Compare two Fp12 elements a and b for equality. */
int Fp12_cmpeq(Fp12_t a, Fp12_t b) {
  return (Fp6_cmpeq(a->a0, b->a0) && Fp6_cmpeq(a->a1, b->a1));
}

/* Negate an Fp12 element. */
void Fp12_neg(Fp12_t c, Fp12_t a) {
	Fp6_neg(c->a0, a->a0);
	Fp6_neg(c->a1, a->a1);
}

/* Add two Fp12 elements coefficient wise. */
void Fp12_add (Fp12_t c, Fp12_t a, Fp12_t b) {
	Fp6_add(c->a0, a->a0, b->a0);
	Fp6_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp12 element b from an Fp12 element a coefficient wise. */
void Fp12_sub (Fp12_t c, Fp12_t a, Fp12_t b) {
	Fp6_sub(c->a0, a->a0, b->a0);
	Fp6_sub(c->a1, a->a1, b->a1);
}

/* Multiply an Fp12 element by an Fp6 element. */
void Fp12_mulFp6(Fp12_t c, Fp12_t a, Fp6_t b) {
	Fp6_mul(c->a0, a->a0, b);
	Fp6_mul(c->a1, a->a1, b);
}

/* Multiply an Fp12 element by an Fp2 element. */
void Fp12_mulFp2(Fp12_t c, Fp12_t a, Fp2_t b) {
	Fp6_mulFp2(c->a0, a->a0, b);
	Fp6_mulFp2(c->a1, a->a1, b);
}

/* Multiply an Fp12 element by an Fp element. */
void Fp12_mulFp(Fp12_t c, Fp12_t a, Fp_t b) {
	Fp6_mulFp(c->a0, a->a0, b);
	Fp6_mulFp(c->a1, a->a1, b);
}

/* Multiply two Fp12 elements. */
void Fp12_mul (Fp12_t c, Fp12_t a, Fp12_t b) {
  Fp12_config.Fp12_mul(c, a, b, Fp12_config.mem);
}

/* Square an Fp12 element. */
void Fp12_squ (Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_squ(c, a, Fp12_config.mem);
}

/* Compute the p-power Frobenius of an Fp12 element. */
void Fp12_ppow(Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_ppow(c, a);
}

/* Compute the p^2-power Frobenius of an Fp12 element. */
void Fp12_p2pow(Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_p2pow(c, a);
}

/* Compute the p^3-power Frobenius of an Fp12 element. */
void Fp12_p3pow(Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_p3pow(c, a);
}

/* Compute the p^6-power Frobenius of an Fp12 element. */
void Fp12_p6pow(Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_p6pow(c, a);
}

/* Invert an Fp12 element. */
void Fp12_inv (Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_inv(c, a, Fp12_config.mem);
}

/* Multiply two Fp12 elements using lazy reduction.*/
void Fp12_mul_lazy(Fp12_t c, Fp12_t a, Fp12_t b) {
  Fp12_config.Fp12_mul_lazy (c, a, b, Fp12_config.mem);
}

/* Square an Fp12 element using lazy reduction.*/
void Fp12_squ_lazy(Fp12_t c, Fp12_t a) {
  Fp12_config.Fp12_squ_lazy (c, a, Fp12_config.mem);
}

/* Multiply two Fp12 elements (the second one sparse) using lazy reduction.*/
void Fp12_mul_sparse_twist_lazy(Fp12_t c, Fp12_t a, Fp2_t b00, Fp2_t b01, Fp2_t b11) {
  Fp12_config.Fp12_mul_sparse_twist_lazy (c, a, b00, b01, b11, Fp12_config.mem);
}

/* Multiply two Fp12 elements (the second one sparse) using lazy reduction.*/
void Fp12_mul_sparse_untwist_lazy(Fp12_t c, Fp12_t a, Fp2_t b01, Fp2_t b11, Fp2_t b12) {
  Fp12_config.Fp12_mul_sparse_untwist_lazy (c, a, b01, b11, b12, Fp12_config.mem);
}

/* Multiply two Fp12 elements. This uses that Fp12 is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 20 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void w2minusv_Fp12_mul(Fp12_t c, Fp12_t a, Fp12_t b, void *mem) {
	Fp6_t t0, t1, t2;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
	t2->a0->a0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  t2->a0->a1->limbs = (uint64_t *) mem + 13*Fp_config.m->n;
	t2->a1->a0->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  t2->a1->a1->limbs = (uint64_t *) mem + 15*Fp_config.m->n;
	t2->a2->a0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  t2->a2->a1->limbs = (uint64_t *) mem + 17*Fp_config.m->n;
	
	Fp6_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp6_add(t1, b->a0, b->a1);	// t1 = b0 + b1
	Fp6_mul(t2, t0, t1);			  // t2 = t0*t1 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1

	Fp6_mul(t0, a->a0, b->a0);		// t0 = a0*b0
	Fp6_mul(t1, a->a1, b->a1);		// t1 = a1*b1
	
	Fp6_sub(t2, t2, t0);					// t2 = t2 - t0 = a1*b0 + a0*b1 + a1*b1
	Fp6_sub(c->a1, t2, t1);				// c1 = t2 - t1 = a1*b0 + a0*b1

	Fp6_mulv(t1, t1);						  // t1 = v*t1 = v*a1*b1
	Fp6_add(c->a0, t0, t1);				// c0 = t0 + t1 = a0*b0 + v*a1*b1
}

/* Square a general Fp12 element. This uses that Fp12 is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 22 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void w2minusv_Fp12_squ(Fp12_t c, Fp12_t a, void *mem) {
	Fp6_t t0, t1;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;

	Fp6_sub(t0, a->a0, a->a1);		// t0 = a0 - a1
	Fp6_mulv(t1, a->a1);					// t1 = v*a1
	Fp6_sub(t1, a->a0, t1);				// t1 = a0 - v*a1
	Fp6_mul(t0, t0, t1);			    // t0 = (a0 - a1)*(a0 - v*a1) = a0^2 + v*a1^2 - a0*a1 - v*a0*a1
	Fp6_mul(t1, a->a0, a->a1);		// t1 = a0*a1
	Fp6_add(t0, t0, t1);					// t0 =  a0^2 + v*a1^2 - v*a0*a1
	Fp6_add(c->a1, t1, t1);			  // c1 = 2*t1 = 2*a0*a1
	Fp6_mulv(t1, t1);						  // t1 = v*t1 = v*a0*a1
	Fp6_add(c->a0, t0, t1);				// c0 = t0 + t1 = a0^2 + v*a1^2
}

/* Compute the p-power Frobenius of an Fp12 element. This uses a constant WPPOW, 
 * such that w^p = WPPOW*w, i.e. WPPOW = w^(p-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and Fp12 is an extension of Fp6 via w^2 - v, and p = 1 (mod 6),
 * we have WPPOW = xi^(p-1)/6 in Fp2.
 */
void w2minusv_v3minusxi_p1mod6_Fp12_ppow(Fp12_t c, Fp12_t a) {
	Fp6_ppow(c->a0, a->a0);
	Fp6_ppow(c->a1, a->a1);
	Fp6_mulFp2(c->a1, c->a1, Fp12_config.wppow);
}

/* Compute the p^2-power Frobenius of an Fp12 element. This uses a constant WP2POW, 
 * such that w^(p^2) = WP2POW*w, i.e. WP2POW = w^(p^2-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and Fp12 is an extension of Fp6 via w^2 - v, and p = 1 (mod 6),
 * we have WPPOW = xi^(p^2-1)/6 in Fp2.
 */
void w2minusv_v3minusxi_p1mod6_Fp12_p2pow(Fp12_t c, Fp12_t a) {
	Fp6_p2pow(c->a0, a->a0);
	Fp6_p2pow(c->a1, a->a1);
	Fp6_mulFp(c->a1, c->a1, Fp12_config.wp2pow);
}

/* Compute the p^3-power Frobenius of an Fp12 element. This uses a constant WP3POW, 
 * such that w^(p^3) = WP3POW*w, i.e. WP3POW = w^(p^3-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and Fp12 is an extension of Fp6 via w^2 - v, and p = 1 (mod 6),
 * we have WP3POW = xi^(p^3-1)/6 in Fp2.
 */
void w2minusv_v3minusxi_p1mod6_Fp12_p3pow(Fp12_t c, Fp12_t a) {
	Fp6_p3pow(c->a0, a->a0);
	Fp6_p3pow(c->a1, a->a1);
	Fp6_mulFp2(c->a1, c->a1, Fp12_config.wp3pow);
}

/* Compute the p^2-power Frobenius of an Fp6 element. This is the same as 
 * conjugation in the quadratic extension over Fp6 if the extension is constructed
 * via a quadratic polynomial x^2 - alpha, alpha in Fp6. 
 */
void w2minusv_Fp12_p6pow(Fp12_t c, Fp12_t a) {
	Fp6_copy(c->a0, a->a0);
	Fp6_neg(c->a1, a->a1);
}

/* Invert an Fp12 element. This uses that Fp12 is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 23 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void w2minusv_Fp12_inv(Fp12_t c, Fp12_t a, void *mem) {
	Fp6_t t0, t1;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;

	Fp6_squ(t0, a->a0);				// t0 = a0^2
	Fp6_squ(t1, a->a1);				// t1 = a1^2
	Fp6_mulv(t1, t1);						// t1 = v*t1 = v*a1^2
	Fp6_sub(t0, t0, t1);					// t0 = a0^2 - v*a1^2 = Norm(a)
	Fp6_inv(t1, t0);				// t1 = t0^(-1) = Norm(a)^(-1)
	Fp6_mul(c->a0, a->a0, t1);				// c0 = a0/(a0^2 - v*a1^2)
	Fp6_neg(t1, t1);						// t1 = -1/(a0^2 - v*a1^2)
	Fp6_mul(c->a1, a->a1, t1);				// c1 = a1*t1 = -a1/(a0^2 - v*a1^2)
}

/* Multiply two Fp12 elements using lazy reduction techniques. 
 * This uses that Fp12 is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. This is Algorithm 4 in Aranha et al. 
 */
void w2minusv_Fp12_mul_lazy(Fp12_t c, Fp12_t a, Fp12_t b, void *mem) {
	Fp6_t t0, t1;
  Fp6l_t T0, T1, T2, T3;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  T0->A0->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T0->A0->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T0->A1->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T0->A1->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
	T0->A2->A0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  T0->A2->A1->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
	T1->A0->A0->limbs = (uint64_t *) mem + 24*Fp_config.m->n;
  T1->A0->A1->limbs = (uint64_t *) mem + 26*Fp_config.m->n;
	T1->A1->A0->limbs = (uint64_t *) mem + 28*Fp_config.m->n;
  T1->A1->A1->limbs = (uint64_t *) mem + 30*Fp_config.m->n;
	T1->A2->A0->limbs = (uint64_t *) mem + 32*Fp_config.m->n;
  T1->A2->A1->limbs = (uint64_t *) mem + 34*Fp_config.m->n;
	T2->A0->A0->limbs = (uint64_t *) mem + 36*Fp_config.m->n;
  T2->A0->A1->limbs = (uint64_t *) mem + 38*Fp_config.m->n;
	T2->A1->A0->limbs = (uint64_t *) mem + 40*Fp_config.m->n;
  T2->A1->A1->limbs = (uint64_t *) mem + 42*Fp_config.m->n;
	T2->A2->A0->limbs = (uint64_t *) mem + 44*Fp_config.m->n;
  T2->A2->A1->limbs = (uint64_t *) mem + 46*Fp_config.m->n;
  T3->A0->A0->limbs = (uint64_t *) mem + 48*Fp_config.m->n;
  T3->A0->A1->limbs = (uint64_t *) mem + 50*Fp_config.m->n;
	T3->A1->A0->limbs = (uint64_t *) mem + 52*Fp_config.m->n;
  T3->A1->A1->limbs = (uint64_t *) mem + 54*Fp_config.m->n;
	T3->A2->A0->limbs = (uint64_t *) mem + 56*Fp_config.m->n;
  T3->A2->A1->limbs = (uint64_t *) mem + 58*Fp_config.m->n;

  Fp6_mul_no_red (T0, a->a0, b->a0);
  Fp6_mul_no_red (T1, a->a1, b->a1);
  Fp6_add (t0, a->a0, a->a1);
  Fp6_add (t1, b->a0, b->a1);
  Fp6_mul_no_red (T2, t0, t1);
  Fp6l_add_o2 (T3, T0, T1);  
  Fp6l_sub_o2 (T2, T2, T3);
  Fp6l_red (c->a1, T2);

  Fpl_sub_o2_ct (T2->A0->A0, T1->A2->A0, T1->A2->A1, Fp_config.mem);
  Fpl_add_o2_ct (T2->A0->A1, T1->A2->A0, T1->A2->A1, Fp_config.mem);
  Fp2l_copy (T2->A1, T1->A0);
  Fp2l_copy (T2->A2, T1->A1);
  Fp6l_add_o2 (T2, T0, T2);
  Fp6l_red (c->a0, T2);
}

/* Multiply two Fp12 elements using lazy reduction techniques, 
 * where the second element b in sparse in the following sense:
 * b = b0 + b1*w, where b0 = b00 + b01*v and b1 = b11*v.
 * This uses that Fp12 is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. 
 */
void w2minusv_Fp12_mul_sparse_twist_lazy(Fp12_t c, Fp12_t a, Fp2_t b00, Fp2_t b01, Fp2_t b11, void *mem) {
	Fp6_t t0, t1;
  Fp6l_t T0, T1, T2, T3;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  T0->A0->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T0->A0->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T0->A1->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T0->A1->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
	T0->A2->A0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  T0->A2->A1->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
	T1->A0->A0->limbs = (uint64_t *) mem + 24*Fp_config.m->n;
  T1->A0->A1->limbs = (uint64_t *) mem + 26*Fp_config.m->n;
	T1->A1->A0->limbs = (uint64_t *) mem + 28*Fp_config.m->n;
  T1->A1->A1->limbs = (uint64_t *) mem + 30*Fp_config.m->n;
	T1->A2->A0->limbs = (uint64_t *) mem + 32*Fp_config.m->n;
  T1->A2->A1->limbs = (uint64_t *) mem + 34*Fp_config.m->n;
	T2->A0->A0->limbs = (uint64_t *) mem + 36*Fp_config.m->n;
  T2->A0->A1->limbs = (uint64_t *) mem + 38*Fp_config.m->n;
	T2->A1->A0->limbs = (uint64_t *) mem + 40*Fp_config.m->n;
  T2->A1->A1->limbs = (uint64_t *) mem + 42*Fp_config.m->n;
	T2->A2->A0->limbs = (uint64_t *) mem + 44*Fp_config.m->n;
  T2->A2->A1->limbs = (uint64_t *) mem + 46*Fp_config.m->n;
  T3->A0->A0->limbs = (uint64_t *) mem + 48*Fp_config.m->n;
  T3->A0->A1->limbs = (uint64_t *) mem + 50*Fp_config.m->n;
	T3->A1->A0->limbs = (uint64_t *) mem + 52*Fp_config.m->n;
  T3->A1->A1->limbs = (uint64_t *) mem + 54*Fp_config.m->n;
	T3->A2->A0->limbs = (uint64_t *) mem + 56*Fp_config.m->n;
  T3->A2->A1->limbs = (uint64_t *) mem + 58*Fp_config.m->n;

  Fp6_mul_sparse01_no_red (T0, a->a0, b00, b01); // b0 = b00 + b01*v
  Fp6_mulv (t0, a->a1);
  Fp6_mulFp2_no_red (T1, t0, b11);
  Fp6_add (t0, a->a0, a->a1);
  Fp2_add (t1->a1, b01, b11);
  Fp6_mul_sparse01_no_red (T2, t0, b00, t1->a1); // t1 = b00 + (b01 + b11)*v
  Fp6l_add_o2 (T3, T0, T1);  
  Fp6l_sub_o2 (T2, T2, T3);
  Fp6l_red (c->a1, T2);

  Fpl_sub_o2_ct (T2->A0->A0, T1->A2->A0, T1->A2->A1, Fp_config.mem);
  Fpl_add_o2_ct (T2->A0->A1, T1->A2->A0, T1->A2->A1, Fp_config.mem);
  Fp2l_copy (T2->A1, T1->A0);
  Fp2l_copy (T2->A2, T1->A1);
  Fp6l_add_o2 (T2, T0, T2);
  Fp6l_red (c->a0, T2);
}

/* Similar as above, but now b = b0 + b1*w, where b0 = b01*v and b1 = b11*v + b12*v^2. */
void w2minusv_Fp12_mul_sparse_untwist_lazy(Fp12_t c, Fp12_t a, Fp2_t b01, Fp2_t b11, Fp2_t b12, void *mem) {
	Fp6_t t0, t1;
  Fp6l_t T0, T1, T2, T3;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  T0->A0->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T0->A0->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T0->A1->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T0->A1->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
	T0->A2->A0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  T0->A2->A1->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
	T1->A0->A0->limbs = (uint64_t *) mem + 24*Fp_config.m->n;
  T1->A0->A1->limbs = (uint64_t *) mem + 26*Fp_config.m->n;
	T1->A1->A0->limbs = (uint64_t *) mem + 28*Fp_config.m->n;
  T1->A1->A1->limbs = (uint64_t *) mem + 30*Fp_config.m->n;
	T1->A2->A0->limbs = (uint64_t *) mem + 32*Fp_config.m->n;
  T1->A2->A1->limbs = (uint64_t *) mem + 34*Fp_config.m->n;
	T2->A0->A0->limbs = (uint64_t *) mem + 36*Fp_config.m->n;
  T2->A0->A1->limbs = (uint64_t *) mem + 38*Fp_config.m->n;
	T2->A1->A0->limbs = (uint64_t *) mem + 40*Fp_config.m->n;
  T2->A1->A1->limbs = (uint64_t *) mem + 42*Fp_config.m->n;
	T2->A2->A0->limbs = (uint64_t *) mem + 44*Fp_config.m->n;
  T2->A2->A1->limbs = (uint64_t *) mem + 46*Fp_config.m->n;
  T3->A0->A0->limbs = (uint64_t *) mem + 48*Fp_config.m->n;
  T3->A0->A1->limbs = (uint64_t *) mem + 50*Fp_config.m->n;
	T3->A1->A0->limbs = (uint64_t *) mem + 52*Fp_config.m->n;
  T3->A1->A1->limbs = (uint64_t *) mem + 54*Fp_config.m->n;
	T3->A2->A0->limbs = (uint64_t *) mem + 56*Fp_config.m->n;
  T3->A2->A1->limbs = (uint64_t *) mem + 58*Fp_config.m->n;


  Fp6_mulv (t0, a->a0);
  Fp6_mulFp2_no_red (T0, t0, b01);
  Fp6_mul_sparse12_no_red (T1, a->a1, b11, b12); 
  Fp6_add (t0, a->a0, a->a1);
  Fp2_add (t1->a1, b01, b11);
  Fp6_mul_sparse12_no_red (T2, t0, t1->a1, b12); // t1 = (b01 + b11)*v + b12*v^2
  Fp6l_add_o2 (T3, T0, T1);  
  Fp6l_sub_o2 (T2, T2, T3);
  Fp6l_red (c->a1, T2);

  Fpl_sub_o2_ct (T2->A0->A0, T1->A2->A0, T1->A2->A1, Fp_config.mem);
  Fpl_add_o2_ct (T2->A0->A1, T1->A2->A0, T1->A2->A1, Fp_config.mem);
  Fp2l_copy (T2->A1, T1->A0);
  Fp2l_copy (T2->A2, T1->A1);
  Fp6l_add_o2 (T2, T0, T2);
  Fp6l_red (c->a0, T2);
}


/* Square an Fp12 element using lazy reduction techniques. 
 * This uses that Fp12 is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. This is Algorithm 5 in Aranha et al. 
 */
void w2minusv_Fp12_squ_lazy(Fp12_t c, Fp12_t a, void *mem) {
	Fp6_t t0, t1;
  Fp6l_t T0;

  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	t0->a2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t0->a2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
	t1->a0->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t1->a1->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	t1->a2->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t1->a2->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  T0->A0->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T0->A0->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T0->A1->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T0->A1->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
	T0->A2->A0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  T0->A2->A1->limbs = (uint64_t *) mem + 22*Fp_config.m->n;

  Fp6_add (t0, a->a0, a->a1);
  Fp_sub (t1->a0->a0, a->a1->a2->a0, a->a1->a2->a1);
  Fp_add (t1->a0->a1, a->a1->a2->a0, a->a1->a2->a1);
  Fp2_copy (t1->a1, a->a1->a0);
  Fp2_copy (t1->a2, a->a1->a1);
  Fp6_add (t1, a->a0, t1);
  Fp6_mul_no_red (T0, a->a0, a->a1);
  Fp6l_red (c->a1, T0);
  
  Fp6_mul_no_red (T0, t0, t1);
  Fp6l_red (t0, T0);
  Fp_sub (t1->a0->a0, c->a1->a2->a0, c->a1->a2->a1);
  Fp_add (t1->a0->a1, c->a1->a2->a0, c->a1->a2->a1);
  Fp2_copy (t1->a1, c->a1->a0);
  Fp2_copy (t1->a2, c->a1->a1);
  Fp6_add (t1, c->a1, t1);
  Fp6_sub (c->a0, t0, t1);
  Fp6_add (c->a1, c->a1, c->a1);
}

int Fp12_initialize_config (void) {
  uint64_t *wpp;

  Fp12_config.mem = NULL;
  Fp12_config.mem = (void *) malloc (60 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp12_config.mem == NULL) return ERR_OUT_OF_MEMORY;

  if (Fp_init (Fp12_config.wp2pow) < 0) printf ("Fp12_config memory error.\n");
  if (Fp_init (Fp12_config.wp2pow3) < 0) printf ("Fp12_config memory error.\n");
  if (Fp2_init (Fp12_config.wppow) < 0) printf ("Fp12_config memory error.\n");
  if (Fp2_init (Fp12_config.wp3pow) < 0) printf ("Fp12_config memory error.\n");
  if (Fp2_init (Fp12_config.wppow3) < 0) printf ("Fp12_config memory error.\n");
  if (Fp2_init (Fp12_config.wppow3inv) < 0) printf ("Fp12_config memory error.\n");
  wpp = (uint64_t *) Fp12_config.mem;

  if (PAIR_CURVE == BN12) {
    wpp[3] = (uint64_t) 0x0000000000000000;
    wpp[2] = (uint64_t) 0x49B3624000000002;
    wpp[1] = (uint64_t) 0x4909000000000006;
    wpp[0] = (uint64_t) 0xCD80000000000008;
    Fp_set (Fp12_config.wp2pow, wpp);

    wpp[3] = (uint64_t) 0x1B377619212E7C8C;
    wpp[2] = (uint64_t) 0xB6499B50A846953F;
    wpp[1] = (uint64_t) 0x850974924D3F77C2;
    wpp[0] = (uint64_t) 0xE17DE6C06F2A6DE9;
    Fp_set (Fp12_config.wppow->a0, wpp);
  
    wpp[3] = (uint64_t) 0x09EBEE691ED18375;
    wpp[2] = (uint64_t) 0x03EAB22F57B96AC8;
    wpp[1] = (uint64_t) 0xDC178B6DB2C08850;
    wpp[0] = (uint64_t) 0xC582193F90D5922A;
    Fp_set (Fp12_config.wppow->a1, wpp);

    wpp[3] = (uint64_t) 0x01439AB09C60B248;
    wpp[2] = (uint64_t) 0xF398C5D77B755F92;
    wpp[1] = (uint64_t) 0xB9EDC5F19D287354;
    wpp[0] = (uint64_t) 0x5BE471151A747E4E;
    Fp_set (Fp12_config.wp3pow->a0, wpp);
  
    wpp[3] = (uint64_t) 0x23DFC9D1A39F4DB8;
    wpp[2] = (uint64_t) 0xC69B87A8848AA075;
    wpp[1] = (uint64_t) 0xA7333A0E62D78CBF;
    wpp[0] = (uint64_t) 0x4B1B8EEAE58B81C5;
    Fp_set (Fp12_config.wp3pow->a1, wpp);

    Fp_mul (Fp12_config.wp2pow3, Fp12_config.wp2pow, Fp12_config.wp2pow);
    Fp_mul (Fp12_config.wp2pow3, Fp12_config.wp2pow3, Fp12_config.wp2pow);
    Fp2_squ (Fp12_config.wppow3, Fp12_config.wppow);
    Fp2_mul (Fp12_config.wppow3, Fp12_config.wppow3, Fp12_config.wppow);
    Fp2_inv (Fp12_config.wppow3inv, Fp12_config.wppow3);

    Fp12_config.Fp12_mul = w2minusv_Fp12_mul;
    Fp12_config.Fp12_squ = w2minusv_Fp12_squ;
    Fp12_config.Fp12_inv = w2minusv_Fp12_inv;
    Fp12_config.Fp12_ppow = w2minusv_v3minusxi_p1mod6_Fp12_ppow;
    Fp12_config.Fp12_p2pow = w2minusv_v3minusxi_p1mod6_Fp12_p2pow;
    Fp12_config.Fp12_p3pow = w2minusv_v3minusxi_p1mod6_Fp12_p3pow;
    Fp12_config.Fp12_p6pow = w2minusv_Fp12_p6pow; 
    Fp12_config.Fp12_mul_lazy = w2minusv_Fp12_mul_lazy;
    Fp12_config.Fp12_squ_lazy = w2minusv_Fp12_squ_lazy;
    Fp12_config.Fp12_mul_sparse_twist_lazy = w2minusv_Fp12_mul_sparse_twist_lazy;
    Fp12_config.Fp12_mul_sparse_untwist_lazy = w2minusv_Fp12_mul_sparse_untwist_lazy;
  } 
  else if (PAIR_CURVE == BN12CP) {
    wpp[3] = (uint64_t) 0x2400005100016456;
    wpp[2] = (uint64_t) 0x09FF9DB1F4ECF9A1;
    wpp[1] = (uint64_t) 0xB461E994A5E5451E;
    wpp[0] = (uint64_t) 0x17903165FFCB801B;
    Fp_set (Fp12_config.wp2pow, wpp);

    wpp[3] = (uint64_t) 0x0AADF027935C8B3D;
    wpp[2] = (uint64_t) 0x4EA45098503276A2;
    wpp[1] = (uint64_t) 0x5E5ABDAA03FF2F8C;
    wpp[0] = (uint64_t) 0x8EBAF6621F84D9C7;
    Fp_set (Fp12_config.wppow->a0, wpp);
  
    wpp[3] = (uint64_t) 0x195210296CA4D919;
    wpp[2] = (uint64_t) 0x035B4D9324BC7755;
    wpp[1] = (uint64_t) 0x17F9F9E68A8C613E;
    wpp[0] = (uint64_t) 0x43C77EC6E054264C;
    Fp_set (Fp12_config.wppow->a1, wpp);

    wpp[3] = (uint64_t) 0x20FAE8ECAA2CAF00;
    wpp[2] = (uint64_t) 0xB5EB7ED253355CA0;
    wpp[1] = (uint64_t) 0x90E06A671E7EB29E;
    wpp[0] = (uint64_t) 0x9F8EFD3F611C0488;
    Fp_set (Fp12_config.wp3pow->a0, wpp);
  
    wpp[3] = (uint64_t) 0x0305176455D4B555;
    wpp[2] = (uint64_t) 0x9C141F5921B99156;
    wpp[1] = (uint64_t) 0xE5744D29700CDE2C;
    wpp[0] = (uint64_t) 0x32F377E99EBCFB8B;
    Fp_set (Fp12_config.wp3pow->a1, wpp);

    Fp_mul (Fp12_config.wp2pow3, Fp12_config.wp2pow, Fp12_config.wp2pow);
    Fp_mul (Fp12_config.wp2pow3, Fp12_config.wp2pow3, Fp12_config.wp2pow);
    Fp2_squ (Fp12_config.wppow3, Fp12_config.wppow);
    Fp2_mul (Fp12_config.wppow3, Fp12_config.wppow3, Fp12_config.wppow);
    Fp2_inv (Fp12_config.wppow3inv, Fp12_config.wppow3);

    Fp12_config.Fp12_mul = w2minusv_Fp12_mul;
    Fp12_config.Fp12_squ = w2minusv_Fp12_squ;
    Fp12_config.Fp12_inv = w2minusv_Fp12_inv;
    Fp12_config.Fp12_ppow = w2minusv_v3minusxi_p1mod6_Fp12_ppow;
    Fp12_config.Fp12_p2pow = w2minusv_v3minusxi_p1mod6_Fp12_p2pow;
    Fp12_config.Fp12_p3pow = w2minusv_v3minusxi_p1mod6_Fp12_p3pow;
    Fp12_config.Fp12_p6pow = w2minusv_Fp12_p6pow; 
    Fp12_config.Fp12_mul_lazy = w2minusv_Fp12_mul_lazy;
    Fp12_config.Fp12_squ_lazy = w2minusv_Fp12_squ_lazy;
    Fp12_config.Fp12_mul_sparse_twist_lazy = w2minusv_Fp12_mul_sparse_twist_lazy;
    Fp12_config.Fp12_mul_sparse_untwist_lazy = w2minusv_Fp12_mul_sparse_untwist_lazy;
  }
  else if (PAIR_CURVE == BN12tiny) {
	  uint64_t wp2pow[1] = { 0xC59AB172BC };
	  uint64_t wppowa0[1] = { 0xF5B82C0AFD83C };
	  uint64_t wppowa1[1] = { 0x66E2DB8BB999F };
	  uint64_t wp3powa0[1] = { 0xCBFF49918E901 };
	  uint64_t wp3powa1[1] = { 0x909BBE05288DA };

	  Fp_set(Fp12_config.wp2pow, wp2pow);
	  Fp_set(Fp12_config.wppow->a0, wppowa0);
	  Fp_set(Fp12_config.wppow->a1, wppowa1);
	  Fp_set(Fp12_config.wp3pow->a0, wp3powa0);
	  Fp_set(Fp12_config.wp3pow->a1, wp3powa1);

	  Fp_mul(Fp12_config.wp2pow3, Fp12_config.wp2pow, Fp12_config.wp2pow);
	  Fp_mul(Fp12_config.wp2pow3, Fp12_config.wp2pow3, Fp12_config.wp2pow);
	  Fp2_squ(Fp12_config.wppow3, Fp12_config.wppow);
	  Fp2_mul(Fp12_config.wppow3, Fp12_config.wppow3, Fp12_config.wppow);
	  Fp2_inv(Fp12_config.wppow3inv, Fp12_config.wppow3);

	  Fp12_config.Fp12_mul = w2minusv_Fp12_mul;
	  Fp12_config.Fp12_squ = w2minusv_Fp12_squ;
	  Fp12_config.Fp12_inv = w2minusv_Fp12_inv;
	  Fp12_config.Fp12_ppow = w2minusv_v3minusxi_p1mod6_Fp12_ppow;
	  Fp12_config.Fp12_p2pow = w2minusv_v3minusxi_p1mod6_Fp12_p2pow;
	  Fp12_config.Fp12_p3pow = w2minusv_v3minusxi_p1mod6_Fp12_p3pow;
	  Fp12_config.Fp12_p6pow = w2minusv_Fp12_p6pow;
	  Fp12_config.Fp12_mul_lazy = w2minusv_Fp12_mul_lazy;
	  Fp12_config.Fp12_squ_lazy = w2minusv_Fp12_squ_lazy;
	  Fp12_config.Fp12_mul_sparse_twist_lazy = w2minusv_Fp12_mul_sparse_twist_lazy;
	  Fp12_config.Fp12_mul_sparse_untwist_lazy = w2minusv_Fp12_mul_sparse_untwist_lazy;
  }
  
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }

  return ERR_SUCCESS;
}

void Fp12_free_config (void) {
  Fp2_free (Fp12_config.wppow);
  Fp2_free (Fp12_config.wppow3);
  Fp_free (Fp12_config.wp2pow);
  Fp_free (Fp12_config.wp2pow3);
  Fp2_free (Fp12_config.wp3pow);
  Fp2_free (Fp12_config.wppow3inv);
  free (Fp12_config.mem);
}

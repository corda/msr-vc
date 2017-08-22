/* Fp6/Fp6overFp3.c
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

#include "Fp6overFp3.h"
#include "../pairing/pairing.h"

Fp6q_config_t Fp6q_config;

/* Initialize an element in Fp6q by initializing two Fp3 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp6q_init(Fp6q_t dst) {
	int err = 0;
	err = Fp3_init(dst->a0);
	if (err < 0) return err;
	return Fp3_init(dst->a1);
}

/* Free the memory allocated in src. */
void Fp6q_free (Fp6q_t src) {
	Fp3_free(src->a0);
	Fp3_free(src->a1);
}

/* Get a random element in Fp6q. */
void Fp6q_rand (Fp6q_t c) {
  Fp3_rand(c->a0);
  Fp3_rand(c->a1);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp6q_print (Fp6q_t a) {
  Fp3_print(a->a0);
  //printf("\n");
  Fp3_print(a->a1);
  //printf("\n");
}

/* Set c to 0. */
void Fp6q_set_zero (Fp6q_t c) {
  Fp3_set_zero (c->a0);
  Fp3_set_zero (c->a1);
}

/* Set c to 1. */
void Fp6q_set_one (Fp6q_t c) {
  Fp3_set_one (c->a0);
  Fp3_set_zero (c->a1);
}

/* Copy an Fp6q element a to an Fp6q element c. */
void Fp6q_copy(Fp6q_t c, Fp6q_t a) {
	Fp3_copy(c->a0, a->a0);
	Fp3_copy(c->a1, a->a1);
}

/* Compare two Fp6q elements a and b for equality. */
int Fp6q_cmpeq(Fp6q_t a, Fp6q_t b) {
  return (Fp3_cmpeq(a->a0, b->a0) && Fp3_cmpeq(a->a1, b->a1));
}

/* Negate an Fp6q element. */
void Fp6q_neg(Fp6q_t c, Fp6q_t a) {
	Fp3_neg(c->a0, a->a0);
	Fp3_neg(c->a1, a->a1);
}

/* Add two Fp6q elements coefficient wise. */
void Fp6q_add (Fp6q_t c, Fp6q_t a, Fp6q_t b) {
	Fp3_add(c->a0, a->a0, b->a0);
	Fp3_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp6q element b from an Fp6q element a coefficient wise. */
void Fp6q_sub (Fp6q_t c, Fp6q_t a, Fp6q_t b) {
	Fp3_sub(c->a0, a->a0, b->a0);
	Fp3_sub(c->a1, a->a1, b->a1);
}

/* Multiply an Fp6q element by an Fp3 element. */
void Fp6q_mulFp3(Fp6q_t c, Fp6q_t a, Fp3_t b) {
	Fp3_mul(c->a0, a->a0, b);
	Fp3_mul(c->a1, a->a1, b);
}

/* Multiply an Fp6q element by an Fp element. */
void Fp6q_mulFp(Fp6q_t c, Fp6q_t a, Fp_t b) {
	Fp3_mulFp(c->a0, a->a0, b);
	Fp3_mulFp(c->a1, a->a1, b);
}

/* Multiply two Fp6q elements. */
void Fp6q_mul (Fp6q_t c, Fp6q_t a, Fp6q_t b) {
  Fp6q_config.Fp6q_mul(c, a, b, Fp6q_config.mem);
}

/* Square an Fp6q element. */
void Fp6q_squ (Fp6q_t c, Fp6q_t a) {
  Fp6q_config.Fp6q_squ(c, a, Fp6q_config.mem);
}

/* Compute the p-power Frobenius of an Fp6q element. */
void Fp6q_ppow(Fp6q_t c, Fp6q_t a) {
  Fp6q_config.Fp6q_ppow(c, a);
}

/* Compute the p^2-power Frobenius of an Fp6q element. */
void Fp6q_p2pow(Fp6q_t c, Fp6q_t a) {
  Fp6q_config.Fp6q_p2pow(c, a);
}

/* Compute the p^3-power Frobenius of an Fp6q element. */
void Fp6q_p3pow(Fp6q_t c, Fp6q_t a) {
  Fp6q_config.Fp6q_p3pow(c, a);
}

/* Invert an Fp6q element. */
void Fp6q_inv (Fp6q_t c, Fp6q_t a) {
  Fp6q_config.Fp6q_inv(c, a, Fp6q_config.mem);
}

/* Multiply two Fp6q elements using lazy reduction.*/
void Fp6q_mul_lazy(Fp6q_t c, Fp6q_t a, Fp6q_t b) {
  Fp6q_config.Fp6q_mul_lazy (c, a, b, Fp6q_config.mem);
}

/* Square an Fp6q element using lazy reduction.*/
void Fp6q_squ_lazy(Fp6q_t c, Fp6q_t a) {
  Fp6q_config.Fp6q_squ_lazy (c, a, Fp6q_config.mem);
}

/* Multiply two Fp6q elements (the second one sparse) using lazy reduction.*/
void Fp6q_mul_sparse_twist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b00, Fp_t b01, Fp_t b11) {
  Fp6q_config.Fp6q_mul_sparse_twist_lazy (c, a, b00, b01, b11, Fp6q_config.mem);
}

/* Multiply two Fp6q elements (the second one sparse) using lazy reduction.*/
void Fp6q_mul_sparse_untwist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b01, Fp_t b11, Fp_t b12) {
  Fp6q_config.Fp6q_mul_sparse_untwist_lazy (c, a, b01, b11, b12, Fp6q_config.mem);
}

/* Multiply two Fp6q elements. This uses that Fp6q is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 20 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void w2minusv_Fp6q_mul(Fp6q_t c, Fp6q_t a, Fp6q_t b, void *mem) {
	Fp3_t t0, t1, t2;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t2->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t2->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	t2->a2->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  
	Fp3_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp3_add(t1, b->a0, b->a1);	// t1 = b0 + b1
	Fp3_mul(t2, t0, t1);			  // t2 = t0*t1 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1

	Fp3_mul(t0, a->a0, b->a0);		// t0 = a0*b0
	Fp3_mul(t1, a->a1, b->a1);		// t1 = a1*b1
	
	Fp3_sub(t2, t2, t0);					// t2 = t2 - t0 = a1*b0 + a0*b1 + a1*b1
	Fp3_sub(c->a1, t2, t1);				// c1 = t2 - t1 = a1*b0 + a0*b1

	Fp3_mulu(t1, t1);						  // t1 = u*t1 = u*a1*b1
	Fp3_add(c->a0, t0, t1);				// c0 = t0 + t1 = a0*b0 + u*a1*b1
}

/* Square a general Fp6q element. This uses that Fp6q is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 22 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void w2minusv_Fp6q_squ(Fp6q_t c, Fp6q_t a, void *mem) {
	Fp3_t t0, t1;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  
	Fp3_sub(t0, a->a0, a->a1);		// t0 = a0 - a1
	Fp3_mulu(t1, a->a1);					// t1 = v*a1
	Fp3_sub(t1, a->a0, t1);				// t1 = a0 - v*a1
	Fp3_mul(t0, t0, t1);			    // t0 = (a0 - a1)*(a0 - v*a1) = a0^2 + v*a1^2 - a0*a1 - v*a0*a1
	Fp3_mul(t1, a->a0, a->a1);		// t1 = a0*a1
	Fp3_add(t0, t0, t1);					// t0 =  a0^2 + v*a1^2 - v*a0*a1
	Fp3_add(c->a1, t1, t1);			  // c1 = 2*t1 = 2*a0*a1
	Fp3_mulu(t1, t1);						  // t1 = v*t1 = v*a0*a1
	Fp3_add(c->a0, t0, t1);				// c0 = t0 + t1 = a0^2 + v*a1^2
}

/* Compute the p-power Frobenius of an Fp6q element. This uses a constant WPPOW, 
 * such that w^p = WPPOW*w, i.e. WPPOW = w^(p-1). In the case, where Fp6 is an 
 * extension of Fp via v^3 - xi and Fp6q is an extension of Fp3 via w^2 - v, and p = 1 (mod 6),
 * we have WPPOW = xi^(p-1)/6 in Fp.
 */
void w2minusv_v3minusxi_p1mod6_Fp6q_ppow(Fp6q_t c, Fp6q_t a) {
	Fp3_ppow(c->a0, a->a0);
	Fp3_ppow(c->a1, a->a1);
	Fp3_mulFp(c->a1, c->a1, Fp6q_config.wppow);
}

/* Compute the p^2-power Frobenius of an Fp6q element. This uses a constant WP2POW, 
 * such that w^(p^2) = WP2POW*w, i.e. WP2POW = w^(p^2-1). In the case, where Fp3 is an 
 * extension of Fp via v^3 - xi and Fp6q is an extension of Fp3 via w^2 - v, and p = 1 (mod 6),
 * we have WP2POW = xi^(p^2-1)/6 in Fp.
 */
void w2minusv_v3minusxi_p1mod6_Fp6q_p2pow(Fp6q_t c, Fp6q_t a) {
	Fp3_ppow(c->a0, a->a0);
	Fp3_ppow(c->a0, c->a0);
	Fp3_ppow(c->a1, a->a1);
	Fp3_ppow(c->a1, c->a1);
	Fp3_mulFp(c->a1, c->a1, Fp6q_config.wp2pow);
}

/* Compute the p^3-power Frobenius of an Fp6q element. This is the same as 
 * conjugation in the quadratic extension over Fp3 if the extension is constructed
 * via a quadratic polynomial x^2 - alpha, alpha in Fp3. 
 */
void w2minusv_Fp6q_p3pow(Fp6q_t c, Fp6q_t a) {
	Fp3_copy(c->a0, a->a0);
	Fp3_neg(c->a1, a->a1);
}

/* Invert an Fp6q element. This uses that Fp6q is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 23 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void w2minusv_Fp6q_inv(Fp6q_t c, Fp6q_t a, void *mem) {
	Fp3_t t0, t1;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;

	Fp3_squ(t0, a->a0);				  // t0 = a0^2
	Fp3_squ(t1, a->a1);				  // t1 = a1^2
	Fp3_mulu(t1, t1);						// t1 = v*t1 = v*a1^2
	Fp3_sub(t0, t0, t1);		  	// t0 = a0^2 - v*a1^2 = Norm(a)
  Fp3_inv(t1, t0);				    // t1 = t0^(-1) = Norm(a)^(-1)
	Fp3_mul(c->a0, a->a0, t1);	// c0 = a0/(a0^2 - v*a1^2)
	Fp3_neg(t1, t1);						// t1 = -1/(a0^2 - v*a1^2)
	Fp3_mul(c->a1, a->a1, t1);	// c1 = a1*t1 = -a1/(a0^2 - v*a1^2)
}


/* Multiply two Fp6q elements using lazy reduction techniques. 
 * This uses that Fp6q is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. This is Algorithm 4 in Aranha et al. 
 */
void w2minusv_Fp6q_mul_lazy(Fp6q_t c, Fp6q_t a, Fp6q_t b, void *mem) {
	Fp3_t t0, t1;
  Fp3l_t T0, T1, T2, T3;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T0->A0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
	T0->A2->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T1->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T1->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T1->A2->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
	T2->A0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
	T2->A2->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 24*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 26*Fp_config.m->n;
	T3->A2->limbs = (uint64_t *) mem + 28*Fp_config.m->n;

  Fp3_mul_no_red (T0, a->a0, b->a0);
  Fp3_mul_no_red (T1, a->a1, b->a1);
  Fp3_add (t0, a->a0, a->a1);
  Fp3_add (t1, b->a0, b->a1);
  Fp3_mul_no_red (T2, t0, t1);
  Fp3l_add_o2 (T3, T0, T1);  
  Fp3l_sub_o2 (T2, T2, T3);
  Fp3l_red (c->a1, T2);

  Fp3l_mulu (T2, T1);
  Fp3l_add_o2 (T2, T0, T2);
  Fp3l_red (c->a0, T2);
}

/* Multiply two Fp6q elements using lazy reduction techniques, 
 * where the second element b is sparse in the following sense:
 * b = b0 + b1*w, where b0 = b00 + b01*v and b1 = b11*v.
 * This uses that Fp6q is constructed over Fp3 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. 
 */
void w2minusv_Fp6q_mul_sparse_twist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b00, Fp_t b01, Fp_t b11, void *mem) {
	Fp3_t t0, t1;
  Fp3l_t T0, T1, T2, T3;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T0->A0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
	T0->A2->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T1->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T1->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T1->A2->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
	T2->A0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
	T2->A2->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 24*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 26*Fp_config.m->n;
	T3->A2->limbs = (uint64_t *) mem + 28*Fp_config.m->n;

  Fp3_mul_sparse01_no_red (T0, a->a0, b00, b01); // b0 = b00 + b01*v
  Fp3_mulu (t0, a->a1);
  Fp3_mulFp_no_red (T1, t0, b11);
  Fp3_add (t0, a->a0, a->a1);
  Fp_add (t1->a1, b01, b11);
  Fp3_mul_sparse01_no_red (T2, t0, b00, t1->a1); // t1 = b00 + (b01 + b11)*v
  Fp3l_add_o2 (T3, T0, T1);  
  Fp3l_sub_o2 (T2, T2, T3);
  Fp3l_red (c->a1, T2);

  Fp3l_mulu (T2, T1);
  Fp3l_add_o2 (T2, T0, T2);
  Fp3l_red (c->a0, T2);
}

/* Similar as above, but now b = b0 + b1*w, where b0 = b01*v and b1 = b11*v + b12*v^2. */
void w2minusv_Fp6q_mul_sparse_untwist_lazy(Fp6q_t c, Fp6q_t a, Fp_t b01, Fp_t b11, Fp_t b12, void *mem) {
	Fp3_t t0, t1;
  Fp3l_t T0, T1, T2, T3;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T0->A0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
	T0->A2->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T1->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T1->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
	T1->A2->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
	T2->A0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
	T2->A2->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 24*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 26*Fp_config.m->n;
	T3->A2->limbs = (uint64_t *) mem + 28*Fp_config.m->n;


  Fp3_mulu (t0, a->a0);
  Fp3_mulFp_no_red (T0, t0, b01);
  Fp3_mul_sparse12_no_red (T1, a->a1, b11, b12); 
  Fp3_add_no_red (t0, a->a0, a->a1);
  Fp_add_no_red (t1->a1, b01, b11);
  Fp3_mul_sparse12_no_red (T2, t0, t1->a1, b12); // t1 = (b01 + b11)*v + b12*v^2
  Fp3l_add_o2 (T3, T0, T1);  
  Fp3l_sub_o2 (T2, T2, T3);
  Fp3l_red (c->a1, T2);

  Fp3l_mulu (T2, T1);
  Fp3l_add_o2 (T2, T0, T2);
  Fp3l_red (c->a0, T2);
}


/* Square an Fp6q element using lazy reduction techniques. 
 * This uses that Fp6q is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. This is Algorithm 5 in Aranha et al. 
 */
void w2minusv_Fp6q_squ_lazy(Fp6q_t c, Fp6q_t a, void *mem) {
	Fp3_t t0, t1;
  Fp3l_t T0;

  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	t0->a2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	t1->a2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T0->A0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
	T0->A2->limbs = (uint64_t *) mem + 10*Fp_config.m->n;

  Fp3_add (t0, a->a0, a->a1);
  Fp3_mulu (t1, a->a1);
  Fp3_add (t1, a->a0, t1);
  Fp3_mul_no_red (T0, a->a0, a->a1);
  Fp3l_red (c->a1, T0);               // c1 = a0*a1
  
  Fp3_mul_no_red (T0, t0, t1);
  Fp3l_red (t0, T0);
  Fp3_mulu (t1, c->a1);
  Fp3_sub (c->a0, t0, t1);
  Fp3_sub (c->a0, c->a0, c->a1);      // c0 = a0^2 + v*a1^2
  Fp3_add (c->a1, c->a1, c->a1);      // c1 = 2*a0*a1
}

int Fp6q_initialize_config (void) {
  uint64_t *wpp;

  Fp6q_config.mem = NULL;
  Fp6q_config.mem = (void *) malloc (60 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp6q_config.mem == NULL) return ERR_OUT_OF_MEMORY;

  if (Fp_init (Fp6q_config.wp2pow) < 0) printf ("Fp6q_config memory error.\n");
  if (Fp_init (Fp6q_config.wp2pow3) < 0) printf ("Fp6q_config memory error.\n");
  if (Fp_init (Fp6q_config.wppow) < 0) printf ("Fp6q_config memory error.\n");
  if (Fp_init (Fp6q_config.wp3pow) < 0) printf ("Fp6q_config memory error.\n");
  if (Fp_init (Fp6q_config.wppow3) < 0) printf ("Fp6q_config memory error.\n");
  if (Fp_init (Fp6q_config.wppow3inv) < 0) printf ("Fp6q_config memory error.\n");
  wpp = (uint64_t *) Fp6q_config.mem;

  if (PAIR_CURVE == CP6) {
    uint64_t wppow[8] =  { 0x5257C718F4B582E2,0x0C1197483FBEF0DD,
                           0x74EDC8603F8B0275,0xBB686D34F710DCBC,
                           0xAFAADB6265FA73C3,0x262A5F326F028EA5,
                           0x33B162CE50A3D98C,0x05100016C8007DD9 };
    uint64_t wp2pow[8] = { 0x5257C718F4B582E1,0x0C1197483FBEF0DD,
                           0x74EDC8603F8B0275,0xBB686D34F710DCBC,
                           0xAFAADB6265FA73C3,0x262A5F326F028EA5,
                           0x33B162CE50A3D98C,0x05100016C8007DD9 };
    Fp_set (Fp6q_config.wppow, wppow);
    Fp_set (Fp6q_config.wp2pow, wp2pow);

    //Fp_mul (Fp6q_config.wp2pow3, Fp6q_config.wp2pow, Fp6q_config.wp2pow);
    //Fp_mul (Fp6q_config.wp2pow3, Fp6q_config.wp2pow3, Fp6q_config.wp2pow);
    //Fp_mul (Fp6q_config.wppow3, Fp6q_config.wppow, Fp6q_config.wppow);
    //Fp_mul (Fp6q_config.wppow3, Fp6q_config.wppow3, Fp6q_config.wppow);
    
    Fp6q_config.Fp6q_mul = w2minusv_Fp6q_mul;
    Fp6q_config.Fp6q_squ = w2minusv_Fp6q_squ;
    Fp6q_config.Fp6q_inv = w2minusv_Fp6q_inv;
    Fp6q_config.Fp6q_ppow = w2minusv_v3minusxi_p1mod6_Fp6q_ppow;
    Fp6q_config.Fp6q_p2pow = w2minusv_v3minusxi_p1mod6_Fp6q_p2pow;
    Fp6q_config.Fp6q_p3pow = w2minusv_Fp6q_p3pow; 
    Fp6q_config.Fp6q_mul_lazy = w2minusv_Fp6q_mul_lazy;
    Fp6q_config.Fp6q_squ_lazy = w2minusv_Fp6q_squ_lazy;
    Fp6q_config.Fp6q_mul_sparse_twist_lazy = w2minusv_Fp6q_mul_sparse_twist_lazy;
    Fp6q_config.Fp6q_mul_sparse_untwist_lazy = w2minusv_Fp6q_mul_sparse_untwist_lazy;
  } else if (PAIR_CURVE == CP6b) {
    uint64_t wppow[8] =  { 0x8C000000000004C4,0x6BF39000000009D9,
                           0x9D35F6980000093D,0xBCFC9A74F8000520,
                           0xA416A709DA2C01D6,0x1833E371C780246F,
                           0xBB99497657F45361,0x1029BFD6DCBA8CD0 };
    uint64_t wp2pow[8] = { 0x8C000000000004C3,0x6BF39000000009D9,
                           0x9D35F6980000093D,0xBCFC9A74F8000520,
                           0xA416A709DA2C01D6,0x1833E371C780246F,
                           0xBB99497657F45361,0x1029BFD6DCBA8CD0 };
    Fp_set (Fp6q_config.wppow, wppow);
    Fp_set (Fp6q_config.wp2pow, wp2pow);

    //Fp_mul (Fp6q_config.wp2pow3, Fp6q_config.wp2pow, Fp6q_config.wp2pow);
    //Fp_mul (Fp6q_config.wp2pow3, Fp6q_config.wp2pow3, Fp6q_config.wp2pow);
    //Fp_mul (Fp6q_config.wppow3, Fp6q_config.wppow, Fp6q_config.wppow);
    //Fp_mul (Fp6q_config.wppow3, Fp6q_config.wppow3, Fp6q_config.wppow);
    
    Fp6q_config.Fp6q_mul = w2minusv_Fp6q_mul;
    Fp6q_config.Fp6q_squ = w2minusv_Fp6q_squ;
    Fp6q_config.Fp6q_inv = w2minusv_Fp6q_inv;
    Fp6q_config.Fp6q_ppow = w2minusv_v3minusxi_p1mod6_Fp6q_ppow;
    Fp6q_config.Fp6q_p2pow = w2minusv_v3minusxi_p1mod6_Fp6q_p2pow;
    Fp6q_config.Fp6q_p3pow = w2minusv_Fp6q_p3pow; 
    Fp6q_config.Fp6q_mul_lazy = w2minusv_Fp6q_mul_lazy;
    Fp6q_config.Fp6q_squ_lazy = w2minusv_Fp6q_squ_lazy;
    Fp6q_config.Fp6q_mul_sparse_twist_lazy = w2minusv_Fp6q_mul_sparse_twist_lazy;
    Fp6q_config.Fp6q_mul_sparse_untwist_lazy = w2minusv_Fp6q_mul_sparse_untwist_lazy;
  }
  else if (PAIR_CURVE == CP6tiny) {
	  uint64_t wppow[2] = { 0xDB70B1732C9457D8, 0x59032611135 };
	  uint64_t wp2pow[2] = { 0xDB70B1732C9457D7, 0x59032611135 };
	  Fp_set(Fp6q_config.wppow, wppow);
	  Fp_set(Fp6q_config.wp2pow, wp2pow);

	  //Fp_mul (Fp6q_config.wp2pow3, Fp6q_config.wp2pow, Fp6q_config.wp2pow);
	  //Fp_mul (Fp6q_config.wp2pow3, Fp6q_config.wp2pow3, Fp6q_config.wp2pow);
	  //Fp_mul (Fp6q_config.wppow3, Fp6q_config.wppow, Fp6q_config.wppow);
	  //Fp_mul (Fp6q_config.wppow3, Fp6q_config.wppow3, Fp6q_config.wppow);

	  Fp6q_config.Fp6q_mul = w2minusv_Fp6q_mul;
	  Fp6q_config.Fp6q_squ = w2minusv_Fp6q_squ;
	  Fp6q_config.Fp6q_inv = w2minusv_Fp6q_inv;
	  Fp6q_config.Fp6q_ppow = w2minusv_v3minusxi_p1mod6_Fp6q_ppow;
	  Fp6q_config.Fp6q_p2pow = w2minusv_v3minusxi_p1mod6_Fp6q_p2pow;
	  Fp6q_config.Fp6q_p3pow = w2minusv_Fp6q_p3pow;
	  Fp6q_config.Fp6q_mul_lazy = w2minusv_Fp6q_mul_lazy;
	  Fp6q_config.Fp6q_squ_lazy = w2minusv_Fp6q_squ_lazy;
	  Fp6q_config.Fp6q_mul_sparse_twist_lazy = w2minusv_Fp6q_mul_sparse_twist_lazy;
	  Fp6q_config.Fp6q_mul_sparse_untwist_lazy = w2minusv_Fp6q_mul_sparse_untwist_lazy;
  }
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }

  return ERR_SUCCESS;
}

void Fp6q_free_config (void) {
  Fp_free (Fp6q_config.wppow);
  Fp_free (Fp6q_config.wppow3);
  Fp_free (Fp6q_config.wp2pow);
  Fp_free (Fp6q_config.wp2pow3);
  Fp_free (Fp6q_config.wp3pow);
  Fp_free (Fp6q_config.wppow3inv);
  free (Fp6q_config.mem);
}

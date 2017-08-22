/* Fp6/Fp6.c
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

#include "Fp6.h"
#include "../pairing/pairing.h"

Fp6_config_t Fp6_config;

/* Initialize an element in Fp6 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp6_init(Fp6_t dst) {
	int err = 0;
	err = Fp2_init(dst->a0);
	if (err < 0) return err;
	err = Fp2_init(dst->a1);
	if (err < 0) return err;
	return Fp2_init(dst->a2);
}

/* Free the memory allocated in src. */
void Fp6_free (Fp6_t src) {
	Fp2_free(src->a0);
	Fp2_free(src->a1);
	Fp2_free(src->a2);
}

/* Get a random element in Fp6. */
void Fp6_rand (Fp6_t c) {
  Fp2_rand(c->a0);
  Fp2_rand(c->a1);
  Fp2_rand(c->a2);
}

/* Set c to 0. */
void Fp6_set_zero (Fp6_t c) {
  Fp2_set_zero (c->a0);
  Fp2_set_zero (c->a1);
  Fp2_set_zero (c->a2);
}

/* Set c to 1. */
void Fp6_set_one (Fp6_t c) {
  Fp2_set_one (c->a0);
  Fp2_set_zero (c->a1);
  Fp2_set_zero (c->a2);
}

/* Copy an Fp6 element a to an Fp6 element c. */
void Fp6_copy(Fp6_t c, Fp6_t a) {
	Fp2_copy(c->a0, a->a0);
	Fp2_copy(c->a1, a->a1);
	Fp2_copy(c->a2, a->a2);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp6_print (Fp6_t a) {
  Fp2_print(a->a0);
  printf("\n");
  Fp2_print(a->a1);
  printf("\n");
  Fp2_print(a->a2);
  printf("\n");
}

/* Compare two Fp6 elements a and b for equality. */
int Fp6_cmpeq(Fp6_t a, Fp6_t b) {
  return (Fp2_cmpeq(a->a0, b->a0) 
       && Fp2_cmpeq(a->a1, b->a1) && Fp2_cmpeq(a->a2, b->a2));
}

/* Negate an Fp6 element. */
void Fp6_neg(Fp6_t c, Fp6_t a) {
	Fp2_neg(c->a0, a->a0);
	Fp2_neg(c->a1, a->a1);
	Fp2_neg(c->a2, a->a2);
}

/* Add two Fp6 elements coefficient wise. */
void Fp6_add (Fp6_t c, Fp6_t a, Fp6_t b) {
	Fp2_add(c->a0, a->a0, b->a0);
	Fp2_add(c->a1, a->a1, b->a1);
	Fp2_add(c->a2, a->a2, b->a2);
}

/* Subtract an Fp6 element b from an Fp6 element a coefficient wise. */
void Fp6_sub (Fp6_t c, Fp6_t a, Fp6_t b) {
	Fp2_sub(c->a0, a->a0, b->a0);
	Fp2_sub(c->a1, a->a1, b->a1);
	Fp2_sub(c->a2, a->a2, b->a2);
}

/* Multiply an Fp6 element by an Fp2 element. */
void Fp6_mulFp2(Fp6_t c, Fp6_t a, Fp2_t b) {
	Fp2_mul(c->a0, a->a0, b);
	Fp2_mul(c->a1, a->a1, b);
	Fp2_mul(c->a2, a->a2, b);
}

/* Multiply an Fp6 element by an Fp element. */
void Fp6_mulFp(Fp6_t c, Fp6_t a, Fp_t b) {
	Fp2_mulFp(c->a0, a->a0, b);
	Fp2_mulFp(c->a1, a->a1, b);
	Fp2_mulFp(c->a2, a->a2, b);
}

/* Multiply two Fp6 elements. */
void Fp6_mul (Fp6_t c, Fp6_t a, Fp6_t b) {
  Fp6_config.Fp6_mul(c, a, b, Fp6_config.mem);
}

/* Square an Fp6 element. */
void Fp6_squ (Fp6_t c, Fp6_t a) {
  Fp6_config.Fp6_squ(c, a, Fp6_config.mem);
}

/* Invert an Fp6 element. */
void Fp6_inv (Fp6_t c, Fp6_t a) {
  Fp6_config.Fp6_inv(c, a, Fp6_config.mem);
}

/* Multiply an fp6 element by the special element v. */
void Fp6_mulv(Fp6_t c, Fp6_t a) {
  Fp6_config.Fp6_mulv(c, a, Fp6_config.mem);
}

/* Compute the p-power Frobenius of an Fp6 element. */
void Fp6_ppow(Fp6_t c, Fp6_t a) {
  Fp6_config.Fp6_ppow(c, a);
}

/* Compute the p^2-power Frobenius of an Fp6 element. */
void Fp6_p2pow(Fp6_t c, Fp6_t a) {
  Fp6_config.Fp6_p2pow(c, a);
}

/* Compute the p^3-power Frobenius of an Fp6 element. */
void Fp6_p3pow(Fp6_t c, Fp6_t a) {
  Fp6_config.Fp6_p3pow(c, a);
}


/* Multiply two Fp6 elements. This uses that Fp6 is constructed over Fp2 via v^3 - xi. 
 * This is Algorithm 13 in Beuchat et al., Pairing 2010.
 */
void v3minusxi_Fp6_mul(Fp6_t c, Fp6_t a, Fp6_t b, void *mem) {
	Fp2_t t0, t1, t2, t3, t4, t5;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t3->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t3->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t4->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t4->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t5->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t5->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;

	Fp2_mul(t0, a->a0, b->a0);    // t0 = a0*b0
	Fp2_mul(t1, a->a1, b->a1);    // t1 = a1*b1
	Fp2_mul(t2, a->a2, b->a2);    // t2 = a2*b2

	Fp2_add(t3, a->a1, a->a2);    // t3 = a1 + a2
	Fp2_add(t4, b->a1, b->a2);    // t4 = b1 + b2
	Fp2_mul(t3, t3, t4);          // t3 = (a1 + a2)*(b1 + b2) = a1*b1 + a2*b1 + a1*b2 + a2*b2
	Fp2_sub(t3, t3, t1);          // t3 = a2*b1 + a1*b2 + a2*b2
	Fp2_sub(t3, t3, t2);          // t3 = a2*b1 + a1*b2 
  Fp2_mulxi(t3, t3);					  // t3 = (a2*b1 + a1*b2)*xi	
                                // t3 now almost has the result for c0, t4 is free.		

	Fp2_add(t4, a->a0, a->a2);    // t4 = a0 + a2
	Fp2_add(t5, b->a0, b->a2);    // t5 = b0 + b2
	Fp2_mul(t4, t4, t5);          // t4 = (a0 + a2)*(b0 + b2) = a0*b0 + a2*b0 + a0*b2 + a2*b2
	Fp2_sub(t4, t4, t0);          // t4 = a2*b0 + a0*b2 + a2*b2
	Fp2_sub(t4, t4, t2);          // t4 = a2*b0 + a0*b2
	Fp2_add(c->a2, t4, t1);				// c2 = a2*b0 + a0*b2 + a1*b1 	
                                // c2 is done, t4 and t5 are free.

	Fp2_add(t4, a->a0, a->a1);    // t4 = a0 + a1
	Fp2_add(t5, b->a0, b->a1);    // t4 = b0 + b1
	Fp2_mul(t4, t4, t5);          // t4 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1
	Fp2_sub(t4, t4, t0);          // t4 = a1*b0 + a0*b1 + a1*b1
	Fp2_sub(t4, t4, t1);          // t4 = a1*b0 + a0*b1
	Fp2_mulxi(t2, t2);            // t2 = a2*b2 * xi
	Fp2_add(c->a1, t4, t2);				// c1 = a1*b0 + a0*b1 + a2*b2*xi	
                                // c1 is done.

	Fp2_add(c->a0, t3, t0);					// c0 is done.
}

/* Square an Fp6 element. This uses that Fp6 is constructed over Fp2 via v^3 - xi. 
 * This is Algorithm 16 in Beuchat et al., Pairing 2010.
 */
void v3minusxi_Fp6_squ(Fp6_t c, Fp6_t a, void *mem) {
	Fp2_t t0, t1, t2, t3, t4;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t3->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t3->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t4->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t4->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
	
	Fp2_mul(t0, a->a0, a->a1);    // t0 = a0*a1  
	Fp2_add(t0, t0, t0);          // t0 = 2*a0*a1
	Fp2_squ(t1, a->a2);           // t1 = a2^2
	Fp2_squ(t2, a->a0);           // t2 = a0^2
	Fp2_sub(t3, a->a0, a->a1);    // t3 = a0 - a1
	Fp2_add(t3, t3, a->a2);       // t3 = a0 - a1 + a2
	Fp2_mul(t4, a->a1, a->a2);    // t4 = a1*a2
	Fp2_add(t4, t4, t4);          // t4 = 2*a1*a2

	Fp2_sub(c->a2, t0, t1);       // c2 = 2*a0*a1 - a2^2
	Fp2_mulxi(t1, t1);            // t1 = xi*a2^2
	Fp2_add(c->a1, t0, t1);       // c1 = 2*a0*a1 + xi*a2^2

	Fp2_squ(t0, t3);              // t0 = (a0 - a1 + a2)^2
	Fp2_add(c->a2, c->a2, t0);    // c2 = a0^2 + a1^2 - 2*a1*a2 + 2*a0*a2
	Fp2_add(c->a2, c->a2, t4);    // c2 = a0^2 + a1^2 + 2*a0*a2
	Fp2_sub(c->a2, c->a2, t2);    // c2 = a1^2 + 2*a0*a2

	Fp2_mulxi(t4, t4);            // t4 = 2*a1*a2*xi
	Fp2_add(c->a0, t4, t2);       // c0 = a0^2 + 2*a1*a2*xi
}

/* Multiplication by the special element v in Fp6, with v^3 = xi in Fp2.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = a2*xi + a0*v + a1*v^2.
 */
void v3minusxi_Fp6_mulv(Fp6_t c, Fp6_t a, void *mem) {
  Fp2_t t0;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	
  Fp2_mulxi(t0, a->a2);
  Fp2_copy(c->a2, a->a1);
	Fp2_copy(c->a1, a->a0);
  Fp2_copy(c->a0, t0);
}

/* Compute the p-power Frobenius of an Fp6 element. This uses a constant VPPOW, 
 * such that v^p = VPPOW*v, i.e. VPPOW = v^(p-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and p = 1 (mod 3), we have VPPOW = xi^(p-1)/3 in Fp2. 
 * VPPOW2 = VPPOW^2.
 */
void p1mod3_v3minusxi_Fp6_ppow(Fp6_t c, Fp6_t a) {
  Fp2_ppow(c->a0, a->a0);
	Fp2_ppow(c->a1, a->a1);
	Fp2_ppow(c->a2, a->a2);
	Fp2_muliFp(c->a1, c->a1, Fp6_config.vppow);
	Fp2_mulFp(c->a2, c->a2, Fp6_config.vppow2);
}

/* Compute the p^2-power Frobenius of an Fp6 element. As above, we use a constant
 * VP2POW such that v^(p^2) = VP2POW*v, i.e. VP2POW = v^(p^2-1).  In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and p = 1 (mod 3), we have VP2POW = xi^(p^2-1)/3 is a third root
 * of unity in Fp. VP2POW2 = VP2POW^2.  
 */
void p1mod3_v3minusxi_Fp6_p2pow(Fp6_t c, Fp6_t a) {
	Fp2_copy(c->a0, a->a0);
	Fp2_mulFp(c->a1, a->a1, Fp6_config.vp2pow);
	Fp2_mulFp(c->a2, a->a2, Fp6_config.vp2pow2);
}

/* Compute the p^3-power Frobenius of an Fp6 element. As above, we use a constant
 * VP3POW such that v^(p^3) = VP3POW*v, i.e. VP3POW = v^(p^3-1).  In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and p = 1 (mod 3), we have VP3POW = xi^(p^3-1)/3 in Fp2. 
 * VP3POW2 = VP3POW^2. In the BN case, VP3POW = i and VP3POW2 = -1. No need for a constant here.
 */
void p1mod3_v3minusxi_Fp6_p3pow(Fp6_t c, Fp6_t a) {
	Fp2_ppow(c->a0, a->a0);
	Fp2_ppow(c->a1, a->a1);
	Fp2_ppow(c->a2, a->a2);
	Fp2_muli(c->a1, c->a1);
	Fp2_neg(c->a2, c->a2);
}

/* Inversion of an Fp6 element. This uses that Fp6 is constructed over Fp2 via v^3 - xi. 
 * This is Algorithm 17 in Beuchat et al., Pairing 2010.
 */
void v3minusxi_Fp6_inv(Fp6_t c, Fp6_t a, void *mem) {
	Fp2_t t0, t1, t2, t3, t4, t5;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t2->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t2->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t3->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t3->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t4->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t4->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t5->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t5->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;

	Fp2_squ(t0, a->a0);         // t0 = a0^2
	Fp2_squ(t1, a->a1);         // t1 = a1^2
	Fp2_squ(t2, a->a2);         // t2 = a2^2

	Fp2_mul(t3, a->a0, a->a1);  // t3 = a0*a1
	Fp2_mul(t4, a->a0, a->a2);  // t4 = a0*a2
	Fp2_mul(t5, a->a1, a->a2);  // t5 = a1*a2

	Fp2_mulxi(t2, t2);          // t2 = xi*a2^2
	Fp2_mulxi(t5, t5);          // t5 = xi*a1*a2

	Fp2_sub(t0, t0, t5);        // t0 = a0^2 - xi*a1*a2
	Fp2_sub(t1, t1, t4);        // t1 = a1^2 - a0*a2
	Fp2_sub(t2, t2, t3);        // t2 = xi*a2^2 - a0*a1

	Fp2_mul(t3, a->a0, t0);     // t3 = a0^3 - xi*a0*a1*a2
	Fp2_mul(t4, a->a2, t2);     // t4 = xi*a2^3 - a0*a1*a2
	Fp2_mulxi(t4, t4);          // t4 = xi^2*a2^3 - xi*a0*a1*a2
	Fp2_add(t3, t3, t4);        // t3 = a0^3 + xi^2*a2^3 - 2*xi*a0*a1*a2 
	Fp2_mul(t4, a->a1, t1);     // t4 = a1^3 - a0*a1*a2
	Fp2_mulxi(t4, t4);          // t4 = xi*a1^3 - xi*a0*a1*a2
	Fp2_add(t3, t3, t4);        // t3 = a0^3 + xi*a1^3 + xi^2*a2^3 - 3*xi*a0*a1*a2

	Fp2_inv(t3, t3);            // t3 = 1/N(a)

	Fp2_mul(c->a0, t0, t3);
	Fp2_mul(c->a1, t2, t3);
	Fp2_mul(c->a2, t1, t3);
}

/************** Fp6 functions for lazy reduction ***************/
/* Functions assume that there is enough space such that no overflow occurs. */

/* Initialize a long Fp6 element. */
int Fp6l_init(Fp6l_t dst) {
	int err = 0;
	if (err = Fp2l_init(dst->A0) < 0) return err;
	if (err = Fp2l_init(dst->A1) < 0) return err;
	return Fp2l_init(dst->A2);
}

/* Free the memory allocated in src. */
void Fp6l_free (Fp6l_t src) {
	Fp2l_free(src->A0);
	Fp2l_free(src->A1);
	Fp2l_free(src->A2);
}

/* Copy an Fp6l element A to an Fp6l element C. */
void Fp6l_copy(Fp6l_t C, Fp6l_t A) {
	Fp2l_copy(C->A0, A->A0);
	Fp2l_copy(C->A1, A->A1);
	Fp2l_copy(C->A2, A->A2);
}

/* Reduction function calling the Fp2 reduction for each coefficient. */
void Fp6l_red(Fp6_t c, Fp6l_t A) {
	Fp2l_red(c->a0, A->A0);
	Fp2l_red(c->a1, A->A1);
	Fp2l_red(c->a2, A->A2);
}

/* Subtraction of Fp6l elements, constant-time version of option 2 in Aranha et al. */
void Fp6l_sub_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b) {
  Fp2l_sub_o2(c->A0, a->A0, b->A0);
  Fp2l_sub_o2(c->A1, a->A1, b->A1);
  Fp2l_sub_o2(c->A2, a->A2, b->A2);
}

/* Addition of Fp6l elements, constant-time version of "option 2" addition in Aranha et al. */
void Fp6l_add_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b) {
  Fp2l_add_o2(c->A0, a->A0, b->A0);
  Fp2l_add_o2(c->A1, a->A1, b->A1);
  Fp2l_add_o2(c->A2, a->A2, b->A2);
}

/* Multiply an Fp6 element by an Fp element. */
void Fp6_mulFp2_no_red (Fp6l_t c, Fp6_t a, Fp2_t b) {
	Fp2_mul_no_red_o2 (c->A0, a->a0, b);
	Fp2_mul_no_red_o2 (c->A1, a->a1, b);
  Fp2_mul_no_red_o2 (c->A2, a->a2, b);
}

/* Multiply two Fp6 elements without reduction. */
void Fp6_mul_no_red (Fp6l_t c, Fp6_t a, Fp6_t b) {
  Fp6_config.Fp6_mul_no_red (c, a, b, Fp6_config.mem);
}

/* Multiply an Fp6 element and a sparse Fp6 element without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void Fp6_mul_sparse01_no_red (Fp6l_t c, Fp6_t a, Fp2_t b0, Fp2_t b1){
  Fp6_config.Fp6_mul_sparse01_no_red (c, a, b0, b1, Fp6_config.mem);
}
void Fp6_mul_sparse12_no_red (Fp6l_t c, Fp6_t a, Fp2_t b1, Fp2_t b2){
  Fp6_config.Fp6_mul_sparse12_no_red (c, a, b1, b2, Fp6_config.mem);
}

/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * This is Algorithm 3 in Aranha et al.
 */
void v3minusxi_Fp6_mul_no_red (Fp6l_t c, Fp6_t a, Fp6_t b, void *mem) {
  Fp2_t t0, t1;
  Fp2l_t T0, T1, T2, T3, T4;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  
  T0->A0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T1->A0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T1->A1->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T2->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  T4->A0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  T4->A1->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  
  Fp2_mul_no_red_o1 (T0, a->a0, b->a0, 2);
  Fp2_mul_no_red_o1 (T1, a->a1, b->a1, 2);
  Fp2_mul_no_red_o1 (T2, a->a2, b->a2, 2);

  Fp2_add_no_red (t0, a->a1, a->a2);
  Fp2_add_no_red (t1, b->a1, b->a2);
  Fp2_mul_no_red_o2 (T3, t0, t1);
  Fp2l_add_no_red (T4, T1, T2);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T4->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T4->A1);
  Fpl_sub_o2_ct (T4->A0, T3->A0, T3->A1, Fp_config.mem);
  Fpl_add_o2_ct (T4->A1, T3->A0, T3->A1, Fp_config.mem);  
  Fp2l_add_o2 (c->A0, T4, T0);

  Fp2_add_no_red (t0, a->a0, a->a1);
  Fp2_add_no_red (t1, b->a0, b->a1);
  Fp2_mul_no_red_o2 (T3, t0, t1);
  Fp2l_add_no_red (T4, T0, T1);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T4->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T4->A1);
  Fpl_sub_o1_ct (T4->A0, T2->A0, T2->A1, 1, Fp_config.mem);
  Fpl_add_no_red (T4->A1, T2->A0, T2->A1);
  Fp2l_add_o2 (c->A1, T3, T4);

  Fp2_add_no_red (t0, a->a0, a->a2);
  Fp2_add_no_red (t1, b->a0, b->a2);
  Fp2_mul_no_red_o2 (T3, t0, t1);
  Fp2l_add_no_red (T4, T0, T2);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T4->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T4->A1);
  Fpl_add_o2_ct (c->A2->A0, T3->A0, T1->A0, Fp_config.mem); 
  Fpl_add_no_red (c->A2->A1, T3->A1, T1->A1);
}


/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * The third coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b0 + b1*v, i.e. b2 = 0.
 */
void v3minusxi_Fp6_mul_sparse01_no_red (Fp6l_t c, Fp6_t a, Fp2_t b0, Fp2_t b1, void *mem) {
  Fp2_t t0, t1;
  Fp2l_t T0, T1, T2, T3;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  
  T0->A0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T1->A0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T1->A1->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T2->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  
  Fp2_mul_no_red_o1 (T0, a->a0, b0, 2);
  Fp2_mul_no_red_o1 (T1, a->a1, b1, 2);

  Fp2_add_no_red (t0, a->a1, a->a2);
  Fp2_mul_no_red_o2 (T3, t0, b1);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T1->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T1->A1);
  Fpl_sub_o2_ct (T2->A0, T3->A0, T3->A1, Fp_config.mem);
  Fpl_add_o2_ct (T2->A1, T3->A0, T3->A1, Fp_config.mem);  
  Fp2l_add_o2 (c->A0, T2, T0);

  Fp2_add_no_red (t0, a->a0, a->a1);
  Fp2_add_no_red (t1, b0, b1);
  Fp2_mul_no_red_o2 (T3, t0, t1);
  Fp2l_add_no_red (T2, T0, T1);
  Fpl_sub_o2_ct (c->A1->A0, T3->A0, T2->A0, Fp_config.mem);
  Fpl_sub_no_red (c->A1->A1, T3->A1, T2->A1);

  Fp2_add_no_red (t0, a->a0, a->a2);
  Fp2_mul_no_red_o2 (T3, t0, b0);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T0->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T0->A1);
  Fpl_add_o2_ct (c->A2->A0, T3->A0, T1->A0, Fp_config.mem); 
  Fpl_add_no_red (c->A2->A1, T3->A1, T1->A1);
}

/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * The first coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b1*v + b2*v^2 = v*(b1 + b2*v).
 */
void v3minusxi_Fp6_mul_sparse12_no_red (Fp6l_t c, Fp6_t a, Fp2_t b1, Fp2_t b2, void *mem) {
  Fp2_t t0, t1;
  Fp2l_t T0, T1, T2, T3;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  
  T0->A0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T0->A1->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T1->A0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T1->A1->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T2->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  
  Fp2_mul_no_red_o1 (T1, a->a1, b1, 2);
  Fp2_mul_no_red_o1 (T2, a->a2, b2, 2);

  Fp2_add_no_red (t0, a->a1, a->a2);
  Fp2_add_no_red (t1, b1, b2);
  Fp2_mul_no_red_o2 (T3, t0, t1);
  Fp2l_add_no_red (T0, T1, T2);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T0->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T0->A1);
  Fpl_sub_o2_ct (c->A0->A0, T3->A0, T3->A1, Fp_config.mem);
  Fpl_add_o2_ct (c->A0->A1, T3->A0, T3->A1, Fp_config.mem);

  Fp2_add_no_red (t0, a->a0, a->a1);
  Fp2_mul_no_red_o2 (T3, t0, b1);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T1->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T1->A1);
  Fpl_sub_o1_ct (T0->A0, T2->A0, T2->A1, 1, Fp_config.mem);
  Fpl_add_no_red (T0->A1, T2->A0, T2->A1);
  Fp2l_add_o2 (c->A1, T3, T0);

  Fp2_add_no_red (t0, a->a0, a->a2);
  Fp2_mul_no_red_o2 (T3, t0, b2);
  Fpl_sub_o2_ct (T3->A0, T3->A0, T2->A0, Fp_config.mem);
  Fpl_sub_no_red (T3->A1, T3->A1, T2->A1);
  Fpl_add_o2_ct (c->A2->A0, T3->A0, T1->A0, Fp_config.mem); 
  Fpl_add_no_red (c->A2->A1, T3->A1, T1->A1);
}

int Fp6_initialize_config (void) {
  uint64_t *vpp;

  Fp6_config.mem = NULL;
  Fp6_config.mem = (void *) malloc (24 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp6_config.mem == NULL) return ERR_OUT_OF_MEMORY;
  
  vpp = (uint64_t *) Fp6_config.mem;
  if (Fp_init (Fp6_config.vppow) < 0) printf ("Fp6_config memory error.\n");
  if (Fp_init (Fp6_config.vppow2) < 0) printf ("Fp6_config memory error.\n");
  if (Fp_init (Fp6_config.vp2pow) < 0) printf ("Fp6_config memory error.\n");
  if (Fp_init (Fp6_config.vp2pow2) < 0) printf ("Fp6_config memory error.\n");
  if (Fp_init (Fp6_config.vppowinv) < 0) printf ("Fp6_config memory error.\n");
  
  if (PAIR_CURVE == BN12) {
    vpp[3] = (uint64_t) 0x2523648240000001;
    vpp[2] = (uint64_t) 0x7080EB4000000006;
    vpp[1] = (uint64_t) 0x181800000000000C;
    vpp[0] = (uint64_t) 0xD98000000000000B;
    Fp_set (Fp6_config.vppow, vpp);
   
    vpp[3] = (uint64_t) 0x0000000000000000;
    vpp[2] = (uint64_t) 0x49B3624000000002;
    vpp[1] = (uint64_t) 0x4909000000000006;
    vpp[0] = (uint64_t) 0xCD80000000000007;
    Fp_set (Fp6_config.vp2pow, vpp);
    
    Fp_mul (Fp6_config.vppow2, Fp6_config.vppow, Fp6_config.vppow);
    Fp_neg (Fp6_config.vppow2, Fp6_config.vppow2);
    Fp_mul (Fp6_config.vp2pow2, Fp6_config.vp2pow, Fp6_config.vp2pow);
    Fp_modinv (Fp6_config.vppowinv, Fp6_config.vppow);
    Fp_neg (Fp6_config.vppowinv, Fp6_config.vppowinv);

    Fp6_config.Fp6_mul = v3minusxi_Fp6_mul;
    Fp6_config.Fp6_squ = v3minusxi_Fp6_squ;
    Fp6_config.Fp6_inv = v3minusxi_Fp6_inv;
    Fp6_config.Fp6_mulv = v3minusxi_Fp6_mulv;
    Fp6_config.Fp6_ppow = p1mod3_v3minusxi_Fp6_ppow;
    Fp6_config.Fp6_p2pow = p1mod3_v3minusxi_Fp6_p2pow;
    Fp6_config.Fp6_p3pow = p1mod3_v3minusxi_Fp6_p3pow;
    Fp6_config.Fp6_mul_no_red = v3minusxi_Fp6_mul_no_red;
    Fp6_config.Fp6_mul_sparse01_no_red = v3minusxi_Fp6_mul_sparse01_no_red;
    Fp6_config.Fp6_mul_sparse12_no_red = v3minusxi_Fp6_mul_sparse12_no_red;
  } 
  else if (PAIR_CURVE == BN12CP) {
    vpp[3] = (uint64_t) 0x0000000000000000;
    vpp[2] = (uint64_t) 0x480000798001F455;
    vpp[1] = (uint64_t) 0xC1F2CDFBE8A64BAC;
    vpp[0] = (uint64_t) 0xBAF243C3000D7FF8;
    Fp_set (Fp6_config.vppow, vpp);
   
    vpp[3] = (uint64_t) 0x2400005100016456;
    vpp[2] = (uint64_t) 0x09FF9DB1F4ECF9A1;
    vpp[1] = (uint64_t) 0xB461E994A5E5451E;
    vpp[0] = (uint64_t) 0x17903165FFCB801A;
    Fp_set (Fp6_config.vp2pow, vpp);
    
    Fp_mul (Fp6_config.vppow2, Fp6_config.vppow, Fp6_config.vppow);
    Fp_neg (Fp6_config.vppow2, Fp6_config.vppow2);
    Fp_mul (Fp6_config.vp2pow2, Fp6_config.vp2pow, Fp6_config.vp2pow);
    Fp_modinv (Fp6_config.vppowinv, Fp6_config.vppow);
    Fp_neg (Fp6_config.vppowinv, Fp6_config.vppowinv);

    Fp6_config.Fp6_mul = v3minusxi_Fp6_mul;
    Fp6_config.Fp6_squ = v3minusxi_Fp6_squ;
    Fp6_config.Fp6_inv = v3minusxi_Fp6_inv;
    Fp6_config.Fp6_mulv = v3minusxi_Fp6_mulv;
    Fp6_config.Fp6_ppow = p1mod3_v3minusxi_Fp6_ppow;
    Fp6_config.Fp6_p2pow = p1mod3_v3minusxi_Fp6_p2pow;
    Fp6_config.Fp6_p3pow = p1mod3_v3minusxi_Fp6_p3pow;
    Fp6_config.Fp6_mul_no_red = v3minusxi_Fp6_mul_no_red;
    Fp6_config.Fp6_mul_sparse01_no_red = v3minusxi_Fp6_mul_sparse01_no_red;
    Fp6_config.Fp6_mul_sparse12_no_red = v3minusxi_Fp6_mul_sparse12_no_red;
  }  
  else if (PAIR_CURVE == BN12tiny) {
	  uint64_t vppow[1] = { 0x15C8EADEB9FF1F };
	  uint64_t vp2pow[1] = { 0xC59AB172BB };
	  Fp_set(Fp6_config.vppow, vppow);
	  Fp_set(Fp6_config.vp2pow, vp2pow);

	  Fp_mul(Fp6_config.vppow2, Fp6_config.vppow, Fp6_config.vppow);
	  Fp_neg(Fp6_config.vppow2, Fp6_config.vppow2);
	  Fp_mul(Fp6_config.vp2pow2, Fp6_config.vp2pow, Fp6_config.vp2pow);
	  Fp_modinv(Fp6_config.vppowinv, Fp6_config.vppow);
	  Fp_neg(Fp6_config.vppowinv, Fp6_config.vppowinv);

	  Fp6_config.Fp6_mul = v3minusxi_Fp6_mul;
	  Fp6_config.Fp6_squ = v3minusxi_Fp6_squ;
	  Fp6_config.Fp6_inv = v3minusxi_Fp6_inv;
	  Fp6_config.Fp6_mulv = v3minusxi_Fp6_mulv;
	  Fp6_config.Fp6_ppow = p1mod3_v3minusxi_Fp6_ppow;
	  Fp6_config.Fp6_p2pow = p1mod3_v3minusxi_Fp6_p2pow;
	  Fp6_config.Fp6_p3pow = p1mod3_v3minusxi_Fp6_p3pow;
	  Fp6_config.Fp6_mul_no_red = v3minusxi_Fp6_mul_no_red;
	  Fp6_config.Fp6_mul_sparse01_no_red = v3minusxi_Fp6_mul_sparse01_no_red;
	  Fp6_config.Fp6_mul_sparse12_no_red = v3minusxi_Fp6_mul_sparse12_no_red;
  }
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }
  
  return ERR_SUCCESS;
}

void Fp6_free_config (void) {
  free (Fp6_config.mem);
  Fp_free (Fp6_config.vppow);
  Fp_free (Fp6_config.vppow2);
  Fp_free (Fp6_config.vp2pow);
  Fp_free (Fp6_config.vp2pow2);
  Fp_free (Fp6_config.vppowinv);
}


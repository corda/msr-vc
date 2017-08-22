/* Fp12/Fp12overFp4.c
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

#include "Fp12overFp4.h"
#include "../pairing/pairing.h"

Fp12c_config_t Fp12c_config;

/* Initialize an element in Fp12 by initializing three Fp4 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp12c_init(Fp12c_t dst) {
	int err = 0;
	err = Fp4_init(dst->a0);
	if (err < 0) return err;
	err = Fp4_init(dst->a1);
	if (err < 0) return err;
	return Fp4_init(dst->a2);
}

/* Free the memory allocated in src. */
void Fp12c_free (Fp12c_t src) {
	Fp4_free(src->a0);
	Fp4_free(src->a1);
	Fp4_free(src->a2);
}

/* Get a random element in Fp12. */
void Fp12c_rand (Fp12c_t c) {
  Fp4_rand(c->a0);
  Fp4_rand(c->a1);
  Fp4_rand(c->a2);
}

/* Set c to 0. */
void Fp12c_set_zero (Fp12c_t c) {
  Fp4_set_zero (c->a0);
  Fp4_set_zero (c->a1);
  Fp4_set_zero (c->a2);
}

/* Set c to 1. */
void Fp12c_set_one (Fp12c_t c) {
  Fp4_set_one (c->a0);
  Fp4_set_zero (c->a1);
  Fp4_set_zero (c->a2);
}

/* Copy an Fp12 element a to an Fp12 element c. */
void Fp12c_copy(Fp12c_t c, Fp12c_t a) {
	Fp4_copy(c->a0, a->a0);
	Fp4_copy(c->a1, a->a1);
	Fp4_copy(c->a2, a->a2);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp12c_print (Fp12c_t a) {
  Fp4_print(a->a0);
  printf("\n");
  Fp4_print(a->a1);
  printf("\n");
  Fp4_print(a->a2);
  printf("\n");
}

/* Compare two Fp12 elements a and b for equality. */
int Fp12c_cmpeq(Fp12c_t a, Fp12c_t b) {
  return (Fp4_cmpeq(a->a0, b->a0) 
       && Fp4_cmpeq(a->a1, b->a1) && Fp4_cmpeq(a->a2, b->a2));
}

/* Negate an Fp12 element. */
void Fp12c_neg(Fp12c_t c, Fp12c_t a) {
	Fp4_neg(c->a0, a->a0);
	Fp4_neg(c->a1, a->a1);
	Fp4_neg(c->a2, a->a2);
}

/* Add two Fp12 elements coefficient wise. */
void Fp12c_add (Fp12c_t c, Fp12c_t a, Fp12c_t b) {
	Fp4_add(c->a0, a->a0, b->a0);
	Fp4_add(c->a1, a->a1, b->a1);
	Fp4_add(c->a2, a->a2, b->a2);
}

/* Subtract an Fp12 element b from an Fp12 element a coefficient wise. */
void Fp12c_sub (Fp12c_t c, Fp12c_t a, Fp12c_t b) {
	Fp4_sub(c->a0, a->a0, b->a0);
	Fp4_sub(c->a1, a->a1, b->a1);
	Fp4_sub(c->a2, a->a2, b->a2);
}

/* Multiply an Fp12 element by an Fp4 element. */
void Fp12c_mulFp4(Fp12c_t c, Fp12c_t a, Fp4_t b) {
	Fp4_mul(c->a0, a->a0, b);
	Fp4_mul(c->a1, a->a1, b);
	Fp4_mul(c->a2, a->a2, b);
}

/* Multiply an Fp12 element by an Fp2 element. */
void Fp12c_mulFp2(Fp12c_t c, Fp12c_t a, Fp2_t b) {
	Fp4_mulFp2(c->a0, a->a0, b);
	Fp4_mulFp2(c->a1, a->a1, b);
	Fp4_mulFp2(c->a2, a->a2, b);
}

/* Multiply an Fp12 element by an Fp element. */
void Fp12c_mulFp(Fp12c_t c, Fp12c_t a, Fp_t b) {
	Fp4_mulFp(c->a0, a->a0, b);
	Fp4_mulFp(c->a1, a->a1, b);
	Fp4_mulFp(c->a2, a->a2, b);
}

/* Multiply two Fp12 elements. */
void Fp12c_mul (Fp12c_t c, Fp12c_t a, Fp12c_t b) {
  Fp12c_config.Fp12c_mul(c, a, b, Fp12c_config.mem);
}

/* Square an Fp12 element. */
void Fp12c_squ (Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_squ(c, a, Fp12c_config.mem);
}

/* Invert an Fp12 element. */
void Fp12c_inv (Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_inv(c, a, Fp12c_config.mem);
}

/* Multiply an Fp12 element by the special element v. */
void Fp12c_mulv(Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_mulv(c, a, Fp12c_config.mem);
}

/* Multiply an Fp12 element by the special element i. */
void Fp12c_muli(Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_muli(c, a);
}

/* Multiply an Fp12 element by the special element s. */
void Fp12c_muls(Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_muls(c, a);
}

/* Compute the p-power Frobenius of an Fp12 element. */
void Fp12c_ppow(Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_ppow(c, a);
}

/* Compute the p^2-power Frobenius of an Fp12 element. */
void Fp12c_p2pow(Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_p2pow(c, a);
}

/* Compute the p^4-power Frobenius of an Fp12 element. */
void Fp12c_p4pow(Fp12c_t c, Fp12c_t a) {
  Fp12c_config.Fp12c_p4pow(c, a);
}

/* Multiply an Fp12c element and a sparse Fp12c element without reduction.
   This function assumes that Fp12 is constructed via v^3 - s.
 */
void Fp12c_mul_sparse01 (Fp12c_t c, Fp12c_t a, Fp4_t b0, Fp4_t b1){
  Fp12c_config.Fp12c_mul_sparse01 (c, a, b0, b1, Fp12c_config.mem);
}
void Fp12c_mul_sparse12 (Fp12c_t c, Fp12c_t a, Fp4_t b1, Fp4_t b2){
  Fp12c_config.Fp12c_mul_sparse12 (c, a, b1, b2, Fp12c_config.mem);
}

/* Multiply two Fp12 elements. This uses that Fp12 is constructed over Fp2 via v^3 - s. 
 * This is Algorithm 13 in Beuchat et al., Pairing 2010.
 */
void v3minuss_Fp12c_mul(Fp12c_t c, Fp12c_t a, Fp12c_t b, void *mem) {
	Fp4_t t0, t1, t2, t3, t4, t5;
  
  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a0->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t1->a1->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t2->a0->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t2->a0->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t2->a1->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t2->a1->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  t3->a0->a0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  t3->a0->a1->limbs = (uint64_t *) mem + 13*Fp_config.m->n;
  t3->a1->a0->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  t3->a1->a1->limbs = (uint64_t *) mem + 15*Fp_config.m->n;
  t4->a0->a0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  t4->a0->a1->limbs = (uint64_t *) mem + 17*Fp_config.m->n;
  t4->a1->a0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  t4->a1->a1->limbs = (uint64_t *) mem + 19*Fp_config.m->n;
  t5->a0->a0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  t5->a0->a1->limbs = (uint64_t *) mem + 21*Fp_config.m->n;
  t5->a1->a0->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  t5->a1->a1->limbs = (uint64_t *) mem + 23*Fp_config.m->n;

	Fp4_mul(t0, a->a0, b->a0);    // t0 = a0*b0
	Fp4_mul(t1, a->a1, b->a1);    // t1 = a1*b1
	Fp4_mul(t2, a->a2, b->a2);    // t2 = a2*b2

	Fp4_add(t3, a->a1, a->a2);    // t3 = a1 + a2
	Fp4_add(t4, b->a1, b->a2);    // t4 = b1 + b2
	Fp4_mul(t3, t3, t4);          // t3 = (a1 + a2)*(b1 + b2) = a1*b1 + a2*b1 + a1*b2 + a2*b2
	Fp4_sub(t3, t3, t1);          // t3 = a2*b1 + a1*b2 + a2*b2
	Fp4_sub(t3, t3, t2);          // t3 = a2*b1 + a1*b2 
  Fp4_muls(t3, t3);					  // t3 = (a2*b1 + a1*b2)*s	
                                // t3 now almost has the result for c0, t4 is free.		

	Fp4_add(t4, a->a0, a->a2);    // t4 = a0 + a2
	Fp4_add(t5, b->a0, b->a2);    // t5 = b0 + b2
	Fp4_mul(t4, t4, t5);          // t4 = (a0 + a2)*(b0 + b2) = a0*b0 + a2*b0 + a0*b2 + a2*b2
	Fp4_sub(t4, t4, t0);          // t4 = a2*b0 + a0*b2 + a2*b2
	Fp4_sub(t4, t4, t2);          // t4 = a2*b0 + a0*b2
	Fp4_add(c->a2, t4, t1);				// c2 = a2*b0 + a0*b2 + a1*b1 	
                                // c2 is done, t4 and t5 are free.

	Fp4_add(t4, a->a0, a->a1);    // t4 = a0 + a1
	Fp4_add(t5, b->a0, b->a1);    // t4 = b0 + b1
	Fp4_mul(t4, t4, t5);          // t4 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1
	Fp4_sub(t4, t4, t0);          // t4 = a1*b0 + a0*b1 + a1*b1
	Fp4_sub(t4, t4, t1);          // t4 = a1*b0 + a0*b1
	Fp4_muls(t2, t2);            // t2 = a2*b2 * s
	Fp4_add(c->a1, t4, t2);				// c1 = a1*b0 + a0*b1 + a2*b2*s	
                                // c1 is done.

	Fp4_add(c->a0, t3, t0);					// c0 is done.
}

/* Square an Fp12 element. This uses that Fp12 is constructed over Fp4 via v^3 - s. 
 * This is Algorithm 16 in Beuchat et al., Pairing 2010.
 */
void v3minuss_Fp12c_squ(Fp12c_t c, Fp12c_t a, void *mem) {
	Fp4_t t0, t1, t2, t3, t4;
  
  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a0->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t1->a1->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t2->a0->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t2->a0->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t2->a1->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t2->a1->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  t3->a0->a0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  t3->a0->a1->limbs = (uint64_t *) mem + 13*Fp_config.m->n;
  t3->a1->a0->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  t3->a1->a1->limbs = (uint64_t *) mem + 15*Fp_config.m->n;
  t4->a0->a0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  t4->a0->a1->limbs = (uint64_t *) mem + 17*Fp_config.m->n;
  t4->a1->a0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  t4->a1->a1->limbs = (uint64_t *) mem + 19*Fp_config.m->n;
	
	Fp4_mul(t0, a->a0, a->a1);    // t0 = a0*a1  
	Fp4_add(t0, t0, t0);          // t0 = 2*a0*a1
	Fp4_squ(t1, a->a2);           // t1 = a2^2
	Fp4_squ(t2, a->a0);           // t2 = a0^2
	Fp4_sub(t3, a->a0, a->a1);    // t3 = a0 - a1
	Fp4_add(t3, t3, a->a2);       // t3 = a0 - a1 + a2
	Fp4_mul(t4, a->a1, a->a2);    // t4 = a1*a2
	Fp4_add(t4, t4, t4);          // t4 = 2*a1*a2

	Fp4_sub(c->a2, t0, t1);       // c2 = 2*a0*a1 - a2^2
	Fp4_muls(t1, t1);             // t1 = s*a2^2
	Fp4_add(c->a1, t0, t1);       // c1 = 2*a0*a1 + s*a2^2

	Fp4_squ(t0, t3);              // t0 = (a0 - a1 + a2)^2
	Fp4_add(c->a2, c->a2, t0);    // c2 = a0^2 + a1^2 - 2*a1*a2 + 2*a0*a2
	Fp4_add(c->a2, c->a2, t4);    // c2 = a0^2 + a1^2 + 2*a0*a2
	Fp4_sub(c->a2, c->a2, t2);    // c2 = a1^2 + 2*a0*a2

	Fp4_muls(t4, t4);             // t4 = 2*a1*a2*s
	Fp4_add(c->a0, t4, t2);       // c0 = a0^2 + 2*a1*a2*s
}

/* Multiplication by the special element v in Fp12, with v^3 = s in Fp4.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = a2*s + a0*v + a1*v^2.
 */
void v3minuss_Fp12c_mulv(Fp12c_t c, Fp12c_t a, void *mem) {
  Fp4_t t0;
  
  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	
  Fp4_muls(t0, a->a2);
  Fp4_copy(c->a2, a->a1);
	Fp4_copy(c->a1, a->a0);
  Fp4_copy(c->a0, t0);
}

/* Multiply an Fp12c element by the special element i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 */
void p3mod4_Fp12c_muli(Fp12c_t c, Fp12c_t a) {
  Fp4_muli (c->a0, a->a0);
  Fp4_muli (c->a1, a->a1);
  Fp4_muli (c->a2, a->a2);
}

/* Multiply an Fp12c element by the special element s. This function assumes that Fp4 is generated
 * over Fp2 by s^2 - i, and that the field extension is as described above.
 */
void s2minusxi_Fp12c_muls(Fp12c_t c, Fp12c_t a) {
  Fp4_muls (c->a0, a->a0);
  Fp4_muls (c->a1, a->a1);
  Fp4_muls (c->a2, a->a2);
}


/* Compute the p-power Frobenius of an Fp12 element. This uses a constant VPPOW, 
 * such that v^p = VPPOW*v, i.e. VPPOW = v^(p-1). In the case, where Fp12 is an 
 * extension of Fp2 via v^3 - s and p = 1 (mod 3), we have VPPOW = s^(p-1)/3 in Fp4. 
 * VPPOW2 = VPPOW^2.
 */
void p1mod3_v3minuss_Fp12c_ppow(Fp12c_t c, Fp12c_t a) {
  Fp4_ppow(c->a0, a->a0);
	Fp4_ppow(c->a1, a->a1);
	Fp4_ppow(c->a2, a->a2);
	Fp4_mulFp2(c->a1, c->a1, Fp12c_config.vppow);
	Fp4_mulFp2(c->a2, c->a2, Fp12c_config.vppow2);
}

/* Compute the p^4-power Frobenius of an Fp12 element. As above, we use a constant
 * VP2POW such that v^(p^4) = VP4POW*v, i.e. VP4POW = v^(p^2-1).  In the case, where Fp12 is an 
 * extension of Fp4 via v^3 - s and p = 1 (mod 3), we have VP4POW = s^(p^4-1)/3 is a third root
 * of unity in Fp. VP4POW2 = VP4POW^2.  
 */
void p1mod3_v3minuss_Fp12c_p4pow(Fp12c_t c, Fp12c_t a) {
	Fp4_copy(c->a0, a->a0);
	Fp4_mulFp(c->a1, a->a1, Fp12c_config.vp4pow);
	Fp4_mulFp(c->a2, a->a2, Fp12c_config.vp4pow2);
}

/* Compute the p^2-power Frobenius of an Fp12 element. As above, we use a constant
 * VP2POW such that v^(p^2) = VP2POW*v, i.e. VP2POW = v^(p^2-1).  In the case, where Fp12 is an 
 * extension of Fp4 via v^3 - s and p = 1 (mod 3), we have VP2POW = s^(p^2-1)/3 in Fp4. 
 * VP2POW2 = VP2POW^2.
 */
void p1mod3_v3minuss_Fp12c_p2pow(Fp12c_t c, Fp12c_t a) {
	Fp4_p2pow(c->a0, a->a0);
	Fp4_p2pow(c->a1, a->a1);
	Fp4_p2pow(c->a2, a->a2);
	Fp4_mulFp(c->a1, c->a1, Fp12c_config.vp2pow);
	Fp4_mulFp(c->a2, c->a2, Fp12c_config.vp2pow2);
}

/* Inversion of an Fp12 element. This uses that Fp12 is constructed over Fp4 via v^3 - s. 
 * This is Algorithm 17 in Beuchat et al., Pairing 2010.
 */
void v3minuss_Fp12c_inv(Fp12c_t c, Fp12c_t a, void *mem) {
	Fp4_t t0, t1, t2, t3, t4, t5;
  
  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a0->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t1->a1->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t2->a0->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t2->a0->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t2->a1->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t2->a1->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  t3->a0->a0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  t3->a0->a1->limbs = (uint64_t *) mem + 13*Fp_config.m->n;
  t3->a1->a0->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  t3->a1->a1->limbs = (uint64_t *) mem + 15*Fp_config.m->n;
  t4->a0->a0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  t4->a0->a1->limbs = (uint64_t *) mem + 17*Fp_config.m->n;
  t4->a1->a0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  t4->a1->a1->limbs = (uint64_t *) mem + 19*Fp_config.m->n;
  t5->a0->a0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  t5->a0->a1->limbs = (uint64_t *) mem + 21*Fp_config.m->n;
  t5->a1->a0->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  t5->a1->a1->limbs = (uint64_t *) mem + 23*Fp_config.m->n;

	Fp4_squ(t0, a->a0);         // t0 = a0^2
	Fp4_squ(t1, a->a1);         // t1 = a1^2
	Fp4_squ(t2, a->a2);         // t2 = a2^2

	Fp4_mul(t3, a->a0, a->a1);  // t3 = a0*a1
	Fp4_mul(t4, a->a0, a->a2);  // t4 = a0*a2
	Fp4_mul(t5, a->a1, a->a2);  // t5 = a1*a2

	Fp4_muls(t2, t2);           // t2 = s*a2^2
	Fp4_muls(t5, t5);           // t5 = s*a1*a2

	Fp4_sub(t0, t0, t5);        // t0 = a0^2 - s*a1*a2
	Fp4_sub(t1, t1, t4);        // t1 = a1^2 - a0*a2
	Fp4_sub(t2, t2, t3);        // t2 = s*a2^2 - a0*a1

	Fp4_mul(t3, a->a0, t0);     // t3 = a0^3 - s*a0*a1*a2
	Fp4_mul(t4, a->a2, t2);     // t4 = s*a2^3 - a0*a1*a2
	Fp4_muls(t4, t4);           // t4 = s^2*a2^3 - s*a0*a1*a2
	Fp4_add(t3, t3, t4);        // t3 = a0^3 + s^2*a2^3 - 2*s*a0*a1*a2 
	Fp4_mul(t4, a->a1, t1);     // t4 = a1^3 - a0*a1*a2
	Fp4_muls(t4, t4);           // t4 = s*a1^3 - s*a0*a1*a2
	Fp4_add(t3, t3, t4);        // t3 = a0^3 + s*a1^3 + s^2*a2^3 - 3*s*a0*a1*a2

	Fp4_inv(t3, t3);            // t3 = 1/N(a)

	Fp4_mul(c->a0, t0, t3);
	Fp4_mul(c->a1, t2, t3);
	Fp4_mul(c->a2, t1, t3);
}


/* Multiply two Fp12c.
 * This function assumes that Fp12c is constructed via v^3 - s.
 * The first coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b1*v + b2*v^2 = v*(b1 + b2*v).
 */
void v3minuss_Fp12c_mul_sparse12 (Fp12c_t c, Fp12c_t a, Fp4_t b1, Fp4_t b2, void *mem) {
  Fp4_t t0, t1, t2, t3, t4, t5;
  
  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a0->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t1->a1->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t2->a0->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t2->a0->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t2->a1->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t2->a1->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  t3->a0->a0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  t3->a0->a1->limbs = (uint64_t *) mem + 13*Fp_config.m->n;
  t3->a1->a0->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  t3->a1->a1->limbs = (uint64_t *) mem + 15*Fp_config.m->n;
  t4->a0->a0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  t4->a0->a1->limbs = (uint64_t *) mem + 17*Fp_config.m->n;
  t4->a1->a0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  t4->a1->a1->limbs = (uint64_t *) mem + 19*Fp_config.m->n;
  t5->a0->a0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  t5->a0->a1->limbs = (uint64_t *) mem + 21*Fp_config.m->n;
  t5->a1->a0->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  t5->a1->a1->limbs = (uint64_t *) mem + 23*Fp_config.m->n;
  
  Fp4_mul (t3, a->a1, b1);
  Fp4_mul (t4, a->a2, b2);

  Fp4_add (t0, a->a1, a->a2);
  Fp4_add (t1, b1, b2);
  Fp4_mul (t5, t0, t1);
  Fp4_add (t2, t3, t4);
  Fp2_sub (t5->a0, t5->a0, t2->a0);
  Fp2_sub (t5->a1, t5->a1, t2->a1);
  Fp2_sub (c->a0->a0, t5->a0, t5->a1);
  Fp2_add (c->a0->a1, t5->a0, t5->a1);

  Fp4_add (t0, a->a0, a->a1);
  Fp4_mul (t5, t0, b1);
  Fp2_sub (t5->a0, t5->a0, t3->a0);
  Fp2_sub (t5->a1, t5->a1, t3->a1);
  Fp2_sub (t2->a0, t4->a0, t4->a1);
  Fp2_add (t2->a1, t4->a0, t4->a1);
  Fp4_add (c->a1, t5, t2);

  Fp4_add (t0, a->a0, a->a2);
  Fp4_mul (t5, t0, b2);
  Fp2_sub (t5->a0, t5->a0, t4->a0);
  Fp2_sub (t5->a1, t5->a1, t4->a1);
  Fp2_add (c->a2->a0, t5->a0, t3->a0); 
  Fp2_add (c->a2->a1, t5->a1, t3->a1);
}

/* Multiply two Fp12c.
 * This function assumes that Fp12c is constructed via v^3 - s.
 * The first coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b0 + b1*v, i.e. b2 = 0.
 */
void v3minuss_Fp12c_mul_sparse01 (Fp12c_t c, Fp12c_t a, Fp4_t b0, Fp4_t b1, void *mem) {
  Fp4_t t0, t1, t2, t3, t4, t5;
  
  t0->a0->a0->limbs = (uint64_t *) mem;
  t0->a0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t0->a1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t0->a1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t1->a0->a0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t1->a0->a1->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t1->a1->a0->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  t1->a1->a1->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
  t2->a0->a0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  t2->a0->a1->limbs = (uint64_t *) mem + 9*Fp_config.m->n;
  t2->a1->a0->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  t2->a1->a1->limbs = (uint64_t *) mem + 11*Fp_config.m->n;
  t3->a0->a0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  t3->a0->a1->limbs = (uint64_t *) mem + 13*Fp_config.m->n;
  t3->a1->a0->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  t3->a1->a1->limbs = (uint64_t *) mem + 15*Fp_config.m->n;
  t4->a0->a0->limbs = (uint64_t *) mem + 16*Fp_config.m->n;
  t4->a0->a1->limbs = (uint64_t *) mem + 17*Fp_config.m->n;
  t4->a1->a0->limbs = (uint64_t *) mem + 18*Fp_config.m->n;
  t4->a1->a1->limbs = (uint64_t *) mem + 19*Fp_config.m->n;
  t5->a0->a0->limbs = (uint64_t *) mem + 20*Fp_config.m->n;
  t5->a0->a1->limbs = (uint64_t *) mem + 21*Fp_config.m->n;
  t5->a1->a0->limbs = (uint64_t *) mem + 22*Fp_config.m->n;
  t5->a1->a1->limbs = (uint64_t *) mem + 23*Fp_config.m->n;
  
  Fp4_mul (t3, a->a0, b0);            // t3 = a0*b0
  Fp4_mul (t4, a->a1, b1);            // t4 = a1*b1

  Fp4_add (t0, a->a1, a->a2);         // t0 = a1 + a2
  Fp4_mul (t5, t0, b1);               // t5 = t0*b1 = a1*b1 + a2*b1
  Fp4_sub (t5, t5, t4);               // t5 = a2*b1
  Fp4_muls (t5, t5);                  // t5 = a2*b1*s
  Fp4_add (t5, t5, t3);               // t5 = a2*b1*s + a0*b0

  Fp4_add (t0, a->a0, a->a1);         // t0 = a0 + a1
  Fp4_add (t1, b0, b1);               // t1 = b0 + b1
  Fp4_mul (c->a1, t0, t1);            // c1 = (a0+a1)(b0+b1)
  Fp4_add (t2, t3, t4);               // t2 = a0*b0 + a1*b1
  Fp4_sub (c->a1, c->a1, t2);         // c1 = a0*b1 + a1*b0

  Fp4_add (t0, a->a0, a->a2);         // t0 = a0 + a2
  Fp4_mul (c->a2, t0, b0);            // c2 = a0*b0 + a2*b0
  Fp4_sub (c->a2, c->a2, t3);         // c2 = a2*b0
  Fp4_add (c->a2, c->a2, t4);         // c2 = a2*b0 + a1*b1
  
  Fp4_copy (c->a0, t5);               // c0 = a2*b1*s + a0*b0
}

int Fp12c_initialize_config (void) {
  uint64_t *vpp;

  Fp12c_config.mem = NULL;
  Fp12c_config.mem = (void *) malloc (72 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp12c_config.mem == NULL) return ERR_OUT_OF_MEMORY;
  
  vpp = (uint64_t *) Fp12c_config.mem;
  if (Fp2_init (Fp12c_config.vppow) < 0) printf ("Fp12c_config memory error.\n");
  if (Fp2_init (Fp12c_config.vppow2) < 0) printf ("Fp12c_config memory error.\n");
  if (Fp_init (Fp12c_config.vp2pow) < 0) printf ("Fp12c_config memory error.\n");
  if (Fp_init (Fp12c_config.vp2pow2) < 0) printf ("Fp12c_config memory error.\n");
  if (Fp_init (Fp12c_config.vp4pow) < 0) printf ("Fp12c_config memory error.\n");
  if (Fp_init (Fp12c_config.vp4pow2) < 0) printf ("Fp12c_config memory error.\n");
  
  
  
  return ERR_SUCCESS;
}

void Fp12c_free_config (void) {
  free (Fp12c_config.mem);
  Fp2_free (Fp12c_config.vppow);
  Fp2_free (Fp12c_config.vppow2);
  Fp_free (Fp12c_config.vp2pow);
  Fp_free (Fp12c_config.vp2pow2);
  Fp_free (Fp12c_config.vp4pow);
  Fp_free (Fp12c_config.vp4pow2);
}


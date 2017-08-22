/* Fp2/Fp2.c
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

#include "Fp2.h"

Fp2_config_t Fp2_config;

/************** Regular Fp2 functions ***************/

/* Initialize an element in Fp2 by initializing two Fp elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp2_init (Fp2_t dst) {
	int err = 0;
	if (err = Fp_init(dst->a0) < 0) return err;
	return Fp_init(dst->a1);
}

/* Free the memory allocated in src. */
void Fp2_free (Fp2_t src) {
	Fp_free(src->a0);
	Fp_free(src->a1);
}

/* Get a random element in Fp2. */
void Fp2_rand (Fp2_t c) {
  Fp_rand(c->a0);
  Fp_rand(c->a1);
}

/* Set c to the value of a0 + a1*i. */
void Fp2_set (Fp2_t c, uint64_t *a0, uint64_t *a1) {
  Fp_set(c->a0, a0);
  Fp_set(c->a1, a1);
}

/* Set c to the value x. */
void Fp2_set_ui (Fp2_t c, uint64_t x) {
  Fp_set_ui (c->a0, x);
  Fp_set_ui (c->a1, 0);
}

/* Set c to 0. */
void Fp2_set_zero (Fp2_t c) {
  Fp_set_ui (c->a0, 0);
  Fp_set_ui (c->a1, 0);
}

/* Set c to 1. */
void Fp2_set_one (Fp2_t c) {
  Fp_set_ui (c->a0, 1);
  Fp_set_ui (c->a1, 0);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp2_print (Fp2_t a) {
  Fp_print(a->a0);
  printf("  ");
  Fp_print(a->a1);
  printf("\n");
}

/* Copy an Fp2 element a to an Fp2 element c. */
void Fp2_copy(Fp2_t c, Fp2_t a) {
	Fp_copy(c->a0, a->a0);
	Fp_copy(c->a1, a->a1);
}

/* Compare two Fp2 elements a and b for equality. */
int Fp2_cmpeq(Fp2_t a, Fp2_t b) {
  return (Fp_cmp(a->a0, b->a0) == 0 && Fp_cmp(a->a1, b->a1) == 0);
}

/* Negate an Fp2 element. */
void Fp2_neg(Fp2_t c, Fp2_t a) {
	Fp_neg(c->a0, a->a0);
	Fp_neg(c->a1, a->a1);
}

/* Add two Fp2 elements coefficient wise. */
void Fp2_add (Fp2_t c, Fp2_t a, Fp2_t b) {
	Fp_add(c->a0, a->a0, b->a0);
	Fp_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp2 element b from an Fp2 element a coefficient wise. */
void Fp2_sub (Fp2_t c, Fp2_t a, Fp2_t b) {
	Fp_sub(c->a0, a->a0, b->a0);
	Fp_sub(c->a1, a->a1, b->a1);
}

/* Divide an Fp2 element a by 2 coefficient wise. */
void Fp2_div2 (Fp2_t c, Fp2_t a) {
	Fp_div2(c->a0, a->a0);
	Fp_div2(c->a1, a->a1);
}

/* Multiply an Fp2 element by 3. */
void Fp2_mul3 (Fp2_t c, Fp2_t a) {
  Fp_mul3 (c->a0, a->a0);
  Fp_mul3 (c->a1, a->a1);
}

/* Multiply an Fp2 element by an Fp element. */
void Fp2_mulFp(Fp2_t c, Fp2_t a, Fp_t b) {
	Fp_mul(c->a0, a->a0, b);
	Fp_mul(c->a1, a->a1, b);
}

/* Multiply two Fp2 elements. */
void Fp2_mul(Fp2_t c, Fp2_t a, Fp2_t b) {
  //Fp2_config.Fp2_mul(c, a, b, Fp2_config.mem);
  Fp2_mul_lazy (c, a, b);
}

/* Square an Fp2 element. */
void Fp2_squ(Fp2_t c, Fp2_t a){
  //Fp2_config.Fp2_squ(c, a, Fp2_config.mem);
  Fp2_squ_lazy(c, a);
}

/* Multiply an Fp2 element by the special element xi.*/
void Fp2_mulxi(Fp2_t c, Fp2_t a) {
  Fp2_config.Fp2_mulxi(c, a, Fp2_config.mem);
}

/* Multiply an Fp2 element by the special element i.*/
void Fp2_muli(Fp2_t c, Fp2_t a) {
  Fp2_config.Fp2_muli(c, a, Fp2_config.mem);
}

/* Multiply an Fp2 element by an element i*b, where b is in Fp. */
void Fp2_muliFp(Fp2_t c, Fp2_t a, Fp_t b) {
  Fp2_config.Fp2_muliFp(c, a, b, Fp2_config.mem);
}

/* Norm of an Fp2 element. */
void Fp2_norm(Fp_t c, Fp2_t a) {
  Fp2_config.Fp2_norm(c, a, Fp2_config.mem);
}

/* Invert an Fp2 element. */
void Fp2_inv(Fp2_t c, Fp2_t a) {
  Fp2_config.Fp2_inv(c, a, Fp2_config.mem);
}

/* Compute the p-power Frobenius of an Fp2 element. */
void Fp2_ppow(Fp2_t c, Fp2_t a) {
	 Fp2_config.Fp2_ppow(c, a);
}

void Fp2l_mulxi (Fp2l_t c, Fp2l_t a) {
  Fp2_config.Fp2l_mulxi(c, a, Fp2_config.mem);
}


/* Compute the p-power Frobenius of an Fp2 element. This function assumes that the field extension
 * is constructed via a binomial i^2 - alpha, alpha in Fp.
 * This is the same as (complex) conjugation, i.e. Fp_ppow(a0 + a1*i) = a0 - a1*i.
 * When we need this, it might be simpler to just negate the second coordinate in place.
 */
void binomial_Fp2_ppow(Fp2_t c, Fp2_t a) {
	Fp_copy(c->a0, a->a0);
	Fp_neg(c->a1, a->a1);
}

/* Multiply two Fp2 elements. This function assumes that p = 3 (mod4)
 * such that the quadratic extension is Fp2 = Fp(i) with i^2 = -1.
 * Then c = c0 + c1*i = a*b = (a0 + a1*i)*(b0 + b1*i) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)*i.
 */
void p3mod4_Fp2_mul(Fp2_t c, Fp2_t a, Fp2_t b, void *mem) {
	Fp_t t0, t1, t2;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
	
	Fp_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp_add(t1, b->a0, b->a1); // t1 = b0 + b1
	Fp_mul(t2, t0, t1);			  // t2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + b0*a1 + a1*b1

	Fp_mul(t0, a->a0, b->a0);	// t0 = a0*b0
	Fp_mul(t1, a->a1, b->a1);	// t1 = a1*b1

	Fp_sub(c->a0, t0, t1);		// c0 = a0*b0 - a1*b1
	Fp_add(t0, t0, t1);			// t0 = a0*b0 + a1*b1
	Fp_sub(c->a1, t2, t0);		// c1 = t2 - t0 = a0*b1 + a1*b0
}

/* Square an Fp2 element. This function assumes that p = 3 (mod4)
 * such that the quadratic extension is as described above and thus
 * c = c0 + c1*i = a^2 = (a0 + a1*i)^2 = (a0^2 - a1^2) + (2*a0*a1)*i.
 */
void p3mod4_Fp2_squ(Fp2_t c, Fp2_t a, void *mem) {
	Fp_t t0, t1, t2;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;

	Fp_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp_sub(t1, a->a0, a->a1);	// t1 = a0 - a1

	Fp_mul(t2, a->a0, a->a1);	// t2 = a0*a1

	Fp_mul(c->a0, t0, t1);		// c0 = (a0 - a1)*(a0 + a1)
	Fp_add(c->a1, t2, t2);		// c1 = 2*a0*a1
}

/* Multiply an Fp2 element by the special element xi = 1+i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 * We have c = c0 + c1*i = (a0 - a1) + (a0 + a1)*i. */
void p3mod4_Fp2_mulxi(Fp2_t c, Fp2_t a, void *mem) {
	Fp_t t0;
  t0->limbs = (uint64_t *) mem;

  Fp_sub(t0, a->a0, a->a1);
	Fp_add(c->a1, a->a0, a->a1);
  Fp_copy(c->a0, t0);
}

/* Multiply an Fp2 element by the special element i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 * We have c = c0 + c1*i = -a1 + a0*i. */
void p3mod4_Fp2_muli(Fp2_t c, Fp2_t a, void *mem) {
	Fp_t t0;
  t0->limbs = (uint64_t *) mem;
  
  Fp_neg(t0, a->a1);
	Fp_copy(c->a1, a->a0);
  Fp_copy(c->a0, t0);
}

/* Multiply an Fp2 element by an element i*b, where b is in Fp. */
void p3mod4_Fp2_muliFp(Fp2_t c, Fp2_t a, Fp_t b, void *mem) {
  Fp_t t0;
  t0->limbs = (uint64_t *) mem;

  Fp_mul(t0, a->a0, b);
  Fp_mul(c->a0, a->a1, b);
  Fp_neg(c->a0, c->a0);
  Fp_copy(c->a1, t0);
}

/* Compute the norm (in Fp) of an Fp2 element. This function assumes that p = 3 (mod 4)
 * such that the quadratic extension is as described above and the norm can be computed as
 * c = N(a) = N(a0 + a1*i) = (a0 + a1*i)*(a0 - a1*i) = a0^2 + a1^2.
 */
void p3mod4_Fp2_norm(Fp_t c, Fp2_t a, void *mem) {
	Fp_t t0; 
  t0->limbs = (uint64_t *) mem;
  
	Fp_mul(c, a->a0, a->a0);	// c = a0^2
	Fp_mul(t0, a->a1, a->a1);	// t0 = a1^2

	Fp_add(c, c, t0);			    // c = a0^2 + a1^2
}

/* Invert an Fp2 element. This function assumes that p = 3 (mod4)
 * such that the quadratic extension is as described above and thus
 * c = c0 + c1*i = a^-1 = (a0 + a1*i)^-1 = (a0 - a1*i)/(a0^2 + a1^2) = a0/N(a) - a1/N(a)*i.
  * mem should have room for 3 Fp elements.
 */
void p3mod4_Fp2_inv(Fp2_t c, Fp2_t a, void *mem) {
	Fp_t t0, t1;
  void *m;

  t0->limbs = (uint64_t *) mem;
  m = (void *) ((uint64_t *) mem + Fp_config.m->n);

	p3mod4_Fp2_norm (t0, a, m); // t0 = a0^2 + a1^2

  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
	Fp_modinv(t1, t0);			 // t1 = (a0^2 + a1^2)^-1

	Fp_mul(c->a0, a->a0, t1);	 // c0 = a0/(a0^2 + a1^2)
	Fp_mul(c->a1, a->a1, t1);	 // c1 = a1/(a0^2 + a1^2)
	Fp_neg(c->a1, c->a1);		 // c1 = -a1/(a0^2 + a1^2)
}

/************** Fp2 functions for lazy reduction ***************/
/* Functions assume that there is enough space such that no overflow occurs. */

/* Initialize a long Fp2 element. */
int Fp2l_init(Fp2l_t dst) {
	int err = 0;
	if (err = Fpl_init(dst->A0) < 0) return err;
	return Fpl_init(dst->A1);
}

/* Free the memory allocated in src. */
void Fp2l_free (Fp2l_t src) {
	Fpl_free(src->A0);
	Fpl_free(src->A1);
}

/* Copy an Fp2l element a to an Fp2l element c. */
void Fp2l_copy(Fp2l_t C, Fp2l_t A) {
	Fpl_copy(C->A0, A->A0);
	Fpl_copy(C->A1, A->A1);
}

/* Reduction function calling the Fp reduction for each coefficient. */
void Fp2l_red(Fp2_t c, Fp2l_t A) {
	Fpl_red(c->a0, A->A0);
	Fpl_red(c->a1, A->A1);
}

/* Add two Fp2 elements coefficient wise without reducing mod p for lazy reduction. */ 
void Fp2_add_no_red (Fp2_t c, Fp2_t a, Fp2_t b) {
	Fp_add_no_red(c->a0, a->a0, b->a0);
	Fp_add_no_red(c->a1, a->a1, b->a1);
}

/* Add two long Fp2 elements coefficient wise without reducing mod p for lazy reduction. */ 
void Fp2l_add_no_red (Fp2l_t C, Fp2l_t A, Fp2l_t B) {
	Fpl_add_no_red(C->A0, A->A0, B->A0);
	Fpl_add_no_red(C->A1, A->A1, B->A1);
}

/* Subtraction of Fp2l elements, option 1 in Aranha et al. */
void Fp2l_sub_o1 (Fp2l_t c, Fp2l_t a, Fp2l_t b, int h) {
  Fpl_sub_o1_ct(c->A0, a->A0, b->A0, h, Fp_config.mem);
  Fpl_sub_o1_ct(c->A1, a->A1, b->A1, h, Fp_config.mem);
}

/* Subtraction of Fp2l elements, constant-time version of option 2 in Aranha et al. */
void Fp2l_sub_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b) {
  Fpl_sub_o2_ct(c->A0, a->A0, b->A0, Fp_config.mem);
  Fpl_sub_o2_ct(c->A1, a->A1, b->A1, Fp_config.mem);
}

/* Addition of Fp2l elements, constant-time version of "option 2" addition in Aranha et al. */
void Fp2l_add_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b) {
  Fpl_add_o2_ct(c->A0, a->A0, b->A0, Fp_config.mem);
  Fpl_add_o2_ct(c->A1, a->A1, b->A1, Fp_config.mem);
}

/* Multiply an Fp2 element by an Fp element. */
void Fp2_mulFp_no_red (Fp2l_t c, Fp2_t a, Fp_t b) {
	Fp_mul_no_red (c->A0, a->a0, b);
	Fp_mul_no_red (c->A1, a->a1, b);
}

/* Multiply two Fp2 elements without reduction. */
void Fp2_mul_no_red_o1(Fp2l_t c, Fp2_t a, Fp2_t b, int h) {
  Fp2_config.Fp2_mul_no_red_o1(c, a, b, h, Fp2_config.mem);
}

void Fp2_mul_no_red_o2(Fp2l_t c, Fp2_t a, Fp2_t b) {
  Fp2_config.Fp2_mul_no_red_o2(c, a, b, Fp2_config.mem);
}

/* Multiply in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void Fp2_mul_lazy (Fp2_t c, Fp2_t a, Fp2_t b) {
  Fp2_config.Fp2_mul_no_red_o2(Fp2_config.T0, a, b, Fp2_config.mem);
  Fp2l_red(c, Fp2_config.T0);
}

/* Square a Fp2 element without reduction. */
void Fp2_squ_no_red(Fp2l_t c, Fp2_t a) {
  Fp2_config.Fp2_squ_no_red(c, a, Fp2_config.mem);
}

/* Square in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void Fp2_squ_lazy (Fp2_t c, Fp2_t a) {
  Fp2_config.Fp2_squ_no_red(Fp2_config.T0, a, Fp2_config.mem);
  Fp2l_red(c, Fp2_config.T0);
}

/* Multiply an Fp2l element by the special element xi = 1+i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 * We have c = c0 + c1*i = (a0 - a1) + (a0 + a1)*i. */
void p3mod4_Fp2l_mulxi(Fp2l_t c, Fp2l_t a, void *mem) {
	Fpl_t T0;
  T0->limbs = (uint64_t *) mem;

  Fpl_sub_o2_ct (T0, a->A0, a->A1, Fp_config.mem);
	Fpl_add_no_red (c->A1, a->A0, a->A1);
  Fpl_copy(c->A0, T0);
}

/* Multiply two Fp2 elements without reduction resulting in a long Fp2l element.
 * This function assumes that p = 3 (mod4). It implements the function in Alg. 2
 * in Aranha et al., Eurocrypt 2011.
 * This function needs 8 n-word Fp elements for temporary memory in mem.
 */
void p3mod4_Fp2_mul_no_red_o1(Fp2l_t C, Fp2_t a, Fp2_t b, int h, void *mem) {
	Fp_t t0, t1;
  Fpl_t T0, T1, T2;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
    
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;

	Fp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
  Fp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
  Fp_mul_no_red(T2, t0, t1);		    // T2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1

	Fp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	Fp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1

	Fpl_sub_o1_ct(C->A0, T0, T1, h, Fp_config.mem);			// C0 = a0*b0 - a1*b1
  Fpl_add_no_red(T0, T0, T1);		    // T0 = a0*b0 + a1*b1
	Fpl_sub_no_red(C->A1, T2, T0);		// C1 = T2 - T0 = a0*b1 + a1*b0
}

void p3mod4_Fp2_mul_no_red_o2(Fp2l_t C, Fp2_t a, Fp2_t b, void *mem) {
	Fp_t t0, t1;
  Fpl_t T0, T1, T2;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
    
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;

	Fp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
	Fp_mul_no_red(T2, t0, t1);		    // T2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1

	Fp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	Fp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1

	Fpl_sub_o2_ct(C->A0, T0, T1, Fp_config.mem);			// C0 = a0*b0 - a1*b1
	Fpl_add_no_red(T0, T0, T1);		    // T0 = a0*b0 + a1*b1
	Fpl_sub_no_red(C->A1, T2, T0);		// C1 = T2 - T0 = a0*b1 + a1*b0
}

/* Squaring an Fp2 element without reduction resulting in a long Fp2l element. 
 * This function assumes that p = 3 (mod 4). It implements the function in Alg. 7
 * in Aranha et al., Eurocrypt 2011.
 */
void p3mod4_Fp2_squ_no_red(Fp2l_t C, Fp2_t a, void *mem) {
	Fp_t t0, t1, t2;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
	
	Fp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp_sub(t1, a->a0, a->a1);	        // t1 = a0 - a1
	Fp_add_no_red(t2, a->a0, a->a0);  // t2 = 2*a0
	Fp_mul_no_red(C->A0, t0, t1);	    // C0 = (a0 + a1)*(a0 - a1) = a0^2 - a1^2
	Fp_mul_no_red(C->A1, t2, a->a1);  // C1 = 2*a0*a1
}


int Fp2_initialize_config (void) {
   
  Fp2_config.mem = NULL;
  Fp2_config.mem = (void *) malloc (8 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp2_config.mem == NULL) return ERR_OUT_OF_MEMORY;
  
  if (Fp2l_init (Fp2_config.T0) < 0) printf ("Fp2_config memory error.\n");

  Fp2_config.Fp2_mul = p3mod4_Fp2_mul;
  Fp2_config.Fp2_squ = p3mod4_Fp2_squ;
  Fp2_config.Fp2_inv = p3mod4_Fp2_inv;
  Fp2_config.Fp2_mulxi = p3mod4_Fp2_mulxi;
  Fp2_config.Fp2_norm = p3mod4_Fp2_norm;
  Fp2_config.Fp2_muli = p3mod4_Fp2_muli;
  Fp2_config.Fp2_muliFp = p3mod4_Fp2_muliFp;
  Fp2_config.Fp2_mul_no_red_o1 = p3mod4_Fp2_mul_no_red_o1;
  Fp2_config.Fp2_mul_no_red_o2 = p3mod4_Fp2_mul_no_red_o2;
  Fp2_config.Fp2_squ_no_red = p3mod4_Fp2_squ_no_red;
  Fp2_config.Fp2_ppow = binomial_Fp2_ppow;
  Fp2_config.Fp2l_mulxi = p3mod4_Fp2l_mulxi;

  return ERR_SUCCESS;
}

void Fp2_free_config (void) {
  free (Fp2_config.mem);
  Fp2l_free(Fp2_config.T0);
}

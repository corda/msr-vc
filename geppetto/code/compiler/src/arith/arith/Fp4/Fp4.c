/* Fp4/Fp4.c
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

#include "Fp4.h"
#include "../pairing/pairing.h"

Fp4_config_t Fp4_config;

/* Initialize an element in Fp4 by initializing two Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp4_init(Fp4_t dst) {
	int err = 0;
	err = Fp2_init(dst->a0);
	if (err < 0) return err;
	return Fp2_init(dst->a1);
}

/* Free the memory allocated in src. */
void Fp4_free (Fp4_t src) {
	Fp2_free(src->a0);
	Fp2_free(src->a1);
}

/* Get a random element in Fp4. */
void Fp4_rand (Fp4_t c) {
  Fp2_rand(c->a0);
  Fp2_rand(c->a1);
}

/* Set c to 0. */
void Fp4_set_zero (Fp4_t c) {
  Fp2_set_zero (c->a0);
  Fp2_set_zero (c->a1);
}

/* Set c to 1. */
void Fp4_set_one (Fp4_t c) {
  Fp2_set_one (c->a0);
  Fp2_set_zero (c->a1);
}

/* Copy an Fp4 element a to an Fp4 element c. */
void Fp4_copy(Fp4_t c, Fp4_t a) {
	Fp2_copy(c->a0, a->a0);
	Fp2_copy(c->a1, a->a1);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp4_print (Fp4_t a) {
  Fp2_print(a->a0);
  printf("\n");
  Fp2_print(a->a1);
  printf("\n");
}

/* Compare two Fp4 elements a and b for equality. */
int Fp4_cmpeq(Fp4_t a, Fp4_t b) {
  return (Fp2_cmpeq(a->a0, b->a0) && Fp2_cmpeq(a->a1, b->a1));
}

/* Negate an Fp4 element. */
void Fp4_neg(Fp4_t c, Fp4_t a) {
	Fp2_neg(c->a0, a->a0);
	Fp2_neg(c->a1, a->a1);
}

/* Add two Fp4 elements coefficient wise. */
void Fp4_add (Fp4_t c, Fp4_t a, Fp4_t b) {
	Fp2_add(c->a0, a->a0, b->a0);
	Fp2_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp4 element b from an Fp4 element a coefficient wise. */
void Fp4_sub (Fp4_t c, Fp4_t a, Fp4_t b) {
	Fp2_sub(c->a0, a->a0, b->a0);
	Fp2_sub(c->a1, a->a1, b->a1);
}

/* Divide an Fp4 element a by 2 coefficient wise. */
void Fp4_div2 (Fp4_t c, Fp4_t a) {
	Fp2_div2(c->a0, a->a0);
	Fp2_div2(c->a1, a->a1);
}

/* Multiply an Fp4 element by an Fp2 element. */
void Fp4_mulFp2(Fp4_t c, Fp4_t a, Fp2_t b) {
	Fp2_mul(c->a0, a->a0, b);
	Fp2_mul(c->a1, a->a1, b);
}

/* Multiply an Fp4 element by an Fp element. */
void Fp4_mulFp(Fp4_t c, Fp4_t a, Fp_t b) {
	Fp2_mulFp(c->a0, a->a0, b);
	Fp2_mulFp(c->a1, a->a1, b);
}

/* Multiply two Fp4 elements. */
void Fp4_mul (Fp4_t c, Fp4_t a, Fp4_t b) {
  //Fp4_config.Fp4_mul(c, a, b, Fp4_config.mem);
  Fp4_mul_lazy (c, a, b);
}

/* Multiply two Fp4 elements using lazy reduction. */
void Fp4_mul_lazy (Fp4_t c, Fp4_t a, Fp4_t b) {
  Fp4_config.Fp4_mul_lazy (c, a, b, Fp4_config.mem);
}


/* Compute 3*a. */
void Fp4_mul3 (Fp4_t c, Fp4_t a) {
  Fp4_t t;
  
  t->a0->a0->limbs = (uint64_t *) Fp4_config.mem;
  t->a0->a1->limbs = (uint64_t *) Fp4_config.mem + 1*Fp_config.m->n;
  t->a1->a0->limbs = (uint64_t *) Fp4_config.mem + 2*Fp_config.m->n;
  t->a1->a1->limbs = (uint64_t *) Fp4_config.mem + 3*Fp_config.m->n;
  
  Fp4_add (t, a, a);
  Fp4_add (c, t, a);
}

/* Square an Fp4 element. */
void Fp4_squ (Fp4_t c, Fp4_t a) {
  //Fp4_config.Fp4_squ(c, a, Fp4_config.mem);
  Fp4_squ_lazy (c, a);
}

/* Square an Fp4 element using lazy reduction. */
void Fp4_squ_lazy (Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_squ_lazy (c, a, Fp4_config.mem);
}


/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void Fp4_squ_sep (Fp2_t c0, Fp2_t c1, Fp2_t a0, Fp2_t a1) {
  Fp4_config.Fp4_squ_sep(c0, c1, a0, a1, Fp4_config.mem);
}

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2. */
void Fp4_muls (Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_muls(c, a, Fp4_config.mem);
}

/* Multiply an Fp4 element by the special element i in Fp2.*/
void Fp4_muli(Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_muli(c, a);
}

/* Multiply an Fp4 element by the special element i+1 in Fp2.*/
void Fp4_mulxi(Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_mulxi(c, a);
}

/* Compute the p^2-power Frobenius of an Fp4 element. */
void Fp4_p2pow(Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_p2pow(c, a);
}

/* Compute the p-power Frobenius of an Fp4 element. */
void Fp4_ppow(Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_ppow(c, a);
}

/* Invert an Fp4 element. */
void Fp4_inv (Fp4_t c, Fp4_t a) {
  Fp4_config.Fp4_inv(c, a, Fp4_config.mem);
}

/* Multiply two Fp4 elements. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 */
void s2minusxi_Fp4_mul(Fp4_t c, Fp4_t a, Fp4_t b, void *mem) {
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

  Fp2_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp2_add(t1, b->a0, b->a1);  // t1 = b0 + b1
	Fp2_mul(t2, t0, t1);			  // t2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + b0*a1 + a1*b1

	Fp2_mul(t0, a->a0, b->a0);	// t0 = a0*b0
	Fp2_mul(t1, a->a1, b->a1);	// t1 = a1*b1

	Fp2_sub(t2, t2, t0);		    // t2 = t2 - t0 = a0*b1 + a1*b0 + a1*b1
  Fp2_sub(c->a1, t2, t1);		    // t2 = t2 - t1 = a0*b1 + a1*b0
  Fp2_mulxi(t1, t1);		      // t1 = a1*b1*xi
	Fp2_add(c->a0, t0, t1);			// t0 = a0*b0 + a1*b1*xi
	
	}

/* Multiply two Fp4 elements using lazy reduction. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 */
void s2minusxi_Fp4_mul_lazy(Fp4_t c, Fp4_t a, Fp4_t b, void *mem) {
	Fp2_t t0, t1;
  Fp2l_t T2, T3, T4;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;

  T2->A0->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->A1->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T3->A0->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T3->A1->limbs = (uint64_t *) mem + 10*Fp_config.m->n;
  T4->A0->limbs = (uint64_t *) mem + 12*Fp_config.m->n;
  T4->A1->limbs = (uint64_t *) mem + 14*Fp_config.m->n;
  

  Fp2_add_no_red (t0, a->a0, a->a1);	// t0 = a0 + a1
	Fp2_add_no_red (t1, b->a0, b->a1);  // t1 = b0 + b1
	Fp2_mul_no_red_o1 (T2, t0, t1, 2); // T2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + b0*a1 + a1*b1

	Fp2_mul_no_red_o1 (T3, a->a0, b->a0, 2);	// T3 = a0*b0
	Fp2_mul_no_red_o1 (T4, a->a1, b->a1, 2);	// T4 = a1*b1

	Fp2l_sub_o2 (T2, T2, T3);		 // T2 = T2 - T3 = a0*b1 + a1*b0 + a1*b1
  Fp2l_sub_o2 (T2, T2, T4);		 // T2 = T2 - T4 = a0*b1 + a1*b0
  Fp2l_red (c->a1, T2);

  Fp2l_mulxi(T4, T4);		            // T4 = a1*b1*xi
	Fp2l_add_no_red (T2, T3, T4);			// T2 = a0*b0 + a1*b1*xi
  Fp2l_red (c->a0, T2);
	}

/* Square an Fp4 element. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 * This is Algorithm 9 in Beuchat et al., Pairing 2010.
 */
void s2minusxi_Fp4_squ(Fp4_t c, Fp4_t a, void *mem) {
	Fp2_t t0, t1;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	
	Fp2_squ (t0, a->a0);            // t0 = a0^2
  Fp2_squ (t1, a->a1);            // t1 = a1^2
  Fp2_add (c->a1, a->a0, a->a1);  // c1 = a0 + a1
  Fp2_mulxi (c->a0, t1);          // c0 = t1*xi = a1^2*xi
  Fp2_add (c->a0, c->a0, t0);     // c0 = a1^2*xi + a0^2
  Fp2_squ (c->a1, c->a1);         // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  Fp2_sub (c->a1, c->a1, t0);     // 2*a0*a1 + a1^2
  Fp2_sub (c->a1, c->a1, t1);     // 2*a0*a1
}

/* Square an Fp4 element using lazy reduction. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 * This is Algorithm 9 in Beuchat et al., Pairing 2010.
 */
void s2minusxi_Fp4_squ_lazy(Fp4_t c, Fp4_t a, void *mem) {
	Fp2_t t0, t1;
  Fp2l_t T0, T1, T2;
  
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
	
	Fp2_squ_no_red (T0, a->a0);     // T0 = a0^2
  Fp2_squ_no_red (T1, a->a1);     // T1 = a1^2
  Fp2_add (t1, a->a0, a->a1);  // t1 = a0 + a1
  Fp2l_mulxi (T2, T1);              // T2 = t1*xi = a1^2*xi
  Fp2l_add_no_red (T2, T2, T0);     // T2 = a1^2*xi + a0^2
  Fp2l_red (c->a0, T2);

  Fp2_squ_no_red (T2, t1);       // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  Fp2l_sub_o2 (T2, T2, T0);     // 2*a0*a1 + a1^2
  Fp2l_sub_o2 (T2, T2, T1);     // 2*a0*a1
  Fp2l_red (c->a1, T2);
}


/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void s2minusxi_Fp4_squ_sep (Fp2_t c0, Fp2_t c1, Fp2_t a0, Fp2_t a1, void *mem) {
	Fp2_t t0, t1;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
	
	Fp2_squ (t0, a0);         // t0 = a0^2
  Fp2_squ (t1, a1);         // t1 = a1^2
  Fp2_add (c1, a0, a1);     // c1 = a0 + a1
  Fp2_mulxi (c0, t1);       // c0 = t1*xi = a1^2*xi
  Fp2_add (c0, c0, t0);     // c0 = a1^2*xi + a0^2
  Fp2_squ (c1, c1);         // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  Fp2_sub (c1, c1, t0);     // 2*a0*a1 + a1^2
  Fp2_sub (c1, c1, t1);     // 2*a0*a1
}

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2.
 * c = c0 + c1*s = (a0 + a1*s)*s = a1*xi + a0*s.
 */
void s2minusxi_Fp4_muls(Fp4_t c, Fp4_t a, void *mem) {
  Fp2_t t0;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
	
  Fp2_mulxi(t0, a->a1);
	Fp2_copy(c->a1, a->a0);
  Fp2_copy(c->a0, t0);
}

/* Multiply an Fp4 element by the special element i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 */
void p3mod4_Fp4_muli(Fp4_t c, Fp4_t a) {
  Fp2_muli (c->a0, a->a0);
  Fp2_muli (c->a1, a->a1);
}

/* Multiply an Fp4 element by the special element xi.
 */
void p3mod4_Fp4_mulxi(Fp4_t c, Fp4_t a) {
  Fp2_mulxi (c->a0, a->a0);
  Fp2_mulxi (c->a1, a->a1);
}

/* Inversion of an Fp4 element. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 */
void s2minusxi_Fp4_inv(Fp4_t c, Fp4_t a, void *mem) {
	Fp2_t t0, t1;
  
  t0->a0->limbs = (uint64_t *) mem;
  t0->a1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t1->a0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t1->a1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
   
	Fp2_squ (t0, a->a0);         // t0 = a0^2
	Fp2_squ (t1, a->a1);         // t1 = a1^2
	Fp2_mulxi (t1, t1);          // t1 = xi*a1^2
  Fp2_sub (t0, t0, t1);        // t0 = a0^2 - xi*a1^2 = N(a)

	Fp2_inv (t0, t0);            // t0 = 1/N(a)

	Fp2_mul (c->a0, a->a0, t0);  // c0 = a0/N(a)
	Fp2_mul (c->a1, a->a1, t0);  
  Fp2_neg (c->a1, c->a1);      // c1 = -a1/N(a)
}


/* Compute the p^2-power Frobenius of an Fp4 element. This function assumes that the field extension
 * is constructed via a binomial s^2 - xi, xi in Fp2.
 * This is the same as (complex) conjugation, i.e. Fp4_p2pow(a0 + a1*s) = a0 - a1*s.
 * When we need this, it might be simpler to just negate the second coordinate in place.
 */
void binomial_Fp4_p2pow(Fp4_t c, Fp4_t a) {
	Fp2_copy(c->a0, a->a0);
	Fp2_neg(c->a1, a->a1);
}

/* Compute the p-power Frobenius of an Fp4 element. This function assumes that the field extension
 * is constructed via a binomial s^2 - xi, xi in Fp2.
 * This is Fp4_ppow(a0 + a1*s) = a0^p + (a1^p*sppow)*s, where sppow = s^(p-1) in Fp2.
 */
void binomial_Fp4_ppow(Fp4_t c, Fp4_t a) {
	Fp2_ppow (c->a0, a->a0);
	Fp2_ppow (c->a1, a->a1);
  Fp2_mul (c->a1, c->a1, Fp4_config.sppow);
}

int Fp4_initialize_config (void) {
  Fp4_config.mem = NULL;
  Fp4_config.mem = (void *) malloc (24 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp4_config.mem == NULL) return ERR_OUT_OF_MEMORY;
  
  if (Fp2_init (Fp4_config.sppow) < 0) printf ("Fp4_config memory error.\n");

  Fp4_config.Fp4_mul = s2minusxi_Fp4_mul;
  Fp4_config.Fp4_squ = s2minusxi_Fp4_squ;
  Fp4_config.Fp4_squ_sep = s2minusxi_Fp4_squ_sep;
  Fp4_config.Fp4_muls = s2minusxi_Fp4_muls;
  Fp4_config.Fp4_inv = s2minusxi_Fp4_inv;
  Fp4_config.Fp4_p2pow = binomial_Fp4_p2pow;
  Fp4_config.Fp4_ppow = binomial_Fp4_ppow;
  Fp4_config.Fp4_muli = p3mod4_Fp4_muli;
  Fp4_config.Fp4_mulxi = p3mod4_Fp4_mulxi;
  Fp4_config.Fp4_mul_lazy = s2minusxi_Fp4_mul_lazy;
  Fp4_config.Fp4_squ_lazy = s2minusxi_Fp4_squ_lazy;

  return ERR_SUCCESS;
}

void Fp4_free_config (void) {
  Fp2_free (Fp4_config.sppow);
  free (Fp4_config.mem);
}


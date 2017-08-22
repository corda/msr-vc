/* Fp3/Fp3.c
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

#include "Fp3.h"
#include "../pairing/pairing.h"

Fp3_config_t Fp3_config;

/* Initialize an element in Fp3 by initializing three Fp elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int Fp3_init(Fp3_t dst) {
	int err = 0;
	err = Fp_init(dst->a0);
	if (err < 0) return err;
	err = Fp_init(dst->a1);
	if (err < 0) return err;
	return Fp_init(dst->a2);
}

/* Free the memory allocated in src. */
void Fp3_free (Fp3_t src) {
	Fp_free(src->a0);
	Fp_free(src->a1);
	Fp_free(src->a2);
}

/* Get a random element in Fp3. */
void Fp3_rand (Fp3_t c) {
  Fp_rand(c->a0);
  Fp_rand(c->a1);
  Fp_rand(c->a2);
}

/* Set c to 0. */
void Fp3_set_zero (Fp3_t c) {
  Fp_set_ui (c->a0, 0);
  Fp_set_ui (c->a1, 0);
  Fp_set_ui (c->a2, 0);
}

/* Set c to 1. */
void Fp3_set_one (Fp3_t c) {
  Fp_set_ui (c->a0, 1);
  Fp_set_ui (c->a1, 0);
  Fp_set_ui (c->a2, 0);
}

/* Copy an Fp3 element a to an Fp3 element c. */
void Fp3_copy(Fp3_t c, Fp3_t a) {
	Fp_copy(c->a0, a->a0);
	Fp_copy(c->a1, a->a1);
	Fp_copy(c->a2, a->a2);
}

/* Print the (regular non-Montgomery form) value of a */
void Fp3_print (Fp3_t a) {
  Fp_print(a->a0);
  printf("\n");
  Fp_print(a->a1);
  printf("\n");
  Fp_print(a->a2);
  printf("\n");
}

/* Compare two Fp3 elements a and b for equality. */
int Fp3_cmpeq(Fp3_t a, Fp3_t b) {
  return (Fp_cmp(a->a0, b->a0) == 0 
       && Fp_cmp(a->a1, b->a1) == 0 && Fp_cmp(a->a2, b->a2) == 0);
}

/* Negate an Fp3 element. */
void Fp3_neg(Fp3_t c, Fp3_t a) {
	Fp_neg(c->a0, a->a0);
	Fp_neg(c->a1, a->a1);
	Fp_neg(c->a2, a->a2);
}

/* Add two Fp3 elements coefficient wise. */
void Fp3_add (Fp3_t c, Fp3_t a, Fp3_t b) {
	Fp_add(c->a0, a->a0, b->a0);
	Fp_add(c->a1, a->a1, b->a1);
	Fp_add(c->a2, a->a2, b->a2);
}

/* Subtract an Fp3 element b from an Fp3 element a coefficient wise. */
void Fp3_sub (Fp3_t c, Fp3_t a, Fp3_t b) {
	Fp_sub(c->a0, a->a0, b->a0);
	Fp_sub(c->a1, a->a1, b->a1);
	Fp_sub(c->a2, a->a2, b->a2);
}

/* Multiply an Fp3 element by an Fp element. */
void Fp3_mulFp(Fp3_t c, Fp3_t a, Fp_t b) {
	Fp_mul(c->a0, a->a0, b);
	Fp_mul(c->a1, a->a1, b);
	Fp_mul(c->a2, a->a2, b);
}

/* Compute 3*a. */
void Fp3_mul3 (Fp3_t c, Fp3_t a) {
  Fp3_t t;
  
  t->a0->limbs = (uint64_t *) Fp3_config.mem;
  t->a1->limbs = (uint64_t *) Fp3_config.mem + 1*Fp_config.m->n;
  t->a2->limbs = (uint64_t *) Fp3_config.mem + 2*Fp_config.m->n;
    
  Fp3_add (t, a, a);
  Fp3_add (c, t, a);
}

/* Multiply two Fp3 elements. */
void Fp3_mul (Fp3_t c, Fp3_t a, Fp3_t b) {
  Fp3_config.Fp3_mul(c, a, b, Fp3_config.mem);
  //Fp3_mul_lazy (c, a, b);
}

/* Square an Fp3 element. */
void Fp3_squ (Fp3_t c, Fp3_t a) {
  Fp3_config.Fp3_squ(c, a, Fp3_config.mem);
  //Fp3_squ_lazy ( c, a);
}

/* Invert an Fp3 element. */
void Fp3_inv (Fp3_t c, Fp3_t a) {
  Fp3_config.Fp3_inv(c, a, Fp3_config.mem);
}

/* Multiply an Fp3 element by the special element u. */
void Fp3_mulu(Fp3_t c, Fp3_t a) {
  Fp3_config.Fp3_mulu(c, a, Fp3_config.mem);
}

/* Compute the p-power Frobenius of an Fp3 element. */
void Fp3_ppow(Fp3_t c, Fp3_t a) {
  Fp3_config.Fp3_ppow(c, a);
}


/* Multiply two Fp3 elements. This uses that Fp3 is constructed over Fp via u^3 + 2. 
 * This is Algorithm 13 in Beuchat et al., Pairing 2010.
 */
void u3plus2_Fp3_mul(Fp3_t c, Fp3_t a, Fp3_t b, void *mem) {
	Fp_t t0, t1, t2, t3, t4, t5;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t3->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t4->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t5->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  
	Fp_mul(t0, a->a0, b->a0);    // t0 = a0*b0
	Fp_mul(t1, a->a1, b->a1);    // t1 = a1*b1
	Fp_mul(t2, a->a2, b->a2);    // t2 = a2*b2

	Fp_add(t3, a->a1, a->a2);    // t3 = a1 + a2
	Fp_add(t4, b->a1, b->a2);    // t4 = b1 + b2
	Fp_mul(t3, t3, t4);          // t3 = (a1 + a2)*(b1 + b2) = a1*b1 + a2*b1 + a1*b2 + a2*b2
	Fp_sub(t3, t3, t1);          // t3 = a2*b1 + a1*b2 + a2*b2
	Fp_sub(t3, t3, t2);          // t3 = a2*b1 + a1*b2 
  Fp_add(t3, t3, t3);					 // t3 = (a2*b1 + a1*b2)*2	
                               // t3 now almost has the result for c0, t4 is free.		

	Fp_add(t4, a->a0, a->a2);    // t4 = a0 + a2
	Fp_add(t5, b->a0, b->a2);    // t5 = b0 + b2
	Fp_mul(t4, t4, t5);          // t4 = (a0 + a2)*(b0 + b2) = a0*b0 + a2*b0 + a0*b2 + a2*b2
	Fp_sub(t4, t4, t0);          // t4 = a2*b0 + a0*b2 + a2*b2
	Fp_sub(t4, t4, t2);          // t4 = a2*b0 + a0*b2
	Fp_add(c->a2, t4, t1);			 // c2 = a2*b0 + a0*b2 + a1*b1 	
                               // c2 is done, t4 and t5 are free.

	Fp_add(t4, a->a0, a->a1);    // t4 = a0 + a1
	Fp_add(t5, b->a0, b->a1);    // t4 = b0 + b1
	Fp_mul(t4, t4, t5);          // t4 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1
	Fp_sub(t4, t4, t0);          // t4 = a1*b0 + a0*b1 + a1*b1
	Fp_sub(t4, t4, t1);          // t4 = a1*b0 + a0*b1
	Fp_add(t2, t2, t2);          // t2 = a2*b2 * 2
  Fp_sub(c->a1, t4, t2);			 // c1 = a1*b0 + a0*b1 - a2*b2*2	
  // c1 is done.

	//Fp_add(c->a0, t3, t0);			 // c0 is done.
  Fp_sub (c->a0, t0, t3);      // c0 is done.
}

/* Multiply two Fp3 elements. This uses that Fp3 is constructed over Fp via u^3 - 3. 
 */
void u3minus3_Fp3_mul(Fp3_t c, Fp3_t a, Fp3_t b, void *mem) {
	Fp_t t0, t1, t2, t3, t4, t5;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t3->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t4->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t5->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  
	Fp_mul(t0, a->a0, b->a0);    // t0 = a0*b0
	Fp_mul(t1, a->a1, b->a1);    // t1 = a1*b1
	Fp_mul(t2, a->a2, b->a2);    // t2 = a2*b2

	Fp_add(t3, a->a1, a->a2);    // t3 = a1 + a2
	Fp_add(t4, b->a1, b->a2);    // t4 = b1 + b2
	Fp_mul(t3, t3, t4);          // t3 = (a1 + a2)*(b1 + b2) = a1*b1 + a2*b1 + a1*b2 + a2*b2
	Fp_sub(t3, t3, t1);          // t3 = a2*b1 + a1*b2 + a2*b2
	Fp_sub(t3, t3, t2);          // t3 = a2*b1 + a1*b2 
  Fp_add(t4, t3, t3);					 // t3 = (a2*b1 + a1*b2)*2
  Fp_add(t3, t4, t3);					 // t3 = (a2*b1 + a1*b2)*3
                               // t3 now almost has the result for c0, t4 is free.		

	Fp_add(t4, a->a0, a->a2);    // t4 = a0 + a2
	Fp_add(t5, b->a0, b->a2);    // t5 = b0 + b2
	Fp_mul(t4, t4, t5);          // t4 = (a0 + a2)*(b0 + b2) = a0*b0 + a2*b0 + a0*b2 + a2*b2
	Fp_sub(t4, t4, t0);          // t4 = a2*b0 + a0*b2 + a2*b2
	Fp_sub(t4, t4, t2);          // t4 = a2*b0 + a0*b2
	Fp_add(c->a2, t4, t1);			 // c2 = a2*b0 + a0*b2 + a1*b1 	
                               // c2 is done, t4 and t5 are free.

	Fp_add(t4, a->a0, a->a1);    // t4 = a0 + a1
	Fp_add(t5, b->a0, b->a1);    // t4 = b0 + b1
	Fp_mul(t4, t4, t5);          // t4 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1
	Fp_sub(t4, t4, t0);          // t4 = a1*b0 + a0*b1 + a1*b1
	Fp_sub(t4, t4, t1);          // t4 = a1*b0 + a0*b1
	Fp_add(c->a1, t2, t2);       // c1 = a2*b2 * 2
  Fp_add(c->a1, c->a1, t2);    // c1 = a2*b2 * 3
  Fp_add(c->a1, t4, c->a1);		 // c1 = a1*b0 + a0*b1 + a2*b2*3	
  // c1 is done.

	Fp_add(c->a0, t3, t0);			 // c0 is done.
}

/* Square an Fp3 element. This uses that Fp3 is constructed over Fp via u^3 + 2. 
 * This is Algorithm 16 in Beuchat et al., Pairing 2010.
 */
void u3plus2_Fp3_squ(Fp3_t c, Fp3_t a, void *mem) {
	Fp_t t0, t1, t2, t3, t4;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t3->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t4->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	
  Fp_mul(t0, a->a0, a->a1);    // t0 = a0*a1  
	Fp_add(t0, t0, t0);          // t0 = 2*a0*a1
	Fp_mul(t1, a->a2, a->a2);    // t1 = a2^2
	Fp_mul(t2, a->a0, a->a0);    // t2 = a0^2
	Fp_add(t3, a->a0, a->a1);    // t3 = a0 + a1
	Fp_add(t3, t3, a->a2);       // t3 = a0 + a1 + a2
	Fp_mul(t4, a->a1, a->a2);    // t4 = a1*a2
	Fp_add(t4, t4, t4);          // t4 = 2*a1*a2

	Fp_add(c->a2, t1, t1);       // c2 = 2*a2^2
  Fp_sub(c->a1, t0, c->a2);    // c1 = 2*a0*a1 - 2*a2^2

	Fp_mul(c->a2, t3, t3);       // t0 = (a0 + a1 + a2)^2
	Fp_sub(c->a2, c->a2, t0);    // c2 = a0^2 + a1^2 + a2^2 + 2*a1*a2 + 2*a0*a2
	Fp_sub(c->a2, c->a2, t4);    // c2 = a0^2 + a1^2 + a2^2 + 2*a0*a2
	Fp_sub(c->a2, c->a2, t2);    // c2 = a1^2 + a2^2 + 2*a0*a2
  Fp_sub(c->a2, c->a2, t1);    // c2 = a1^2 + 2*a0*a2

	Fp_add(t4, t4, t4);          // t4 = 4*a1*a2
  Fp_sub(c->a0, t2, t4);       // c0 = a0^2 - 4*a1*a2
}

/* Square an Fp3 element. This uses that Fp3 is constructed over Fp via u^3 - 3. 
 */
void u3minus3_Fp3_squ(Fp3_t c, Fp3_t a, void *mem) {
	Fp_t t0, t1, t2, t3, t4;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t3->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t4->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
	
  Fp_mul(t0, a->a0, a->a1);    // t0 = a0*a1  
	Fp_add(t0, t0, t0);          // t0 = 2*a0*a1
	Fp_mul(t1, a->a2, a->a2);    // t1 = a2^2
	Fp_mul(t2, a->a0, a->a0);    // t2 = a0^2
	Fp_add(t3, a->a0, a->a1);    // t3 = a0 + a1
	Fp_add(t3, t3, a->a2);       // t3 = a0 + a1 + a2
	Fp_mul(t4, a->a1, a->a2);    // t4 = a1*a2
	Fp_add(t4, t4, t4);          // t4 = 2*a1*a2

	Fp_add(c->a2, t1, t1);       // c2 = 2*a2^2
	Fp_add(c->a2, c->a2, t1);    // c2 = 3*a2^2
  Fp_add(c->a1, t0, c->a2);    // c1 = 2*a0*a1 + 3*a2^2

	Fp_mul(c->a2, t3, t3);       // t0 = (a0 + a1 + a2)^2
	Fp_sub(c->a2, c->a2, t0);    // c2 = a0^2 + a1^2 + a2^2 + 2*a1*a2 + 2*a0*a2
	Fp_sub(c->a2, c->a2, t4);    // c2 = a0^2 + a1^2 + a2^2 + 2*a0*a2
	Fp_sub(c->a2, c->a2, t2);    // c2 = a1^2 + a2^2 + 2*a0*a2
  Fp_sub(c->a2, c->a2, t1);    // c2 = a1^2 + 2*a0*a2

	Fp_add(t3, t4, t4);          // t4 = 4*a1*a2
  Fp_add(t4, t4, t3);          // t4 = 6*a1*a2
  Fp_add(c->a0, t2, t4);       // c0 = a0^2 + 6*a1*a2
}


/* Multiplication by the special element v in Fp3, with v^3 = -2 in Fp.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = -2*a2 + a0*v + a1*v^2.
 */
void u3plus2_Fp3_mulu(Fp3_t c, Fp3_t a, void *mem) {
  Fp_t t0;
  
  t0->limbs = (uint64_t *) mem;
  
  Fp_add(t0, a->a2, a->a2);
  Fp_copy(c->a2, a->a1);
	Fp_copy(c->a1, a->a0);
  Fp_neg(c->a0, t0);
}

/* Multiplication by the special element v in Fp3, with v^3 = 3 in Fp.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = 3*a2 + a0*v + a1*v^2.
 */
void u3minus3_Fp3_mulu(Fp3_t c, Fp3_t a, void *mem) {
  Fp_t t0;
  
  t0->limbs = (uint64_t *) mem;
  
  Fp_add(t0, a->a2, a->a2);
  Fp_add(t0, t0, a->a2);
  Fp_copy(c->a2, a->a1);
	Fp_copy(c->a1, a->a0);
  Fp_copy(c->a0, t0);
}


/* Compute the p-power Frobenius of an Fp3 element. This uses a constant VPPOW, 
 * such that v^p = VPPOW*v, i.e. VPPOW = v^(p-1). In the case, where Fp3 is an 
 * extension of Fp via v^3 - xi and p = 1 (mod 3), we have VPPOW = xi^(p-1)/3 in Fp. 
 * VPPOW2 = VPPOW^2.
 */
void p1mod3_Fp3_ppow(Fp3_t c, Fp3_t a) {
  Fp_copy(c->a0, a->a0);
	Fp_mul(c->a1, a->a1, Fp3_config.uppow);
	Fp_mul(c->a2, a->a2, Fp3_config.uppow2);
}


/* Inversion of an Fp3 element. This uses that Fp3 is constructed over Fp via v^3 - xi. 
 * This is Algorithm 17 in Beuchat et al., Pairing 2010.
 */
void u3plus2_Fp3_inv(Fp3_t c, Fp3_t a, void *mem) {
	Fp_t t0, t1, t2, t3, t4, t5;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t3->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t4->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t5->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  
	Fp_mul(t0, a->a0, a->a0);         // t0 = a0^2
	Fp_mul(t1, a->a1, a->a1);         // t1 = a1^2
	Fp_mul(t2, a->a2, a->a2);         // t2 = a2^2

	Fp_mul(t3, a->a0, a->a1);  // t3 = a0*a1
	Fp_mul(t4, a->a0, a->a2);  // t4 = a0*a2
	Fp_mul(t5, a->a1, a->a2);  // t5 = a1*a2

	Fp_add(t2, t2, t2);          // t2 = 2*a2^2
	Fp_add(t5, t5, t5);          // t5 = 2*a1*a2

	//Fp_sub(t0, t0, t5);        // t0 = a0^2 - xi*a1*a2
  Fp_add(t0, t0, t5);        // t0 = a0^2 - xi*a1*a2
	Fp_sub(t1, t1, t4);        // t1 = a1^2 - a0*a2
	//Fp_sub(t2, t2, t3);        // t2 = xi*a2^2 - a0*a1
  Fp_add(t2, t2, t3);        // t2 = xi*a2^2 - a0*a1
  Fp_neg (t2, t2);

	Fp_mul(t3, a->a0, t0);     // t3 = a0^3 - xi*a0*a1*a2
	Fp_mul(t4, a->a2, t2);     // t4 = xi*a2^3 - a0*a1*a2
	Fp_add(t4, t4, t4);        // t4 = xi^2*a2^3 - xi*a0*a1*a2
	//Fp_add(t3, t3, t4);        // t3 = a0^3 + xi^2*a2^3 - 2*xi*a0*a1*a2 
  Fp_sub(t3, t3, t4);        // t3 = a0^3 + xi^2*a2^3 - 2*xi*a0*a1*a2 
	Fp_mul(t4, a->a1, t1);     // t4 = a1^3 - a0*a1*a2
	Fp_add(t4, t4, t4);          // t4 = xi*a1^3 - xi*a0*a1*a2
	//Fp_add(t3, t3, t4);        // t3 = a0^3 + xi*a1^3 + xi^2*a2^3 - 3*xi*a0*a1*a2
  Fp_sub(t3, t3, t4);        // t3 = a0^3 + xi*a1^3 + xi^2*a2^3 - 3*xi*a0*a1*a2

	Fp_modinv(t3, t3);            // t3 = 1/N(a)

	Fp_mul(c->a0, t0, t3);
	Fp_mul(c->a1, t2, t3);
	Fp_mul(c->a2, t1, t3);
}

/* Inversion of an Fp3 element. This uses that Fp3 is constructed over Fp via v^3 - xi. 
 */
void u3minus3_Fp3_inv(Fp3_t c, Fp3_t a, void *mem) {
	Fp_t t0, t1, t2, t3, t4, t5, t6;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  t2->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  t3->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  t4->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  t5->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  t6->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  
	Fp_mul(t0, a->a0, a->a0);         // t0 = a0^2
	Fp_mul(t1, a->a1, a->a1);         // t1 = a1^2
	Fp_mul(t2, a->a2, a->a2);         // t2 = a2^2

	Fp_mul(t3, a->a0, a->a1);  // t3 = a0*a1
	Fp_mul(t4, a->a0, a->a2);  // t4 = a0*a2
	Fp_mul(t5, a->a1, a->a2);  // t5 = a1*a2

	Fp_add(t6, t2, t2);          // t2 = 2*a2^2
	Fp_add(t2, t6, t2);          // t2 = 3*a2^2
	Fp_add(t6, t5, t5);          // t5 = 2*a1*a2
  Fp_add(t5, t6, t5);          // t5 = 3*a1*a2

	Fp_sub(t0, t0, t5);        // t0 = a0^2 - 3*a1*a2
	Fp_sub(t1, t1, t4);        // t1 = a1^2 - a0*a2
	Fp_sub(t2, t2, t3);        // t2 = 3*a2^2 - a0*a1

	Fp_mul(t3, a->a0, t0);     // t3 = a0^3 - 3*a0*a1*a2
	Fp_mul(t4, a->a2, t2);     // t4 = 3*a2^3 - a0*a1*a2
	Fp_add(t6, t4, t4);        // t4 = 6*a2^3 - 2*a0*a1*a2
	Fp_add(t4, t6, t4);        // t4 = 9*a2^3 - 3*a0*a1*a2
  Fp_add(t3, t3, t4);        // t3 = a0^3 + 9*a2^3 - 6*a0*a1*a2 
	Fp_mul(t4, a->a1, t1);     // t4 = a1^3 - a0*a1*a2
	Fp_add(t6, t4, t4);        // t6 = 2*a1^3 - 2*a0*a1*a2
  Fp_add(t4, t6, t4);        // t4 = 3*a1^3 - 3*a0*a1*a2
	Fp_add(t3, t3, t4);        // t3 = a0^3 + 3*a1^3 + 9*a2^3 - 9*a0*a1*a2

	Fp_modinv(t3, t3);         // t3 = 1/N(a)

	Fp_mul(c->a0, t0, t3);
	Fp_mul(c->a1, t2, t3);
	Fp_mul(c->a2, t1, t3);
}

/************** Fp3 functions for lazy reduction ***************/
/* Functions assume that there is enough space such that no overflow occurs. */

/* Initialize a long Fp3 element. */
int Fp3l_init(Fp3l_t dst) {
	int err = 0;
	if (err = Fpl_init(dst->A0) < 0) return err;
	if (err = Fpl_init(dst->A1) < 0) return err;
	return Fpl_init(dst->A2);
}

/* Free the memory allocated in src. */
void Fp3l_free (Fp3l_t src) {
	Fpl_free(src->A0);
	Fpl_free(src->A1);
  Fpl_free(src->A2);
}

/* Copy an Fp3l element a to an Fp3l element c. */
void Fp3l_copy(Fp3l_t C, Fp3l_t A) {
	Fpl_copy(C->A0, A->A0);
	Fpl_copy(C->A1, A->A1);
  Fpl_copy(C->A2, A->A2);
}

/* Reduction function calling the Fp reduction for each coefficient. */
void Fp3l_red(Fp3_t c, Fp3l_t A) {
	Fpl_red(c->a0, A->A0);
	Fpl_red(c->a1, A->A1);
  Fpl_red(c->a2, A->A2);
}

/* Add two Fp3 elements coefficient wise without reducing mod p for lazy reduction. */ 
void Fp3_add_no_red (Fp3_t c, Fp3_t a, Fp3_t b) {
	Fp_add_no_red(c->a0, a->a0, b->a0);
	Fp_add_no_red(c->a1, a->a1, b->a1);
  Fp_add_no_red(c->a2, a->a2, b->a2);
}

/* Add two long Fp3 elements coefficient wise without reducing mod p for lazy reduction. */ 
void Fp3l_add_no_red (Fp3l_t C, Fp3l_t A, Fp3l_t B) {
	Fpl_add_no_red(C->A0, A->A0, B->A0);
	Fpl_add_no_red(C->A1, A->A1, B->A1);
	Fpl_add_no_red(C->A2, A->A2, B->A2);
}

/* Subtraction of Fp3l elements, option 1 in Aranha et al. */
void Fp3l_sub_o1 (Fp3l_t c, Fp3l_t a, Fp3l_t b, int h) {
  Fpl_sub_o1_ct(c->A0, a->A0, b->A0, h, Fp_config.mem);
  Fpl_sub_o1_ct(c->A1, a->A1, b->A1, h, Fp_config.mem);
  Fpl_sub_o1_ct(c->A2, a->A2, b->A2, h, Fp_config.mem);
}

/* Subtraction of Fp3l elements, constant-time version of option 2 in Aranha et al. */
void Fp3l_sub_o2 (Fp3l_t c, Fp3l_t a, Fp3l_t b) {
  Fpl_sub_o2_ct(c->A0, a->A0, b->A0, Fp_config.mem);
  Fpl_sub_o2_ct(c->A1, a->A1, b->A1, Fp_config.mem);
  Fpl_sub_o2_ct(c->A2, a->A2, b->A2, Fp_config.mem);
}

/* Addition of Fp3l elements, constant-time version of "option 2" addition in Aranha et al. */
void Fp3l_add_o2 (Fp3l_t c, Fp3l_t a, Fp3l_t b) {
  Fpl_add_o2_ct(c->A0, a->A0, b->A0, Fp_config.mem);
  Fpl_add_o2_ct(c->A1, a->A1, b->A1, Fp_config.mem);
  Fpl_add_o2_ct(c->A2, a->A2, b->A2, Fp_config.mem);
}

/* Multiply an Fp3 element by an Fp element. */
void Fp3_mulFp_no_red (Fp3l_t c, Fp3_t a, Fp_t b) {
	Fp_mul_no_red (c->A0, a->a0, b);
	Fp_mul_no_red (c->A1, a->a1, b);
	Fp_mul_no_red (c->A2, a->a2, b);
}

/* Multiplication by the special element v in Fp3, with v^3 = 3 in Fp.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = 3*a2 + a0*v + a1*v^2.
 */
void Fp3l_mulu (Fp3l_t c, Fp3l_t a) {
  Fp3_config.Fp3l_mulu(c, a, Fp3_config.mem);
}

/* Multiply two Fp3 elements without reduction. */
void Fp3_mul_no_red_o1(Fp3l_t c, Fp3_t a, Fp3_t b, int h) {
  Fp3_config.Fp3_mul_no_red_o1(c, a, b, h, Fp3_config.mem);
}

void Fp3_mul_no_red(Fp3l_t c, Fp3_t a, Fp3_t b) {
  Fp3_config.Fp3_mul_no_red(c, a, b, Fp3_config.mem);
}

void Fp3_mul_no_red_o2(Fp3l_t c, Fp3_t a, Fp3_t b) {
  Fp3_config.Fp3_mul_no_red_o2(c, a, b, Fp3_config.mem);
}

/* Multiply in Fp3 using the lazy reduction multiplications and 3 reductions only. */
void Fp3_mul_lazy (Fp3_t c, Fp3_t a, Fp3_t b) {
  Fp3_config.Fp3_mul_no_red_o2(Fp3_config.T0, a, b, Fp3_config.mem);
  Fp3l_red(c, Fp3_config.T0);
}

/* Square an Fp3 element without reduction. */
void Fp3_squ_no_red_o1(Fp3l_t c, Fp3_t a, int h) {
  Fp3_config.Fp3_squ_no_red_o1(c, a, h, Fp3_config.mem);
}

void Fp3_squ_no_red_o2(Fp3l_t c, Fp3_t a) {
  Fp3_config.Fp3_squ_no_red_o2(c, a, Fp3_config.mem);
}

void Fp3_squ_no_red(Fp3l_t c, Fp3_t a) {
  Fp3_config.Fp3_squ_no_red(c, a, Fp3_config.mem);
}


/* Square in Fp3 using the lazy reduction multiplications and 3 reductions only. */
void Fp3_squ_lazy (Fp3_t c, Fp3_t a) {
  Fp3_config.Fp3_squ_no_red_o2(Fp3_config.T0, a, Fp3_config.mem);
  Fp3l_red(c, Fp3_config.T0);
}

/* Multiply an Fp3 element and a sparse Fp3 element without reduction.
   This function assumes that Fp3 is constructed via v^3 - xi.
 */
void Fp3_mul_sparse01_no_red (Fp3l_t c, Fp3_t a, Fp_t b0, Fp_t b1){
  Fp3_config.Fp3_mul_sparse01_no_red (c, a, b0, b1, Fp3_config.mem);
}
void Fp3_mul_sparse12_no_red (Fp3l_t c, Fp3_t a, Fp_t b1, Fp_t b2){
  Fp3_config.Fp3_mul_sparse12_no_red (c, a, b1, b2, Fp_config.mem);
}

/* Multiplication by the special element v in Fp3, with v^3 = 3 in Fp.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = 3*a2 + a0*v + a1*v^2.
 */
void u3minus3_Fp3l_mulu(Fp3l_t c, Fp3l_t a, void *mem) {
  Fpl_t T0;
  
  T0->limbs = (uint64_t *) mem;
  
  Fpl_add_no_red(T0, a->A2, a->A2);
  Fpl_add_no_red(T0, T0, a->A2);
  Fpl_copy(c->A2, a->A1);
	Fpl_copy(c->A1, a->A0);
  Fpl_copy(c->A0, T0);
}

/* Multiply two Fp3 elements without reduction resulting in a long Fp3l element.
 * This function assumes that Fp3 = Fp(u) where u^3 = -2. 
 * This function needs 12 n-word Fp elements for temporary memory in mem.
 */
void u3plus2_Fp3_mul_no_red_o1(Fp3l_t C, Fp3_t a, Fp3_t b, int h, void *mem) {
	Fp_t t0, t1;
  Fpl_t T0, T1, T2, T3, T4;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
    
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T4->limbs = (uint64_t *) mem + 10*Fp_config.m->n;

  Fp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	Fp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1
  Fp_mul_no_red(T2, a->a2, b->a2);	// T2 = a2*b2

	Fp_add_no_red(t0, a->a1, a->a2);	// t0 = a1 + a2
  Fp_add_no_red(t1, b->a1, b->a2);  // t1 = b1 + b2
  Fp_mul_no_red(T3, t0, t1);		    // T3 = (a1 + a2)*(b1 + b2) = a1*b1 + a1*b2 + a2*b1 + a2*b2
  Fpl_add_no_red(T4, T1, T2);		    // T4 = a1*b1 + a2*b2
	Fpl_sub_no_red(T3, T3, T4);		    // T3 = T3 - T4 = a1*b2 + a2*b1
  Fpl_add_no_red(T3, T3, T3);       // T3 = 2*(a1*b2 + a2*b1)
	
  Fp_add_no_red(t0, a->a0, a->a2);	// t0 = a0 + a2
  Fp_add_no_red(t1, b->a0, b->a2);  // t1 = b0 + b2
  Fp_mul_no_red(T4, t0, t1);		    // T4 = (a0 + a2)*(b0 + b2) = a0*b0 + a0*b2 + a2*b0 + a2*b2
  Fpl_sub_no_red(T4, T4, T0);
  Fpl_sub_no_red(T4, T4, T2);       // T4 = a0*b2 + a2*b0
  Fpl_add_no_red(C->A2, T4, T1);    // C2 = a0*b2 + a2*b0 + a1*b1
  
  Fp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
  Fp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
  Fp_mul_no_red(T4, t0, t1);		    // T4 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1
  Fpl_sub_no_red(T4, T4, T0);
  Fpl_sub_no_red(T4, T4, T1);       // T4 = a0*b1 + a1*b0
  Fpl_add_no_red(T2, T2, T2);       // T2 = 2*a2*b2
  Fpl_sub_o1_ct(C->A1, T4, T2, h, Fp_config.mem); // C1 = a0*b1 + a1*b0 - 2*a2*b2

  Fpl_sub_o1_ct(C->A0, T0, T3, h, Fp_config.mem); // C0 = a0*b0 - 2*(a1*b2 + a2*b1)
}

void u3plus2_Fp3_mul_no_red_o2(Fp3l_t C, Fp3_t a, Fp3_t b, void *mem) {
	Fp_t t0, t1;
  Fpl_t T0, T1, T2, T3, T4;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
    
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T4->limbs = (uint64_t *) mem + 10*Fp_config.m->n;

  Fp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	Fp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1
  Fp_mul_no_red(T2, a->a2, b->a2);	// T2 = a2*b2

	Fp_add_no_red(t0, a->a1, a->a2);	// t0 = a1 + a2
  Fp_add_no_red(t1, b->a1, b->a2);  // t1 = b1 + b2
  Fp_mul_no_red(T3, t0, t1);		    // T3 = (a1 + a2)*(b1 + b2) = a1*b1 + a1*b2 + a2*b1 + a2*b2
  Fpl_add_no_red(T4, T1, T2);		    // T4 = a1*b1 + a2*b2
	Fpl_sub_no_red(T3, T3, T4);		    // T3 = T3 - T4 = a1*b2 + a2*b1
  Fpl_add_no_red(T3, T3, T3);       // T3 = 2*(a1*b2 + a2*b1)
	
  Fp_add_no_red(t0, a->a0, a->a2);	// t0 = a0 + a2
  Fp_add_no_red(t1, b->a0, b->a2);  // t1 = b0 + b2
  Fp_mul_no_red(T4, t0, t1);		    // T4 = (a0 + a2)*(b0 + b2) = a0*b0 + a0*b2 + a2*b0 + a2*b2
  Fpl_sub_no_red(T4, T4, T0);
  Fpl_sub_no_red(T4, T4, T2);       // T4 = a0*b2 + a2*b0
  Fpl_add_no_red(C->A2, T4, T1);    // C2 = a0*b2 + a2*b0 + a1*b1
  
  Fp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
  Fp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
  Fp_mul_no_red(T4, t0, t1);		    // T4 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1
  Fpl_sub_no_red(T4, T4, T0);
  Fpl_sub_no_red(T4, T4, T1);       // T4 = a0*b1 + a1*b0
  Fpl_add_no_red(T2, T2, T2);       // T2 = 2*a2*b2
  Fpl_sub_o2_ct(C->A1, T4, T2, Fp_config.mem); // C1 = a0*b1 + a1*b0 - 2*a2*b2

  Fpl_sub_o2_ct(C->A0, T0, T3, Fp_config.mem); // C0 = a0*b0 - 2*(a1*b2 + a2*b1)
}

/* Multiply two Fp3 elements without reduction resulting in a long Fp3l element.
 * This function assumes that Fp3 = Fp(u) where u^3 = 3. 
 * This function needs 12 n-word Fp elements for temporary memory in mem.
 */
void u3minus3_Fp3_mul_no_red(Fp3l_t C, Fp3_t a, Fp3_t b, void *mem) {
	Fp_t t0, t1;
  Fpl_t T0, T1, T2, T3, T4;

  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + Fp_config.m->n;
    
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  T4->limbs = (uint64_t *) mem + 10*Fp_config.m->n;

  Fp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	Fp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1
  Fp_mul_no_red(T2, a->a2, b->a2);	// T2 = a2*b2

	Fp_add_no_red(t0, a->a1, a->a2);	// t0 = a1 + a2
  Fp_add_no_red(t1, b->a1, b->a2);  // t1 = b1 + b2
  Fp_mul_no_red(T3, t0, t1);		    // T3 = (a1 + a2)*(b1 + b2) = a1*b1 + a1*b2 + a2*b1 + a2*b2
  Fpl_add_no_red(T4, T1, T2);		    // T4 = a1*b1 + a2*b2
	Fpl_sub_no_red(T3, T3, T4);		    // T3 = T3 - T4 = a1*b2 + a2*b1
  Fpl_add_no_red(T4, T3, T3);       // T4 = 2*(a1*b2 + a2*b1)
  Fpl_add_no_red(T3, T4, T3);       // T3 = 3*(a1*b2 + a2*b1)
	
  Fp_add_no_red(t0, a->a0, a->a2);	// t0 = a0 + a2
  Fp_add_no_red(t1, b->a0, b->a2);  // t1 = b0 + b2
  Fp_mul_no_red(T4, t0, t1);		    // T4 = (a0 + a2)*(b0 + b2) = a0*b0 + a0*b2 + a2*b0 + a2*b2
  Fpl_sub_no_red(T4, T4, T0);
  Fpl_sub_no_red(T4, T4, T2);       // T4 = a0*b2 + a2*b0
  Fpl_add_no_red(C->A2, T4, T1);    // C2 = a0*b2 + a2*b0 + a1*b1
  
  Fp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
  Fp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
  Fp_mul_no_red(T4, t0, t1);		    // T4 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1
  Fpl_sub_no_red(T4, T4, T0);
  Fpl_sub_no_red(T4, T4, T1);       // T4 = a0*b1 + a1*b0
  Fpl_add_no_red(C->A1, T2, T2);    // C1 = 2*a2*b2
  Fpl_add_no_red(T2, C->A1, T2);    // C1 = 3*a2*b2
  Fpl_add_no_red(C->A1, T4, T2);    // C1 = a0*b1 + a1*b0 + 3*a2*b2

  Fpl_add_no_red(C->A0, T0, T3); // C0 = a0*b0 + 3*(a1*b2 + a2*b1)
}

/* Squaring an Fp3 element without reduction resulting in a long Fp3l element. 
 * This function assumes that Fp3 = Fp(u) with u^3 = -2. 
 */
void u3plus2_Fp3_squ_no_red_o1(Fp3l_t C, Fp3_t a, int h, void *mem) {
	Fp_t t0;
  Fpl_t T0, T1, T2, T3;

  t0->limbs = (uint64_t *) mem;
  T0->limbs = (uint64_t *) mem + Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	
	Fp_mul_no_red(T0, a->a0, a->a1);    // T0 = a0*a1  
	Fpl_add_no_red(T0, T0, T0);         // T0 = 2*a0*a1
	Fp_mul_no_red(T1, a->a2, a->a2);    // T1 = a2^2
	Fp_mul_no_red(T2, a->a0, a->a0);    // T2 = a0^2
	Fp_add_no_red(t0, a->a0, a->a1);    // t0 = a0 + a1
	Fp_add_no_red(t0, t0, a->a2);       // t0 = a0 + a1 + a2
	Fp_mul_no_red(T3, a->a1, a->a2);    // T3 = a1*a2
	Fpl_add_no_red(T3, T3, T3);         // T3 = 2*a1*a2

	Fpl_add_no_red(C->A2, T1, T1);      // C2 = 2*a2^2
  Fpl_sub_o1_ct(C->A1, T0, C->A2, h, Fp_config.mem); // C1 = 2*a0*a1 - 2*a2^2

	Fp_mul_no_red(C->A2, t0, t0);       // C2 = (a0 + a1 + a2)^2
	Fpl_sub_no_red(C->A2, C->A2, T0);   // C2 = a0^2 + a1^2 + a2^2 + 2*a1*a2 + 2*a0*a2
	Fpl_sub_no_red(C->A2, C->A2, T3);   // C2 = a0^2 + a1^2 + a2^2 + 2*a0*a2
	Fpl_sub_no_red(C->A2, C->A2, T2);   // C2 = a1^2 + a2^2 + 2*a0*a2
  Fpl_sub_no_red(C->A2, C->A2, T1);   // C2 = a1^2 + 2*a0*a2

	Fpl_add_no_red(T3, T3, T3);         // T3 = 4*a1*a2
  Fpl_sub_o1_ct(C->A0, T2, T3, h, Fp_config.mem); // C0 = a0^2 - 4*a1*a2
}


void u3plus2_Fp3_squ_no_red_o2(Fp3l_t C, Fp3_t a, void *mem) {
	Fp_t t0;
  Fpl_t T0, T1, T2, T3;

  t0->limbs = (uint64_t *) mem;
  T0->limbs = (uint64_t *) mem + Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	
	Fp_mul_no_red(T0, a->a0, a->a1);    // T0 = a0*a1  
	Fpl_add_no_red(T0, T0, T0);         // T0 = 2*a0*a1
	Fp_mul_no_red(T1, a->a2, a->a2);    // T1 = a2^2
	Fp_mul_no_red(T2, a->a0, a->a0);    // T2 = a0^2
	Fp_add_no_red(t0, a->a0, a->a1);    // t0 = a0 + a1
	Fp_add_no_red(t0, t0, a->a2);       // t0 = a0 + a1 + a2
	Fp_mul_no_red(T3, a->a1, a->a2);    // T3 = a1*a2
	Fpl_add_no_red(T3, T3, T3);         // T3 = 2*a1*a2

	Fpl_add_no_red(C->A2, T1, T1);      // C2 = 2*a2^2
  Fpl_sub_o2_ct(C->A1, T0, C->A2, Fp_config.mem); // C1 = 2*a0*a1 - 2*a2^2

	Fp_mul_no_red(C->A2, t0, t0);       // C2 = (a0 + a1 + a2)^2
	Fpl_sub_no_red(C->A2, C->A2, T0);   // C2 = a0^2 + a1^2 + a2^2 + 2*a1*a2 + 2*a0*a2
	Fpl_sub_no_red(C->A2, C->A2, T3);   // C2 = a0^2 + a1^2 + a2^2 + 2*a0*a2
	Fpl_sub_no_red(C->A2, C->A2, T2);   // C2 = a1^2 + a2^2 + 2*a0*a2
  Fpl_sub_no_red(C->A2, C->A2, T1);   // C2 = a1^2 + 2*a0*a2

	Fpl_add_no_red(T3, T3, T3);          // T3 = 4*a1*a2
  Fpl_sub_o2_ct(C->A0, T2, T3, Fp_config.mem); // C0 = a0^2 - 4*a1*a2
}

/* Squaring an Fp3 element without reduction resulting in a long Fp3l element. 
 * This function assumes that Fp3 = Fp(u) with u^3 = 3. 
 */
void u3minus3_Fp3_squ_no_red(Fp3l_t C, Fp3_t a, void *mem) {
	Fp_t t0;
  Fpl_t T0, T1, T2, T3;

  t0->limbs = (uint64_t *) mem;
  T0->limbs = (uint64_t *) mem + Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 3*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 5*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 7*Fp_config.m->n;
	
	Fp_mul_no_red(T0, a->a0, a->a1);    // T0 = a0*a1  
	Fpl_add_no_red(T0, T0, T0);         // T0 = 2*a0*a1
	Fp_mul_no_red(T1, a->a2, a->a2);    // T1 = a2^2
	Fp_mul_no_red(T2, a->a0, a->a0);    // T2 = a0^2
	Fp_add_no_red(t0, a->a0, a->a1);    // t0 = a0 + a1
	Fp_add_no_red(t0, t0, a->a2);       // t0 = a0 + a1 + a2
	Fp_mul_no_red(T3, a->a1, a->a2);    // T3 = a1*a2
	Fpl_add_no_red(T3, T3, T3);         // T3 = 2*a1*a2

	Fpl_add_no_red(C->A2, T1, T1);      // C2 = 2*a2^2
  Fpl_add_no_red(C->A2, C->A2, T1);   // C2 = 3*a2^2
  Fpl_add_no_red(C->A1, T0, C->A2);   // C1 = 2*a0*a1 + 3*a2^2

	Fp_mul_no_red(C->A2, t0, t0);       // C2 = (a0 + a1 + a2)^2
	Fpl_sub_no_red(C->A2, C->A2, T0);   // C2 = a0^2 + a1^2 + a2^2 + 2*a1*a2 + 2*a0*a2
	Fpl_sub_no_red(C->A2, C->A2, T3);   // C2 = a0^2 + a1^2 + a2^2 + 2*a0*a2
	Fpl_sub_no_red(C->A2, C->A2, T2);   // C2 = a1^2 + a2^2 + 2*a0*a2
  Fpl_sub_no_red(C->A2, C->A2, T1);   // C2 = a1^2 + 2*a0*a2

	Fpl_add_no_red(C->A0, T3, T3);      // C0 = 4*a1*a2
  Fpl_add_no_red(C->A0, C->A0, T3);   // C0 = 6*a1*a2
  Fpl_add_no_red(C->A0, T2, C->A0);   // C0 = a0^2 + 6*a1*a2
}

/* Multiply two Fp3 elements without reduction.
 * This function assumes that Fp3 is constructed via v^3 - xi.
 * The third coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b0 + b1*v, i.e. b2 = 0.
 */
void v3minusxi_Fp3_mul_sparse01_no_red (Fp3l_t c, Fp3_t a, Fp_t b0, Fp_t b1, void *mem) {
  Fp_t t0, t1;
  Fpl_t T0, T1, T2, T3;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  
  Fp_mul_no_red (T0, a->a0, b0);
  Fp_mul_no_red (T1, a->a1, b1);

  Fp_mul_no_red (T3, a->a2, b1);
  Fpl_add_no_red (T2, T3, T3);  
  Fpl_add_no_red (T2, T2, T3);
  Fpl_add_no_red (c->A0, T2, T0);

  Fp_add_no_red (t0, a->a0, a->a1);
  Fp_add_no_red (t1, b0, b1);
  Fp_mul_no_red (T3, t0, t1);
  Fpl_add_no_red (T2, T0, T1);
  Fpl_sub_no_red (c->A1, T3, T2);
  
  Fp_mul_no_red (T3, a->a2, b0);
  Fpl_add_no_red (c->A2, T3, T1);
}

/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * The first coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b1*v + b2*v^2 = v*(b1 + b2*v).
 */
void v3minusxi_Fp3_mul_sparse12_no_red (Fp3l_t c, Fp3_t a, Fp_t b1, Fp_t b2, void *mem) {
  Fp_t t0, t1;
  Fpl_t T0, T1, T2, T3;
  
  t0->limbs = (uint64_t *) mem;
  t1->limbs = (uint64_t *) mem + 1*Fp_config.m->n;
  
  T0->limbs = (uint64_t *) mem + 2*Fp_config.m->n;
  T1->limbs = (uint64_t *) mem + 4*Fp_config.m->n;
  T2->limbs = (uint64_t *) mem + 6*Fp_config.m->n;
  T3->limbs = (uint64_t *) mem + 8*Fp_config.m->n;
  
  Fp_mul_no_red (T0, a->a2, b2);
  Fp_mul_no_red (T1, a->a1, b1);

  Fp_mul_no_red (T3, a->a0, b1);
  Fpl_add_no_red (T2, T0, T0);  
  Fpl_add_no_red (T2, T2, T0);
  Fpl_add_no_red (c->A1, T3, T2);

  Fp_add_no_red (t0, a->a1, a->a2);
  Fp_add_no_red (t1, b1, b2);
  Fp_mul_no_red (T3, t0, t1);
  Fpl_add_no_red (T2, T0, T1);
  Fpl_sub_no_red (T3, T3, T2);
  Fpl_add_no_red (c->A0, T3, T3);
  Fpl_add_no_red (c->A0, c->A0, T3);

  Fp_mul_no_red (T3, a->a0, b2);
  Fpl_add_no_red (c->A2, T3, T1);
}



int Fp3_initialize_config (void) {
  uint64_t *vpp;

  Fp3_config.mem = NULL;
  Fp3_config.mem = (void *) malloc (12 * Fp_config.m->n * sizeof (uint64_t));
  if (Fp3_config.mem == NULL) return ERR_OUT_OF_MEMORY;
  
  vpp = (uint64_t *) Fp3_config.mem;
  if (Fp_init (Fp3_config.uppow) < 0) printf ("Fp3_config memory error.\n");
  if (Fp_init (Fp3_config.uppow2) < 0) printf ("Fp3_config memory error.\n");
  if (Fp3l_init (Fp3_config.T0) < 0) printf ("Fp3_config memory error.\n");
  
 if (PAIR_CURVE == CP6) {
    uint64_t upp[8] = { 0x5257C718F4B582E1,0x0C1197483FBEF0DD,
                        0x74EDC8603F8B0275,0xBB686D34F710DCBC,
                        0xAFAADB6265FA73C3,0x262A5F326F028EA5,
                        0x33B162CE50A3D98C,0x05100016C8007DD9 };
    
    Fp_set (Fp3_config.uppow, upp);
    Fp_mul (Fp3_config.uppow2, Fp3_config.uppow, Fp3_config.uppow);
  
    Fp3_config.Fp3_mul = u3minus3_Fp3_mul;
    Fp3_config.Fp3_squ = u3minus3_Fp3_squ;
    Fp3_config.Fp3_inv = u3minus3_Fp3_inv;
    Fp3_config.Fp3_mulu = u3minus3_Fp3_mulu;
    Fp3_config.Fp3l_mulu = u3minus3_Fp3l_mulu;
    Fp3_config.Fp3_ppow = p1mod3_Fp3_ppow;
    Fp3_config.Fp3_mul_no_red = u3minus3_Fp3_mul_no_red;
    Fp3_config.Fp3_squ_no_red = u3minus3_Fp3_squ_no_red;
    Fp3_config.Fp3_mul_sparse01_no_red = v3minusxi_Fp3_mul_sparse01_no_red;
    Fp3_config.Fp3_mul_sparse12_no_red = v3minusxi_Fp3_mul_sparse12_no_red;
  } else if (PAIR_CURVE == CP6b) {
    uint64_t upp[8] = { 0x8C000000000004C3,0x6BF39000000009D9,
                        0x9D35F6980000093D,0xBCFC9A74F8000520,
                        0xA416A709DA2C01D6,0x1833E371C780246F,
                        0xBB99497657F45361,0x1029BFD6DCBA8CD0 };
    
    Fp_set (Fp3_config.uppow, upp);
    Fp_mul (Fp3_config.uppow2, Fp3_config.uppow, Fp3_config.uppow);
  
    Fp3_config.Fp3_mul = u3minus3_Fp3_mul;
    Fp3_config.Fp3_squ = u3minus3_Fp3_squ;
    Fp3_config.Fp3_inv = u3minus3_Fp3_inv;
    Fp3_config.Fp3_mulu = u3minus3_Fp3_mulu;
    Fp3_config.Fp3l_mulu = u3minus3_Fp3l_mulu;
    Fp3_config.Fp3_ppow = p1mod3_Fp3_ppow;
    Fp3_config.Fp3_mul_no_red = u3minus3_Fp3_mul_no_red;
    Fp3_config.Fp3_squ_no_red = u3minus3_Fp3_squ_no_red;
    Fp3_config.Fp3_mul_sparse01_no_red = v3minusxi_Fp3_mul_sparse01_no_red;
    Fp3_config.Fp3_mul_sparse12_no_red = v3minusxi_Fp3_mul_sparse12_no_red;
  }
  else if (PAIR_CURVE == CP6tiny) {
	  uint64_t upp[2] = { 0xDB70B1732C9457D7, 0x59032611135 };

	  Fp_set(Fp3_config.uppow, upp);
	  Fp_mul(Fp3_config.uppow2, Fp3_config.uppow, Fp3_config.uppow);

	  Fp3_config.Fp3_mul = u3minus3_Fp3_mul;
	  Fp3_config.Fp3_squ = u3minus3_Fp3_squ;
	  Fp3_config.Fp3_inv = u3minus3_Fp3_inv;
	  Fp3_config.Fp3_mulu = u3minus3_Fp3_mulu;
	  Fp3_config.Fp3l_mulu = u3minus3_Fp3l_mulu;
	  Fp3_config.Fp3_ppow = p1mod3_Fp3_ppow;
	  Fp3_config.Fp3_mul_no_red = u3minus3_Fp3_mul_no_red;
	  Fp3_config.Fp3_squ_no_red = u3minus3_Fp3_squ_no_red;
	  Fp3_config.Fp3_mul_sparse01_no_red = v3minusxi_Fp3_mul_sparse01_no_red;
	  Fp3_config.Fp3_mul_sparse12_no_red = v3minusxi_Fp3_mul_sparse12_no_red;
  }
  else if (PAIR_CURVE == CP3) {
	  uint64_t upp[16] = { 0x6500000004D749F9, 0x7432800013F55DF7,
		  0xFC0F14902752C06B, 0x518B2D58741C2675,
		  0x4657568B85098E0F, 0xEFF11F2E1B0F0BEF,
		  0x433121B22515C924, 0x108A2330A035F429,
		  0x1BAFC6FCC865F0C0, 0xCFF1F71321CDB3C6,
		  0x81951022D2BADD18, 0x448A0FAF16AA71EE,
		  0xB995C1D855434EF8, 0x254C47AA246EEF4E,
		  0xF3182284379A111C, 0x3F03E5315731FFF2 };

	  Fp_set(Fp3_config.uppow, upp);
	  Fp_mul(Fp3_config.uppow2, Fp3_config.uppow, Fp3_config.uppow);

	  Fp3_config.Fp3_mul = u3minus3_Fp3_mul;
	  Fp3_config.Fp3_squ = u3minus3_Fp3_squ;
	  Fp3_config.Fp3_inv = u3minus3_Fp3_inv;
	  Fp3_config.Fp3_mulu = u3minus3_Fp3_mulu;
	  Fp3_config.Fp3l_mulu = u3minus3_Fp3l_mulu;
	  Fp3_config.Fp3_ppow = p1mod3_Fp3_ppow;
	  Fp3_config.Fp3_mul_no_red = u3minus3_Fp3_mul_no_red;
	  Fp3_config.Fp3_squ_no_red = u3minus3_Fp3_squ_no_red;
	  Fp3_config.Fp3_mul_sparse01_no_red = v3minusxi_Fp3_mul_sparse01_no_red;
	  Fp3_config.Fp3_mul_sparse12_no_red = v3minusxi_Fp3_mul_sparse12_no_red;
  }
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }
  
  return ERR_SUCCESS;
}

void Fp3_free_config (void) {
  free (Fp3_config.mem);
  Fp_free (Fp3_config.uppow);
  Fp_free (Fp3_config.uppow2);
  Fp3l_free(Fp3_config.T0);
}


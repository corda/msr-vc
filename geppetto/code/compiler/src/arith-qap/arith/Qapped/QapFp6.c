#include "QapFp6.h"
#include "Qappairing.h"
#include "QapTempMem.h"

Fp6_config_t Fp6_config;

/* Initialize an element in Fp6 by initializing three Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp6_init(QapFp6_t dst) {
	int err = 0;
	err = QapFp2_init(dst->a0);
	if (err < 0) return err;
	err = QapFp2_init(dst->a1);
	if (err < 0) return err;
	return QapFp2_init(dst->a2);
}

/* Free the memory allocated in src. */
void QapFp6_free (QapFp6_t src) {
	QapFp2_free(src->a0);
	QapFp2_free(src->a1);
	QapFp2_free(src->a2);
}

/* Get a random element in Fp6. */
void QapFp6_rand (QapFp6_t c) {
  QapFp2_rand(c->a0);
  QapFp2_rand(c->a1);
  QapFp2_rand(c->a2);
}

/* Set c to 0. */
void QapFp6_set_zero (QapFp6_t c) {
  QapFp2_set_zero (c->a0);
  QapFp2_set_zero (c->a1);
  QapFp2_set_zero (c->a2);
}

/* Set c to 1. */
void QapFp6_set_one (QapFp6_t c) {
  QapFp2_set_one (c->a0);
  QapFp2_set_zero (c->a1);
  QapFp2_set_zero (c->a2);
}

/* Copy an Fp6 element a to an Fp6 element c. */
void QapFp6_copy(QapFp6_t c, QapFp6_t a) {
	QapFp2_copy(c->a0, a->a0);
	QapFp2_copy(c->a1, a->a1);
	QapFp2_copy(c->a2, a->a2);
}

void QapFp6_select(QapFp6_t c, int bit, QapFp6_t a, QapFp6_t b) {
	QapFp2_select(c->a0, bit, a->a0, b->a0);
	QapFp2_select(c->a1, bit, a->a1, b->a1);
	QapFp2_select(c->a2, bit, a->a2, b->a2);
}

/* Print the (regular non-Montgomery form) value of a */
void QapFp6_print (QapFp6_t a) {
  QapFp2_print(a->a0);
  printf("\n");
  QapFp2_print(a->a1);
  printf("\n");
  QapFp2_print(a->a2);
  printf("\n");
}

/* Compare two Fp6 elements a and b for equality. */
int QapFp6_cmpeq(QapFp6_t a, QapFp6_t b) {
#if MQAP
	return (QapFp2_cmpeq(a->a0, b->a0) * QapFp2_cmpeq(a->a1, b->a1) * QapFp2_cmpeq(a->a2, b->a2));
#else
	return (QapFp2_cmpeq(a->a0, b->a0) && QapFp2_cmpeq(a->a1, b->a1) && QapFp2_cmpeq(a->a2, b->a2));
#endif
}

/* Negate an Fp6 element. */
void QapFp6_neg(QapFp6_t c, QapFp6_t a) {
	QapFp2_neg(c->a0, a->a0);
	QapFp2_neg(c->a1, a->a1);
	QapFp2_neg(c->a2, a->a2);
}

/* Add two Fp6 elements coefficient wise. */
void QapFp6_add (QapFp6_t c, QapFp6_t a, QapFp6_t b) {
	QapFp2_add(c->a0, a->a0, b->a0);
	QapFp2_add(c->a1, a->a1, b->a1);
	QapFp2_add(c->a2, a->a2, b->a2);
}

/* Subtract an Fp6 element b from an Fp6 element a coefficient wise. */
void QapFp6_sub (QapFp6_t c, QapFp6_t a, QapFp6_t b) {
	QapFp2_sub(c->a0, a->a0, b->a0);
	QapFp2_sub(c->a1, a->a1, b->a1);
	QapFp2_sub(c->a2, a->a2, b->a2);
}

/* Multiply an Fp6 element by an Fp2 element. */
void QapFp6_mulFp2(QapFp6_t c, QapFp6_t a, QapFp2_t b) {
	QapFp2_mul(c->a0, a->a0, b);
	QapFp2_mul(c->a1, a->a1, b);
	QapFp2_mul(c->a2, a->a2, b);
}

/* Multiply an Fp6 element by an Fp element. */
void QapFp6_mulFp(QapFp6_t c, QapFp6_t a, QapFp_t b) {
	QapFp2_mulFp(c->a0, a->a0, b);
	QapFp2_mulFp(c->a1, a->a1, b);
	QapFp2_mulFp(c->a2, a->a2, b);
}

/* Multiply two Fp6 elements. */
void QapFp6_mul (QapFp6_t c, QapFp6_t a, QapFp6_t b) {
  Qapv3minusxi_Fp6_mul(c, a, b, TEMP_PASS_MEM6());
}

/* Square an Fp6 element. */
void QapFp6_squ (QapFp6_t c, QapFp6_t a) {
  Qapv3minusxi_Fp6_squ(c, a, TEMP_PASS_MEM6());
}

/* Invert an Fp6 element. */
void QapFp6_inv (QapFp6_t c, QapFp6_t a) {
  Qapv3minusxi_Fp6_inv(c, a, TEMP_PASS_MEM6());
}

/* Multiply an fp6 element by the special element v. */
void QapFp6_mulv(QapFp6_t c, QapFp6_t a) {
  Qapv3minusxi_Fp6_mulv(c, a, TEMP_PASS_MEM6());
}

/* Compute the p-power Frobenius of an Fp6 element. */
void QapFp6_ppow(QapFp6_t c, QapFp6_t a) {
  Qapp1mod3_v3minusxi_Fp6_ppow(c, a);
}

/* Compute the p^2-power Frobenius of an Fp6 element. */
void QapFp6_p2pow(QapFp6_t c, QapFp6_t a) {
  Qapp1mod3_v3minusxi_Fp6_p2pow(c, a);
}

/* Compute the p^3-power Frobenius of an Fp6 element. */
void QapFp6_p3pow(QapFp6_t c, QapFp6_t a) {
  Qapp1mod3_v3minusxi_Fp6_p3pow(c, a);
}


/* Multiply two Fp6 elements. This uses that Fp6 is constructed over Fp2 via v^3 - xi. 
 * This is Algorithm 13 in Beuchat et al., Pairing 2010.
 */
void Qapv3minusxi_Fp6_mul(QapFp6_t c, QapFp6_t a, QapFp6_t b, TEMP_FORMAL) {
	QapFp2_t t0, t1, t2, t3, t4, t5;
	TEMP_ALLOC_FP2(t0, 0);
	TEMP_ALLOC_FP2(t1, 2);
	TEMP_ALLOC_FP2(t2, 4);
	TEMP_ALLOC_FP2(t3, 6);
	TEMP_ALLOC_FP2(t4, 8);
	TEMP_ALLOC_FP2(t5, 10);

	QapFp2_mul(t0, a->a0, b->a0);    // t0 = a0*b0
	QapFp2_mul(t1, a->a1, b->a1);    // t1 = a1*b1
	QapFp2_mul(t2, a->a2, b->a2);    // t2 = a2*b2

	QapFp2_add(t3, a->a1, a->a2);    // t3 = a1 + a2
	QapFp2_add(t4, b->a1, b->a2);    // t4 = b1 + b2
	QapFp2_mul(t3, t3, t4);          // t3 = (a1 + a2)*(b1 + b2) = a1*b1 + a2*b1 + a1*b2 + a2*b2
	QapFp2_sub(t3, t3, t1);          // t3 = a2*b1 + a1*b2 + a2*b2
	QapFp2_sub(t3, t3, t2);          // t3 = a2*b1 + a1*b2 
	QapFp2_mulxi(t3, t3);					  // t3 = (a2*b1 + a1*b2)*xi	
	// t3 now almost has the result for c0, t4 is free.		

	QapFp2_add(t4, a->a0, a->a2);    // t4 = a0 + a2
	QapFp2_add(t5, b->a0, b->a2);    // t5 = b0 + b2
	QapFp2_mul(t4, t4, t5);          // t4 = (a0 + a2)*(b0 + b2) = a0*b0 + a2*b0 + a0*b2 + a2*b2
	QapFp2_sub(t4, t4, t0);          // t4 = a2*b0 + a0*b2 + a2*b2
	QapFp2_sub(t4, t4, t2);          // t4 = a2*b0 + a0*b2
	QapFp2_add(c->a2, t4, t1);				// c2 = a2*b0 + a0*b2 + a1*b1 	
	// c2 is done, t4 and t5 are free.

	QapFp2_add(t4, a->a0, a->a1);    // t4 = a0 + a1
	QapFp2_add(t5, b->a0, b->a1);    // t4 = b0 + b1
	QapFp2_mul(t4, t4, t5);          // t4 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1
	QapFp2_sub(t4, t4, t0);          // t4 = a1*b0 + a0*b1 + a1*b1
	QapFp2_sub(t4, t4, t1);          // t4 = a1*b0 + a0*b1
	QapFp2_mulxi(t2, t2);            // t2 = a2*b2 * xi
	QapFp2_add(c->a1, t4, t2);				// c1 = a1*b0 + a0*b1 + a2*b2*xi	
	// c1 is done.

	QapFp2_add(c->a0, t3, t0);					// c0 is done.

	TEMP_FREE_FP2(t0);
	TEMP_FREE_FP2(t1);
	TEMP_FREE_FP2(t2);
	TEMP_FREE_FP2(t3);
	TEMP_FREE_FP2(t4);
	TEMP_FREE_FP2(t5);
}

/* Square an Fp6 element. This uses that Fp6 is constructed over Fp2 via v^3 - xi. 
 * This is Algorithm 16 in Beuchat et al., Pairing 2010.
 */
void Qapv3minusxi_Fp6_squ(QapFp6_t c, QapFp6_t a, TEMP_FORMAL) {
	QapFp2_t t0, t1, t2, t3, t4;
	TEMP_ALLOC_FP2(t0, 0);
	TEMP_ALLOC_FP2(t1, 2);
	TEMP_ALLOC_FP2(t2, 4);
	TEMP_ALLOC_FP2(t3, 6);
	TEMP_ALLOC_FP2(t4, 8);

	QapFp2_mul(t0, a->a0, a->a1);    // t0 = a0*a1  
	QapFp2_add(t0, t0, t0);          // t0 = 2*a0*a1
	QapFp2_squ(t1, a->a2);           // t1 = a2^2
	QapFp2_squ(t2, a->a0);           // t2 = a0^2
	QapFp2_sub(t3, a->a0, a->a1);    // t3 = a0 - a1
	QapFp2_add(t3, t3, a->a2);       // t3 = a0 - a1 + a2
	QapFp2_mul(t4, a->a1, a->a2);    // t4 = a1*a2
	QapFp2_add(t4, t4, t4);          // t4 = 2*a1*a2

	QapFp2_sub(c->a2, t0, t1);       // c2 = 2*a0*a1 - a2^2
	QapFp2_mulxi(t1, t1);            // t1 = xi*a2^2
	QapFp2_add(c->a1, t0, t1);       // c1 = 2*a0*a1 + xi*a2^2

	QapFp2_squ(t0, t3);              // t0 = (a0 - a1 + a2)^2
	QapFp2_add(c->a2, c->a2, t0);    // c2 = a0^2 + a1^2 - 2*a1*a2 + 2*a0*a2
	QapFp2_add(c->a2, c->a2, t4);    // c2 = a0^2 + a1^2 + 2*a0*a2
	QapFp2_sub(c->a2, c->a2, t2);    // c2 = a1^2 + 2*a0*a2

	QapFp2_mulxi(t4, t4);            // t4 = 2*a1*a2*xi
	QapFp2_add(c->a0, t4, t2);       // c0 = a0^2 + 2*a1*a2*xi

	TEMP_FREE_FP2(t0);
	TEMP_FREE_FP2(t1);
	TEMP_FREE_FP2(t2);
	TEMP_FREE_FP2(t3);
	TEMP_FREE_FP2(t4);
}

/* Multiplication by the special element v in Fp6, with v^3 = xi in Fp2.
 * c = c0 + c1*v + c2*v^2 = (a0 + a1*v + a2*v^2)*v = a2*xi + a0*v + a1*v^2.
 */
void Qapv3minusxi_Fp6_mulv(QapFp6_t c, QapFp6_t a, TEMP_FORMAL) {
  QapFp2_t t0;
  TEMP_ALLOC_FP2(t0, 0);
 	
  QapFp2_mulxi(t0, a->a2);
  QapFp2_copy(c->a2, a->a1);
	QapFp2_copy(c->a1, a->a0);
  QapFp2_copy(c->a0, t0);
}

/* Compute the p-power Frobenius of an Fp6 element. This uses a constant VPPOW, 
 * such that v^p = VPPOW*v, i.e. VPPOW = v^(p-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and p = 1 (mod 3), we have VPPOW = xi^(p-1)/3 in Fp2. 
 * VPPOW2 = VPPOW^2.
 */
void Qapp1mod3_v3minusxi_Fp6_ppow(QapFp6_t c, QapFp6_t a) {
  QapFp2_ppow(c->a0, a->a0);
	QapFp2_ppow(c->a1, a->a1);
	QapFp2_ppow(c->a2, a->a2);
	QapFp2_muliFp(c->a1, c->a1, Fp6_config.vppow);
	QapFp2_mulFp(c->a2, c->a2, Fp6_config.vppow2);
}

/* Compute the p^2-power Frobenius of an Fp6 element. As above, we use a constant
 * VP2POW such that v^(p^2) = VP2POW*v, i.e. VP2POW = v^(p^2-1).  In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and p = 1 (mod 3), we have VP2POW = xi^(p^2-1)/3 is a third root
 * of unity in Fp. VP2POW2 = VP2POW^2.  
 */
void Qapp1mod3_v3minusxi_Fp6_p2pow(QapFp6_t c, QapFp6_t a) {
	QapFp2_copy(c->a0, a->a0);
	QapFp2_mulFp(c->a1, a->a1, Fp6_config.vp2pow);
	QapFp2_mulFp(c->a2, a->a2, Fp6_config.vp2pow2);
}

/* Compute the p^3-power Frobenius of an Fp6 element. As above, we use a constant
 * VP3POW such that v^(p^3) = VP3POW*v, i.e. VP3POW = v^(p^3-1).  In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and p = 1 (mod 3), we have VP3POW = xi^(p^3-1)/3 in Fp2. 
 * VP3POW2 = VP3POW^2. In the BN case, VP3POW = i and VP3POW2 = -1. No need for a constant here.
 */
void Qapp1mod3_v3minusxi_Fp6_p3pow(QapFp6_t c, QapFp6_t a) {
	QapFp2_ppow(c->a0, a->a0);
	QapFp2_ppow(c->a1, a->a1);
	QapFp2_ppow(c->a2, a->a2);
	QapFp2_muli(c->a1, c->a1);
	QapFp2_neg(c->a2, c->a2);
}

/* Inversion of an Fp6 element. This uses that Fp6 is constructed over Fp2 via v^3 - xi. 
 * This is Algorithm 17 in Beuchat et al., Pairing 2010.
 */
void Qapv3minusxi_Fp6_inv(QapFp6_t c, QapFp6_t a, TEMP_FORMAL) {
	QapFp2_t t0, t1, t2, t3, t4, t5;
	TEMP_ALLOC_FP2(t0, 0);
	TEMP_ALLOC_FP2(t1, 2);
	TEMP_ALLOC_FP2(t2, 4);
	TEMP_ALLOC_FP2(t3, 6);
	TEMP_ALLOC_FP2(t4, 8);
	TEMP_ALLOC_FP2(t5, 10);

	QapFp2_squ(t0, a->a0);         // t0 = a0^2
	QapFp2_squ(t1, a->a1);         // t1 = a1^2
	QapFp2_squ(t2, a->a2);         // t2 = a2^2

	QapFp2_mul(t3, a->a0, a->a1);  // t3 = a0*a1
	QapFp2_mul(t4, a->a0, a->a2);  // t4 = a0*a2
	QapFp2_mul(t5, a->a1, a->a2);  // t5 = a1*a2

	QapFp2_mulxi(t2, t2);          // t2 = xi*a2^2
	QapFp2_mulxi(t5, t5);          // t5 = xi*a1*a2

	QapFp2_sub(t0, t0, t5);        // t0 = a0^2 - xi*a1*a2
	QapFp2_sub(t1, t1, t4);        // t1 = a1^2 - a0*a2
	QapFp2_sub(t2, t2, t3);        // t2 = xi*a2^2 - a0*a1

	QapFp2_mul(t3, a->a0, t0);     // t3 = a0^3 - xi*a0*a1*a2
	QapFp2_mul(t4, a->a2, t2);     // t4 = xi*a2^3 - a0*a1*a2
	QapFp2_mulxi(t4, t4);          // t4 = xi^2*a2^3 - xi*a0*a1*a2
	QapFp2_add(t3, t3, t4);        // t3 = a0^3 + xi^2*a2^3 - 2*xi*a0*a1*a2 
	QapFp2_mul(t4, a->a1, t1);     // t4 = a1^3 - a0*a1*a2
	QapFp2_mulxi(t4, t4);          // t4 = xi*a1^3 - xi*a0*a1*a2
	QapFp2_add(t3, t3, t4);        // t3 = a0^3 + xi*a1^3 + xi^2*a2^3 - 3*xi*a0*a1*a2

	QapFp2_inv(t3, t3);            // t3 = 1/N(a)

	QapFp2_mul(c->a0, t0, t3);
	QapFp2_mul(c->a1, t2, t3);
	QapFp2_mul(c->a2, t1, t3);

	TEMP_FREE_FP2(t0);
	TEMP_FREE_FP2(t1);
	TEMP_FREE_FP2(t2);
	TEMP_FREE_FP2(t3);
	TEMP_FREE_FP2(t4);
	TEMP_FREE_FP2(t5);
}

/************** Fp6 functions for lazy reduction ***************/
/* Functions assume that there is enough space such that no overflow occurs. */

/* Initialize a long Fp6 element. */
int QapFp6l_init(Fp6l_t dst) {
	int err = 0;
	if ((err = QapFp2l_init(dst->A0)) < 0) return err;
	if ((err = QapFp2l_init(dst->A1)) < 0) return err;
	return QapFp2l_init(dst->A2);
}

/* Free the memory allocated in src. */
void QapFp6l_free (Fp6l_t src) {
	QapFp2l_free(src->A0);
	QapFp2l_free(src->A1);
	QapFp2l_free(src->A2);
}

/* Copy an Fp6l element A to an Fp6l element C. */
void QapFp6l_copy(Fp6l_t C, Fp6l_t A) {
	QapFp2l_copy(C->A0, A->A0);
	QapFp2l_copy(C->A1, A->A1);
	QapFp2l_copy(C->A2, A->A2);
}

/* Reduction function calling the Fp2 reduction for each coefficient. */
void QapFp6l_red(QapFp6_t c, Fp6l_t A) {
	QapFp2l_red(c->a0, A->A0);
	QapFp2l_red(c->a1, A->A1);
	QapFp2l_red(c->a2, A->A2);
}

/* Subtraction of Fp6l elements, constant-time version of option 2 in Aranha et al. */
void QapFp6l_sub_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b) {
  QapFp2l_sub_o2(c->A0, a->A0, b->A0);
  QapFp2l_sub_o2(c->A1, a->A1, b->A1);
  QapFp2l_sub_o2(c->A2, a->A2, b->A2);
}

/* Addition of Fp6l elements, constant-time version of "option 2" addition in Aranha et al. */
void QapFp6l_add_o2 (Fp6l_t c, Fp6l_t a, Fp6l_t b) {
  QapFp2l_add_o2(c->A0, a->A0, b->A0);
  QapFp2l_add_o2(c->A1, a->A1, b->A1);
  QapFp2l_add_o2(c->A2, a->A2, b->A2);
}

/* Multiply an Fp6 element by an Fp element. */
void QapFp6_mulFp2_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b) {
	QapFp2_mul_no_red_o2 (c->A0, a->a0, b);
	QapFp2_mul_no_red_o2 (c->A1, a->a1, b);
  QapFp2_mul_no_red_o2 (c->A2, a->a2, b);
}

/* Multiply two Fp6 elements without reduction. */
void QapFp6_mul_no_red (Fp6l_t c, QapFp6_t a, QapFp6_t b) {
  Qapv3minusxi_Fp6_mul_no_red(c, a, b, TEMP_PASS_MEM6());
}

/* Multiply an Fp6 element and a sparse Fp6 element without reduction.
   This function assumes that Fp6 is constructed via v^3 - xi.
 */
void QapFp6_mul_sparse01_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b0, QapFp2_t b1){
  Qapv3minusxi_Fp6_mul_sparse01_no_red(c, a, b0, b1, TEMP_PASS_MEM6());
}
void QapFp6_mul_sparse12_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b1, QapFp2_t b2){
  Qapv3minusxi_Fp6_mul_sparse12_no_red(c, a, b1, b2, TEMP_PASS_MEM6());
}

/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * This is Algorithm 3 in Aranha et al.
 */
void Qapv3minusxi_Fp6_mul_no_red (Fp6l_t c, QapFp6_t a, QapFp6_t b, TEMP_FORMAL) {
  QapFp2_t t0, t1;
  Fp2l_t T0, T1, T2, T3, T4;

  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
  TEMP_ALLOC_FP2L(T0, 4);
  TEMP_ALLOC_FP2L(T1, 8);
  TEMP_ALLOC_FP2L(T2, 12);
  TEMP_ALLOC_FP2L(T3, 16);
  TEMP_ALLOC_FP2L(T4, 20);

  QapFp2_mul_no_red_o1 (T0, a->a0, b->a0, 2);
  QapFp2_mul_no_red_o1 (T1, a->a1, b->a1, 2);
  QapFp2_mul_no_red_o1 (T2, a->a2, b->a2, 2);

  QapFp2_add_no_red (t0, a->a1, a->a2);
  QapFp2_add_no_red (t1, b->a1, b->a2);
  QapFp2_mul_no_red_o2 (T3, t0, t1);
  QapFp2l_add_no_red (T4, T1, T2);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T4->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T4->A1);
  QapFpl_sub_o2_ct (T4->A0, T3->A0, T3->A1, TEMP_PASS_MEM());
  QapFpl_add_o2_ct (T4->A1, T3->A0, T3->A1, TEMP_PASS_MEM());  
  QapFp2l_add_o2 (c->A0, T4, T0);

  QapFp2_add_no_red (t0, a->a0, a->a1);
  QapFp2_add_no_red (t1, b->a0, b->a1);
  QapFp2_mul_no_red_o2 (T3, t0, t1);
  QapFp2l_add_no_red (T4, T0, T1);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T4->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T4->A1);
  QapFpl_sub_o1_ct (T4->A0, T2->A0, T2->A1, 1, TEMP_PASS_MEM());
  QapFpl_add_no_red (T4->A1, T2->A0, T2->A1);
  QapFp2l_add_o2 (c->A1, T3, T4);

  QapFp2_add_no_red (t0, a->a0, a->a2);
  QapFp2_add_no_red (t1, b->a0, b->a2);
  QapFp2_mul_no_red_o2 (T3, t0, t1);
  QapFp2l_add_no_red (T4, T0, T2);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T4->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T4->A1);
  QapFpl_add_o2_ct (c->A2->A0, T3->A0, T1->A0, TEMP_PASS_MEM()); 
  QapFpl_add_no_red (c->A2->A1, T3->A1, T1->A1);

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
  TEMP_FREE_FP2L(T0);
  TEMP_FREE_FP2L(T1);
  TEMP_FREE_FP2L(T2);
  TEMP_FREE_FP2L(T3);
  TEMP_FREE_FP2L(T4);
}


/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * The third coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b0 + b1*v, i.e. b2 = 0.
 */
void Qapv3minusxi_Fp6_mul_sparse01_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b0, QapFp2_t b1, TEMP_FORMAL) {
  QapFp2_t t0, t1;
  Fp2l_t T0, T1, T2, T3;
  
  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
  TEMP_ALLOC_FP2L(T0, 4);
  TEMP_ALLOC_FP2L(T1, 8);
  TEMP_ALLOC_FP2L(T2, 12);
  TEMP_ALLOC_FP2L(T3, 16);
  
  QapFp2_mul_no_red_o1 (T0, a->a0, b0, 2);
  QapFp2_mul_no_red_o1 (T1, a->a1, b1, 2);

  QapFp2_add_no_red (t0, a->a1, a->a2);
  QapFp2_mul_no_red_o2 (T3, t0, b1);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T1->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T1->A1);
  QapFpl_sub_o2_ct (T2->A0, T3->A0, T3->A1, TEMP_PASS_MEM());
  QapFpl_add_o2_ct (T2->A1, T3->A0, T3->A1, TEMP_PASS_MEM());  
  QapFp2l_add_o2 (c->A0, T2, T0);

  QapFp2_add_no_red (t0, a->a0, a->a1);
  QapFp2_add_no_red (t1, b0, b1);
  QapFp2_mul_no_red_o2 (T3, t0, t1);
  QapFp2l_add_no_red (T2, T0, T1);
  QapFpl_sub_o2_ct (c->A1->A0, T3->A0, T2->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (c->A1->A1, T3->A1, T2->A1);

  QapFp2_add_no_red (t0, a->a0, a->a2);
  QapFp2_mul_no_red_o2 (T3, t0, b0);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T0->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T0->A1);
  QapFpl_add_o2_ct (c->A2->A0, T3->A0, T1->A0, TEMP_PASS_MEM()); 
  QapFpl_add_no_red (c->A2->A1, T3->A1, T1->A1);

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
  TEMP_FREE_FP2L(T0);
  TEMP_FREE_FP2L(T1);
  TEMP_FREE_FP2L(T2);
  TEMP_FREE_FP2L(T3);
}

/* Multiply two Fp6 elements without reduction.
 * This function assumes that Fp6 is constructed via v^3 - xi.
 * The first coefficient of b is 0, so the element b is sparse in the following way: 
 * b = b1*v + b2*v^2 = v*(b1 + b2*v).
 */
void Qapv3minusxi_Fp6_mul_sparse12_no_red (Fp6l_t c, QapFp6_t a, QapFp2_t b1, QapFp2_t b2, TEMP_FORMAL) {
  QapFp2_t t0, t1;
  Fp2l_t T0, T1, T2, T3;

  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
  TEMP_ALLOC_FP2L(T0, 4);
  TEMP_ALLOC_FP2L(T1, 8);
  TEMP_ALLOC_FP2L(T2, 12);
  TEMP_ALLOC_FP2L(T3, 16);
  
  QapFp2_mul_no_red_o1 (T1, a->a1, b1, 2);
  QapFp2_mul_no_red_o1 (T2, a->a2, b2, 2);

  QapFp2_add_no_red (t0, a->a1, a->a2);
  QapFp2_add_no_red (t1, b1, b2);
  QapFp2_mul_no_red_o2 (T3, t0, t1);
  QapFp2l_add_no_red (T0, T1, T2);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T0->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T0->A1);
  QapFpl_sub_o2_ct (c->A0->A0, T3->A0, T3->A1, TEMP_PASS_MEM());
  QapFpl_add_o2_ct (c->A0->A1, T3->A0, T3->A1, TEMP_PASS_MEM());

  QapFp2_add_no_red (t0, a->a0, a->a1);
  QapFp2_mul_no_red_o2 (T3, t0, b1);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T1->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T1->A1);
  QapFpl_sub_o1_ct (T0->A0, T2->A0, T2->A1, 1, TEMP_PASS_MEM());
  QapFpl_add_no_red (T0->A1, T2->A0, T2->A1);
  QapFp2l_add_o2 (c->A1, T3, T0);

  QapFp2_add_no_red (t0, a->a0, a->a2);
  QapFp2_mul_no_red_o2 (T3, t0, b2);
  QapFpl_sub_o2_ct (T3->A0, T3->A0, T2->A0, TEMP_PASS_MEM());
  QapFpl_sub_no_red (T3->A1, T3->A1, T2->A1);
  QapFpl_add_o2_ct (c->A2->A0, T3->A0, T1->A0, TEMP_PASS_MEM()); 
  QapFpl_add_no_red (c->A2->A1, T3->A1, T1->A1);

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
  TEMP_FREE_FP2L(T0);
  TEMP_FREE_FP2L(T1);
  TEMP_FREE_FP2L(T2);
  TEMP_FREE_FP2L(T3);
}

int QapFp6_initialize_config (void) {
  if (QapFp_init (Fp6_config.vppow) < 0) printf ("Fp6_config memory error.\n");
  if (QapFp_init (Fp6_config.vppow2) < 0) printf ("Fp6_config memory error.\n");
  if (QapFp_init (Fp6_config.vp2pow) < 0) printf ("Fp6_config memory error.\n");
  if (QapFp_init (Fp6_config.vp2pow2) < 0) printf ("Fp6_config memory error.\n");
  if (QapFp_init (Fp6_config.vppowinv) < 0) printf ("Fp6_config memory error.\n");

  if (PAIR_CURVE == BN12) {
	  QapFp_set_u32(Fp6_config.vppow, 
		0x25236482, 0x40000001, 0x7080EB40, 0x00000006,
		0x18180000, 0x0000000C, 0xD9800000, 0x0000000B);
   
	  QapFp_set_u32(Fp6_config.vp2pow, 
		0x00000000, 0x00000000, 0x49B36240, 0x00000002,
		0x49090000, 0x00000006, 0xCD800000, 0x00000007);
    
    QapFp_mul (Fp6_config.vppow2, Fp6_config.vppow, Fp6_config.vppow);
    QapFp_neg (Fp6_config.vppow2, Fp6_config.vppow2);
    QapFp_mul (Fp6_config.vp2pow2, Fp6_config.vp2pow, Fp6_config.vp2pow);
    QapFp_modinv (Fp6_config.vppowinv, Fp6_config.vppow);
    QapFp_neg (Fp6_config.vppowinv, Fp6_config.vppowinv);

#if 0
    Fp6_config.QapFp6_mul = Qapv3minusxi_Fp6_mul;
    Fp6_config.QapFp6_squ = Qapv3minusxi_Fp6_squ;
    Fp6_config.QapFp6_inv = Qapv3minusxi_Fp6_inv;
    Fp6_config.QapFp6_mulv = Qapv3minusxi_Fp6_mulv;
    Fp6_config.QapFp6_ppow = Qapp1mod3_v3minusxi_Fp6_ppow;
    Fp6_config.QapFp6_p2pow = Qapp1mod3_v3minusxi_Fp6_p2pow;
    Fp6_config.QapFp6_p3pow = Qapp1mod3_v3minusxi_Fp6_p3pow;
    Fp6_config.QapFp6_mul_no_red = Qapv3minusxi_Fp6_mul_no_red;
    Fp6_config.QapFp6_mul_sparse01_no_red = Qapv3minusxi_Fp6_mul_sparse01_no_red;
    Fp6_config.QapFp6_mul_sparse12_no_red = Qapv3minusxi_Fp6_mul_sparse12_no_red;
#endif
  } else if (PAIR_CURVE == BN12tiny) {
    QapFp_set_u32(Fp6_config.vppow,
      0, 0, 0, 0,
      0, 0, 0x15C8EA, 0xDEB9FF1F);

    QapFp_set_u32(Fp6_config.vp2pow,
      0, 0, 0, 0,
      0, 0, 0xc5, 0x9AB172BB);

    QapFp_mul(Fp6_config.vppow2, Fp6_config.vppow, Fp6_config.vppow);
    QapFp_neg(Fp6_config.vppow2, Fp6_config.vppow2);
    QapFp_mul(Fp6_config.vp2pow2, Fp6_config.vp2pow, Fp6_config.vp2pow);
    QapFp_modinv(Fp6_config.vppowinv, Fp6_config.vppow);
    QapFp_neg(Fp6_config.vppowinv, Fp6_config.vppowinv);
  }
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }


  return ERR_SUCCESS;
}

void QapFp6_free_config (void) {
  QapFp_free (Fp6_config.vppow);
  QapFp_free (Fp6_config.vppow2);
  QapFp_free (Fp6_config.vp2pow);
  QapFp_free (Fp6_config.vp2pow2);
  QapFp_free (Fp6_config.vppowinv);
}


#include "QapFp12.h"
#include "Qappairing.h"
#include "QapTempMem.h"

Fp12_config_t Fp12_config;

/* Initialize an element in Fp12 by initializing two Fp6 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp12_init(QapFp12_t dst) {
	int err = 0;
	err = QapFp6_init(dst->a0);
	if (err < 0) return err;
	return QapFp6_init(dst->a1);
}

/* Free the memory allocated in src. */
void QapFp12_free (QapFp12_t src) {
	QapFp6_free(src->a0);
	QapFp6_free(src->a1);
}

/* Get a random element in Fp12. */
void QapFp12_rand (QapFp12_t c) {
  QapFp6_rand(c->a0);
  QapFp6_rand(c->a1);
}

/* Print the (regular non-Montgomery form) value of a */
void QapFp12_print (QapFp12_t a) {
  QapFp6_print(a->a0);
  //printf("\n");
  QapFp6_print(a->a1);
  //printf("\n");
}

/* Set c to 0. */
void QapFp12_set_zero (QapFp12_t c) {
  QapFp6_set_zero (c->a0);
  QapFp6_set_zero (c->a1);
}

/* Set c to 1. */
void QapFp12_set_one (QapFp12_t c) {
  QapFp6_set_one (c->a0);
  QapFp6_set_zero (c->a1);
}

/* Copy an Fp12 element a to an Fp12 element c. */
void QapFp12_copy(QapFp12_t c, QapFp12_t a) {
	QapFp6_copy(c->a0, a->a0);
	QapFp6_copy(c->a1, a->a1);
}

void QapFp12_select(QapFp12_t c, int bit, QapFp12_t a, QapFp12_t b) {
	QapFp6_select(c->a0, bit, a->a0, b->a0);
	QapFp6_select(c->a1, bit, a->a1, b->a1);
}

/* Compare two Fp12 elements a and b for equality. */
int QapFp12_cmpeq(QapFp12_t a, QapFp12_t b) {
#if MQAP
	return (QapFp6_cmpeq(a->a0, b->a0)  * QapFp6_cmpeq(a->a1, b->a1));
#else
	return (QapFp6_cmpeq(a->a0, b->a0) && QapFp6_cmpeq(a->a1, b->a1));
#endif
}

/* Negate an Fp12 element. */
void QapFp12_neg(QapFp12_t c, QapFp12_t a) {
	QapFp6_neg(c->a0, a->a0);
	QapFp6_neg(c->a1, a->a1);
}

/* Add two Fp12 elements coefficient wise. */
void QapFp12_add (QapFp12_t c, QapFp12_t a, QapFp12_t b) {
	QapFp6_add(c->a0, a->a0, b->a0);
	QapFp6_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp12 element b from an Fp12 element a coefficient wise. */
void QapFp12_sub (QapFp12_t c, QapFp12_t a, QapFp12_t b) {
	QapFp6_sub(c->a0, a->a0, b->a0);
	QapFp6_sub(c->a1, a->a1, b->a1);
}

/* Multiply an Fp12 element by an Fp6 element. */
void QapFp12_mulFp6(QapFp12_t c, QapFp12_t a, QapFp6_t b) {
	QapFp6_mul(c->a0, a->a0, b);
	QapFp6_mul(c->a1, a->a1, b);
}

/* Multiply an Fp12 element by an Fp2 element. */
void QapFp12_mulFp2(QapFp12_t c, QapFp12_t a, QapFp2_t b) {
	QapFp6_mulFp2(c->a0, a->a0, b);
	QapFp6_mulFp2(c->a1, a->a1, b);
}

/* Multiply an Fp12 element by an Fp element. */
void QapFp12_mulFp(QapFp12_t c, QapFp12_t a, QapFp_t b) {
	QapFp6_mulFp(c->a0, a->a0, b);
	QapFp6_mulFp(c->a1, a->a1, b);
}

/* Multiply two Fp12 elements. */
void QapFp12_mul (QapFp12_t c, QapFp12_t a, QapFp12_t b) {
  Qapw2minusv_Fp12_mul(c, a, b, TEMP_PASS_MEM12());
}

/* Square an Fp12 element. */
void QapFp12_squ (QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_Fp12_squ(c, a, TEMP_PASS_MEM12());
}

/* Compute the p-power Frobenius of an Fp12 element. */
void QapFp12_ppow(QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_v3minusxi_p1mod6_Fp12_ppow(c, a);
}

/* Compute the p^2-power Frobenius of an Fp12 element. */
void QapFp12_p2pow(QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_v3minusxi_p1mod6_Fp12_p2pow(c, a);
}

/* Compute the p^3-power Frobenius of an Fp12 element. */
void QapFp12_p3pow(QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_v3minusxi_p1mod6_Fp12_p3pow(c, a);
}

/* Compute the p^6-power Frobenius of an Fp12 element. */
void QapFp12_p6pow(QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_Fp12_p6pow(c, a);
}

/* Invert an Fp12 element. */
void QapFp12_inv (QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_Fp12_inv(c, a, TEMP_PASS_MEM12());
}

/* Multiply two Fp12 elements using lazy reduction.*/
void QapFp12_mul_lazy(QapFp12_t c, QapFp12_t a, QapFp12_t b) {
  Qapw2minusv_Fp12_mul_lazy(c, a, b, TEMP_PASS_MEM12());
}

/* Square an Fp12 element using lazy reduction.*/
void QapFp12_squ_lazy(QapFp12_t c, QapFp12_t a) {
  Qapw2minusv_Fp12_squ_lazy(c, a, TEMP_PASS_MEM12());
}

/* Multiply two Fp12 elements (the second one sparse) using lazy reduction.*/
void QapFp12_mul_sparse_twist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b00, QapFp2_t b01, QapFp2_t b11) {
  Qapw2minusv_Fp12_mul_sparse_twist_lazy(c, a, b00, b01, b11, TEMP_PASS_MEM12());
}

/* Multiply two Fp12 elements (the second one sparse) using lazy reduction.*/
void QapFp12_mul_sparse_untwist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b01, QapFp2_t b11, QapFp2_t b12) {
  Qapw2minusv_Fp12_mul_sparse_untwist_lazy(c, a, b01, b11, b12, TEMP_PASS_MEM12());
}

/* Multiply two Fp12 elements. This uses that Fp12 is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 20 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void Qapw2minusv_Fp12_mul(QapFp12_t c, QapFp12_t a, QapFp12_t b, TEMP_FORMAL) {
	QapFp6_t t0, t1, t2;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);
	TEMP_ALLOC_FP6(t2, 12);
	
	QapFp6_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp6_add(t1, b->a0, b->a1);	// t1 = b0 + b1
	QapFp6_mul(t2, t0, t1);			  // t2 = t0*t1 = (a0 + a1)*(b0 + b1) = a0*b0 + a1*b0 + a0*b1 + a1*b1

	QapFp6_mul(t0, a->a0, b->a0);		// t0 = a0*b0
	QapFp6_mul(t1, a->a1, b->a1);		// t1 = a1*b1
	
	QapFp6_sub(t2, t2, t0);					// t2 = t2 - t0 = a1*b0 + a0*b1 + a1*b1
	QapFp6_sub(c->a1, t2, t1);				// c1 = t2 - t1 = a1*b0 + a0*b1

	QapFp6_mulv(t1, t1);						  // t1 = v*t1 = v*a1*b1
	QapFp6_add(c->a0, t0, t1);				// c0 = t0 + t1 = a0*b0 + v*a1*b1

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
	TEMP_FREE_FP6(t2);
}

/* Square a general Fp12 element. This uses that Fp12 is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 22 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void Qapw2minusv_Fp12_squ(QapFp12_t c, QapFp12_t a, TEMP_FORMAL) {
	QapFp6_t t0, t1;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);
  
	QapFp6_sub(t0, a->a0, a->a1);		// t0 = a0 - a1
	QapFp6_mulv(t1, a->a1);					// t1 = v*a1
	QapFp6_sub(t1, a->a0, t1);				// t1 = a0 - v*a1
	QapFp6_mul(t0, t0, t1);			    // t0 = (a0 - a1)*(a0 - v*a1) = a0^2 + v*a1^2 - a0*a1 - v*a0*a1
	QapFp6_mul(t1, a->a0, a->a1);		// t1 = a0*a1
	QapFp6_add(t0, t0, t1);					// t0 =  a0^2 + v*a1^2 - v*a0*a1
	QapFp6_add(c->a1, t1, t1);			  // c1 = 2*t1 = 2*a0*a1
	QapFp6_mulv(t1, t1);						  // t1 = v*t1 = v*a0*a1
	QapFp6_add(c->a0, t0, t1);				// c0 = t0 + t1 = a0^2 + v*a1^2

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
}

/* Compute the p-power Frobenius of an Fp12 element. This uses a constant WPPOW, 
 * such that w^p = WPPOW*w, i.e. WPPOW = w^(p-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and Fp12 is an extension of Fp6 via w^2 - v, and p = 1 (mod 6),
 * we have WPPOW = xi^(p-1)/6 in Fp2.
 */
void Qapw2minusv_v3minusxi_p1mod6_Fp12_ppow(QapFp12_t c, QapFp12_t a) {
	QapFp6_ppow(c->a0, a->a0);
	QapFp6_ppow(c->a1, a->a1);
	QapFp6_mulFp2(c->a1, c->a1, Fp12_config.wppow);
}

/* Compute the p^2-power Frobenius of an Fp12 element. This uses a constant WP2POW, 
 * such that w^(p^2) = WP2POW*w, i.e. WP2POW = w^(p^2-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and Fp12 is an extension of Fp6 via w^2 - v, and p = 1 (mod 6),
 * we have WPPOW = xi^(p^2-1)/6 in Fp2.
 */
void Qapw2minusv_v3minusxi_p1mod6_Fp12_p2pow(QapFp12_t c, QapFp12_t a) {
	QapFp6_p2pow(c->a0, a->a0);
	QapFp6_p2pow(c->a1, a->a1);
	QapFp6_mulFp(c->a1, c->a1, Fp12_config.wp2pow);
}

/* Compute the p^3-power Frobenius of an Fp12 element. This uses a constant WP3POW, 
 * such that w^(p^3) = WP3POW*w, i.e. WP3POW = w^(p^3-1). In the case, where Fp6 is an 
 * extension of Fp2 via v^3 - xi and Fp12 is an extension of Fp6 via w^2 - v, and p = 1 (mod 6),
 * we have WP3POW = xi^(p^3-1)/6 in Fp2.
 */
void Qapw2minusv_v3minusxi_p1mod6_Fp12_p3pow(QapFp12_t c, QapFp12_t a) {
	QapFp6_p3pow(c->a0, a->a0);
	QapFp6_p3pow(c->a1, a->a1);
	QapFp6_mulFp2(c->a1, c->a1, Fp12_config.wp3pow);
}

/* Compute the p^2-power Frobenius of an Fp6 element. This is the same as 
 * conjugation in the quadratic extension over Fp6 if the extension is constructed
 * via a quadratic polynomial x^2 - alpha, alpha in Fp6. 
 */
void Qapw2minusv_Fp12_p6pow(QapFp12_t c, QapFp12_t a) {
	QapFp6_copy(c->a0, a->a0);
	QapFp6_neg(c->a1, a->a1);
}

/* Invert an Fp12 element. This uses that Fp12 is constructed over Fp6 via w^2 - v. 
 * This is Algorithm 23 in the ePrint version of Beuchat et al., Pairing 2010.
 */
void Qapw2minusv_Fp12_inv(QapFp12_t c, QapFp12_t a, TEMP_FORMAL) {
	QapFp6_t t0, t1;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);

	QapFp6_squ(t0, a->a0);				// t0 = a0^2
	QapFp6_squ(t1, a->a1);				// t1 = a1^2
	QapFp6_mulv(t1, t1);						// t1 = v*t1 = v*a1^2
	QapFp6_sub(t0, t0, t1);					// t0 = a0^2 - v*a1^2 = Norm(a)
	QapFp6_inv(t1, t0);				// t1 = t0^(-1) = Norm(a)^(-1)
	QapFp6_mul(c->a0, a->a0, t1);				// c0 = a0/(a0^2 - v*a1^2)
	QapFp6_neg(t1, t1);						// t1 = -1/(a0^2 - v*a1^2)
	QapFp6_mul(c->a1, a->a1, t1);				// c1 = a1*t1 = -a1/(a0^2 - v*a1^2)

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
}

/* Multiply two Fp12 elements using lazy reduction techniques. 
 * This uses that Fp12 is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. This is Algorithm 4 in Aranha et al. 
 */
void Qapw2minusv_Fp12_mul_lazy(QapFp12_t c, QapFp12_t a, QapFp12_t b, TEMP_FORMAL) {
	QapFp6_t t0, t1;
  Fp6l_t T0, T1, T2, T3;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);
	TEMP_ALLOC_FP6L(T0, 12);
	TEMP_ALLOC_FP6L(T1, 24);
	TEMP_ALLOC_FP6L(T2, 36);
	TEMP_ALLOC_FP6L(T3, 48);

  QapFp6_mul_no_red (T0, a->a0, b->a0);
  QapFp6_mul_no_red (T1, a->a1, b->a1);
  QapFp6_add (t0, a->a0, a->a1);
  QapFp6_add (t1, b->a0, b->a1);
  QapFp6_mul_no_red (T2, t0, t1);
  QapFp6l_add_o2 (T3, T0, T1);  
  QapFp6l_sub_o2 (T2, T2, T3);
  QapFp6l_red (c->a1, T2);

  QapFpl_sub_o2_ct (T2->A0->A0, T1->A2->A0, T1->A2->A1, TEMP_PASS_MEM());
  QapFpl_add_o2_ct (T2->A0->A1, T1->A2->A0, T1->A2->A1, TEMP_PASS_MEM());
  QapFp2l_copy (T2->A1, T1->A0);
  QapFp2l_copy (T2->A2, T1->A1);
  QapFp6l_add_o2 (T2, T0, T2);
  QapFp6l_red (c->a0, T2);

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
	TEMP_FREE_FP6L(T0);
	TEMP_FREE_FP6L(T1);
	TEMP_FREE_FP6L(T2);
	TEMP_FREE_FP6L(T3);
}

/* Multiply two Fp12 elements using lazy reduction techniques, 
 * where the second element b in sparse in the following sense:
 * b = b0 + b1*w, where b0 = b00 + b01*v and b1 = b11*v.
 * This uses that Fp12 is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. 
 */
void Qapw2minusv_Fp12_mul_sparse_twist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b00, QapFp2_t b01, QapFp2_t b11, TEMP_FORMAL) {
	QapFp6_t t0, t1;
  Fp6l_t T0, T1, T2, T3;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);
	TEMP_ALLOC_FP6L(T0, 12);
	TEMP_ALLOC_FP6L(T1, 24);
	TEMP_ALLOC_FP6L(T2, 36);
	TEMP_ALLOC_FP6L(T3, 48);

  QapFp6_mul_sparse01_no_red (T0, a->a0, b00, b01); // b0 = b00 + b01*v
  QapFp6_mulv (t0, a->a1);
  QapFp6_mulFp2_no_red (T1, t0, b11);
  QapFp6_add (t0, a->a0, a->a1);
  QapFp2_add (t1->a1, b01, b11);
  QapFp6_mul_sparse01_no_red (T2, t0, b00, t1->a1); // t1 = b00 + (b01 + b11)*v
  QapFp6l_add_o2 (T3, T0, T1);  
  QapFp6l_sub_o2 (T2, T2, T3);
  QapFp6l_red (c->a1, T2);

  QapFpl_sub_o2_ct (T2->A0->A0, T1->A2->A0, T1->A2->A1, TEMP_PASS_MEM());
  QapFpl_add_o2_ct (T2->A0->A1, T1->A2->A0, T1->A2->A1, TEMP_PASS_MEM());
  QapFp2l_copy (T2->A1, T1->A0);
  QapFp2l_copy (T2->A2, T1->A1);
  QapFp6l_add_o2 (T2, T0, T2);
  QapFp6l_red (c->a0, T2);

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
	TEMP_FREE_FP6L(T0);
	TEMP_FREE_FP6L(T1);
	TEMP_FREE_FP6L(T2);
	TEMP_FREE_FP6L(T3);
}

/* Similar as above, but now b = b0 + b1*w, where b0 = b01*v and b1 = b11*v + b12*v^2. */
void Qapw2minusv_Fp12_mul_sparse_untwist_lazy(QapFp12_t c, QapFp12_t a, QapFp2_t b01, QapFp2_t b11, QapFp2_t b12, TEMP_FORMAL) {
	QapFp6_t t0, t1;
  Fp6l_t T0, T1, T2, T3;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);
	TEMP_ALLOC_FP6L(T0, 12);
	TEMP_ALLOC_FP6L(T1, 24);
	TEMP_ALLOC_FP6L(T2, 36);
	TEMP_ALLOC_FP6L(T3, 48);

  QapFp6_mulv (t0, a->a0);
  QapFp6_mulFp2_no_red (T0, t0, b01);
  QapFp6_mul_sparse12_no_red (T1, a->a1, b11, b12); 
  QapFp6_add (t0, a->a0, a->a1);
  QapFp2_add (t1->a1, b01, b11);
  QapFp6_mul_sparse12_no_red (T2, t0, t1->a1, b12); // t1 = (b01 + b11)*v + b12*v^2
  QapFp6l_add_o2 (T3, T0, T1);  
  QapFp6l_sub_o2 (T2, T2, T3);
  QapFp6l_red (c->a1, T2);

  QapFpl_sub_o2_ct (T2->A0->A0, T1->A2->A0, T1->A2->A1, TEMP_PASS_MEM());
  QapFpl_add_o2_ct (T2->A0->A1, T1->A2->A0, T1->A2->A1, TEMP_PASS_MEM());
  QapFp2l_copy (T2->A1, T1->A0);
  QapFp2l_copy (T2->A2, T1->A1);
  QapFp6l_add_o2 (T2, T0, T2);
  QapFp6l_red (c->a0, T2);

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
	TEMP_FREE_FP6L(T0);
	TEMP_FREE_FP6L(T1);
	TEMP_FREE_FP6L(T2);
	TEMP_FREE_FP6L(T3);
}


/* Square an Fp12 element using lazy reduction techniques. 
 * This uses that Fp12 is constructed over Fp6 via w^2 - v and 
 * that the prime p is a few bits shorter than a multiple of the word size
 * and leaves some space for lazy reduction. This is Algorithm 5 in Aranha et al. 
 */
void Qapw2minusv_Fp12_squ_lazy(QapFp12_t c, QapFp12_t a, TEMP_FORMAL) {
	QapFp6_t t0, t1;
  Fp6l_t T0;

	TEMP_ALLOC_FP6(t0, 0);
	TEMP_ALLOC_FP6(t1, 6);
	TEMP_ALLOC_FP6L(T0, 12);

  QapFp6_add (t0, a->a0, a->a1);
  QapFp_sub (t1->a0->a0, a->a1->a2->a0, a->a1->a2->a1);
  QapFp_add (t1->a0->a1, a->a1->a2->a0, a->a1->a2->a1);
  QapFp2_copy (t1->a1, a->a1->a0);
  QapFp2_copy (t1->a2, a->a1->a1);
  QapFp6_add (t1, a->a0, t1);
  QapFp6_mul_no_red (T0, a->a0, a->a1);
  QapFp6l_red (c->a1, T0);
  
  QapFp6_mul_no_red (T0, t0, t1);
  QapFp6l_red (t0, T0);
  QapFp_sub (t1->a0->a0, c->a1->a2->a0, c->a1->a2->a1);
  QapFp_add (t1->a0->a1, c->a1->a2->a0, c->a1->a2->a1);
  QapFp2_copy (t1->a1, c->a1->a0);
  QapFp2_copy (t1->a2, c->a1->a1);
  QapFp6_add (t1, c->a1, t1);
  QapFp6_sub (c->a0, t0, t1);
  QapFp6_add (c->a1, c->a1, c->a1);

	TEMP_FREE_FP6(t0);
	TEMP_FREE_FP6(t1);
	TEMP_FREE_FP6L(T0);
}

int QapFp12_initialize_config (void) {

  if (QapFp_init (Fp12_config.wp2pow) < 0) printf ("Fp12_config memory error.\n");
  if (QapFp_init (Fp12_config.wp2pow3) < 0) printf ("Fp12_config memory error.\n");
  if (QapFp2_init (Fp12_config.wppow) < 0) printf ("Fp12_config memory error.\n");
  if (QapFp2_init (Fp12_config.wp3pow) < 0) printf ("Fp12_config memory error.\n");
  if (QapFp2_init (Fp12_config.wppow3) < 0) printf ("Fp12_config memory error.\n");
  if (QapFp2_init (Fp12_config.wppow3inv) < 0) printf ("Fp12_config memory error.\n");

  if (PAIR_CURVE == BN12) {
  	QapFp_set_u32(Fp12_config.wp2pow,
		0x00000000,0x00000000, 0x49B36240,0x00000002,
		0x49090000,0x00000006, 0xCD800000,0x00000008);

    QapFp_set_u32 (Fp12_config.wppow->a0,
		0x1B377619,0x212E7C8C, 0xB6499B50,0xA846953F,
		0x85097492,0x4D3F77C2, 0xE17DE6C0,0x6F2A6DE9);
  
    QapFp_set_u32 (Fp12_config.wppow->a1,
		0x09EBEE69,0x1ED18375, 0x03EAB22F,0x57B96AC8,
		0xDC178B6D,0xB2C08850, 0xC582193F,0x90D5922A);

    QapFp_set_u32 (Fp12_config.wp3pow->a0,
		0x01439AB0,0x9C60B248, 0xF398C5D7,0x7B755F92,
		0xB9EDC5F1,0x9D287354, 0x5BE47115,0x1A747E4E);
  
    QapFp_set_u32 (Fp12_config.wp3pow->a1,
		0x23DFC9D1,0xA39F4DB8, 0xC69B87A8,0x848AA075,
		0xA7333A0E,0x62D78CBF, 0x4B1B8EEA,0xE58B81C5);

    QapFp_mul (Fp12_config.wp2pow3, Fp12_config.wp2pow, Fp12_config.wp2pow);
    QapFp_mul (Fp12_config.wp2pow3, Fp12_config.wp2pow3, Fp12_config.wp2pow);
    QapFp2_squ (Fp12_config.wppow3, Fp12_config.wppow);
    QapFp2_mul (Fp12_config.wppow3, Fp12_config.wppow3, Fp12_config.wppow);
    QapFp2_inv (Fp12_config.wppow3inv, Fp12_config.wppow3);

#if 0
    Fp12_config.QapFp12_mul = Qapw2minusv_Fp12_mul;
    Fp12_config.QapFp12_squ = Qapw2minusv_Fp12_squ;
    Fp12_config.QapFp12_inv = Qapw2minusv_Fp12_inv;
    Fp12_config.QapFp12_ppow = Qapw2minusv_v3minusxi_p1mod6_Fp12_ppow;
    Fp12_config.QapFp12_p2pow = Qapw2minusv_v3minusxi_p1mod6_Fp12_p2pow;
    Fp12_config.QapFp12_p3pow = Qapw2minusv_v3minusxi_p1mod6_Fp12_p3pow;
    Fp12_config.QapFp12_p6pow = Qapw2minusv_Fp12_p6pow; 
    Fp12_config.QapFp12_mul_lazy = Qapw2minusv_Fp12_mul_lazy;
    Fp12_config.QapFp12_squ_lazy = Qapw2minusv_Fp12_squ_lazy;
    Fp12_config.QapFp12_mul_sparse_twist_lazy = Qapw2minusv_Fp12_mul_sparse_twist_lazy;
    Fp12_config.QapFp12_mul_sparse_untwist_lazy = Qapw2minusv_Fp12_mul_sparse_untwist_lazy;
#endif
  } else if (PAIR_CURVE == BN12tiny) {
    QapFp_set_u32(Fp12_config.wp2pow,
      0, 0, 0, 0,
      0, 0, 0xc5, 0x9AB172BC);

    QapFp_set_u32(Fp12_config.wppow->a0,
      0, 0, 0, 0,
      0, 0, 0xF5B82, 0xC0AFD83C);

    QapFp_set_u32(Fp12_config.wppow->a1,
      0, 0, 0, 0,
      0, 0, 0x66E2D, 0xB8BB999F);

    QapFp_set_u32(Fp12_config.wp3pow->a0,
      0, 0, 0, 0,
      0, 0, 0xCBFF4, 0x9918E901);

    QapFp_set_u32(Fp12_config.wp3pow->a1,
      0, 0, 0, 0,
      0, 0, 0x909BB, 0xE05288DA);

    QapFp_mul(Fp12_config.wp2pow3, Fp12_config.wp2pow, Fp12_config.wp2pow);
    QapFp_mul(Fp12_config.wp2pow3, Fp12_config.wp2pow3, Fp12_config.wp2pow);
    QapFp2_squ(Fp12_config.wppow3, Fp12_config.wppow);
    QapFp2_mul(Fp12_config.wppow3, Fp12_config.wppow3, Fp12_config.wppow);
    QapFp2_inv(Fp12_config.wppow3inv, Fp12_config.wppow3);
  }
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }

  return ERR_SUCCESS;
}

void QapFp12_free_config (void) {
  QapFp2_free (Fp12_config.wppow);
  QapFp2_free (Fp12_config.wppow3);
  QapFp_free (Fp12_config.wp2pow);
  QapFp_free (Fp12_config.wp2pow3);
  QapFp2_free (Fp12_config.wp3pow);
  QapFp2_free (Fp12_config.wppow3inv);
}

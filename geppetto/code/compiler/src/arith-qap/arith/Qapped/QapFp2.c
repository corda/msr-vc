#include "QapFp2.h"
#include "QapTempMem.h"

Fp2_config_t Fp2_config;

/************** Regular Fp2 functions ***************/

/* Initialize an element in Fp2 by initializing two Fp elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp2_init (QapFp2_t dst) {
	int err = 0;
	if ((err = QapFp_init(dst->a0)) < 0) return err;
	return QapFp_init(dst->a1);
}

/* Free the memory allocated in src. */
void QapFp2_free (QapFp2_t src) {
	QapFp_free(src->a0);
	QapFp_free(src->a1);
}

/* Get a random element in Fp2. */
void QapFp2_rand (QapFp2_t c) {
  QapFp_rand(c->a0);
  QapFp_rand(c->a1);
}

/* Set c to the value of a0 + a1*i. */
void QapFp2_set (QapFp2_t c, QapFp_t a0, QapFp_t a1) {
  QapFp_set(c->a0, a0);
  QapFp_set(c->a1, a1);
}

/* Set c to the value x. */
void QapFp2_set_ui (QapFp2_t c, uint64_t x) {
  QapFp_set_ui (c->a0, x);
  QapFp_set_ui (c->a1, 0);
}

/* Set c to 0. */
void QapFp2_set_zero (QapFp2_t c) {
  QapFp_set_ui (c->a0, 0);
  QapFp_set_ui (c->a1, 0);
}

/* Set c to 1. */
void QapFp2_set_one (QapFp2_t c) {
  QapFp_set_ui (c->a0, 1);
  QapFp_set_ui (c->a1, 0);
}

/* Print the (regular non-Montgomery form) value of a */
void QapFp2_print (QapFp2_t a) {
  QapFp_print(a->a0);
  printf("  ");
  QapFp_print(a->a1);
  //printf("\n");
}

/* Copy an Fp2 element a to an Fp2 element c. */
void QapFp2_copy(QapFp2_t c, QapFp2_t a) {
	QapFp_copy(c->a0, a->a0);
	QapFp_copy(c->a1, a->a1);
}

void QapFp2_select(QapFp2_t c, int bit, QapFp2_t a, QapFp2_t b) {
	QapFp_select(c->a0, bit, a->a0, b->a0);
	QapFp_select(c->a1, bit, a->a1, b->a1);
}

/* Compare two Fp2 elements a and b for equality. */
int QapFp2_cmpeq(QapFp2_t a, QapFp2_t b) {
#if MQAP
	// forcing evaluation of both tests; not great if one of the sides is known at compile-time.
	return (QapFp_cmpeq(a->a0, b->a0)  * QapFp_cmpeq(a->a1, b->a1));
#else
	return (QapFp_cmpeq(a->a0, b->a0) && QapFp_cmpeq(a->a1, b->a1));
#endif
	// was: return (QapFp_cmp(a->a0, b->a0) == 0 && QapFp_cmp(a->a1, b->a1) == 0);
}

/* Negate an Fp2 element. */
void QapFp2_neg(QapFp2_t c, QapFp2_t a) {
	QapFp_neg(c->a0, a->a0);
	QapFp_neg(c->a1, a->a1);
}

/* Add two Fp2 elements coefficient wise. */
void QapFp2_add (QapFp2_t c, QapFp2_t a, QapFp2_t b) {
	QapFp_add(c->a0, a->a0, b->a0);
	QapFp_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp2 element b from an Fp2 element a coefficient wise. */
void QapFp2_sub (QapFp2_t c, QapFp2_t a, QapFp2_t b) {
	QapFp_sub(c->a0, a->a0, b->a0);
	QapFp_sub(c->a1, a->a1, b->a1);
}

/* Divide an Fp2 element a by 2 coefficient wise. */
void QapFp2_div2 (QapFp2_t c, QapFp2_t a) {
	QapFp_div2(c->a0, a->a0);
	QapFp_div2(c->a1, a->a1);
}

/* Multiply an Fp2 element by 3. */
void QapFp2_mul3 (QapFp2_t c, QapFp2_t a) {
  QapFp_mul3 (c->a0, a->a0);
  QapFp_mul3 (c->a1, a->a1);
}

/* Multiply an Fp2 element by an Fp element. */
void QapFp2_mulFp(QapFp2_t c, QapFp2_t a, QapFp_t b) {
	QapFp_mul(c->a0, a->a0, b);
	QapFp_mul(c->a1, a->a1, b);
}

/* Multiply two Fp2 elements. */
void QapFp2_mul(QapFp2_t c, QapFp2_t a, QapFp2_t b) {
  //Fp2_config.QapFp2_mul(c, a, b, Fp2_config.mem);
  QapFp2_mul_lazy (c, a, b);
}

/* Square an Fp2 element. */
void QapFp2_squ(QapFp2_t c, QapFp2_t a){
  QapFp2_squ_lazy(c, a);
}

/* Multiply an Fp2 element by the special element xi.*/
void QapFp2_mulxi(QapFp2_t c, QapFp2_t a) {
  Qapp3mod4_Fp2_mulxi(c, a, Fp2_config.mem);
}

/* Multiply an Fp2 element by the special element i.*/
void QapFp2_muli(QapFp2_t c, QapFp2_t a) {
  Qapp3mod4_Fp2_muli(c, a, Fp2_config.mem);
}

/* Multiply an Fp2 element by an element i*b, where b is in Fp. */
void QapFp2_muliFp(QapFp2_t c, QapFp2_t a, QapFp_t b) {
  Qapp3mod4_Fp2_muliFp(c, a, b, Fp2_config.mem);
}

/* Norm of an Fp2 element. */
void QapFp2_norm(QapFp_t c, QapFp2_t a) {
  Qapp3mod4_Fp2_norm(c, a, Fp2_config.mem);
}

/* Invert an Fp2 element. */
void QapFp2_inv(QapFp2_t c, QapFp2_t a) {
  Qapp3mod4_Fp2_inv(c, a, Fp2_config.mem);
}

/* Compute the p-power Frobenius of an Fp2 element. */
void QapFp2_ppow(QapFp2_t c, QapFp2_t a) {
	 Qapbinomial_Fp2_ppow(c, a);
}

void QapFp2l_mulxi (Fp2l_t c, Fp2l_t a) {
  Qapp3mod4_Fp2l_mulxi(c, a, Fp2_config.mem);
}


/* Compute the p-power Frobenius of an Fp2 element. This function assumes that the field extension
 * is constructed via a binomial i^2 - alpha, alpha in Fp.
 * This is the same as (complex) conjugation, i.e. Fp_ppow(a0 + a1*i) = a0 - a1*i.
 * When we need this, it might be simpler to just negate the second coordinate in place.
 */
void Qapbinomial_Fp2_ppow(QapFp2_t c, QapFp2_t a) {
	QapFp_copy(c->a0, a->a0);
	QapFp_neg(c->a1, a->a1);
}

/* Multiply two Fp2 elements. This function assumes that p = 3 (mod4)
 * such that the quadratic extension is Fp2 = Fp(i) with i^2 = -1.
 * Then c = c0 + c1*i = a*b = (a0 + a1*i)*(b0 + b1*i) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)*i.
 */
void Qapp3mod4_Fp2_mul(QapFp2_t c, QapFp2_t a, QapFp2_t b, TEMP_FORMAL) {
	QapFp_t t0, t1, t2;
	TEMP_ALLOC_FP(t0, 0);
	TEMP_ALLOC_FP(t1, 1);
	TEMP_ALLOC_FP(t2, 2);

	QapFp_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp_add(t1, b->a0, b->a1); // t1 = b0 + b1
	QapFp_mul(t2, t0, t1);			  // t2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + b0*a1 + a1*b1

	QapFp_mul(t0, a->a0, b->a0);	// t0 = a0*b0
	QapFp_mul(t1, a->a1, b->a1);	// t1 = a1*b1

	QapFp_sub(c->a0, t0, t1);		// c0 = a0*b0 - a1*b1
	QapFp_add(t0, t0, t1);			// t0 = a0*b0 + a1*b1
	QapFp_sub(c->a1, t2, t0);		// c1 = t2 - t0 = a0*b1 + a1*b0

	TEMP_FREE_FP(t0);
	TEMP_FREE_FP(t1);
	TEMP_FREE_FP(t2);
}

/* Square an Fp2 element. This function assumes that p = 3 (mod4)
 * such that the quadratic extension is as described above and thus
 * c = c0 + c1*i = a^2 = (a0 + a1*i)^2 = (a0^2 - a1^2) + (2*a0*a1)*i.
 */
void Qapp3mod4_Fp2_squ(QapFp2_t c, QapFp2_t a, TEMP_FORMAL) {
	QapFp_t t0, t1, t2;
	TEMP_ALLOC_FP(t0, 0);
	TEMP_ALLOC_FP(t1, 1);
	TEMP_ALLOC_FP(t2, 2);
  
	QapFp_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp_sub(t1, a->a0, a->a1);	// t1 = a0 - a1

	QapFp_mul(t2, a->a0, a->a1);	// t2 = a0*a1

	QapFp_mul(c->a0, t0, t1);		// c0 = (a0 - a1)*(a0 + a1)
	QapFp_add(c->a1, t2, t2);		// c1 = 2*a0*a1

	TEMP_FREE_FP(t0);
	TEMP_FREE_FP(t1);
	TEMP_FREE_FP(t2);
}

/* Multiply an Fp2 element by the special element xi = 1+i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 * We have c = c0 + c1*i = (a0 - a1) + (a0 + a1)*i. */
void Qapp3mod4_Fp2_mulxi(QapFp2_t c, QapFp2_t a, TEMP_FORMAL) {
	QapFp_t t0;
	TEMP_ALLOC_FP(t0, 0);

  QapFp_sub(t0, a->a0, a->a1);
	QapFp_add(c->a1, a->a0, a->a1);
  QapFp_copy(c->a0, t0);

	TEMP_FREE_FP(t0);
}

/* Multiply an Fp2 element by the special element i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 * We have c = c0 + c1*i = -a1 + a0*i. */
void Qapp3mod4_Fp2_muli(QapFp2_t c, QapFp2_t a, TEMP_FORMAL) {
	QapFp_t t0;
	TEMP_ALLOC_FP(t0, 0);
  
  QapFp_neg(t0, a->a1);
	QapFp_copy(c->a1, a->a0);
  QapFp_copy(c->a0, t0);

	TEMP_FREE_FP(t0);
}

/* Multiply an Fp2 element by an element i*b, where b is in Fp. */
void Qapp3mod4_Fp2_muliFp(QapFp2_t c, QapFp2_t a, QapFp_t b, TEMP_FORMAL) {
  QapFp_t t0;
  TEMP_ALLOC_FP(t0, 0);

  QapFp_mul(t0, a->a0, b);
  QapFp_mul(c->a0, a->a1, b);
  QapFp_neg(c->a0, c->a0);
  QapFp_copy(c->a1, t0);
}

/* Compute the norm (in Fp) of an Fp2 element. This function assumes that p = 3 (mod 4)
 * such that the quadratic extension is as described above and the norm can be computed as
 * c = N(a) = N(a0 + a1*i) = (a0 + a1*i)*(a0 - a1*i) = a0^2 + a1^2.
 */
void Qapp3mod4_Fp2_norm(QapFp_t c, QapFp2_t a, TEMP_FORMAL) {
	QapFp_t t0;
	TEMP_ALLOC_FP(t0, 0);
  
	QapFp_mul(c, a->a0, a->a0);	// c = a0^2
	QapFp_mul(t0, a->a1, a->a1);	// t0 = a1^2

	QapFp_add(c, c, t0);			    // c = a0^2 + a1^2
}

/* Invert an Fp2 element. This function assumes that p = 3 (mod4)
 * such that the quadratic extension is as described above and thus
 * c = c0 + c1*i = a^-1 = (a0 + a1*i)^-1 = (a0 - a1*i)/(a0^2 + a1^2) = a0/N(a) - a1/N(a)*i.
  * mem should have room for 3 Fp elements.
 */
void Qapp3mod4_Fp2_inv(QapFp2_t c, QapFp2_t a, TEMP_FORMAL) {
	QapFp_t t0, t1;
  TEMP_ALLOC_FP(t0, 0);
  TEMP_ALLOC_FP(t1, 1);

  Qapp3mod4_Fp2_norm (t0, a, TEMP_PASS(1) ); // t0 = a0^2 + a1^2

	QapFp_modinv(t1, t0);			 // t1 = (a0^2 + a1^2)^-1

	QapFp_mul(c->a0, a->a0, t1);	 // c0 = a0/(a0^2 + a1^2)
	QapFp_mul(c->a1, a->a1, t1);	 // c1 = a1/(a0^2 + a1^2)
	QapFp_neg(c->a1, c->a1);		 // c1 = -a1/(a0^2 + a1^2)

	TEMP_FREE_FP(t0);
	TEMP_FREE_FP(t1);
}

/************** Fp2 functions for lazy reduction ***************/
/* Functions assume that there is enough space such that no overflow occurs. */

/* Initialize a long Fp2 element. */
int QapFp2l_init(Fp2l_t dst) {
	int err = 0;
	if ((err = QapFpl_init(dst->A0)) < 0) return err;
	return QapFpl_init(dst->A1);
}

/* Free the memory allocated in src. */
void QapFp2l_free (Fp2l_t src) {
	QapFpl_free(src->A0);
	QapFpl_free(src->A1);
}

/* Copy an Fp2l element a to an Fp2l element c. */
void QapFp2l_copy(Fp2l_t C, Fp2l_t A) {
	QapFpl_copy(C->A0, A->A0);
	QapFpl_copy(C->A1, A->A1);
}

/* Reduction function calling the Fp reduction for each coefficient. */
void QapFp2l_red(QapFp2_t c, Fp2l_t A) {
	QapFpl_red(c->a0, A->A0);
	QapFpl_red(c->a1, A->A1);
}

/* Add two Fp2 elements coefficient wise without reducing mod p for lazy reduction. */ 
void QapFp2_add_no_red (QapFp2_t c, QapFp2_t a, QapFp2_t b) {
	QapFp_add_no_red(c->a0, a->a0, b->a0);
	QapFp_add_no_red(c->a1, a->a1, b->a1);
}

/* Add two long Fp2 elements coefficient wise without reducing mod p for lazy reduction. */ 
void QapFp2l_add_no_red (Fp2l_t C, Fp2l_t A, Fp2l_t B) {
	QapFpl_add_no_red(C->A0, A->A0, B->A0);
	QapFpl_add_no_red(C->A1, A->A1, B->A1);
}

/* Subtraction of Fp2l elements, option 1 in Aranha et al. */
void QapFp2l_sub_o1 (Fp2l_t c, Fp2l_t a, Fp2l_t b, int h) {
  QapFpl_sub_o1_ct(c->A0, a->A0, b->A0, h, TEMP_PASS_MEM());
  QapFpl_sub_o1_ct(c->A1, a->A1, b->A1, h, TEMP_PASS_MEM());
}

/* Subtraction of Fp2l elements, constant-time version of option 2 in Aranha et al. */
void QapFp2l_sub_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b) {
  QapFpl_sub_o2_ct(c->A0, a->A0, b->A0, TEMP_PASS_MEM());
  QapFpl_sub_o2_ct(c->A1, a->A1, b->A1, TEMP_PASS_MEM());
}

/* Addition of Fp2l elements, constant-time version of "option 2" addition in Aranha et al. */
void QapFp2l_add_o2 (Fp2l_t c, Fp2l_t a, Fp2l_t b) {
  QapFpl_add_o2_ct(c->A0, a->A0, b->A0, TEMP_PASS_MEM());
  QapFpl_add_o2_ct(c->A1, a->A1, b->A1, TEMP_PASS_MEM());
}

/* Multiply an Fp2 element by an Fp element. */
void QapFp2_mulFp_no_red (Fp2l_t c, QapFp2_t a, QapFp_t b) {
	QapFp_mul_no_red (c->A0, a->a0, b);
	QapFp_mul_no_red (c->A1, a->a1, b);
}

/* Multiply two Fp2 elements without reduction. */
void QapFp2_mul_no_red_o1(Fp2l_t c, QapFp2_t a, QapFp2_t b, int h) {
  Qapp3mod4_Fp2_mul_no_red_o1(c, a, b, h, TEMP_PASS_MEM2());
}

void QapFp2_mul_no_red_o2(Fp2l_t c, QapFp2_t a, QapFp2_t b) {
  Qapp3mod4_Fp2_mul_no_red_o2(c, a, b, TEMP_PASS_MEM2());
}

/* Multiply in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void QapFp2_mul_lazy (QapFp2_t c, QapFp2_t a, QapFp2_t b) {
  Qapp3mod4_Fp2_mul_no_red_o2(Fp2_config.T0, a, b, TEMP_PASS_MEM2());
  QapFp2l_red(c, Fp2_config.T0);
}

/* Square a Fp2 element without reduction. */
void QapFp2_squ_no_red(Fp2l_t c, QapFp2_t a) {
  Qapp3mod4_Fp2_squ_no_red(c, a, TEMP_PASS_MEM2());
}

/* Square in Fp2 using the lazy reduction multiplications and 2 reductions only. */
void QapFp2_squ_lazy (QapFp2_t c, QapFp2_t a) {
  Qapp3mod4_Fp2_squ_no_red(Fp2_config.T0, a, TEMP_PASS_MEM2());
  QapFp2l_red(c, Fp2_config.T0);
}

/* Multiply an Fp2l element by the special element xi = 1+i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 * We have c = c0 + c1*i = (a0 - a1) + (a0 + a1)*i. */
void Qapp3mod4_Fp2l_mulxi(Fp2l_t c, Fp2l_t a, TEMP_FORMAL) {
	QapFpl_t T0;
	TEMP_ALLOC_FPL(T0, 0);

  QapFpl_sub_o2_ct (T0, a->A0, a->A1, TEMP_PASS_MEM());
	QapFpl_add_no_red (c->A1, a->A0, a->A1);
  QapFpl_copy(c->A0, T0);

	TEMP_FREE_FPL(T0);
}

/* Multiply two Fp2 elements without reduction resulting in a long Fp2l element.
 * This function assumes that p = 3 (mod4). It implements the function in Alg. 2
 * in Aranha et al., Eurocrypt 2011.
 * This function needs 8 n-word Fp elements for temporary memory in mem.
 */
void Qapp3mod4_Fp2_mul_no_red_o1(Fp2l_t C, QapFp2_t a, QapFp2_t b, int h, TEMP_FORMAL) {
	QapFp_t t0, t1;
  QapFpl_t T0, T1, T2;
  TEMP_ALLOC_FP(t0, 0);
  TEMP_ALLOC_FP(t1, 1);
  TEMP_ALLOC_FPL(T0, 2);
  TEMP_ALLOC_FPL(T1, 4);
  TEMP_ALLOC_FPL(T2, 6);

  QapFp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
  QapFp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
  QapFp_mul_no_red(T2, t0, t1);		    // T2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1

	QapFp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	QapFp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1

	QapFpl_sub_o1_ct(C->A0, T0, T1, h, TEMP_PASS_MEM());			// C0 = a0*b0 - a1*b1
  QapFpl_add_no_red(T0, T0, T1);		    // T0 = a0*b0 + a1*b1
	QapFpl_sub_no_red(C->A1, T2, T0);		// C1 = T2 - T0 = a0*b1 + a1*b0

  TEMP_FREE_FP(t0);
  TEMP_FREE_FP(t1);
  TEMP_FREE_FPL(T0);
  TEMP_FREE_FPL(T1);
  TEMP_FREE_FPL(T2);
}

void Qapp3mod4_Fp2_mul_no_red_o2(Fp2l_t C, QapFp2_t a, QapFp2_t b, TEMP_FORMAL) {
	QapFp_t t0, t1;
  QapFpl_t T0, T1, T2;
  TEMP_ALLOC_FP(t0, 0);
  TEMP_ALLOC_FP(t1, 1);
  TEMP_ALLOC_FPL(T0, 2);
  TEMP_ALLOC_FPL(T1, 4);
  TEMP_ALLOC_FPL(T2, 6);

	QapFp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp_add_no_red(t1, b->a0, b->a1);  // t1 = b0 + b1
	QapFp_mul_no_red(T2, t0, t1);		    // T2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + a1*b0 + a1*b1

	QapFp_mul_no_red(T0, a->a0, b->a0);	// T0 = a0*b0
	QapFp_mul_no_red(T1, a->a1, b->a1);	// T1 = a1*b1

	QapFpl_sub_o2_ct(C->A0, T0, T1, TEMP_PASS_MEM());			// C0 = a0*b0 - a1*b1
	QapFpl_add_no_red(T0, T0, T1);		    // T0 = a0*b0 + a1*b1
	QapFpl_sub_no_red(C->A1, T2, T0);		// C1 = T2 - T0 = a0*b1 + a1*b0

  TEMP_FREE_FP(t0);
  TEMP_FREE_FP(t1);
  TEMP_FREE_FPL(T0);
  TEMP_FREE_FPL(T1);
  TEMP_FREE_FPL(T2);
}

/* Squaring an Fp2 element without reduction resulting in a long Fp2l element. 
 * This function assumes that p = 3 (mod 4). It implements the function in Alg. 7
 * in Aranha et al., Eurocrypt 2011.
 */
void Qapp3mod4_Fp2_squ_no_red(Fp2l_t C, QapFp2_t a, TEMP_FORMAL) {
	QapFp_t t0, t1, t2;
  TEMP_ALLOC_FP(t0, 0);
  TEMP_ALLOC_FP(t1, 1);
  TEMP_ALLOC_FP(t2, 2);
  
	QapFp_add_no_red(t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp_sub(t1, a->a0, a->a1);	        // t1 = a0 - a1
	QapFp_add_no_red(t2, a->a0, a->a0);  // t2 = 2*a0
	QapFp_mul_no_red(C->A0, t0, t1);	    // C0 = (a0 + a1)*(a0 - a1) = a0^2 - a1^2
	QapFp_mul_no_red(C->A1, t2, a->a1);  // C1 = 2*a0*a1

  TEMP_FREE_FP(t0);
  TEMP_FREE_FP(t1);
  TEMP_FREE_FP(t2);
}


int QapFp2_initialize_config (void) {
   
#if 0
  Fp2_config.mem = NULL;
  Fp2_config.mem = (void *) malloc (8 * FP_CONFIG_N() * sizeof (uint64_t));
  if (Fp2_config.mem == NULL) return ERR_OUT_OF_MEMORY;
#endif
  
  if (QapFp2l_init (Fp2_config.T0) < 0) printf ("Fp2_config memory error.\n");

  return ERR_SUCCESS;
}

void QapFp2_free_config (void) {
#if 0
  free (Fp2_config.mem);
#endif
}

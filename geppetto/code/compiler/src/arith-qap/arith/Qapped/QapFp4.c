#include "QapFp4.h"
#include "Qappairing.h"
#include "QapTempMem.h"

Fp4_config_t Fp4_config;

/* Initialize an element in Fp4 by initializing two Fp2 elements.
 * Return > 0 on success and
 * return < 0 on error. 
 */
int QapFp4_init(QapFp4_t dst) {
	int err = 0;
	err = QapFp2_init(dst->a0);
	if (err < 0) return err;
	return QapFp2_init(dst->a1);
}

/* Free the memory allocated in src. */
void QapFp4_free (QapFp4_t src) {
	QapFp2_free(src->a0);
	QapFp2_free(src->a1);
}

/* Get a random element in Fp4. */
void QapFp4_rand (QapFp4_t c) {
  QapFp2_rand(c->a0);
  QapFp2_rand(c->a1);
}

/* Set c to 0. */
void QapFp4_set_zero (QapFp4_t c) {
  QapFp2_set_zero (c->a0);
  QapFp2_set_zero (c->a1);
}

/* Set c to 1. */
void QapFp4_set_one (QapFp4_t c) {
  QapFp2_set_one (c->a0);
  QapFp2_set_zero (c->a1);
}

/* Copy an Fp4 element a to an Fp4 element c. */
void QapFp4_copy(QapFp4_t c, QapFp4_t a) {
	QapFp2_copy(c->a0, a->a0);
	QapFp2_copy(c->a1, a->a1);
}

/* Print the (regular non-Montgomery form) value of a */
void QapFp4_print (QapFp4_t a) {
  QapFp2_print(a->a0);
  printf("\n");
  QapFp2_print(a->a1);
  printf("\n");
}

/* Compare two Fp4 elements a and b for equality. */
int QapFp4_cmpeq(QapFp4_t a, QapFp4_t b) {
  return (QapFp2_cmpeq(a->a0, b->a0) && QapFp2_cmpeq(a->a1, b->a1));
}

/* Negate an Fp4 element. */
void QapFp4_neg(QapFp4_t c, QapFp4_t a) {
	QapFp2_neg(c->a0, a->a0);
	QapFp2_neg(c->a1, a->a1);
}

/* Add two Fp4 elements coefficient wise. */
void QapFp4_add (QapFp4_t c, QapFp4_t a, QapFp4_t b) {
	QapFp2_add(c->a0, a->a0, b->a0);
	QapFp2_add(c->a1, a->a1, b->a1);
}

/* Subtract an Fp4 element b from an Fp4 element a coefficient wise. */
void QapFp4_sub (QapFp4_t c, QapFp4_t a, QapFp4_t b) {
	QapFp2_sub(c->a0, a->a0, b->a0);
	QapFp2_sub(c->a1, a->a1, b->a1);
}

/* Divide an Fp4 element a by 2 coefficient wise. */
void QapFp4_div2 (QapFp4_t c, QapFp4_t a) {
	QapFp2_div2(c->a0, a->a0);
	QapFp2_div2(c->a1, a->a1);
}

/* Multiply an Fp4 element by an Fp2 element. */
void QapFp4_mulFp2(QapFp4_t c, QapFp4_t a, QapFp2_t b) {
	QapFp2_mul(c->a0, a->a0, b);
	QapFp2_mul(c->a1, a->a1, b);
}

/* Multiply an Fp4 element by an Fp element. */
void QapFp4_mulFp(QapFp4_t c, QapFp4_t a, QapFp_t b) {
	QapFp2_mulFp(c->a0, a->a0, b);
	QapFp2_mulFp(c->a1, a->a1, b);
}

/* Multiply two Fp4 elements. */
void QapFp4_mul (QapFp4_t c, QapFp4_t a, QapFp4_t b) {
  //Qaps2minusxi_Fp4_mul(c, a, b, Fp4_config.mem);
  QapFp4_mul_lazy (c, a, b);
}

/* Multiply two Fp4 elements using lazy reduction. */
void QapFp4_mul_lazy (QapFp4_t c, QapFp4_t a, QapFp4_t b) {
  Qaps2minusxi_Fp4_mul_lazy(c, a, b, TEMP_PASS_MEM4());
}


/* Compute 3*a. */
void QapFp4_mul3 (QapFp4_t c, QapFp4_t a) {
  QapFp4_t t;
  TEMP_ALLOC_MEMFP4(t, 0);
  
  QapFp4_add (t, a, a);
  QapFp4_add (c, t, a);
  TEMP_FREE_MEMFP4(t);
}

/* Square an Fp4 element. */
void QapFp4_squ (QapFp4_t c, QapFp4_t a) {
  //Qaps2minusxi_Fp4_squ(c, a, Fp4_config.mem);
  QapFp4_squ_lazy (c, a);
}

/* Square an Fp4 element using lazy reduction. */
void QapFp4_squ_lazy (QapFp4_t c, QapFp4_t a) {
  Qaps2minusxi_Fp4_squ_lazy(c, a, TEMP_PASS_MEM4());
}


/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void QapFp4_squ_sep (QapFp2_t c0, QapFp2_t c1, QapFp2_t a0, QapFp2_t a1) {
  Qaps2minusxi_Fp4_squ_sep(c0, c1, a0, a1, TEMP_PASS_MEM4());
}

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2. */
void QapFp4_muls (QapFp4_t c, QapFp4_t a) {
  Qaps2minusxi_Fp4_muls(c, a, TEMP_PASS_MEM4());
}

/* Multiply an Fp4 element by the special element i in Fp2.*/
void QapFp4_muli(QapFp4_t c, QapFp4_t a) {
  Qapp3mod4_Fp4_muli(c, a);
}

/* Multiply an Fp4 element by the special element i+1 in Fp2.*/
void QapFp4_mulxi(QapFp4_t c, QapFp4_t a) {
  Qapp3mod4_Fp4_mulxi(c, a);
}

/* Compute the p^2-power Frobenius of an Fp4 element. */
void QapFp4_p2pow(QapFp4_t c, QapFp4_t a) {
  Qapbinomial_Fp4_p2pow(c, a);
}

/* Compute the p-power Frobenius of an Fp4 element. */
void QapFp4_ppow(QapFp4_t c, QapFp4_t a) {
  Qapbinomial_Fp4_ppow(c, a);
}

/* Invert an Fp4 element. */
void QapFp4_inv (QapFp4_t c, QapFp4_t a) {
  Qaps2minusxi_Fp4_inv(c, a, TEMP_PASS_MEM4());
}

/* Multiply two Fp4 elements. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 */
void Qaps2minusxi_Fp4_mul(QapFp4_t c, QapFp4_t a, QapFp4_t b, TEMP_FORMAL) {
	QapFp2_t t0, t1, t2, t3, t4, t5;

  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
  TEMP_ALLOC_FP2(t2, 4);
  TEMP_ALLOC_FP2(t3, 6);
  TEMP_ALLOC_FP2(t4, 8);
  TEMP_ALLOC_FP2(t5, 10);

  QapFp2_add(t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp2_add(t1, b->a0, b->a1);  // t1 = b0 + b1
	QapFp2_mul(t2, t0, t1);			  // t2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + b0*a1 + a1*b1

	QapFp2_mul(t0, a->a0, b->a0);	// t0 = a0*b0
	QapFp2_mul(t1, a->a1, b->a1);	// t1 = a1*b1

	QapFp2_sub(t2, t2, t0);		    // t2 = t2 - t0 = a0*b1 + a1*b0 + a1*b1
  QapFp2_sub(c->a1, t2, t1);		    // t2 = t2 - t1 = a0*b1 + a1*b0
  QapFp2_mulxi(t1, t1);		      // t1 = a1*b1*xi
	QapFp2_add(c->a0, t0, t1);			// t0 = a0*b0 + a1*b1*xi
	
  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
  TEMP_FREE_FP2(t2);
  TEMP_FREE_FP2(t3);
  TEMP_FREE_FP2(t4);
  TEMP_FREE_FP2(t5);
}

/* Multiply two Fp4 elements using lazy reduction. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 */
void Qaps2minusxi_Fp4_mul_lazy(QapFp4_t c, QapFp4_t a, QapFp4_t b, TEMP_FORMAL) {
	QapFp2_t t0, t1;
  Fp2l_t T2, T3, T4;
  
  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
  TEMP_ALLOC_FP2L(T2, 4);
  TEMP_ALLOC_FP2L(T3, 8);
  TEMP_ALLOC_FP2L(T4, 12);  

  QapFp2_add_no_red (t0, a->a0, a->a1);	// t0 = a0 + a1
	QapFp2_add_no_red (t1, b->a0, b->a1);  // t1 = b0 + b1
	QapFp2_mul_no_red_o1 (T2, t0, t1, 2); // T2 = (a0 + a1)*(b0 + b1) = a0*b0 + a0*b1 + b0*a1 + a1*b1

	QapFp2_mul_no_red_o1 (T3, a->a0, b->a0, 2);	// T3 = a0*b0
	QapFp2_mul_no_red_o1 (T4, a->a1, b->a1, 2);	// T4 = a1*b1

	QapFp2l_sub_o2 (T2, T2, T3);		 // T2 = T2 - T3 = a0*b1 + a1*b0 + a1*b1
  QapFp2l_sub_o2 (T2, T2, T4);		 // T2 = T2 - T4 = a0*b1 + a1*b0
  QapFp2l_red (c->a1, T2);

  QapFp2l_mulxi(T4, T4);		            // T4 = a1*b1*xi
	QapFp2l_add_no_red (T2, T3, T4);			// T2 = a0*b0 + a1*b1*xi
  QapFp2l_red (c->a0, T2);

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
  TEMP_FREE_FP2L(T2);
  TEMP_FREE_FP2L(T3);
  TEMP_FREE_FP2L(T4);  
}

/* Square an Fp4 element. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 * This is Algorithm 9 in Beuchat et al., Pairing 2010.
 */
void Qaps2minusxi_Fp4_squ(QapFp4_t c, QapFp4_t a, TEMP_FORMAL) {
	QapFp2_t t0, t1;
  
  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
	
	QapFp2_squ (t0, a->a0);            // t0 = a0^2
  QapFp2_squ (t1, a->a1);            // t1 = a1^2
  QapFp2_add (c->a1, a->a0, a->a1);  // c1 = a0 + a1
  QapFp2_mulxi (c->a0, t1);          // c0 = t1*xi = a1^2*xi
  QapFp2_add (c->a0, c->a0, t0);     // c0 = a1^2*xi + a0^2
  QapFp2_squ (c->a1, c->a1);         // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  QapFp2_sub (c->a1, c->a1, t0);     // 2*a0*a1 + a1^2
  QapFp2_sub (c->a1, c->a1, t1);     // 2*a0*a1

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
}

/* Square an Fp4 element using lazy reduction. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 * This is Algorithm 9 in Beuchat et al., Pairing 2010.
 */
void Qaps2minusxi_Fp4_squ_lazy(QapFp4_t c, QapFp4_t a, TEMP_FORMAL) {
	QapFp2_t t1;
  Fp2l_t T0, T1, T2;

  TEMP_ALLOC_FP2(t1, 2);
  TEMP_ALLOC_FP2L(T0, 4);
  TEMP_ALLOC_FP2L(T1, 8);
  TEMP_ALLOC_FP2L(T2, 12);  

  QapFp2_squ_no_red (T0, a->a0);     // T0 = a0^2
  QapFp2_squ_no_red (T1, a->a1);     // T1 = a1^2
  QapFp2_add (t1, a->a0, a->a1);  // t1 = a0 + a1
  QapFp2l_mulxi (T2, T1);              // T2 = t1*xi = a1^2*xi
  QapFp2l_add_no_red (T2, T2, T0);     // T2 = a1^2*xi + a0^2
  QapFp2l_red (c->a0, T2);

  QapFp2_squ_no_red (T2, t1);       // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  QapFp2l_sub_o2 (T2, T2, T0);     // 2*a0*a1 + a1^2
  QapFp2l_sub_o2 (T2, T2, T1);     // 2*a0*a1
  QapFp2l_red (c->a1, T2);

  TEMP_FREE_FP2(t1);
  TEMP_FREE_FP2L(T0);
  TEMP_FREE_FP2L(T1);
  TEMP_FREE_FP2L(T2);  
}


/* Square an Fp4 element given by coefficients separately, return coefficients separately. */
void Qaps2minusxi_Fp4_squ_sep (QapFp2_t c0, QapFp2_t c1, QapFp2_t a0, QapFp2_t a1, TEMP_FORMAL) {
	QapFp2_t t0, t1;
  
  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
	
  QapFp2_squ (t0, a0);         // t0 = a0^2
  QapFp2_squ (t1, a1);         // t1 = a1^2
  QapFp2_add (c1, a0, a1);     // c1 = a0 + a1
  QapFp2_mulxi (c0, t1);       // c0 = t1*xi = a1^2*xi
  QapFp2_add (c0, c0, t0);     // c0 = a1^2*xi + a0^2
  QapFp2_squ (c1, c1);         // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  QapFp2_sub (c1, c1, t0);     // 2*a0*a1 + a1^2
  QapFp2_sub (c1, c1, t1);     // 2*a0*a1

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
}

/* Multiplication by the special element s in Fp4, with s^2 = xi in Fp2.
 * c = c0 + c1*s = (a0 + a1*s)*s = a1*xi + a0*s.
 */
void Qaps2minusxi_Fp4_muls(QapFp4_t c, QapFp4_t a, TEMP_FORMAL) {
  QapFp2_t t0;
  TEMP_ALLOC_FP2(t0, 0);

  QapFp2_mulxi(t0, a->a1);
	QapFp2_copy(c->a1, a->a0);
  QapFp2_copy(c->a0, t0);

  TEMP_FREE_FP2(t0);
}

/* Multiply an Fp4 element by the special element i. This function assumes that p = 3 (mod4)
 * and that the field extension is as described above.
 */
void Qapp3mod4_Fp4_muli(QapFp4_t c, QapFp4_t a) {
  QapFp2_muli (c->a0, a->a0);
  QapFp2_muli (c->a1, a->a1);
}

/* Multiply an Fp4 element by the special element xi.
 */
void Qapp3mod4_Fp4_mulxi(QapFp4_t c, QapFp4_t a) {
  QapFp2_mulxi (c->a0, a->a0);
  QapFp2_mulxi (c->a1, a->a1);
}

/* Inversion of an Fp4 element. This uses that Fp4 is constructed over Fp2 via s^2 - xi. 
 */
void Qaps2minusxi_Fp4_inv(QapFp4_t c, QapFp4_t a, TEMP_FORMAL) {
	QapFp2_t t0, t1;

  TEMP_ALLOC_FP2(t0, 0);
  TEMP_ALLOC_FP2(t1, 2);
   
	QapFp2_squ (t0, a->a0);         // t0 = a0^2
	QapFp2_squ (t1, a->a1);         // t1 = a1^2
	QapFp2_mulxi (t1, t1);          // t1 = xi*a1^2
  QapFp2_sub (t0, t0, t1);        // t0 = a0^2 - xi*a1^2 = N(a)

	QapFp2_inv (t0, t0);            // t0 = 1/N(a)

	QapFp2_mul (c->a0, a->a0, t0);  // c0 = a0/N(a)
	QapFp2_mul (c->a1, a->a1, t0);  
  QapFp2_neg (c->a1, c->a1);      // c1 = -a1/N(a)

  TEMP_FREE_FP2(t0);
  TEMP_FREE_FP2(t1);
}


/* Compute the p^2-power Frobenius of an Fp4 element. This function assumes that the field extension
 * is constructed via a binomial s^2 - xi, xi in Fp2.
 * This is the same as (complex) conjugation, i.e. QapFp4_p2pow(a0 + a1*s) = a0 - a1*s.
 * When we need this, it might be simpler to just negate the second coordinate in place.
 */
void Qapbinomial_Fp4_p2pow(QapFp4_t c, QapFp4_t a) {
	QapFp2_copy(c->a0, a->a0);
	QapFp2_neg(c->a1, a->a1);
}

/* Compute the p-power Frobenius of an Fp4 element. This function assumes that the field extension
 * is constructed via a binomial s^2 - xi, xi in Fp2.
 * This is QapFp4_ppow(a0 + a1*s) = a0^p + (a1^p*sppow)*s, where sppow = s^(p-1) in Fp2.
 */
void Qapbinomial_Fp4_ppow(QapFp4_t c, QapFp4_t a) {
	QapFp2_ppow (c->a0, a->a0);
	QapFp2_ppow (c->a1, a->a1);
  QapFp2_mul (c->a1, c->a1, Fp4_config.sppow);
}

int QapFp4_initialize_config (void) {
  if (QapFp2_init (Fp4_config.sppow) < 0) printf ("Fp4_config memory error.\n");

  return ERR_SUCCESS;
}

void QapFp4_free_config (void) {
  QapFp2_free (Fp4_config.sppow);
}


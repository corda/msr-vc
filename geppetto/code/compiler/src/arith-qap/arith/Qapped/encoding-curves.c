
#undef weierstrass_add
#undef weierstrass_dbl
#undef coordinate_eq
#undef wp_print

#if SIDENUM == SIDEL
#define weierstrass_add QapFp_weierstrass_affadd
#define weierstrass_dbl QapFp_weierstrass_affdbl
#define coordinate_eq QapFp_cmpeq
#define wp_print QapFp_print
#elif SIDENUM == SIDER
#define weierstrass_add QapFp2_weierstrass_affadd
#define weierstrass_dbl QapFp2_weierstrass_affdbl
#define coordinate_eq QapFp2_cmpeq
#define wp_print QapFp2_print
#else
#error Uhoh
#endif

#include "QapFp.h"
#include "..\..\..\..\input\pinocchio.h"


struct SIDE(s_EncodedElt) {
#if SIDENUM == SIDEL
	QapPoint_wp_t e;
#elif SIDENUM == SIDER
	QapPoint_wp2_t e;
#else
#error Uhoh
#endif
};
#define EncodedElt struct SIDE(s_EncodedElt)

#define true (1)
#define false (0)

//////////////////////////////////////////////////////////////////////////////
// declarations

void SIDE(_encoding_init)(EncodedElt* e);
void SIDE(_encoding_encode)(QapFp_t in, EncodedElt* out);
bool SIDE(_encoding_equal)(EncodedElt* a, EncodedElt* b);
void SIDE(_encoding_add)(EncodedElt* a, EncodedElt* b, EncodedElt* r);
void SIDE(_encoding_add_dd)(EncodedElt* a, EncodedElt* b, EncodedElt* r, bool definitely_distinct);
void SIDE(_encoding_sub)(EncodedElt* a, EncodedElt* b, EncodedElt* r);
int num_exponent_bits();
void SIDE(_encoding_mul)(EncodedElt* g, QapFp_t c, EncodedElt* r);
void SIDE(testOnCurve)(EncodedElt* a);
bool SIDE(isZero)(EncodedElt* elt);
bool SIDE(equal)(EncodedElt* a, EncodedElt* b);
void SIDE(set_zero)(EncodedElt* a);
void SIDE(copy)(EncodedElt* src, EncodedElt* dst);

//////////////////////////////////////////////////////////////////////////////
// definitions

// debug

void SIDE(enc_print_raw)(const char* msg, EncodedElt* p) {
  printf("%s:\n(", msg);
  wp_print(p->e->X);
  printf(",\n ");
  wp_print(p->e->Y);
  printf(")\n");
}

void SIDE(enc_print_rawraw)(EncodedElt* p) {
  wp_print(p->e->X);
  wp_print(p->e->Y);
}

void SIDE(enc_print)(const char* msg, EncodedElt *orig_p) {
	if (SIDE(isZero)(orig_p))
	{
		printf("%s: ZERO\n", msg);
	}
	else
	{
		SIDE(enc_print_raw)(msg, orig_p);
	}
}


void SIDE(_encoding_init)(EncodedElt* e)
{
#if SIDENUM == SIDEL
	QapFp_wp_init(e->e);
#elif SIDENUM == SIDER
	QapFp2_wp_init(e->e);
#else
#error ouch
#endif
}

EncodedElt SIDE(e_gen), SIDE(e_zero);

void SIDE(encoding_static_init)() {
	SIDE(_encoding_init)(&SIDE(e_gen));
	SIDE(_encoding_init)(&SIDE(e_zero));
#if SIDENUM == SIDEL
	QapFp_wp_copy(SIDE(e_gen).e, qap_ifc.Lg);
	QapFp_wp_copy(SIDE(e_zero).e, qap_ifc.Lzero);
#elif SIDENUM == SIDER
	QapFp2_wp_copy(SIDE(e_gen).e, qap_ifc.Rg);
	QapFp2_wp_copy(SIDE(e_zero).e, qap_ifc.Rzero);
#else
#error sad
#endif
}

// selecting between (possibly abstract) points; could be more modular
void SIDE(_encoding_select)(int bit, EncodedElt* a, EncodedElt* b, EncodedElt* r)
{
#if SIDENUM == SIDEL
	QapFp_select(r->e->X, bit, a->e->X, b->e->X);
	QapFp_select(r->e->Y, bit, a->e->Y, b->e->Y);
#elif SIDENUM == SIDER
	QapFp_select(r->e->X->a0, bit, a->e->X->a0, b->e->X->a0);
	QapFp_select(r->e->X->a1, bit, a->e->X->a1, b->e->X->a1);
	QapFp_select(r->e->Y->a0, bit, a->e->Y->a0, b->e->Y->a0);
	QapFp_select(r->e->Y->a1, bit, a->e->Y->a1, b->e->Y->a1);
#else
	never
#endif
}

void SIDE(_encoding_encode)(QapFp_t in, EncodedElt* out)
{
	SIDE(_encoding_mul)(&SIDE(e_gen), in, out);
}

bool SIDE(_encoding_equal)(EncodedElt* a, EncodedElt* b)
{
	return SIDE(equal)(a, b);
//#if SIDENUM == SIDEL
//	return QapFp_cmpeq(a, b);
//#elif SIDENUM == SIDER
//	return QapFp2_cmpeq(a, b);
//#else
//#error ouch
//#endif
}

void SIDE(copy)(EncodedElt* src, EncodedElt* dst)
{
	// Note argument order reversal as we cross interfaces.
#if SIDENUM == SIDEL
	QapFp_wp_copy(dst->e, src->e);
#elif SIDENUM == SIDER
	QapFp2_wp_copy(dst->e, src->e);
#else
#error ouch
#endif
}

void SIDE(_encoding_add)(EncodedElt* a, EncodedElt* b, EncodedElt* r)
{
	// dangerous with AFF ?
	SIDE(_encoding_add_dd)(a, b, r, 0);
}

// This variant of addition works on *all* cases for input variables
// but that's grotesque -- better to exclude those cases before calling, or to "patch" the divider with a test result.
void SIDE(_encoding_add_raw)(EncodedElt* a, EncodedElt* b, EncodedElt* r)
{
	bool a0 = SIDE(isZero)(a);
	bool b0 = SIDE(isZero)(b);

#ifndef MQAP
	// sanity check
	EncodedElt spec;
	SIDE(_encoding_init)(&spec);
	SIDE(_encoding_add)(a, b, &spec);
#endif
	EncodedElt a_nonzero;
	SIDE(_encoding_init)(&a_nonzero);
	SIDE(_encoding_select)(a0, &SIDE(e_gen), a, &a_nonzero);  // a' := (a == 0 ? 1 : a)

	EncodedElt b_nonzero;
	SIDE(_encoding_init)(&b_nonzero);
	SIDE(_encoding_select)(b0, &SIDE(e_gen), b, &b_nonzero); // b' := (b == 0 ? 1 : b)  

	EncodedElt r_dbl;
	SIDE(_encoding_init)(&r_dbl);
	weierstrass_dbl(r_dbl.e, a_nonzero.e); // r_dbl := 2a' 

	bool ab = SIDE(equal)(&a_nonzero, &b_nonzero);
	SIDE(_encoding_select)(ab, &r_dbl, &b_nonzero, &b_nonzero);	// if ab, then b' becomes 2a', and a' b' are distinct non-zero.

	EncodedElt r_add;
	SIDE(_encoding_init)(&r_add);
	weierstrass_add(r_add.e, a_nonzero.e, b_nonzero.e); // r_add := a' + b' 
	SIDE(_encoding_select)(ab, &r_dbl, &r_add, &r_add);
	SIDE(_encoding_select)(b0, a, &r_add, &r_add);
	SIDE(_encoding_select)(a0, b, &r_add, r);


#ifndef MQAP
	assert(SIDE(equal)(r, &spec));
#endif

}

void SIDE(_encoding_add_dd)(EncodedElt* a, EncodedElt* b, EncodedElt* r, bool definitely_distinct)
{
	SIDE(testOnCurve)(a);
	SIDE(testOnCurve)(b);

	if (SIDE(isZero)(a)) {
		// TODO this shouldn't work, because SIDE()==L will always be true (== is numerical),
		// so we should be on this branch even when r->e is a QapFp2, but we're passing it to QapFp_wp_copy(QapFp1). Broken.
		//CF ???
		SIDE(copy)(b, r);
	}
	else if (SIDE(isZero)(b)) {
		SIDE(copy)(a, r);
	}
	else if (!definitely_distinct && SIDE(equal)(a, b)) {
		// Library doesn't support addition of a point to itself
		weierstrass_dbl(r->e, a->e);
	}
	else {
		weierstrass_add(r->e, a->e, b->e);
	}
	SIDE(testOnCurve)(r);
}

void SIDE(_encoding_sub)(EncodedElt* a, EncodedElt* b, EncodedElt* r)
{
	EncodedElt minusb;
	SIDE(_encoding_init)(&minusb);
#if SIDENUM == SIDEL
	QapFp_wp_neg(minusb.e, b->e);
#elif SIDENUM == SIDER
	QapFp2_wp_neg(minusb.e, b->e);
#else
#error sad
#endif
//	enc_print("b:", b);
//	enc_print("minusb:", &minusb);
//	enc_print("a->e:", a);
	SIDE(_encoding_add)(a, &minusb, r);
//	enc_print("output:", r);
//	QapFp_wp_free(&minusb);	// TODO Leak
}

void SIDE(_encoding_sub_raw)(EncodedElt* a, EncodedElt* b, EncodedElt* r)
{
	EncodedElt minusb;
	SIDE(_encoding_init)(&minusb);
#if SIDENUM == SIDEL
	QapFp_wp_neg(minusb.e, b->e);
#elif SIDENUM == SIDER
	QapFp2_wp_neg(minusb.e, b->e);
#else
#error sad
#endif
	weierstrass_add(r->e, minusb.e, a->e);
}

void SIDE(_encoding_new_prod)();

#ifndef TINY
int SIDE(num_exponent_bits)() { return 254; }
#else
int SIDE(num_exponent_bits)() { return 53; }
#endif


//#define TICKER
int nRoot(void);
#ifdef TICKER
#define TICK log0 = log; log = nRoot() ; if (log != log0) { print(10000000 + log - log0); }
#else 
#define TICK {}
#endif


void SIDE(_encoding_mul)(EncodedElt* base, QapFp_t c, EncodedElt* r)
{
	// Naive square and multiply
	// SIDE(enc_print_rawraw)(&SIDE(e_gen).e);
	// printf("encoding_mul:\n"); QapFp_print(c); SIDE(enc_print_raw)(" * ", base->e);
	SIDE(testOnCurve)(base);

	int log;
	log = nRoot();

	EncodedElt result;
	SIDE(_encoding_init)(&result);

#if SELECTMUL
	// we'd like to exclude those cases, but they are needed with input keys
	// printf("Commitments to zero not supported in this mode.\n");
	zeroAssert(elem_eq_zero(c->e));

	// instead of excluding this one, we select after the computation
	// zeroAssert(SIDE(isZero)(base));
	//  if (SIDE(isZero)(base)) SIDE(set_zero)(&result); 
	//  else
#else
	if (SIDE(isZero)(base) || elem_eq_zero(c->e)) {
		printf("excluded cases %d %d\n", SIDE(isZero)(base), elem_eq_zero(c->e));
		SIDE(set_zero)(&result);
	}
	else
#endif
	{
		EncodedElt sum;
		SIDE(_encoding_init)(&sum);
		// we start with 1*g instead of 0 to avoid special cases of add/dbl
		// and we finally subtract 2^num_exponent_bits()*g 
		SIDE(copy)(&SIDE(e_gen), &result);

		// TODO: This should only be done once at compile time
		EncodedElt g254;
		SIDE(_encoding_init)(&g254);
		SIDE(copy)(&SIDE(e_gen), &g254);
		for (int bit_index = SIDE(num_exponent_bits)(); bit_index >= 0; bit_index--) weierstrass_dbl(g254.e, g254.e);

		for (int bit_index = SIDE(num_exponent_bits)(); bit_index >= 1; bit_index--) {
			weierstrass_dbl(result.e, result.e);
			int bitset = elem_get_bit(c->e, bit_index);

			weierstrass_add(sum.e, base->e, result.e);

			// result := bitset ? sum : result 
			SIDE(_encoding_select)(bitset, &sum, &result, &result);
			// printf(" index %3d, bitset %1d", bit_index, bitset);
		}
		// SIDE(enc_print_raw)("before last bit: ", result.e);

		// Need to be careful with the last dbl and add, since we may have overflowed the order of the group 
		// and arrived at 0 or the inverse of base
#if SELECTMUL
		// these two conditions could be turned into selects.
		zeroAssert(SIDE(isZero)(&result));
#else
		if (!SIDE(isZero)(&result))
#endif
		{
			weierstrass_dbl(result.e, result.e);
		}
#if SELECTMUL
		zeroAssert(coordinate_eq(result.e->X, base->e->X));
#else
		if (coordinate_eq(result.e->X, base->e->X) != 0)
#endif
		{  // Not an inverse, so safe to add
			int bitset = elem_get_bit(c->e, 0);
			// SIDE(enc_print_raw)("base: ", base->e);
			// SIDE(enc_print_raw)("result: ", result.e);
			weierstrass_add(sum.e, base->e, result.e);
			// SIDE(enc_print_raw)("during last bit: ", sum.e);
			SIDE(_encoding_select)(bitset, &sum, &result, &result);
		}
		// SIDE(enc_print_raw)("after last bit: ", result.e);

		// Subtract off 2^254
		SIDE(_encoding_sub_raw)(&result, &g254, &result);
		// SIDE(enc_print_raw)("after subtraction: ", result.e);
	} // ! special cases

	// was, without select
	// SIDE(copy)(&result, r);

	// with select:
	EncodedElt zero;
	SIDE(_encoding_init)(&zero);
	SIDE(set_zero)(&zero);
	SIDE(_encoding_select)(SIDE(isZero)(base), &zero, &result, r);

	SIDE(testOnCurve)(r);
	// for comparative debugging:
	// SIDE(enc_print_rawraw)(&result);
}

void SIDE(testOnCurve)(EncodedElt* a)
{
#ifdef _DEBUG
  #if SIDENUM == SIDEL
  QapFp_weierstrass_oncurve_aff(a->e);
#elif SIDENUM == SIDER
  QapFp2_weierstrass_oncurve_aff(a->e);
#else
#error Uhoh
#endif
#endif
}

bool SIDE(isZero)(EncodedElt* elt)
{
	// assert(elt != 0);
	SIDE(testOnCurve)(elt);
	bool ret;

#if SIDENUM == SIDEL
	ret = QapFp_cmpeq(elt->e->Y, SIDE(e_zero).e->Y);
#elif SIDENUM == SIDER
	ret = QapFp2_cmpeq(elt->e->Y, SIDE(e_zero).e->Y);
#else
#error Uhoh
#endif

	return ret;
}

bool SIDE(equal)(EncodedElt* a, EncodedElt* b)
{
	assert(a != 0 && b != 0);
	SIDE(testOnCurve)(a);
	SIDE(testOnCurve)(b);
	
	bool ret = false;

#if SIDENUM == SIDEL
	ret = QapFp_cmpeq(a->e->X, b->e->X) * QapFp_cmpeq(a->e->Y, b->e->Y);
	//ret = Fp_weierstrass_equal_aff(aff_a, aff_b);

#elif SIDENUM == SIDER
	ret = QapFp2_cmpeq(a->e->X, b->e->X) * QapFp2_cmpeq(a->e->Y, b->e->Y);
#else
#error Uhoh
#endif

#if 0 // only when non-affine? 
	// we need an abstract select. 
	// return a0 ? b0 : (b0? false : ret);
	return select(a0, b0, select(b0, false, ret));
#else
	return ret;
#endif
}

void SIDE(set_zero)(EncodedElt* a) { 
#if SIDENUM == SIDEL
	QapFp_wp_copy(a->e, SIDE(e_zero).e);	// Use our precomputed value	
#elif SIDENUM == SIDER
	QapFp2_wp_copy(a->e, SIDE(e_zero).e);	// Use our precomputed value	
#else
#error Uhoh
#endif
}

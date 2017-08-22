#include "PrimitiveIfc.h"
#include "Qaparith_error.h"
#include "QapFp.h"

#if MQAP
#else
//#include <assert.h>
#endif

QapFp_config_t QapFp_config;

int QapFp_init (QapFp_t a) {
  elem_init(a->e);
  return ERR_SUCCESS;
}

int QapFpl_init (QapFpl_t a) {
  elem_init(a->e);
  return ERR_SUCCESS;
}

void QapFp_free (QapFp_t a) {
	elem_free(a->e);
}

void QapFp_free_config (void) {
}

void QapFp_modinv (QapFp_t dst, QapFp_t src) {
	elem_modinv(dst->e, src->e);
}

void QapFp_mul (QapFp_t c, QapFp_t a, QapFp_t b) {
	elem_mul(c->e, a->e, b->e);
}

void QapFp_add (QapFp_t c, QapFp_t a, QapFp_t b) {
	elem_add(c->e, a->e, b->e);
}

// QAP-specific.
// the first variant is more efficient (but fails unless we use temp variables!)
void QapFp_select(QapFp_t c, int bit, QapFp_t a, QapFp_t b) {
	QapFp_t bitElt;
	QapFp_init(bitElt);
	QapFp_set_ui(bitElt, bit);
	QapFp_t c0, c1;
	QapFp_init(c0);
	QapFp_init(c1);
#if 1
	QapFp_sub(c0, a, b);
	QapFp_mul(c1, c0, bitElt);
	QapFp_add(c, c1, b);
#else
	QapFp_mul(c1, a, bitElt);
	QapFp_set_ui(bitElt, 1 - bit);
	QapFp_mul(c0, b, bitElt);
	QapFp_add(c, c1, c0);
#endif
//	QapFp_print(c);
}

/* Addition without reduction. No check for overflow! */
void QapFp_add_no_red (QapFp_t c, QapFp_t a, QapFp_t b) {
	QapFp_add(c,a,b);
}

/* Addition without reduction on long elements. */
void QapFpl_add_no_red (QapFpl_t c, QapFpl_t a, QapFpl_t b) {
	QapFp_add(c,a,b);
}

void QapFp_sub (QapFp_t c, QapFp_t a, QapFp_t b) {
	elem_sub(c->e, a->e, b->e);
}

/* Addition without reduction. No check for overflow! */
void QapFp_sub_no_red (QapFp_t c, QapFp_t a, QapFp_t b) {
	QapFp_sub(c,a,b);
}

void QapFp_div2 (QapFp_t c, QapFp_t a) {
	elem_div2(c->e, a->e);
}

void QapFpl_sub_no_red (QapFpl_t c, QapFpl_t a, QapFpl_t b) {
  QapFp_sub(c, a, b);
}

void QapFpl_sub_o1_ct (QapFpl_t c, QapFpl_t a, QapFpl_t b, int h, void *mem) {
  QapFp_sub(c, a, b);
}

void QapFpl_sub_o2_ct (QapFpl_t c, QapFpl_t a, QapFpl_t b, void *mem) {
  QapFp_sub(c, a, b);
}

void QapFpl_add_o2_ct (QapFpl_t c, QapFpl_t a, QapFpl_t b, void *mem) {
  QapFp_add(c,a,b);
}

int QapFp_cmpeq (QapFp_t a, QapFp_t b) {
	QapFp_t c;
	QapFp_init(c);
	QapFp_sub(c, a, b);
	return elem_eq_zero(c->e);
}
int QapFp_cmp (QapFp_t a, QapFp_t b) {
	return elem_cmp(a->e, b->e);
}

void QapFp_mul_no_red (QapFpl_t c, QapFp_t a, QapFp_t b) {
	QapFp_mul(c,a,b);
}

void QapFpl_red (QapFp_t c, QapFpl_t a) {
	// reduced is how we roll
	QapFp_copy(c,a);
}

void QapFp_mul3 (QapFp_t c, QapFp_t a) {
	elem_mul3(c->e, a->e);
}

void QapFp_initialize_config () {
	elem_static_init();
}

int QapFp_get_bit(QapFp_t a, int bit) {
  return elem_get_bit(a->e, bit);
}

void QapFp_print (QapFp_t a) {
	elem_dbg_print(a->e);
}

void QapFp_set(QapFp_t c, QapFp_t a)
{
	QapFp_copy(c,a);
}

void QapFp_set_ui (QapFp_t c, uint64_t x)
{
	elem_set_ui(c->e, x);
}

void QapFp_set_u32 (QapFp_t c,
	uint32_t a7, uint32_t a6, uint32_t a5, uint32_t a4,
	uint32_t a3, uint32_t a2, uint32_t a1, uint32_t a0)
{
	elem_set_u32(c->e, a7, a6, a5, a4, a3, a2, a1, a0);
}

void QapFp_copy (QapFp_t c, QapFp_t a) {
	elem_copy(c->e, a->e);
}

void QapFpl_copy (QapFpl_t c, QapFpl_t a) {
	QapFp_copy(c,a);
}

void QapFp_neg (QapFp_t c, QapFp_t a) {
	elem_neg(c->e, a->e);
}

void QapFp_rand (QapFp_t c) {
	elem_rand(c->e);
}

void QapDbgAssert(int condition)
{
	elem_assert(condition);
}

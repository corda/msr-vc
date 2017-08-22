#include "QapFp2_weierstrass.h"

typedef struct {
  QapFp2_t b;
} curve_wp2_t;

typedef curve_wp2_t Curve_wp2_t[1];


typedef struct {
  QapFp2_t t0;
  QapFp2_t t1;
  QapFp2_t t2;
  QapFp2_t t3;
  QapFp2_t t4;
  QapFp2_t t5;
  QapFp2_t t6;
  QapFp2_t t7;
  QapFp2_t t8;
} aux_wp2_t;

typedef aux_wp2_t Aux_wp2_t[1];

static Aux_wp2_t Fp2_w_aux;
static Curve_wp2_t Fp2_w_curve;


#define POINT QapPoint_wp2_t
#define MEM   Aux_wp2_t
#define CURVE Curve_wp2_t

// TODO: Fix Mul2, Mul4, Mul8
#define MUL(c,a,b) QapFp2_mul (c,a,b)
#define MUL2(c,a)  QapFp2_add (c,a,a)
#define MUL3(c,a)  QapFp2_mul3 (c,a)
#define MUL4(c,a)  MUL2(c,a); MUL2(c,c)
#define MUL8(c,a)  MUL4(c,a); MUL2 (c,c)
#define SQR(c,a)   QapFp2_squ (c,a)
#define ADD(c,a,b) QapFp2_add (c,a,b)
#define SUB(c,a,b) QapFp2_sub (c,a,b)
#define INV(c,a)   QapFp2_inv (c,a)
#define SET(c,a)   QapFp2_copy (c,a)
#define SETONE(c)  QapFp2_set_ui (c, 1)
#define CMPEQ(a,b) QapFp2_cmpeq(a,b)
#define VARIANT(n)	w2##n

#include "Qapweierstrass.h"

// undo the defines so we can be compiled in the same unit as other users of Qapweierstrass
#undef POINT
#undef MEM
#undef CURVE
#undef MUL
#undef MUL2
#undef MUL3
#undef MUL4
#undef MUL8
#undef SQR
#undef ADD
#undef SUB
#undef INV
#undef SET
#undef SETONE
#undef CMPEQ
#undef VARIANT

int QapFp2_wp_init (QapPoint_wp2_t p) {
  int ret;

  if ((ret=QapFp2_init (p->X)) < 0) return ret;
  if ((ret=QapFp2_init (p->Y)) < 0) return ret;
  
  return ERR_SUCCESS;
}

void QapFp2_wp_free (QapPoint_wp2_t p) {
  QapFp2_free (p->X);
  QapFp2_free (p->Y);
}


int QapFp2_initialize_weierstrass (QapFp2_t b) {
  int ret;

  if ((ret=QapFp2_init (Fp2_w_aux->t0)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t1)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t2)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t3)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t4)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t5)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t6)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t7)) < 0) return ret;
  if ((ret=QapFp2_init (Fp2_w_aux->t8)) < 0) return ret;
  
  if ((ret=QapFp2_init (Fp2_w_curve->b)) < 0) return ret;
  QapFp2_copy (Fp2_w_curve->b, b);
  return ERR_SUCCESS;
}

void QapFp2_free_weierstrass (void) {
  QapFp2_free (Fp2_w_aux->t0);
  QapFp2_free (Fp2_w_aux->t1);
  QapFp2_free (Fp2_w_aux->t2);
  QapFp2_free (Fp2_w_aux->t3);
  QapFp2_free (Fp2_w_aux->t4);
  QapFp2_free (Fp2_w_aux->t5);
  QapFp2_free (Fp2_w_aux->t6);
  QapFp2_free (Fp2_w_aux->t7);
  QapFp2_free (Fp2_w_aux->t8);

  QapFp2_free (Fp2_w_curve->b);
}

void QapFp2_wp_copy (QapPoint_wp2_t c, QapPoint_wp2_t a) {
  QapFp2_copy (c->X, a->X);
  QapFp2_copy (c->Y, a->Y);
}

void QapFp2_wp_neg (QapPoint_wp2_t c, QapPoint_wp2_t a) {
  QapFp2_copy (c->X, a->X);
  QapFp2_neg (c->Y, a->Y);
}

void QapFp2_weierstrass_affdbl(QapPoint_wp2_t c, QapPoint_wp2_t a) {
  w2weierstrass_affdbl_template(c, a, Fp2_w_aux);
}

void QapFp2_weierstrass_affadd(QapPoint_wp2_t c, QapPoint_wp2_t a, QapPoint_wp2_t b) {
  w2weierstrass_affadd_template(c, a, b, Fp2_w_aux);
}

int QapFp2_weierstrass_oncurve_aff (QapPoint_wp2_t x) {
  return w2weierstrass_oncurve_aff_template (x, Fp2_w_aux, Fp2_w_curve);
}

int QapFp2_weierstrass_equal_aff (QapPoint_wp2_t x, QapPoint_wp2_t y) {
  return w2weierstrass_equal_aff_template (x, y, Fp2_w_curve);
}

#include "QapFp_weierstrass.h"

typedef struct {
  QapFp_t B[2];
  QapFp_t b;
} curve_wp_t;

typedef curve_wp_t Curve_wp_t[1];

typedef struct {
  QapFp_t t0;
  QapFp_t t1;
  QapFp_t t2;
  QapFp_t t3;
  QapFp_t t4;
  QapFp_t t5;
  QapFp_t t6;
  QapFp_t t7;
  QapFp_t t8;
} aux_wp_t;

typedef aux_wp_t Aux_wp_t[1];

static Aux_wp_t Fp_w_aux;
static Curve_wp_t Fp_w_curve;


#define POINT QapPoint_wp_t
#define MEM   Aux_wp_t
#define CURVE Curve_wp_t

// TODO: Fix Mul2, Mul4, Mul8
#define MUL(c,a,b) QapFp_mul (c,a,b)
#define MUL2(c,a)  QapFp_add (c,a,a)
#define MUL3(c,a)  QapFp_mul3 (c,a)
#define MUL4(c,a)  MUL2(c,a); MUL2(c,c)
#define MUL8(c,a)  MUL4(c,a); MUL2 (c,c)
#define SQR(c,a)   QapFp_sqr (c,a)
#define ADD(c,a,b) QapFp_add (c,a,b)
#define SUB(c,a,b) QapFp_sub (c,a,b)
#define INV(c,a)   QapFp_modinv (c,a)
#define SET(c,a)   QapFp_copy (c,a)
#define SETONE(c)  QapFp_set_ui(c, 1)
#define CMPEQ(a,b) (QapFp_cmpeq(a,b))
#define VARIANT(n)	w1##n

#include "Qaparith_error.h"
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

int QapFp_wp_init (QapPoint_wp_t p) {
  int ret;

  if ((ret=QapFp_init (p->X)) < 0) return ret;
  if ((ret=QapFp_init (p->Y)) < 0) return ret;
  
  return ERR_SUCCESS;
}

void QapFp_wp_free (QapPoint_wp_t p) {
  QapFp_free (p->X);
  QapFp_free (p->Y);
}

void QapFp_wp_copy (QapPoint_wp_t c, QapPoint_wp_t a) {
  QapFp_copy (c->X, a->X);
  QapFp_copy (c->Y, a->Y);
}


void QapFp_wp_neg (QapPoint_wp_t c, QapPoint_wp_t a) {
  QapFp_copy (c->X, a->X);
  QapFp_neg (c->Y, a->Y);
}

int QapFp_initialize_weierstrass (QapFp_t *b, int num) {
  int ret;

  if ((ret=QapFp_init (Fp_w_aux->t0)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t1)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t2)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t3)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t4)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t5)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t6)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t7)) < 0) return ret;
  if ((ret=QapFp_init (Fp_w_aux->t8)) < 0) return ret;
   
   if ((ret=QapFp_init (Fp_w_curve->B[0])) < 0) return ret;
   QapFp_copy (Fp_w_curve->B[0], b[0]);
   if ((ret=QapFp_init (Fp_w_curve->b)) < 0) return ret;
   QapFp_copy (Fp_w_curve->b, b[0]);
   if (num == 2) {
     if ((ret=QapFp_init (Fp_w_curve->B[1])) < 0) return ret;
     QapFp_copy (Fp_w_curve->B[1], b[1]);
   } else if (num > 2)
   {
     printf ("ERROR: unsupported number of curves!"); /*getchar (); */
   }

  return ERR_SUCCESS;
}

void QapFp_free_weierstrass (void) {
  QapFp_free (Fp_w_aux->t0);
  QapFp_free (Fp_w_aux->t1);
  QapFp_free (Fp_w_aux->t2);
  QapFp_free (Fp_w_aux->t3);
  QapFp_free (Fp_w_aux->t4);
  QapFp_free (Fp_w_aux->t5);
  QapFp_free (Fp_w_aux->t6);
  QapFp_free (Fp_w_aux->t7);
  QapFp_free (Fp_w_aux->t8);

  QapFp_free (Fp_w_curve->B[0]);
  QapFp_free (Fp_w_curve->B[1]);
  QapFp_free (Fp_w_curve->b);
}

void QapFp_weierstrass_affdbl(QapPoint_wp_t c, QapPoint_wp_t a) {
  w1weierstrass_affdbl_template(c, a, Fp_w_aux);  
}

void QapFp_weierstrass_affadd(QapPoint_wp_t c, QapPoint_wp_t a, QapPoint_wp_t b) {
  w1weierstrass_affadd_template(c, a, b, Fp_w_aux);
}

int QapFp_weierstrass_oncurve_aff (QapPoint_wp_t x) {
  return w1weierstrass_oncurve_aff_template (x, Fp_w_aux, Fp_w_curve);
}

int QapFp_weierstrass_oncurve_aff_select_b (QapPoint_wp_t x, int num) {
  QapFp_copy (Fp_w_curve->b, Fp_w_curve->B[num]);
  return w1weierstrass_oncurve_aff_template (x, Fp_w_aux, Fp_w_curve);
}

int QapFp_weierstrass_equal_aff(QapPoint_wp_t x, QapPoint_wp_t y) {
  return w1weierstrass_equal_aff_template(x, y, Fp_w_curve);
}
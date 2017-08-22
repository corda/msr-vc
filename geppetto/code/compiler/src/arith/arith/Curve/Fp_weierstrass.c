/* curves/Fp_weierstrass.c
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

#include "Fp_weierstrass.h"

typedef struct {
  Fp_t B[2];
  Fp_t b;
} curve_wp_t;

typedef curve_wp_t Curve_wp_t[1];

typedef struct {
  Fp_t t0;
  Fp_t t1;
  Fp_t t2;
  Fp_t t3;
  Fp_t t4;
  Fp_t t5;
  Fp_t t6;
  Fp_t t7;
  Fp_t t8;
} aux_wp_t;

typedef aux_wp_t Aux_wp_t[1];

static Aux_wp_t Fp_w_aux;
static Curve_wp_t Fp_w_curve;


#define POINT Point_wp_t
#define POINTAFF Point_wpaff_t
#define MEM   Aux_wp_t
#define CURVE Curve_wp_t

// TODO: Fix Mul2, Mul4, Mul8
#define MUL(c,a,b) Fp_mul (c,a,b)
#define MUL2(c,a)  Fp_add (c,a,a)
#define MUL3(c,a)  Fp_mul3 (c,a)
#define MUL4(c,a)  MUL2(c,a); MUL2(c,c)
#define MUL8(c,a)  MUL4(c,a); MUL2 (c,c)
#define SQR(c,a)   Fp_sqr (c,a)
#define ADD(c,a,b) Fp_add (c,a,b)
#define SUB(c,a,b) Fp_sub (c,a,b)
#define INV(c,a)   Fp_modinv (c,a)
#define SET(c,a)   Fp_copy (c,a)
#define SETONE(c)  Fp_set_ui (c, 1)
#define CMPEQ(a,b) (Fp_cmp(a,b) == 0)

#include "weierstrass.h"

int Fp_wp_init (Point_wp_t p) {
  int ret;

  if ((ret=Fp_init (p->X)) < 0) return ret;
  if ((ret=Fp_init (p->Y)) < 0) return ret;
  if ((ret=Fp_init (p->Z)) < 0) return ret;
  
  return ERR_SUCCESS;
}

int Fp_wpaff_init(Point_wpaff_t p) {
	int ret;

	if ((ret = Fp_init(p->X)) < 0) return ret;
	if ((ret = Fp_init(p->Y)) < 0) return ret;
	
	return ERR_SUCCESS;
}

void Fp_wp_free (Point_wp_t p) {
  Fp_free (p->X);
  Fp_free (p->Y);
  Fp_free (p->Z);
}

void Fp_wpaff_free(Point_wpaff_t p) {
	Fp_free(p->X);
	Fp_free(p->Y);
}


void Fp_wp_copy (Point_wp_t c, Point_wp_t a) {
  Fp_copy (c->X, a->X);
  Fp_copy (c->Y, a->Y);
  Fp_copy (c->Z, a->Z);
}

void Fp_wpaff_copy(Point_wpaff_t c, Point_wpaff_t a) {
	Fp_copy(c->X, a->X);
	Fp_copy(c->Y, a->Y);
}

void Fp_wp_neg (Point_wp_t c, Point_wp_t a) {
  Fp_copy (c->X, a->X);
  Fp_neg (c->Y, a->Y);
  Fp_copy (c->Z, a->Z);
}

void Fp_wpaff_neg(Point_wpaff_t c, Point_wpaff_t a) {
	Fp_copy(c->X, a->X);
	Fp_neg(c->Y, a->Y);
}

int Fp_initialize_weierstrass (Fp_t *b, int num) {
  int ret;

  if ((ret=Fp_init (Fp_w_aux->t0)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t1)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t2)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t3)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t4)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t5)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t6)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t7)) < 0) return ret;
  if ((ret=Fp_init (Fp_w_aux->t8)) < 0) return ret;
   
   if ((ret=Fp_init (Fp_w_curve->B[0])) < 0) return ret;
   Fp_copy (Fp_w_curve->B[0], b[0]);
   if ((ret=Fp_init (Fp_w_curve->b)) < 0) return ret;
   Fp_copy (Fp_w_curve->b, b[0]);
   if (num == 2) {
     if ((ret=Fp_init (Fp_w_curve->B[1])) < 0) return ret;
     Fp_copy (Fp_w_curve->B[1], b[1]);
   } else if (num > 2)
   {
     printf ("ERROR: unsupported number of curves!"); getchar ();
   }

  return ERR_SUCCESS;
}

void Fp_free_weierstrass (void) {
  Fp_free (Fp_w_aux->t0);
  Fp_free (Fp_w_aux->t1);
  Fp_free (Fp_w_aux->t2);
  Fp_free (Fp_w_aux->t3);
  Fp_free (Fp_w_aux->t4);
  Fp_free (Fp_w_aux->t5);
  Fp_free (Fp_w_aux->t6);
  Fp_free (Fp_w_aux->t7);
  Fp_free (Fp_w_aux->t8);

  Fp_free (Fp_w_curve->B[0]);
  Fp_free (Fp_w_curve->B[1]);
  Fp_free (Fp_w_curve->b);
}

void Fp_weierstrass_dbl (Point_wp_t c, Point_wp_t a) {
  weierstrass_dbl_template (c, a, Fp_w_aux);
}

void Fp_weierstrass_add (Point_wp_t c, Point_wp_t a, Point_wp_t b) {
  weierstrass_add_template (c, a, b, Fp_w_aux);
}

void Fp_weierstrass_mix (Point_wp_t c, Point_wp_t a, Point_wp_t b) {
 weierstrass_mix_template (c, a, b, Fp_w_aux);
}

void Fp_weierstrass_aff (Point_wp_t c, Point_wp_t a, Point_wp_t b) {
  weierstrass_aff_template (c, a, b, Fp_w_aux);
}

void Fp_weierstrass_affdbl(Point_wp_t c, Point_wp_t a) {
	weierstrass_affdbl_template(c, a, Fp_w_aux);
}

void Fp_weierstrass_affadd(Point_wp_t c, Point_wp_t a, Point_wp_t b) {
	weierstrass_affadd_template(c, a, b, Fp_w_aux);
}

int Fp_weierstrass_oncurve_jacproj (Point_wp_t x) {
  return weierstrass_oncurve_jacproj_template (x, Fp_w_aux, Fp_w_curve);
}

int Fp_weierstrass_oncurve_aff (Point_wp_t x) {
  return weierstrass_oncurve_aff_template (x, Fp_w_aux, Fp_w_curve);
}

int Fp_weierstrass_oncurve_jacproj_select_b (Point_wp_t x, int num) {
  Fp_copy (Fp_w_curve->b, Fp_w_curve->B[num]);
  return weierstrass_oncurve_jacproj_template (x, Fp_w_aux, Fp_w_curve);
}

int Fp_weierstrass_oncurve_aff_select_b (Point_wp_t x, int num) {
  Fp_copy (Fp_w_curve->b, Fp_w_curve->B[num]);
  return weierstrass_oncurve_aff_template (x, Fp_w_aux, Fp_w_curve);
}

void Fp_weierstrass_aff2proj (Point_wp_t c, Point_wp_t a) {
  weierstrass_aff2proj_template (c, a);
}

void Fp_weierstrass_jacproj2aff (Point_wp_t c, Point_wp_t a) {
  weierstrass_jacproj2aff_template (c, a, Fp_w_aux);
}

int Fp_weierstrass_equal_aff(Point_wp_t a, Point_wp_t b) {
	return weierstrass_equal_aff_template(a, b);
}

int Fp_weierstrass_equal_jacproj(Point_wp_t a, Point_wp_t b) {
	return weierstrass_equal_jacproj_template(a, b, Fp_w_aux);
}

/*
 * Test vectors:
 * p = 23011393939926379162305986076771341882121011402274813698163587504247481406579
 *   = 0x32DFFCC760EFFAE264E023891D54BB6F6C440805F58068C31478B893C9159473
 * x = 16799001282946594873539054437307522853065069661512796544241769834773293970758
 *   = 0x2523E5D4D7533E02FA05B3E58ACD569202951D998AEFD912CF2C2C094B123546
 * y = 20703541998597624734657103243505617270365165626643744749675339538629924006806
 *   = 0x2DC5CA1D38DDF17AE3268541AD518D0BCEEDFD4DA8D6B033CFFEB8DD96F0AB96
 * b = 19630023706425094729206547151154687245152993473663944785960441199559622259465
 *   = 0x2B66331EB33AF7F7DAFF5FAAE325C96D7CA737E7409DF1B694D849B599F7F309
 */

#if 0
int main () {
  uint64_t p[4], b[4];
  Fp_t B;
  Point_wp_t P, Q;

  p[3] = 0x32DFFCC760EFFAE2;
  p[2] = 0x64E023891D54BB6F;
  p[1] = 0x6C440805F58068C3;
  p[0] = 0x1478B893C9159473;
  
  b[3] = 0x2B66331EB33AF7F7;
  b[2] = 0xDAFF5FAAE325C96D;
  b[1] = 0x7CA737E7409DF1B6;
  b[0] = 0x94D849B599F7F309;

  Fp_initialize_config (p, 4);
  
  Fp_init (B);
  Fp_set (B, b);
  Fp_initialize_weierstrass (&B, 1);

  Fp_wp_init (P);
  Fp_wp_init (Q);

  b[3] = 0x2523E5D4D7533E02;
  b[2] = 0xFA05B3E58ACD5692;
  b[1] = 0x02951D998AEFD912;
  b[0] = 0xCF2C2C094B123546;
  Fp_set (P->X, b);

  b[3] = 0x2DC5CA1D38DDF17A;
  b[2] = 0xE3268541AD518D0B;
  b[1] = 0xCEEDFD4DA8D6B033;
  b[0] = 0xCFFEB8DD96F0AB96;
  Fp_set (P->Y, b);

  printf ("Affine point on curve?: %d\n", Fp_weierstrass_oncurve_aff (P));
  
  Fp_weierstrass_aff2proj (P, P);
  printf ("Projective point on curve?: %d\n", Fp_weierstrass_oncurve_jacproj (P));

  Fp_weierstrass_dbl (Q, P);
  Fp_weierstrass_add (Q, Q, P);
  Fp_weierstrass_add (Q, Q, P);
  Fp_weierstrass_add (Q, Q, P);
  Fp_weierstrass_dbl (Q, Q);

  printf ("Projective point on curve?: %d\n", Fp_weierstrass_oncurve_jacproj (Q));

  Fp_weierstrass_jacproj2aff (P, Q);
  
  printf ("Affine point on curve?: %d\n", Fp_weierstrass_oncurve_aff (P));

  Fp_wp_free (P);
  Fp_wp_free (Q);

  Fp_free_weierstrass ();
  Fp_free_config ();
  getchar ();
}
#endif

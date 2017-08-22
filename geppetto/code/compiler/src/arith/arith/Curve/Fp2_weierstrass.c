/* curves/Fp2_weierstrass.c
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

#include "Fp2_weierstrass.h"

typedef struct {
  Fp2_t b;
} curve_wp2_t;

typedef curve_wp2_t Curve_wp2_t[1];


typedef struct {
  Fp2_t t0;
  Fp2_t t1;
  Fp2_t t2;
  Fp2_t t3;
  Fp2_t t4;
  Fp2_t t5;
  Fp2_t t6;
  Fp2_t t7;
  Fp2_t t8;
} aux_wp2_t;

typedef aux_wp2_t Aux_wp2_t[1];

static Aux_wp2_t Fp2_w_aux;
static Curve_wp2_t Fp2_w_curve;


#define POINT Point_wp2_t
#define POINTAFF Point_wp2aff_t
#define MEM   Aux_wp2_t
#define CURVE Curve_wp2_t

// TODO: Fix Mul2, Mul4, Mul8
#define MUL(c,a,b) Fp2_mul (c,a,b)
#define MUL2(c,a)  Fp2_add (c,a,a)
#define MUL3(c,a)  Fp2_mul3 (c,a)
#define MUL4(c,a)  MUL2(c,a); MUL2(c,c)
#define MUL8(c,a)  MUL4(c,a); MUL2 (c,c)
#define SQR(c,a)   Fp2_squ (c,a)
#define ADD(c,a,b) Fp2_add (c,a,b)
#define SUB(c,a,b) Fp2_sub (c,a,b)
#define INV(c,a)   Fp2_inv (c,a)
#define SET(c,a)   Fp2_copy (c,a)
#define SETONE(c)  Fp2_set_ui (c, 1)
#define CMPEQ(a,b) Fp2_cmpeq(a,b)

#include "weierstrass.h"

int Fp2_wp_init (Point_wp2_t p) {
  int ret;

  if ((ret=Fp2_init (p->X)) < 0) return ret;
  if ((ret=Fp2_init (p->Y)) < 0) return ret;
  if ((ret=Fp2_init (p->Z)) < 0) return ret;
  
  return ERR_SUCCESS;
}

int Fp2_wpaff_init(Point_wp2aff_t p) {
	int ret;

	if ((ret = Fp2_init(p->X)) < 0) return ret;
	if ((ret = Fp2_init(p->Y)) < 0) return ret;

	return ERR_SUCCESS;
}

void Fp2_wp_free (Point_wp2_t p) {
  Fp2_free (p->X);
  Fp2_free (p->Y);
  Fp2_free (p->Z);
}

void Fp2_wpaff_free(Point_wp2aff_t p) {
	Fp2_free(p->X);
	Fp2_free(p->Y);
}


int Fp2_initialize_weierstrass (Fp2_t b) {
  int ret;

  if ((ret=Fp2_init (Fp2_w_aux->t0)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t1)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t2)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t3)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t4)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t5)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t6)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t7)) < 0) return ret;
  if ((ret=Fp2_init (Fp2_w_aux->t8)) < 0) return ret;
  
  if ((ret=Fp2_init (Fp2_w_curve->b)) < 0) return ret;
  Fp2_copy (Fp2_w_curve->b, b);
  return ERR_SUCCESS;
}

void Fp2_free_weierstrass (void) {
  Fp2_free (Fp2_w_aux->t0);
  Fp2_free (Fp2_w_aux->t1);
  Fp2_free (Fp2_w_aux->t2);
  Fp2_free (Fp2_w_aux->t3);
  Fp2_free (Fp2_w_aux->t4);
  Fp2_free (Fp2_w_aux->t5);
  Fp2_free (Fp2_w_aux->t6);
  Fp2_free (Fp2_w_aux->t7);
  Fp2_free (Fp2_w_aux->t8);

  Fp2_free (Fp2_w_curve->b);
}

void Fp2_wp_copy (Point_wp2_t c, Point_wp2_t a) {
  Fp2_copy (c->X, a->X);
  Fp2_copy (c->Y, a->Y);
  Fp2_copy (c->Z, a->Z);
}

void Fp2_wpaff_copy(Point_wp2aff_t c, Point_wp2aff_t a) {
	Fp2_copy(c->X, a->X);
	Fp2_copy(c->Y, a->Y);
}

void Fp2_wp_neg (Point_wp2_t c, Point_wp2_t a) {
  Fp2_copy (c->X, a->X);
  Fp2_neg (c->Y, a->Y);
  Fp2_copy (c->Z, a->Z);
}

void Fp2_wpaff_neg(Point_wp2aff_t c, Point_wp2aff_t a) {
	Fp2_copy(c->X, a->X);
	Fp2_neg(c->Y, a->Y);
}

void Fp2_weierstrass_dbl (Point_wp2_t c, Point_wp2_t a) {
  weierstrass_dbl_template (c, a, Fp2_w_aux);
}

void Fp2_weierstrass_add (Point_wp2_t c, Point_wp2_t a, Point_wp2_t b) {
  weierstrass_add_template (c, a, b, Fp2_w_aux);
}

void Fp2_weierstrass_mix (Point_wp2_t c, Point_wp2_t a, Point_wp2_t b) {
 weierstrass_mix_template (c, a, b, Fp2_w_aux);
}

void Fp2_weierstrass_aff (Point_wp2_t c, Point_wp2_t a, Point_wp2_t b) {
  weierstrass_aff_template (c, a, b, Fp2_w_aux);
}

void Fp2_weierstrass_affdbl(Point_wp2_t c, Point_wp2_t a) {
	weierstrass_affdbl_template(c, a, Fp2_w_aux);
}

void Fp2_weierstrass_affadd(Point_wp2_t c, Point_wp2_t a, Point_wp2_t b) {
	weierstrass_affadd_template(c, a, b, Fp2_w_aux);
}


int Fp2_weierstrass_oncurve_jacproj (Point_wp2_t x) {
  return weierstrass_oncurve_jacproj_template (x, Fp2_w_aux, Fp2_w_curve);
}

int Fp2_weierstrass_oncurve_homproj (Point_wp2_t x) {
  return weierstrass_oncurve_homproj_template (x, Fp2_w_aux, Fp2_w_curve);
}

int Fp2_weierstrass_oncurve_aff (Point_wp2_t x) {
  return weierstrass_oncurve_aff_template (x, Fp2_w_aux, Fp2_w_curve);
}

int Fp2_weierstrass_equal_aff (Point_wp2_t x, Point_wp2_t y) {
  return weierstrass_equal_aff_template (x, y);
}

int Fp2_weierstrass_equal_jacproj(Point_wp2_t a, Point_wp2_t b) {
	return weierstrass_equal_jacproj_template(a, b, Fp2_w_aux);
}

void Fp2_weierstrass_aff2proj (Point_wp2_t c, Point_wp2_t a) {
  weierstrass_aff2proj_template (c, a);
}

void Fp2_weierstrass_jacproj2aff (Point_wp2_t c, Point_wp2_t a) {
  weierstrass_jacproj2aff_template (c, a, Fp2_w_aux);
}

void Fp2_weierstrass_homproj2aff (Point_wp2_t c, Point_wp2_t a) {
  weierstrass_homproj2aff_template (c, a, Fp2_w_aux);
}


#if 0
int main () {
  uint64_t p[4], b0[4], b1[4];
  Fp2_t B;
  Point_wp2_t P, Q;

  p[3] = 0x2523648240000001;
  p[2] = 0xBA344D8000000008;
  p[1] = 0x6121000000000013;
  p[0] = 0xA700000000000013;
  
  Fp_initialize_config (p, 4);
  Fp2_initialize_config ();
  Fp2_init (B);

  b1[3] = 0x2523648240000001;
  b1[2] = 0xBA344D8000000008;
  b1[1] = 0x6121000000000013;
  b1[0] = 0xA700000000000012;
  
  b0[3] = 0x0000000000000000;
  b0[2] = 0x0000000000000000;
  b0[1] = 0x0000000000000000;
  b0[0] = 0x0000000000000001;
  Fp2_set (B, b0, b1);
  
  Fp2_initialize_weierstrass (B);

  Fp2_wp_init (P);
  Fp2_wp_init (Q);

  b1[3] = 0x1C9EF62920C85936;
  b1[2] = 0x371AED1A1EB30C83;
  b1[1] = 0x006A310CBA6D1E26;
  b1[0] = 0xE429E3CBDB314928;
  
  b0[3] = 0x1ABE36D06DD07358;
  b0[2] = 0xB7CBB532CF6961DF;
  b0[1] = 0x743BD07961A287A3;
  b0[0] = 0x59B0D4C307B5C914;
  Fp2_set (P->X, b0, b1);
  

  b1[3] = 0x085B1A1EF635AD65;
  b1[2] = 0x29F489BAF930310E;
  b1[1] = 0x70B25864B6178492;
  b1[0] = 0x5E66EEAB664F1C4B;
  
  b0[3] = 0x079A4686E85D1054;
  b0[2] = 0xFEEAE258FB2103C8;
  b0[1] = 0xD6A34D4C2FDD968F;
  b0[0] = 0x6526DE28C4CC4EF3;
  Fp2_set (P->Y, b0, b1);
  
  //printf("X = "); Fp2_print(P->X);
  //printf("Y = "); Fp2_print(P->Y);
  
  printf ("Affine point on curve?: %d\n", Fp2_weierstrass_oncurve_aff (P));
  
  Fp2_weierstrass_aff2proj (P, P);
  printf ("Projective point on curve?: %d\n", Fp2_weierstrass_oncurve_jacproj (P));

  Fp2_weierstrass_dbl (Q, P);
  Fp2_weierstrass_add (Q, Q, P);
  Fp2_weierstrass_add (Q, Q, P);
  Fp2_weierstrass_add (Q, Q, P);
  Fp2_weierstrass_dbl (Q, Q);

  printf ("Projective point on curve?: %d\n", Fp2_weierstrass_oncurve_jacproj (Q));

  Fp2_weierstrass_jacproj2aff (P, Q);
  
  printf ("Affine point on curve?: %d\n", Fp2_weierstrass_oncurve_aff (P));

  Fp2_wp_free (P);
  Fp2_wp_free (Q);

  Fp2_free_weierstrass ();
  Fp2_free_config (); 
  Fp_free_config ();
  getchar ();
}
#endif

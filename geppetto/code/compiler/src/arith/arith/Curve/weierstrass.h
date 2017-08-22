/* curve/weierstrass.h
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

/* weierstrass.c
 * This is a template file for Weierstrass curve arithmetic. 
 * The user should define the data-types:
 * POINT
 * MEM
 * CURVE
 *
 * and the arithmetic routines:
 * MUL
 * MUL2
 * MUL3
 * MUL4
 * MUL8
 * SQR
 * ADD
 * SUB
 * INV
 * SET
 * SETONE
 * CMP
 */

//2M + 5S + 6add + 1*8 + 3*2 + 1*3
static void weierstrass_dbl_template (POINT c, POINT a, MEM mem) {
  MUL  (c->Z, a->Y, a->Z);
  MUL2 (c->Z, c->Z);
  SQR  (c->Y, a->Y);
  ADD  (mem->t0, a->X, c->Y);
  SQR  (c->X, a->X);
  SQR  (c->Y, c->Y);
  SQR  (mem->t0, mem->t0);
  SUB  (mem->t0, mem->t0, c->X);
  SUB  (mem->t0, mem->t0, c->Y);
  MUL2 (mem->t0, mem->t0);
  MUL2 (mem->t1, mem->t0);
  MUL8 (c->Y, c->Y); 
  MUL3 (mem->t2, c->X);
  SQR  (mem->t3, mem->t2);
  SUB  (c->X, mem->t3, mem->t1);
  SUB  (mem->t0, mem->t0, c->X);
  MUL  (mem->t0, mem->t2, mem->t0);
  SUB  (c->Y, mem->t0, c->Y);
}

//11M+5S+9a+4m2
static void weierstrass_add_template (POINT c, POINT a, POINT b, MEM mem) {
  SQR (mem->t5, a->Z);
  SQR (mem->t6, b->Z);
  MUL (mem->t7, a->X, mem->t6);
  MUL (mem->t8, b->X, mem->t5);
  MUL (mem->t0, b->Z, mem->t6);
  MUL (mem->t0, a->Y, mem->t0);
  MUL (mem->t1, a->Z, mem->t5);
  MUL (mem->t1, b->Y, mem->t1);
  SUB (mem->t8, mem->t8, mem->t7);
  ADD (c->Z, a->Z, b->Z);
  SQR (c->Z, c->Z);
  SUB (c->Z, c->Z, mem->t5);
  SUB (c->Z, c->Z, mem->t6);
  MUL (c->Z, c->Z, mem->t8);
  MUL2 (mem->t2, mem->t8);
  SQR (mem->t2, mem->t2);
  SUB (mem->t3, mem->t1, mem->t0);
  MUL (mem->t4, mem->t8, mem->t2);
  MUL2 (mem->t3, mem->t3);
  MUL (c->Y, mem->t7, mem->t2);
  SQR (c->X, mem->t3);
  SUB (c->X, c->X, mem->t4);
  MUL2 (mem->t2, c->Y);
  SUB (c->X, c->X, mem->t2);
  SUB (c->Y, c->Y, c->X);
  MUL (mem->t0, mem->t0, mem->t4);
  MUL2 (mem->t0, mem->t0);
  MUL (c->Y, mem->t3, c->Y);
  SUB (c->Y, c->Y, mem->t0);
}

// 7M + 4S + 9add + 1*4 + 3*2
static void weierstrass_mix_template (POINT c, POINT a, POINT b, MEM mem) {
  SQR (mem->t4, a->Z);
  MUL (mem->t0, a->Z, mem->t4);
  MUL (mem->t7, b->X, mem->t4);
  MUL (mem->t1, b->Y, mem->t0);
  SUB (mem->t7, mem->t7, a->X);
  SQR (mem->t6, mem->t7);
  ADD (c->Z, a->Z, mem->t7);
  SQR (c->Z, c->Z);
  SUB (c->Z, c->Z, mem->t4);
  SUB (c->Z, c->Z, mem->t6);
  MUL4 (mem->t4, mem->t6);
  MUL (mem->t7, mem->t7, mem->t4);
  SUB (mem->t1, mem->t1, a->Y);
  MUL2 (mem->t5, mem->t1);
  MUL (mem->t4, a->X, mem->t4);
  SQR (mem->t2, mem->t5);
  MUL2 (mem->t3, mem->t4);
  SUB (mem->t2, mem->t2, mem->t7);
  SUB (c->X, mem->t2, mem->t3);
  SUB (mem->t1, mem->t4, c->X);
  MUL (mem->t0, a->Y, mem->t7);
  MUL (c->Y, mem->t5, mem->t1);
  MUL2 (mem->t0, mem->t0);
  SUB (c->Y, c->Y, mem->t0);
}

//4M + 2S + 6add + 1*4 + 4*2
static void weierstrass_aff_template (POINT c, POINT a, POINT b, MEM mem) {
  SUB (c->Z, b->X, a->X);
  SUB (mem->t0, b->Y, a->Y);
  SQR (mem->t6, c->Z);
  MUL4 (mem->t6, mem->t6);
  MUL (mem->t4, c->Z, mem->t6);
  MUL2 (mem->t5, mem->t0);
  MUL (mem->t0, a->X, mem->t6);
  SQR (mem->t1, mem->t5);
  MUL2 (mem->t2, mem->t0);
  SUB (mem->t3, mem->t1, mem->t4);
  SUB (c->X, mem->t3, mem->t2);
  SUB (mem->t2, mem->t0, c->X);
  MUL (mem->t3, a->Y, mem->t4);
  MUL2 (mem->t3, mem->t3);
  MUL (mem->t2, mem->t5, mem->t2);
  SUB (c->Y, mem->t2, mem->t3);
  MUL2 (c->Z, c->Z);
}

static void weierstrass_affdbl_template (POINT c, POINT a, MEM mem) {
	//lambda: = x1 ^ 2;
	SQR(mem->t0, a->X);
	//c: = lambda + lambda;
	MUL2(mem->t1, mem->t0);
	//lambda: = c + lambda;
	ADD(mem->t0, mem->t0, mem->t1);
	//c: = 2 * y1;
	MUL2(mem->t1, a->Y);
	//c: = 1 / c;
	INV(mem->t1, mem->t1);
	//lambda: = lambda*c;
	MUL(mem->t0, mem->t0, mem->t1);
	//c: = lambda ^ 2;
	SQR(mem->t1, mem->t0);
	//c: = c - x1;
	SUB(mem->t1, mem->t1, a->X);
	//x3: = c - x1;
	SUB(mem->t2, mem->t1, a->X);
	//y3: = lambda*x3;
	MUL(mem->t3, mem->t0, mem->t2);
	//lambda: = lambda*x1;
	MUL(mem->t0, mem->t0, a->X);
	//y3: = y3 + y1;
	ADD(c->Y, mem->t3, a->Y);
	//y3: = lambda - y3;
    SUB(c->Y, mem->t0, c->Y);
	SET(c->X, mem->t2);
}

static void weierstrass_affadd_template(POINT c, POINT a, POINT b, MEM mem) {
	//lambda: = y2 - y1;
	SUB(mem->t0, b->Y, a->Y);
	//c: = x2 - x1;
	SUB(mem->t1, b->X, a->X);
	//c: = 1 / c;
	INV(mem->t1, mem->t1);
	//lambda: = lambda*c;
	MUL(mem->t0, mem->t0, mem->t1);
	//c: = lambda ^ 2;
	SQR(mem->t1, mem->t0);
	//x3: = c - x1;
	SUB(mem->t2, mem->t1, a->X);
	//x3: = x3 - x2;
	SUB(mem->t2, mem->t2, b->X);
	//y3: = lambda*x3;
	MUL(mem->t3, mem->t0, mem->t2);
	//lambda: = x1*lambda;
	MUL(mem->t0, mem->t0, a->X);
	//y3: = y3 + y1;
	ADD(c->Y, mem->t3, a->Y);
	//y3: = lambda - y3;
	SUB(c->Y, mem->t0, c->Y);
	SET(c->X, mem->t2);
}


static int weierstrass_oncurve_jacproj_template (POINT x, MEM mem, CURVE curve) {
  SQR (mem->t2, x->Y); // left = t2
  SQR (mem->t3, x->X); // right = t3
	MUL (mem->t3, mem->t3, x->X);
	SQR (mem->t0, x->Z);
	SQR (mem->t1, mem->t0);
	MUL (mem->t0, mem->t0, mem->t1);
	MUL (mem->t0, curve->b, mem->t0);
	ADD (mem->t3, mem->t3, mem->t0);
  return (CMPEQ(mem->t2, mem->t3));
}

static int weierstrass_oncurve_homproj_template (POINT x, MEM mem, CURVE curve) {
  SQR (mem->t2, x->Y); // left = t2
  MUL (mem->t2, mem->t2, x->Z);
  SQR (mem->t3, x->X); // right = t3
	MUL (mem->t3, mem->t3, x->X);
	SQR (mem->t0, x->Z);
	MUL (mem->t0, mem->t0, x->Z);
	MUL (mem->t0, curve->b, mem->t0);
	ADD (mem->t3, mem->t3, mem->t0);
  return (CMPEQ(mem->t2, mem->t3));
}

static int weierstrass_oncurve_aff_template (POINT x, MEM mem, CURVE curve) {
  SQR (mem->t0, x->Y); // left = t0
  SQR (mem->t1, x->X); // right = t1
  MUL (mem->t1, mem->t1, x->X);
  ADD (mem->t1, mem->t1, curve->b);
  return (CMPEQ (mem->t0, mem->t1));
}

static int weierstrass_equal_aff_template (POINT x, POINT y) {
  return ((CMPEQ (x->X, y->X)) && (CMPEQ (x->Y, y->Y) ));
}

static int weierstrass_equal_jacproj_template(POINT a, POINT b, MEM mem) {
	SQR(mem->t0, a->Z);
	SQR(mem->t1, b->Z);
	MUL(mem->t2, mem->t0, b->X);
	MUL(mem->t3, mem->t1, a->X);
	MUL(mem->t0, mem->t0, a->Z);
	MUL(mem->t1, mem->t1, b->Z);
	MUL(mem->t0, mem->t0, b->Y);
	MUL(mem->t1, mem->t1, a->Y);
	return ((CMPEQ(mem->t2, mem->t3)) && (CMPEQ(mem->t0, mem->t1)));
}


static  void weierstrass_aff2proj_template (POINT c, POINT a) {
  SET (c->X, a->X);
  SET (c->Y, a->Y);
  SETONE (c->Z);
}

static void weierstrass_jacproj2aff_template (POINT c, POINT a, MEM mem) {
  SQR (mem->t0, a->Z);
  MUL (mem->t0, mem->t0, a->Z);
  INV (mem->t0, mem->t0);
  MUL (c->Y, a->Y, mem->t0);
  MUL (c->X, a->X, mem->t0);
  MUL (c->X, c->X, a->Z);
  SETONE (c->Z);
}

static void weierstrass_homproj2aff_template (POINT c, POINT a, MEM mem) {
  INV (mem->t0, a->Z);
  MUL (c->Y, a->Y, mem->t0);
  MUL (c->X, a->X, mem->t0);
  SETONE (c->Z);
}


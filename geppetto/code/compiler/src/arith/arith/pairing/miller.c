/* pairing/miller.c
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

#include "pairing.h"

#define BIT(a,i) (((a)[i/64] >> (uint64_t) (i%64))&1)

int PAIR_CURVE = 0;

bn_config_t bn_config;
cp6_config_t cp6_config;
cp3_config_t cp3_config;


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
  Fp2_t l00;
  Fp2_t l01;
  Fp2_t l10;
  Fp_t negxP;
  Fp_t negyP;
  Fp2_t negyP_Fp2;
  Point_wp2_t R;
  Point_wp2_t negQ;
  Point_wp2_t Qp;
  Point_wp2_t Qp2;
  Fp12_t f;
  uint64_t *premem;
} aux_miller12_t;

typedef aux_miller12_t Aux_miller12_t[1];

static Aux_miller12_t miller12_aux;


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
  Fp_t l00;
  Fp_t l01;
  Fp_t l10;
  Fp_t negxP;
  Fp_t negyP;
  Point_wp_t R;
  Point_wp_t U;
  Point_wp_t V;
  Point_wp_t negQ;
  Point_wp_t negU;
  Point_wp_t negV;
  Fp6q_t f;
  Fp6q_t f1;
  Fp6q_t f1inv;
  Fp6q_t F;
} aux_miller6_t;

typedef aux_miller6_t Aux_miller6_t[1];

static Aux_miller6_t miller6_aux;

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
	Fp_t l00;
	Fp_t l01;
	Fp_t l10;
	Fp_t l20;
	Fp_t xP2;
	Fp_t negxP;
	Fp_t negyP;
	Point_wp_t R;
	Point_wp_t Ptmp;
	Point_wp_t negQ;
	Fp3_t f;
	Fp3_t f1;
	Fp3_t f1inv;
	Fp3_t F;
} aux_miller3_t;

typedef aux_miller3_t Aux_miller3_t[1];

static Aux_miller3_t miller3_aux;


/* Line function for a doubling step on BN curves.
 * The coefficients returned by this function are the coefficients of the line
 * l := L00+L01*yQ+L10*xQ, note that this assumes that yQ:=-yQ has been precomputed.
 * Also note that b3 = 3*bt, where the twist curve is y^2 = x^3 + bt.
 */
static void bn_bls12_line_dbl_mem (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Fp2_t b3, Aux_miller12_t mem) {
    Fp2_squ (mem->t0, R->X);              // t0 = A:=X1^2;
    Fp2_squ (mem->t1, R->Y);              // t1 = B:=Y1^2; 
	Fp2_squ (mem->t2, R->Z);              // t2 = C:=Z1^2;
	Fp2_mul (mem->t3, mem->t2, b3);       // t3 = D:=Mulb3(C,b);                  		
	Fp2_add (mem->t4, R->X, R->Y);        // t4 = E:=X1+Y1;
	Fp2_squ (mem->t4, mem->t4);           // t4 = E:=E^2;    
	Fp2_sub (mem->t4, mem->t4, mem->t0);  // t4 = E:=E-A;
	Fp2_sub (mem->t4, mem->t4, mem->t1);  // t4 = E:=E-B;
    Fp2_add (l01, R->Y, R->Z);            // L01:=Y1+Z1;
	Fp2_squ (l01, l01);                   // L01:=L01^2;
	Fp2_sub (l01, l01, mem->t1);          // L01:=L01-B;
	Fp2_sub (l01, l01, mem->t2);          // L01:=L01-C;       		
	Fp2_add (mem->t5, mem->t3, mem->t3);  // t5 = H:=2*D;
	Fp2_add (mem->t6, mem->t5, mem->t3);  // t6 = G:=H+D;
	Fp2_sub (R->X, mem->t1, mem->t6);     // X3:=B-G;
	Fp2_mul (R->X, mem->t4, R->X);        // X3:=E*X3; 
	Fp2_add (R->Y, mem->t1, mem->t6);     // Y3:=B+G;
	Fp2_squ (R->Y, R->Y);                 // Y3:=Y3^2;  
	Fp2_squ (mem->t5, mem->t5);           // H:=H^2;  
	Fp2_add (mem->t4, mem->t5, mem->t5);  // E:=2*H;
	Fp2_add (mem->t4, mem->t4, mem->t5);  // E:=E+H; 		
	Fp2_sub (R->Y, R->Y, mem->t4);        // Y3:=Y3-E;   
	Fp2_mul (R->Z, mem->t1, l01);         // Z3:=B*L01;
	Fp2_add (R->Z, R->Z, R->Z);           // Z3:=2*Z3;
	Fp2_add (R->Z, R->Z, R->Z);           // Z3:=2*Z3;   		
	Fp2_add (l10, mem->t0, mem->t0);      // L10:=2*A;
	Fp2_add (l10, l10, mem->t0);          // L10:=L10+A;
	Fp2_sub (l00, mem->t3, mem->t1);      // L00:=D-B;
}

void bn_bls12_line_dbl (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Fp2_t b3) {
  bn_bls12_line_dbl_mem (l00, l01, l10, R, b3, miller12_aux);
}

/* Line function for a doubling step on BN curves using affine coordinates.
 * The coefficients returned by this function are the coefficients of the line
 * l: = -yQ*v + (lambda*xQ*v*w + nu*v ^ 2 * w), note that this assumes that yQ:=-yQ has been precomputed.
 */
static void bn_bls12_line_affdbl_mem(Fp2_t lambda, Fp2_t nu, Point_wp2_t R, Aux_miller12_t mem) {
	Fp2_squ(lambda, R->X); //lambda: = x1 ^ 2;
	Fp2_add(mem->t1, lambda, lambda); //c: = lambda + lambda;
	Fp2_add(lambda, lambda, mem->t1); //lambda: = c + lambda;
	Fp2_add(mem->t1, R->Y, R->Y); //c: = 2 * y1;
	Fp2_inv(mem->t1, mem->t1); //c: = 1 / c;
	Fp2_mul(lambda, lambda, mem->t1); //lambda: = lambda*c;
	Fp2_mul(mem->t1, lambda, R->X); //c: = lambda*x1;
	Fp2_sub(nu, R->Y, mem->t1);//nu: = y1 - c;
	Fp2_squ(mem->t0, lambda); //sqr: = lambda ^ 2;
	Fp2_add(mem->t1, R->X, R->X); //c: = 2 * x1;
	Fp2_sub(R->X, mem->t0, mem->t1);//x3: = sqr - c;
	Fp2_mul(R->Y, R->X, lambda); //y3: = lambda*x3;
	Fp2_add(R->Y, R->Y, nu); //y3: = y3 + nu;
	Fp2_neg(R->Y, R->Y); //y3: = -y3;
}

void bn_bls12_line_affdbl(Fp2_t lambda, Fp2_t nu, Point_wp2_t R) {
	bn_bls12_line_affdbl_mem(lambda, nu, R, miller12_aux);
}


/* Line function for an addition step on BN curves.
 * The coefficients returned by this function are the coefficients of the line
 * l := L00+L01*yQ+L10*xQ, not that this assumes that xQ:=-xQ has been precomputed.
 */
static void bn_bls12_line_add_mem (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Point_wp2_t P, Aux_miller12_t mem) {
    Fp2_mul (mem->t1, R->Z, P->X);    //t1:=Z1*X2; 
	Fp2_sub (mem->t1, R->X, mem->t1); //t1:=X1-t1; 
	Fp2_mul (mem->t2, R->Z, P->Y);    //t2:=Z1*Y2; 
	Fp2_sub (mem->t2, R->Y, mem->t2); //t2:=Y1-t2; 
	Fp2_mul (l00, P->X, mem->t2);     //L00:=X2*t2;
	Fp2_mul (mem->t3, mem->t1, P->Y); //t3:=t1*Y2;
	Fp2_sub (l00, l00, mem->t3);      //L00:=L00-t3;	
	Fp2_copy (l10, mem->t2);          //L10:=t2;
	Fp2_copy (l01, mem->t1);          //L01:=t1;
	Fp2_squ (mem->t3, mem->t1);       //t3:=t1^2;  
	Fp2_mul (R->X, mem->t3, R->X);    //X3:=t3*X1; 
	Fp2_mul (mem->t3, mem->t1, mem->t3);//t3:=t1*t3; 
	Fp2_squ (mem->t4, mem->t2);       //t4:=t2^2; 
	Fp2_mul (mem->t4, mem->t4, R->Z); //t4:=t4*Z1; 
	Fp2_add (mem->t4, mem->t3, mem->t4);//t4:=t3+t4; 
	Fp2_sub (mem->t4, mem->t4, R->X); //t4:=t4-X3; 
	Fp2_sub (mem->t4, mem->t4, R->X); //t4:=t4-X3; 
	Fp2_sub (R->X, R->X, mem->t4);    //X3:=X3-t4; 
	Fp2_mul (mem->t2, mem->t2, R->X); //t2:=t2*X3; 
	Fp2_mul (R->Y, mem->t3, R->Y);    //Y3:=t3*Y1; 
	Fp2_sub (R->Y, mem->t2, R->Y);    //Y3:=t2-Y3; 	
	Fp2_mul (R->X, mem->t1, mem->t4); //X3:=t1*t4; 
    Fp2_mul (R->Z, R->Z, mem->t3);    //Z3:=Z1*t3; 
}

void bn_bls12_line_add (Fp2_t l00, Fp2_t l01, Fp2_t l10, Point_wp2_t R, Point_wp2_t P) {
  bn_bls12_line_add_mem (l00, l01, l10, R, P, miller12_aux);
}

/* Line function for an addition step on BN curves using affine coordinates.
* The coefficients returned by this function are the coefficients of the line
* l: = -yQ*v + (lambda*xQ*v*w + nu*v ^ 2 * w), not that this assumes that xQ:=-xQ has been precomputed.
*/
static void bn_bls12_line_affadd_mem(Fp2_t lambda, Fp2_t nu, Point_wp2_t R, Point_wp2_t P, Aux_miller12_t mem) {
	Fp2_sub(lambda, P->Y, R->Y); //lambda: = y2 - y1;
	Fp2_sub(mem->t0, P->X, R->X); //c: = x2 - x1;
	Fp2_inv(mem->t0, mem->t0); //c: = 1 / c;
	Fp2_mul(lambda, lambda, mem->t0); //lambda: = lambda*c;
	Fp2_squ(mem->t1, lambda); //sqr: = lambda ^ 2;
	Fp2_mul(nu, R->X, lambda); //nu: = x1*lambda;
	Fp2_sub(R->X, mem->t1, R->X); //x3: = sqr - x1;
	Fp2_sub(R->X, R->X, P->X); //x3: = x3 - x2;
	Fp2_sub(nu, R->Y, nu); //nu: = y1 - nu;
	Fp2_mul(R->Y, lambda, R->X); //y3: = lambda*x3;
	Fp2_add(R->Y, R->Y, nu); //y3: = y3 + nu;
	Fp2_neg(R->Y, R->Y); //y3: = -y3;
}

void bn_bls12_line_affadd(Fp2_t lambda, Fp2_t nu, Point_wp2_t R, Point_wp2_t P) {
	bn_bls12_line_affadd_mem(lambda, nu, R, P, miller12_aux);
}


static void bn_Fp2_wp_ppow_untwist (Point_wp2_t Qp, Point_wp2_t Q) {
  Fp2_ppow (Qp->X, Q->X);
  Fp2_ppow (Qp->Y, Q->Y);
  Fp2_ppow (Qp->Z, Q->Z);
  Fp2_mul (Qp->Y, Qp->Y, Fp12_config.wppow3);
  Fp2_muli (Qp->X, Qp->X);
  Fp2_mulFp (Qp->X, Qp->X, Fp6_config.vppow);
}

static void bn_Fp2_wpaff_ppow_untwist(Point_wp2_t Qp, Point_wp2_t Q) {
	Fp2_ppow(Qp->X, Q->X);
	Fp2_ppow(Qp->Y, Q->Y);
	Fp2_mul(Qp->Y, Qp->Y, Fp12_config.wppow3);
	Fp2_muli(Qp->X, Qp->X);
	Fp2_mulFp(Qp->X, Qp->X, Fp6_config.vppow);
}

static void bn_Fp2_wp_p2pow_untwist (Point_wp2_t Qp, Point_wp2_t Q) {
  Fp2_mulFp (Qp->X, Q->X, Fp6_config.vp2pow);
  Fp2_mulFp (Qp->Y, Q->Y, Fp12_config.wp2pow3);
}

static void bn_optimal_ate_miller_mem (Fp12_t f, Point_wp2_t Q, Point_wp_t P, Aux_miller12_t mem) {
  int i;
  
  Fp12_set_zero(f);
  Fp_neg (mem->negxP, P->X);
  Fp_neg (mem->negyP, P->Y);
  Fp2_wp_copy (mem->R, Q);
  Fp2_wp_neg (mem->negQ, Q);

  bn_bls12_line_dbl (f->a1->a2, mem->l01, mem->l10, mem->R, bn_config.bt3);
  Fp2_mulFp (f->a0->a1, mem->l01, mem->negyP);
  Fp2_mulFp (f->a1->a1, mem->l10, P->X);
  if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
    bn_bls12_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
    Fp2_mulFp (mem->l01, mem->l01, P->Y);
    Fp2_mulFp (mem->l10, mem->l10, mem->negxP);
    Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00); 
  }
  if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
    bn_bls12_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->negQ);
    Fp2_mulFp (mem->l01, mem->l01, P->Y);
    Fp2_mulFp (mem->l10, mem->l10, mem->negxP);
    Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00); 
  }
  //printf ("1: f = "); Fp12_print (f);
  for (i=bn_config.len6u2 - 3; i>=0; i--) {
    Fp12_squ_lazy (f, f);
    //if (i>62) {printf ("1squ.%d: f = ", i); Fp12_print (f);}
    bn_bls12_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, bn_config.bt3);
    Fp2_mulFp (mem->l01, mem->l01, mem->negyP);
    Fp2_mulFp (mem->l10, mem->l10, P->X);
    Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00);
    //if (i>62) {printf ("1squstep.%d: f = ", i); Fp12_print (f);}
    if (bn_config.naf6u2[i] == 1) {
      bn_bls12_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
      Fp2_mulFp (mem->l01, mem->l01, P->Y);
      Fp2_mulFp (mem->l10, mem->l10, mem->negxP);
      Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00);
    }
    if (bn_config.naf6u2[i] == -1) {
      bn_bls12_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->negQ);
      Fp2_mulFp (mem->l01, mem->l01, P->Y);
      Fp2_mulFp (mem->l10, mem->l10, mem->negxP);
      Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00);
    }
    //if (i>62) {printf ("i = %d: f = ", i); Fp12_print (f);}  
  }
  //printf ("2: f = "); Fp12_print (f);
  bn_Fp2_wp_ppow_untwist (mem->Qp, Q);
  bn_Fp2_wp_p2pow_untwist (mem->Qp2, Q);
  Fp2_wp_neg (mem->Qp2, mem->Qp2);
  
  if (bn_config.signu == -1) {
    Fp2_wp_neg (mem->R, mem->R);
    Fp12_p6pow (f, f);
  }
  
  bn_bls12_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->Qp);
  Fp2_mulFp (mem->l01, mem->l01, P->Y);
  Fp2_mulFp (mem->l10, mem->l10, mem->negxP);
  Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00);

  bn_bls12_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->Qp2);
  Fp2_mulFp (mem->l01, mem->l01, P->Y);
  Fp2_mulFp (mem->l10, mem->l10, mem->negxP);
  Fp12_mul_sparse_untwist_lazy (f, f, mem->l01, mem->l10, mem->l00);
  //printf ("3: f = "); Fp12_print (f);
}

void bn_optimal_ate_miller (Fp12_t f, Point_wp2_t Q, Point_wp_t P) {
  bn_optimal_ate_miller_mem (f, Q, P, miller12_aux);
}


static void bn_optimal_ate_miller_aff_mem(Fp12_t f, Point_wp2_t Q, Point_wp_t P, Aux_miller12_t mem) {
	int i;

	Fp12_set_zero(f);
	Fp_neg(mem->negyP, P->Y);
	Fp_set_ui(mem->negyP_Fp2->a1, 0);
	Fp_copy(mem->negyP_Fp2->a0, mem->negyP);
	Fp2_wp_copy(mem->R, Q);
	Fp2_wp_neg(mem->negQ, Q);

	bn_bls12_line_affdbl(mem->l01, mem->l10, mem->R);
	Fp_set_ui(f->a0->a1->a1, 0);
	Fp_copy (f->a0->a1->a0, mem->negyP);
	Fp2_mulFp(f->a1->a1, mem->l01, P->X);
	Fp2_copy(f->a1->a2, mem->l10);
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
		bn_bls12_line_affadd(mem->l01, mem->l10, mem->R, Q);
		Fp2_mulFp(mem->l01, mem->l01, P->X);
		Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
	}
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
		bn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->negQ);
		Fp2_mulFp(mem->l01, mem->l01, P->X);
		Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
	}
	//printf ("1: f = "); Fp12_print (f);
	for (i = bn_config.len6u2 - 3; i >= 0; i--) {
		Fp12_squ_lazy(f, f);
		//if (i>62) {printf ("1squ.%d: f = ", i); Fp12_print (f);}
		bn_bls12_line_affdbl(mem->l01, mem->l10, mem->R);
		Fp2_mulFp(mem->l01, mem->l01, P->X);
		Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
		//if (i>62) {printf ("1squstep.%d: f = ", i); Fp12_print (f);}
		if (bn_config.naf6u2[i] == 1) {
			bn_bls12_line_affadd(mem->l01, mem->l10, mem->R, Q);
			Fp2_mulFp(mem->l01, mem->l01, P->X);
			Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
		}
		if (bn_config.naf6u2[i] == -1) {
			bn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->negQ);
			Fp2_mulFp(mem->l01, mem->l01, P->X);
			Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
		}
		//if (i>62) {printf ("i = %d: f = ", i); Fp12_print (f);}  
	}
	//printf ("2: f = "); Fp12_print (f);
	bn_Fp2_wpaff_ppow_untwist(mem->Qp, Q);
	bn_Fp2_wp_p2pow_untwist(mem->Qp2, Q);
	Fp2_wp_neg(mem->Qp2, mem->Qp2);

	if (bn_config.signu == -1) {
		Fp2_wp_neg(mem->R, mem->R);
		Fp12_p6pow(f, f);
	}

	bn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->Qp);
	Fp2_mulFp(mem->l01, mem->l01, P->X);
	Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);

	bn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->Qp2);
	Fp2_mulFp(mem->l01, mem->l01, P->X);
	Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
	//printf ("3: f = "); Fp12_print (f);
}

void bn_optimal_ate_miller_aff(Fp12_t f, Point_wp2_t Q, Point_wp_t P) {
	bn_optimal_ate_miller_aff_mem(f, Q, P, miller12_aux);
}

static void bn_optimal_ate_miller_precompute_aff_mem(Point_wp2_t Q, uint64_t *premem, Aux_miller12_t mem) {
	int i, j;
	Fp2_t l01, l10;

	Fp2_wp_copy(mem->R, Q);
	Fp2_wp_neg(mem->negQ, Q);

	j = 0;
	l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
	l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
	l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
	l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
	bn_bls12_line_affdbl(l01, l10, mem->R);
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
		j++;
		l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
		l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
		l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
		l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
		bn_bls12_line_affadd(l01, l10, mem->R, Q);
	}
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
		j++;
		l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
		l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
		l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
		l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
		bn_bls12_line_affadd(l01, l10, mem->R, mem->negQ);
	}
	for (i = bn_config.len6u2 - 3; i >= 0; i--) {
		j++;
		l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
		l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
		l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
		l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
		bn_bls12_line_affdbl(l01, l10, mem->R);
		if (bn_config.naf6u2[i] == 1) {
			j++;
			l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
			l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
			l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
			l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
			bn_bls12_line_affadd(l01, l10, mem->R, Q);
		}
		if (bn_config.naf6u2[i] == -1) {
			j++;
			l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
			l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
			l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
			l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
			bn_bls12_line_affadd(l01, l10, mem->R, mem->negQ);
		}
	}
	bn_Fp2_wpaff_ppow_untwist(mem->Qp, Q);
	bn_Fp2_wp_p2pow_untwist(mem->Qp2, Q);
	Fp2_wp_neg(mem->Qp2, mem->Qp2);

	if (bn_config.signu == -1) {
		Fp2_wp_neg(mem->R, mem->R);
	}

	j++;
	l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
	l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
	l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
	l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
	bn_bls12_line_affadd(l01, l10, mem->R, mem->Qp);

	j++;
	l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
	l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
	l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
	l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
	bn_bls12_line_affadd(l01, l10, mem->R, mem->Qp2);
}

void bn_optimal_ate_miller_precompute_aff(Point_wp2_t Q, uint64_t *premem) {
	bn_optimal_ate_miller_precompute_aff_mem(Q, premem, miller12_aux);
}

static void bn_optimal_ate_miller_useprecomputed_aff_mem(Fp12_t f, Point_wp_t P, uint64_t *premem, Aux_miller12_t mem) {
	int i, j;
	Fp2_t l01, l10;

	Fp12_set_zero(f);
	Fp_neg(mem->negyP, P->Y);
	Fp_set_ui(mem->negyP_Fp2->a1, 0);
	Fp_copy(mem->negyP_Fp2->a0, mem->negyP);

	j = 0;
	l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
	l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
	l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
	l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;

	Fp_set_ui(f->a0->a1->a1, 0);
	Fp_copy(f->a0->a1->a0, mem->negyP);
	Fp2_mulFp(f->a1->a1, l01, P->X);
	Fp2_copy(f->a1->a2, l10);
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
		j++;
		l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
		l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
		l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
		l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
		Fp2_mulFp(l01, l01, P->X);
		Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, l01, l10);
	}
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
		j++;
		l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
		l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
		l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
		l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
		Fp2_mulFp(mem->l01, l01, P->X);
		Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, l10);
	}
	for (i = bn_config.len6u2 - 3; i >= 0; i--) {
		Fp12_squ_lazy(f, f);
		j++;
		l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
		l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
		l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
		l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
		Fp2_mulFp(mem->l01, l01, P->X);
		Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, l10);
		if (bn_config.naf6u2[i] == 1) {
			j++;
			l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
			l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
			l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
			l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
			Fp2_mulFp(mem->l01, l01, P->X);
			Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, l10);
		}
		if (bn_config.naf6u2[i] == -1) {
			j++;
			l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
			l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
			l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
			l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
			Fp2_mulFp(mem->l01, l01, P->X);
			Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, l10);
		}
	}

	if (bn_config.signu == -1) {
		Fp12_p6pow(f, f);
	}
	
	j++;
	l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
	l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
	l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
	l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
	Fp2_mulFp(mem->l01, l01, P->X);
	Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, l10);

	j++;
	l01->a0->limbs = premem + 4 * j * Fp_config.m->n;
	l01->a1->limbs = premem + (4 * j + 1) * Fp_config.m->n;
	l10->a0->limbs = premem + (4 * j + 2) * Fp_config.m->n;
	l10->a1->limbs = premem + (4 * j + 3) * Fp_config.m->n;
	Fp2_mulFp(mem->l01, l01, P->X);
	Fp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, l10);
}

void bn_optimal_ate_miller_useprecomputed_aff(Fp12_t f, Point_wp_t P, uint64_t *premem) {
	bn_optimal_ate_miller_useprecomputed_aff_mem(f, P, premem, miller12_aux);
}


/********************************************************************************************/

/* Line function for a doubling step on the CP6 curve.
 * The coefficients returned by this function are the coefficients of the line
 * l := L00+L01*yQ+L10*xQ, note that this assumes that yQ:=-yQ has been precomputed.
 * Also note that b3 = 3*bt, where the twist curve is y^2 = x^3 + bt.
 */
static void cp6_line_dbl_mem (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Fp_t b3, Aux_miller6_t mem) {
  Fp_mul (mem->t0, R->X, R->X);        // t0 = A:=X1^2;
  Fp_mul (mem->t1, R->Y, R->Y);        // t1 = B:=Y1^2; 
	Fp_mul (mem->t2, R->Z, R->Z);        // t2 = C:=Z1^2;
	Fp_mul (mem->t3, mem->t2, b3);       // t3 = D:=Mulb3(C,b);                  		
	Fp_add (mem->t4, R->X, R->Y);        // t4 = E:=X1+Y1;
	Fp_mul (mem->t4, mem->t4, mem->t4);  // t4 = E:=E^2;    
	Fp_sub (mem->t4, mem->t4, mem->t0);  // t4 = E:=E-A;
	Fp_sub (mem->t4, mem->t4, mem->t1);  // t4 = E:=E-B;
  Fp_add (l01, R->Y, R->Z);            // L01:=Y1+Z1;
	Fp_mul (l01, l01, l01);              // L01:=L01^2;
	Fp_sub (l01, l01, mem->t1);          // L01:=L01-B;
	Fp_sub (l01, l01, mem->t2);          // L01:=L01-C;       		
	Fp_add (mem->t5, mem->t3, mem->t3);  // t5 = H:=2*D;
	Fp_add (mem->t6, mem->t5, mem->t3);  // t6 = G:=H+D;
	Fp_sub (R->X, mem->t1, mem->t6);     // X3:=B-G;
	Fp_mul (R->X, mem->t4, R->X);        // X3:=E*X3; 
	Fp_add (R->Y, mem->t1, mem->t6);     // Y3:=B+G;
	Fp_mul (R->Y, R->Y, R->Y);           // Y3:=Y3^2;  
	Fp_mul (mem->t5, mem->t5, mem->t5);  // H:=H^2;  
	Fp_add (mem->t4, mem->t5, mem->t5);  // E:=2*H;
	Fp_add (mem->t4, mem->t4, mem->t5);  // E:=E+H; 		
	Fp_sub (R->Y, R->Y, mem->t4);        // Y3:=Y3-E;   
	Fp_mul (R->Z, mem->t1, l01);         // Z3:=B*L01;
	Fp_add (R->Z, R->Z, R->Z);           // Z3:=2*Z3;
	Fp_add (R->Z, R->Z, R->Z);           // Z3:=2*Z3;   		
	Fp_add (l10, mem->t0, mem->t0);      // L10:=2*A;
	Fp_add (l10, l10, mem->t0);          // L10:=L10+A;
	Fp_sub (l00, mem->t3, mem->t1);      // L00:=D-B;
}

void cp6_line_dbl (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Fp_t b3) {
  cp6_line_dbl_mem (l00, l01, l10, R, b3, miller6_aux);
}

/* Line function for an addition step on BN curves.
 * The coefficients returned by this function are the coefficients of the line
 * l := L00+L01*yQ+L10*xQ, not that this assumes that xQ:=-xQ has been precomputed.
 */
static void cp6_line_add_mem (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Point_wp_t P, Aux_miller6_t mem) {
  Fp_mul (mem->t1, R->Z, P->X);    //t1:=Z1*X2; 
	Fp_sub (mem->t1, R->X, mem->t1); //t1:=X1-t1; 
	Fp_mul (mem->t2, R->Z, P->Y);    //t2:=Z1*Y2; 
	Fp_sub (mem->t2, R->Y, mem->t2); //t2:=Y1-t2; 
	Fp_mul (l00, P->X, mem->t2);     //L00:=X2*t2;
	Fp_mul (mem->t3, mem->t1, P->Y); //t3:=t1*Y2;
	Fp_sub (l00, l00, mem->t3);      //L00:=L00-t3;	
	Fp_copy (l10, mem->t2);          //L10:=t2;
	Fp_copy (l01, mem->t1);          //L01:=t1;
	Fp_mul (mem->t3, mem->t1, mem->t1);//t3:=t1^2;  
	Fp_mul (R->X, mem->t3, R->X);    //X3:=t3*X1; 
	Fp_mul (mem->t3, mem->t1, mem->t3);//t3:=t1*t3; 
	Fp_mul (mem->t4, mem->t2, mem->t2);//t4:=t2^2; 
	Fp_mul (mem->t4, mem->t4, R->Z); //t4:=t4*Z1; 
	Fp_add (mem->t4, mem->t3, mem->t4);//t4:=t3+t4; 
	Fp_sub (mem->t4, mem->t4, R->X); //t4:=t4-X3; 
	Fp_sub (mem->t4, mem->t4, R->X); //t4:=t4-X3; 
	Fp_sub (R->X, R->X, mem->t4);    //X3:=X3-t4; 
	Fp_mul (mem->t2, mem->t2, R->X); //t2:=t2*X3; 
	Fp_mul (R->Y, mem->t3, R->Y);    //Y3:=t3*Y1; 
	Fp_sub (R->Y, mem->t2, R->Y);    //Y3:=t2-Y3; 	
	Fp_mul (R->X, mem->t1, mem->t4); //X3:=t1*t4; 
  Fp_mul (R->Z, R->Z, mem->t3);    //Z3:=Z1*t3; 
}

void cp6_line_add (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Point_wp_t P) {
  cp6_line_add_mem (l00, l01, l10, R, P, miller6_aux);
}

/* Line function for an addition step on BN curves.
 * The coefficients returned by this function are the coefficients of the line
 * l := L00+L01*yQ+L10*xQ, not that this assumes that xQ:=-xQ has been precomputed.
 */
static void cp6_line_add_proj_mem (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Point_wp_t P, Aux_miller6_t mem) {
  Fp_mul (mem->t1, R->Z, P->X);        //t1:=Z1*X2; 
	Fp_mul (mem->t2, P->Z, R->X);        //t2:=X1*Z2; 
	Fp_sub (mem->t1, mem->t2, mem->t1);  //t1:=t2-t1;
	Fp_mul (mem->t2, R->Z, P->Y);        //t2:=Z1*Y2;
  Fp_mul (l00, R->Y, P->Z);            //L00:=Y1*Z2
  Fp_sub (mem->t2, l00, mem->t2);      //t2:=L00-t2
  Fp_mul (l00, P->X, mem->t2);         //L00:=X2*t2;
	Fp_mul (mem->t3, mem->t1, P->Y);     //t3:=t1*Y2;	
 	Fp_sub (l00, l00, mem->t3);          //L00:=L00-t3;	
	Fp_mul (l10, mem->t2, P->Z);         //L10:=t2*Z2;
  Fp_mul (l01, mem->t1, P->Z);         //L01:=t1*Z2;
	Fp_mul (mem->t3, mem->t1, mem->t1);  //t3:=t1^2;  
  Fp_mul (R->X, mem->t3, R->X);        //X3:=t3*X1; 
	Fp_mul (mem->t3, mem->t1, mem->t3);  //t3:=t1*t3; 
	Fp_mul (mem->t4, mem->t2, mem->t2);  //t4:=t2^2; 
	Fp_mul (mem->t4, mem->t4, R->Z);     //t4:=t4*Z1; 
	Fp_mul (mem->t5, mem->t4, P->Z);     //t5:=t4*Z2;
  Fp_add (mem->t4, mem->t3, mem->t5);  //t4:=t3+t5;
  Fp_mul (R->X, R->X, P->Z);           //X3:=X3*Z2;
	Fp_sub (mem->t4, mem->t4, R->X);     //t4:=t4-X3; 
	Fp_sub (mem->t4, mem->t4, R->X);     //t4:=t4-X3; 
  Fp_sub (R->X, R->X, mem->t4);        //X3:=X3-t4; 
	Fp_mul (mem->t2, mem->t2, R->X);     //t2:=t2*X3; 
	Fp_mul (R->Y, mem->t3, R->Y);        //Y3:=t3*Y1; 
 	Fp_mul (mem->t5, R->Y, P->Z);        //t5:=Y3*Z2;
  Fp_sub (R->Y, mem->t2, mem->t5);     //Y3:=t2-t5; 	
	Fp_mul (R->X, mem->t1, mem->t4);     //X3:=t1*t4; 
  Fp_mul (R->Z, R->Z, mem->t3);        //Z3:=Z1*t3; 
  Fp_mul (R->Z, R->Z, P->Z);           //Z3:=Z3*Z2; 
}

void cp6_line_add_proj (Fp_t l00, Fp_t l01, Fp_t l10, Point_wp_t R, Point_wp_t P) {
  cp6_line_add_proj_mem (l00, l01, l10, R, P, miller6_aux);
}

/* Line function for a doubling step on CP6 curves using affine coordinates.
* The coefficients returned by this function are the coefficients of the line
* l: = -yQ*v*w + lambda*xQ*v + nu.
*/
static void cp6_line_affdbl_mem(Fp_t lambda, Fp_t nu, Point_wp_t R, Aux_miller6_t mem) {
	Fp_mul(lambda, R->X, R->X); //lambda: = x1 ^ 2;
	Fp_add(mem->t1, lambda, lambda); //c: = lambda + lambda;
	Fp_add(lambda, lambda, mem->t1); //lambda: = c + lambda;
	Fp_add(mem->t1, R->Y, R->Y); //c: = 2 * y1;
	Fp_modinv(mem->t1, mem->t1); //c: = 1 / c;
	Fp_mul(lambda, lambda, mem->t1); //lambda: = lambda*c;
	Fp_mul(mem->t1, lambda, R->X); //c: = lambda*x1;
	Fp_sub(nu, R->Y, mem->t1);//nu: = y1 - c;
	Fp_mul(mem->t0, lambda, lambda); //sqr: = lambda ^ 2;
	Fp_add(mem->t1, R->X, R->X); //c: = 2 * x1;
	Fp_sub(R->X, mem->t0, mem->t1);//x3: = sqr - c;
	Fp_mul(R->Y, R->X, lambda); //y3: = lambda*x3;
	Fp_add(R->Y, R->Y, nu); //y3: = y3 + nu;
	Fp_neg(R->Y, R->Y); //y3: = -y3;
}

void cp6_line_affdbl(Fp_t lambda, Fp_t nu, Point_wp_t R) {
	cp6_line_affdbl_mem(lambda, nu, R, miller6_aux);
}

/* Line function for an addition step on CP6 curves using affine coordinates.
* The coefficients returned by this function are the coefficients of the line
* l: = -yQ*v*w + lambda*xQ*v + nu.
*/
static void cp6_line_affadd_mem(Fp_t lambda, Fp_t nu, Point_wp_t R, Point_wp_t P, Aux_miller6_t mem) {
	Fp_sub(lambda, P->Y, R->Y); //lambda: = y2 - y1;
	Fp_sub(mem->t0, P->X, R->X); //c: = x2 - x1;
	Fp_modinv(mem->t0, mem->t0); //c: = 1 / c;
	Fp_mul(lambda, lambda, mem->t0); //lambda: = lambda*c;
	Fp_mul(mem->t1, lambda, lambda); //sqr: = lambda ^ 2;
	Fp_mul(nu, R->X, lambda); //nu: = x1*lambda;
	Fp_sub(R->X, mem->t1, R->X); //x3: = sqr - x1;
	Fp_sub(R->X, R->X, P->X); //x3: = x3 - x2;
	Fp_sub(nu, R->Y, nu); //nu: = y1 - nu;
	Fp_mul(R->Y, lambda, R->X); //y3: = lambda*x3;
	Fp_add(R->Y, R->Y, nu); //y3: = y3 + nu;
	Fp_neg(R->Y, R->Y); //y3: = -y3;
}

void cp6_line_affadd(Fp_t lambda, Fp_t nu, Point_wp_t R, Point_wp_t P) {
	cp6_line_affadd_mem(lambda, nu, R, P, miller6_aux);
}


static void cp6_optimal_ate_miller_mem (Fp6q_t f, Point_wp_t Q, Point_wp_t P, Aux_miller6_t mem) {
  int i;
  
  Fp6q_set_zero(f);
  Fp_neg (mem->negxP, P->X);
  Fp_neg (mem->negyP, P->Y);
  Fp_wp_copy (mem->R, Q);
  Fp_wp_neg (mem->negQ, Q);

  /**************************/
  cp6_line_dbl (f->a0->a0, mem->l01, mem->l10, mem->R, cp6_config.bt3);
  Fp_mul (f->a1->a1, mem->l01, mem->negyP);
  Fp_mul (f->a0->a1, mem->l10, P->X);
  /*if (cp6_config.nafm[cp6_config.lenm - 2] == 1) {
    cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
    Fp_mul (mem->l01, mem->l01, P->Y);
    Fp_mul (mem->l10, mem->l10, mem->negxP);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01); 
  }
  if (cp6_config.nafm[cp6_config.lenm - 2] == -1) {
    cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->negQ);
    Fp_mul (mem->l01, mem->l01, P->Y);
    Fp_mul (mem->l10, mem->l10, mem->negxP);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01); 
  }*/
  for (i=cp6_config.lenm - 3; i>=0; i--) {
    Fp6q_squ (f, f);
    cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, cp6_config.bt3);
    Fp_mul (mem->l01, mem->l01, mem->negyP);
    Fp_mul (mem->l10, mem->l10, P->X);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    if (cp6_config.nafm[i] == 1) {
      cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    }
    if (cp6_config.nafm[i] == -1) {
      cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->negQ);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    }
  }
  Fp6q_squ (f, f);
  cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, cp6_config.bt3);
  Fp_mul (mem->l01, mem->l01, mem->negyP);
  Fp_mul (mem->l10, mem->l10, P->X);
  Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
  
  // We have now computed f_{2*x0,Q}(P)
  /**************************/
  Fp6q_copy (mem->f1, f);
  Fp6q_p3pow (mem->f1inv, mem->f1);
  Fp_wp_copy (mem->U, mem->R);
  Fp_wp_neg (mem->negU, mem->U); 
 
  for (i=cp6_config.lenm - 2; i>=0; i--) {
    Fp6q_squ (f, f);
    cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, cp6_config.bt3);
    Fp_mul (mem->l01, mem->l01, mem->negyP);
    Fp_mul (mem->l10, mem->l10, P->X);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    if (cp6_config.nafm[i] == 1) {
      cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->U);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
      Fp6q_mul (f, f, mem->f1);
            //printf ("f ="); Fp6q_print (f);
            //printf ("R3 ="); Fp_print (mem->R->Z);
    }
    if (cp6_config.nafm[i] == -1) {
      cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->negU);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
      Fp6q_mul (f, f, mem->f1inv);
            //printf ("f ="); Fp6q_print (f);
            //printf ("R3 ="); Fp_print (mem->R->Z);
    }
    //if (i == cp6_config.lenm-2) {printf ("f ="); Fp6q_print (f);}
  }    


  // We have now computed f_{2*x0^2,Q}(P)
  /**************************/
  Fp_wp_copy (mem->V, mem->R);
  
  Fp6q_squ (mem->F, f);
  cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->V, cp6_config.bt3);
  Fp_mul (mem->l01, mem->l01, mem->negyP);
  Fp_mul (mem->l10, mem->l10, P->X);
  Fp6q_mul_sparse_twist_lazy (mem->F, mem->F, mem->l00, mem->l10, mem->l01);
  Fp6q_mul (mem->F, mem->F, f);
  
  cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->V);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (mem->F, mem->F, mem->l00, mem->l10, mem->l01);
  // F = f_{6x^2,Q}(P)

  Fp6q_mul (mem->F, mem->F, mem->f1);
  cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->U);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (mem->F, mem->F, mem->l00, mem->l10, mem->l01);

  cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (mem->F, mem->F, mem->l00, mem->l10, mem->l01);
  //f_{6x^2+2x+1,Q}

  Fp6q_ppow (f, mem->F);
  Fp6q_mul (f, f, mem->f1);
            
  cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->U);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);

}

void cp6_optimal_ate_miller (Fp6q_t f, Point_wp_t Q, Point_wp_t P) {
  cp6_optimal_ate_miller_mem (f, Q, P, miller6_aux);
}



static void cp6b_optimal_ate_miller_mem (Fp6q_t f, Point_wp_t Q, Point_wp_t P, Aux_miller6_t mem) {
  int i;
  
  Fp6q_set_zero(f);
  Fp_neg (mem->negxP, P->X);
  Fp_neg (mem->negyP, P->Y);
  Fp_wp_copy (mem->R, Q);
  Fp_wp_neg (mem->negQ, Q);

  /**************************/
  cp6_line_dbl (f->a0->a0, mem->l01, mem->l10, mem->R, cp6_config.bt3);
  Fp_mul (f->a1->a1, mem->l01, mem->negyP);
  Fp_mul (f->a0->a1, mem->l10, P->X);
  /*if (cp6_config.nafm[cp6_config.lenm - 2] == 1) {
    cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
    Fp_mul (mem->l01, mem->l01, P->Y);
    Fp_mul (mem->l10, mem->l10, mem->negxP);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01); 
  }
  if (cp6_config.nafm[cp6_config.lenm - 2] == -1) {
    cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->negQ);
    Fp_mul (mem->l01, mem->l01, P->Y);
    Fp_mul (mem->l10, mem->l10, mem->negxP);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01); 
  }*/
  for (i=cp6_config.lenm - 3; i>=0; i--) {
    Fp6q_squ (f, f);
    cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, cp6_config.bt3);
    Fp_mul (mem->l01, mem->l01, mem->negyP);
    Fp_mul (mem->l10, mem->l10, P->X);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    if (cp6_config.nafm[i] == 1) {
      cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    }
    if (cp6_config.nafm[i] == -1) {
      cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, mem->negQ);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    }
  }
  Fp6q_squ (f, f);
  cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, cp6_config.bt3);
  Fp_mul (mem->l01, mem->l01, mem->negyP);
  Fp_mul (mem->l10, mem->l10, P->X);
  Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
  
  // We have now computed f_{2*x0,Q}(P)
  /**************************/
  Fp6q_copy (mem->f1, f);
  Fp6q_p3pow (mem->f1inv, mem->f1);
  Fp_wp_copy (mem->U, mem->R);
  Fp_wp_neg (mem->negU, mem->U); 
 
  for (i=cp6_config.len3m - 2; i>=0; i--) {
    Fp6q_squ (f, f);
    cp6_line_dbl (mem->l00, mem->l01, mem->l10, mem->R, cp6_config.bt3);
    Fp_mul (mem->l01, mem->l01, mem->negyP);
    Fp_mul (mem->l10, mem->l10, P->X);
    Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
    if (cp6_config.naf3m[i] == 1) {
      cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->U);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
      Fp6q_mul (f, f, mem->f1);
            //printf ("f ="); Fp6q_print (f);
            //printf ("R3 ="); Fp_print (mem->R->Z);
    }
    if (cp6_config.naf3m[i] == -1) {
      cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->negU);
      Fp_mul (mem->l01, mem->l01, P->Y);
      Fp_mul (mem->l10, mem->l10, mem->negxP);
      Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
      Fp6q_mul (f, f, mem->f1inv);
            //printf ("f ="); Fp6q_print (f);
            //printf ("R3 ="); Fp_print (mem->R->Z);
    }
    //if (i == cp6_config.lenm-2) {printf ("f ="); Fp6q_print (f);}
  }    

  // We have now computed f_{6*x0^2,Q}(P)
  /**************************/
  
  Fp6q_mul (f, f, mem->f1inv);
  
  cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->negU);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);

  cp6_line_add (mem->l00, mem->l01, mem->l10, mem->R, Q);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
 
  Fp6q_ppow (f, f);
  cp6_line_add_proj (mem->l00, mem->l01, mem->l10, mem->R, mem->negU);
  Fp_mul (mem->l01, mem->l01, P->Y);
  Fp_mul (mem->l10, mem->l10, mem->negxP);
  Fp6q_mul_sparse_twist_lazy (f, f, mem->l00, mem->l10, mem->l01);
  Fp6q_mul (f, f, mem->f1inv);
}

void cp6b_optimal_ate_miller (Fp6q_t f, Point_wp_t Q, Point_wp_t P) {
  cp6b_optimal_ate_miller_mem (f, Q, P, miller6_aux);
}

static void cp6b_optimal_ate_miller_precompute_mem(Point_wp_t Q, uint64_t *premem, Aux_miller6_t mem) {
	int i, j;
	Fp_t l00, l01, l10;

	Fp_wp_copy(mem->R, Q);
	Fp_wp_neg(mem->negQ, Q);

	/**************************/
	j = 0;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	cp6_line_dbl(l00, l01, l10, mem->R, cp6_config.bt3);
	for (i = cp6_config.lenm - 3; i >= 0; i--) {
		j++;
		l00->limbs = premem + 3 * j * Fp_config.m->n;
		l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
		l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
		cp6_line_dbl(l00, l01, l10, mem->R, cp6_config.bt3);
		if (cp6_config.nafm[i] == 1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			cp6_line_add(l00, l01, l10, mem->R, Q);
		}
		if (cp6_config.nafm[i] == -1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			cp6_line_add(l00, l01, l10, mem->R, mem->negQ);
		}
	}
	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	cp6_line_dbl(l00, l01, l10, mem->R, cp6_config.bt3);

	// We have now computed the line coefficients up to f_{2*x0,Q}(P).
	/**************************/
	Fp_wp_copy(mem->U, mem->R);
	Fp_wp_neg(mem->negU, mem->U);

	for (i = cp6_config.len3m - 2; i >= 0; i--) {
		j++;
		l00->limbs = premem + 3 * j * Fp_config.m->n;
		l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
		l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
		cp6_line_dbl(l00, l01, l10, mem->R, cp6_config.bt3);
		if (cp6_config.naf3m[i] == 1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			cp6_line_add_proj(l00, l01, l10, mem->R, mem->U);
		}
		if (cp6_config.naf3m[i] == -1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			cp6_line_add_proj(l00, l01, l10, mem->R, mem->negU);
		}
	}

	// We have now computed the line coefficients up to f_{6*x0^2,Q}(P).
	/**************************/
	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	cp6_line_add_proj(l00, l01, l10, mem->R, mem->negU);

	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	cp6_line_add(l00, l01, l10, mem->R, Q);

	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	cp6_line_add_proj(l00, l01, l10, mem->R, mem->negU);
}

void cp6b_optimal_ate_miller_precompute(Point_wp_t Q, uint64_t *premem) {
	cp6b_optimal_ate_miller_precompute_mem(Q, premem, miller6_aux);
}


static void cp6b_optimal_ate_miller_useprecomputed_mem(Fp6q_t f, Point_wp_t P, uint64_t *premem, Aux_miller6_t mem) {
	int i, j;
	Fp_t l00, l01, l10;

	Fp6q_set_zero(f);
	Fp_neg(mem->negxP, P->X);
	Fp_neg(mem->negyP, P->Y);

	/**************************/
	j = 0;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	Fp_copy(f->a0->a0, l00);
	Fp_mul(f->a1->a1, l01, mem->negyP);
	Fp_mul(f->a0->a1, l10, P->X);
	for (i = cp6_config.lenm - 3; i >= 0; i--) {
		Fp6q_squ(f, f);
		j++;
		l00->limbs = premem + 3 * j * Fp_config.m->n;
		l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
		l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
		Fp_mul(l01, l01, mem->negyP);
		Fp_mul(l10, l10, P->X);
		Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
		if (cp6_config.nafm[i] == 1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			Fp_mul(l01, l01, P->Y);
			Fp_mul(l10, l10, mem->negxP);
			Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
		}
		if (cp6_config.nafm[i] == -1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			Fp_mul(l01, l01, P->Y);
			Fp_mul(l10, l10, mem->negxP);
			Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
		}
	}
	Fp6q_squ(f, f);
	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	Fp_mul(l01, l01, mem->negyP);
	Fp_mul(l10, l10, P->X);
	Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);

	// We have now computed f_{2*x0,Q}(P).
	/**************************/
	Fp6q_copy(mem->f1, f);
	Fp6q_p3pow(mem->f1inv, mem->f1);
	Fp_wp_copy(mem->U, mem->R);
	Fp_wp_neg(mem->negU, mem->U);

	for (i = cp6_config.len3m - 2; i >= 0; i--) {
		Fp6q_squ(f, f);
		j++;
		l00->limbs = premem + 3 * j * Fp_config.m->n;
		l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
		l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
		Fp_mul(l01, l01, mem->negyP);
		Fp_mul(l10, l10, P->X);
		Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
		if (cp6_config.naf3m[i] == 1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			Fp_mul(l01, l01, P->Y);
			Fp_mul(l10, l10, mem->negxP);
			Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
			Fp6q_mul(f, f, mem->f1);
		}
		if (cp6_config.naf3m[i] == -1) {
			j++;
			l00->limbs = premem + 3 * j * Fp_config.m->n;
			l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
			l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
			Fp_mul(l01, l01, P->Y);
			Fp_mul(l10, l10, mem->negxP);
			Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
			Fp6q_mul(f, f, mem->f1inv);
		}
	}

	// We have now computed f_{6*x0^2,Q}(P).
	/**************************/

	Fp6q_mul(f, f, mem->f1inv);

	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	Fp_mul(l01, l01, P->Y);
	Fp_mul(l10, l10, mem->negxP);
	Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);

	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	Fp_mul(l01, l01, P->Y);
	Fp_mul(l10, l10, mem->negxP);
	Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);

	Fp6q_ppow(f, f);
	j++;
	l00->limbs = premem + 3 * j * Fp_config.m->n;
	l01->limbs = premem + (3 * j + 1) * Fp_config.m->n;
	l10->limbs = premem + (3 * j + 2) * Fp_config.m->n;
	Fp_mul(l01, l01, P->Y);
	Fp_mul(l10, l10, mem->negxP);
	Fp6q_mul_sparse_twist_lazy(f, f, l00, l10, l01);
	Fp6q_mul(f, f, mem->f1inv);
}

void cp6b_optimal_ate_miller_useprecomputed(Fp6q_t f, Point_wp_t P, uint64_t *premem) {
	cp6b_optimal_ate_miller_useprecomputed_mem(f, P, premem, miller6_aux);
}


static void cp6b_optimal_ate_miller_aff_mem(Fp6q_t f, Point_wp_t Q, Point_wp_t P, Aux_miller6_t mem) {
	int i;

	Fp6q_set_one(f);
	//Fp_neg(mem->negxP, P->X);
	Fp_neg(mem->negyP, P->Y);
	Fp_wp_copy(mem->R, Q);
	Fp_wp_neg(mem->negQ, Q);

	///**************************/
	cp6_line_affdbl(mem->l00, f->a0->a0, mem->R);
	Fp_mul(f->a0->a1, mem->l00, P->X);
	Fp_copy(f->a1->a1, mem->negyP);
	for (i = cp6_config.lenm - 3; i >= 0; i--) {
		Fp6q_squ(f, f);
		cp6_line_affdbl(mem->l00, mem->l10, mem->R);
		Fp_mul(mem->l00, mem->l00, P->X);
		Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
		if (cp6_config.nafm[i] == 1) {
			cp6_line_affadd(mem->l00, mem->l10, mem->R, Q);
			Fp_mul(mem->l00, mem->l00, P->X);
			Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
		}
		if (cp6_config.nafm[i] == -1) {
			cp6_line_affadd(mem->l00, mem->l10, mem->R, mem->negQ);
			Fp_mul(mem->l00, mem->l00, P->X);
			Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
		}
	}
	Fp6q_squ(f, f);
	cp6_line_affdbl(mem->l00, mem->l10, mem->R);
	Fp_mul(mem->l00, mem->l00, P->X);
	Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);

	// We have now computed f_{2*x0,Q}(P)
	/**************************/
	Fp6q_copy(mem->f1, f);
	Fp6q_p3pow(mem->f1inv, mem->f1);
	Fp_wp_copy(mem->U, mem->R);
	Fp_wp_neg(mem->negU, mem->U);

	for (i = cp6_config.len3m - 2; i >= 0; i--) {
		Fp6q_squ(f, f);
		cp6_line_affdbl(mem->l00, mem->l10, mem->R);
		Fp_mul(mem->l00, mem->l00, P->X);
		Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
		if (cp6_config.naf3m[i] == 1) {
			cp6_line_affadd(mem->l00, mem->l10, mem->R, mem->U);
			Fp_mul(mem->l00, mem->l00, P->X);
			Fp6q_mul(f, f, mem->f1);
			Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
		}
		if (cp6_config.naf3m[i] == -1) {
			cp6_line_affadd(mem->l00, mem->l10, mem->R, mem->negU);
			Fp_mul(mem->l00, mem->l00, P->X);
			Fp6q_mul(f, f, mem->f1inv);
			Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
		}
		//if (i == cp6_config.lenm-2) {printf ("f ="); Fp6q_print (f);}
	}

	// We have now computed f_{6*x0^2,Q}(P)
	/**************************/

	Fp6q_mul(f, f, mem->f1inv);

	cp6_line_affadd(mem->l00, mem->l10, mem->R, mem->negU);
	Fp_mul(mem->l00, mem->l00, P->X);
	Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);

	cp6_line_affadd(mem->l00, mem->l10, mem->R, Q);
	Fp_mul(mem->l00, mem->l00, P->X);
	Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);

	Fp6q_ppow(f, f);
	cp6_line_affadd(mem->l00, mem->l10, mem->R, mem->negU);
	Fp_mul(mem->l00, mem->l00, P->X);
	Fp6q_mul_sparse_twist_lazy(f, f, mem->l10, mem->l00, mem->negyP);
	Fp6q_mul(f, f, mem->f1inv);
}

void cp6b_optimal_ate_miller_aff(Fp6q_t f, Point_wp_t Q, Point_wp_t P) {
	cp6b_optimal_ate_miller_aff_mem(f, Q, P, miller6_aux);
}

/********************************************************************************************/



/* Line function for a doubling step on the CP3 curve.
* The coefficients returned by this function are the coefficients L00, L01, L10, L20 of the line
* l := L00 + L01*yP + L10*xP + L20*xP^2.
*/
static void cp3_line_dbl_mem(Fp_t l00, Fp_t l01, Fp_t l10, Fp_t l20, Point_wp_t R, Fp_t b, Aux_miller3_t mem) {
	Fp_mul(mem->t0, R->X, R->X);        // t0 = A:=X1^2;
	Fp_mul(mem->t1, R->Y, R->Y);        // t1 = B:=Y1^2; 
	Fp_mul(mem->t2, R->Z, R->Z);        // t2 = C:=Z1^2;
	Fp_mul(mem->t3, mem->t2, b);        // t3 = D:=C*b;                  		
	Fp_add(mem->t4, mem->t3, mem->t3);
	Fp_add(mem->t4, mem->t4, mem->t3);  // t4 = E:=3*b*Z1^2
	Fp_mul(mem->t5, R->X, R->Y);        // t5 = F:=X1*Y1;
	Fp_add(mem->t5, mem->t5, mem->t5);  // t5 = F:=2*X1*Y1;
	Fp_mul(mem->t6, R->Y, R->Z);        // t6 = G:=Y1*Z1;
	Fp_add(mem->t6, mem->t6, mem->t6);  // t6 = G:=2*Y1*Z1;
	Fp_add(mem->t2, mem->t4, mem->t4);  // t2 = C:=2*E;
	Fp_add(mem->t7, mem->t2, mem->t4);  // t7 = H:=C+E;
	Fp_sub(R->X, mem->t1, mem->t7);     // X3:=B-H;
	Fp_mul(R->X, mem->t5, R->X);        // X3:=F*X3; 
	Fp_add(R->Y, mem->t1, mem->t7);     // Y3:=B+H;
	Fp_mul(R->Y, R->Y, R->Y);           // Y3:=Y3^2;  
	Fp_mul(mem->t2, mem->t2, mem->t2);  // t2 = C:=C^2;
	Fp_add(mem->t4, mem->t2, mem->t2);  // t4 = E:=2*C;
	Fp_add(mem->t2, mem->t2, mem->t4);  // t2 = C:=C+E;
	Fp_sub(R->Y, R->Y, mem->t2);        // Y3:=Y3-C;  
	Fp_mul(R->Z, mem->t1, mem->t6);     // Z3:=B*G;
	Fp_add(R->Z, R->Z, R->Z);           // Z3:=2*Z3;
	Fp_add(R->Z, R->Z, R->Z);           // Z3:=2*Z3; 
	Fp_sub(l10, mem->t1, mem->t7);      // L10:=B-H;
	Fp_mul(l10, l10, mem->t0);          // L10:=L10*A;
	Fp_mul(l20, mem->t5, mem->t6);      // L20:=F*G;
	Fp_mul(l01, mem->t0, mem->t5);      // L01:=A*F;
	Fp_add(mem->t4, l01, l01);          // t4 = E:=2*L01;
	Fp_add(l01, l01, mem->t4);          // L01:=L01+E;
	//Fp_neg(l01, l01);                   // L01:=-L01;  leave out, instead negate yP in the Miller function
	Fp_sub(mem->t4, mem->t1, mem->t3);  // t4 = E:=B-D;
	Fp_add(l00, mem->t1, mem->t7);      // L00:=B+H;
	Fp_mul(l00, l00, mem->t4);          // L00:= L00*E;
}

void cp3_line_dbl(Fp_t l00, Fp_t l01, Fp_t l10, Fp_t l20, Point_wp_t R, Fp_t b) {
	cp3_line_dbl_mem(l00, l01, l10, l20, R, b, miller3_aux);
}

/* Line function for an addition step on the CP3 curve.
* The coefficients returned by this function are the coefficients L00, L01, L10, L20 of the line
* l := L00 + L01*yP + L10*xP + L20*xP^2.
*/
static void cp3_line_add_mem(Fp_t l00, Fp_t l01, Fp_t l10, Fp_t l20, Point_wp_t R, Point_wp_t Q, Aux_miller3_t mem) {
	Fp_mul(mem->t0, R->Z, Q->X);       // t0 = D:=Z1*X2; 
	Fp_sub(mem->t0, mem->t0, R->X);    // t0 = D:=D-X1; 
	Fp_mul(mem->t1, R->Z, Q->Y);       // t1 = E:=Z1*Y2; 
	Fp_sub(mem->t1, R->Y, mem->t1);    // t1 = E:=Y1-E; 
	Fp_mul(mem->t2, mem->t0, mem->t0); // t2 = F:=D^2;
	Fp_mul(mem->t3, mem->t1, mem->t1); // t3 = G:=E^2;
	Fp_mul(mem->t4, mem->t0, mem->t2);
	Fp_neg(mem->t4, mem->t4);          // t4 = H:=-D*F;
	Fp_mul(mem->t5, mem->t2, R->X);    // t5 = I:=F*X1;
	Fp_mul(mem->t6, mem->t3, R->Z);    // t6 = J:=Z1*G;
	Fp_add(mem->t6, mem->t6, mem->t4); // t6 = J:=J+H;
	Fp_sub(mem->t6, mem->t6, mem->t5); // t6 = J:=J-I;
	Fp_sub(mem->t6, mem->t5, mem->t6); // t6 = J:=I-J;
	Fp_mul(mem->t7, mem->t2, R->Z);    // t7 = K:=Z1*F;
	Fp_mul(mem->t7, mem->t7, mem->t1); // t7 = K:=K*E;
	Fp_mul(R->X, mem->t0, mem->t6);    // X3:=D*J;
	Fp_mul(R->Y, R->Y, mem->t4);       // Y3:=Y1*H;
	Fp_add(mem->t5, mem->t5, mem->t6); // t5 = I:=I+J;
	Fp_mul(mem->t5, mem->t1, mem->t5); // t5 = I:=E*I;
	Fp_sub(R->Y, mem->t5, R->Y);       // Y3:=I-Y3;
	Fp_mul(R->Z, R->Z, mem->t4);       // Z3:=Z1*H;
	Fp_mul(mem->t0, R->X, R->X);       // t0 = L:=X3^2;
	Fp_mul(mem->t1, R->Z, R->Z);       // t1 = M:=Z3^2;
	Fp_mul(l10, R->X, R->Z);           // L10:=X3*Z3;
	Fp_add(l10, l10, l10);             // L10:=2*L10;
	Fp_add(l20, mem->t1, mem->t1);     // L20;=2*M;
	Fp_mul(l00, mem->t7, R->Y);        // L00:=K*Y3;
	Fp_add(l00, l00, mem->t0);         // L00:=L00+L;
	Fp_add(l00, l00, l00);             // L00:=2*L00;
	Fp_mul(l01, mem->t7, R->Z);        // L01:=K*Z3;
	Fp_add(l01, l01, l01);             // L01:=2*L01;
	//Fp_neg(l01, l01);                  // L01:=-L01; leave out, instead negate yP in the Miller function
}

void cp3_line_add(Fp_t l00, Fp_t l01, Fp_t l10, Fp_t l20, Point_wp_t R, Point_wp_t Q) {
	cp3_line_add_mem(l00, l01, l10, l20, R, Q, miller3_aux);
}

static void cp3_ate_miller_mem(Fp3_t f, Point_wp_t Q, Point_wp_t P, Aux_miller3_t mem) {
	int i;

	Fp3_set_one(f);         // f = 1
	Fp_wp_copy(mem->R, Q);  // R = Q

	Fp_mul(mem->Ptmp->X, P->X, cp3_config.c1);
	Fp_neg(mem->negyP, P->Y);
	Fp_mul(mem->Ptmp->Y, mem->negyP, cp3_config.s11);
	Fp_mul(mem->xP2, mem->Ptmp->X, mem->Ptmp->X);

	//printf("\n\nPX = "); Fp_print(mem->Ptmp->X);
	//printf("\n\nPY = "); Fp_print(mem->Ptmp->Y);
	//printf("\n\n");

	//printf("\nPoint R before: "); Fp_wp_print(mem->R);

	cp3_line_dbl(mem->l00, mem->l01, mem->l10, mem->l20, mem->R, cp3_config.bt);
	Fp_mul(mem->f->a0, mem->l01, mem->Ptmp->Y);
	Fp_add(mem->f->a0, mem->f->a0, mem->l00);
	Fp_mul(mem->f->a1, mem->l10, mem->Ptmp->X);
	Fp_mul(mem->f->a2, mem->l20, mem->xP2);
	Fp3_copy(f, mem->f);

	//printf("\n\nf_1 = "); Fp3_print(f); printf("\n\n");
	//printf("\nPoint R: "); Fp_wp_print(mem->R);


	/**************************/
	for (i = cp3_config.lenT - 3; i >= 0; i--) {
		Fp3_squ(f, f);
		cp3_line_dbl(mem->l00, mem->l01, mem->l10, mem->l20, mem->R, cp3_config.bt);
		Fp_mul(mem->f->a0, mem->l01, mem->Ptmp->Y);
		Fp_add(mem->f->a0, mem->f->a0, mem->l00);
		Fp_mul(mem->f->a1, mem->l10, mem->Ptmp->X);
		Fp_mul(mem->f->a2, mem->l20, mem->xP2);
		Fp3_mul(f, f, mem->f);
		if (BIT(cp3_config.T,i) == 1) {
			cp3_line_add(mem->l00, mem->l01, mem->l10, mem->l20, mem->R, Q);
			Fp_mul(mem->f->a0, mem->l01, mem->Ptmp->Y);
			Fp_add(mem->f->a0, mem->f->a0, mem->l00);
			Fp_mul(mem->f->a1, mem->l10, mem->Ptmp->X);
			Fp_mul(mem->f->a2, mem->l20, mem->xP2);
			Fp3_mul(f, f, mem->f);
		}
	}
}

void cp3_ate_miller(Fp3_t f, Point_wp_t Q, Point_wp_t P) {
	cp3_ate_miller_mem(f, Q, P, miller3_aux);
}

/********************************************************************************************/

int bn_bls12_miller_initialize_config (Fp2_t bt) {
  int i, ret;

  if ((ret=Fp2_init (miller12_aux->t0)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t1)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t2)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t3)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t4)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t5)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t6)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t7)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->t8)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->l00)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->l01)) < 0) return ret;
  if ((ret=Fp2_init (miller12_aux->l10)) < 0) return ret;
  if ((ret=Fp_init (miller12_aux->negxP)) < 0) return ret;
  if ((ret=Fp_init (miller12_aux->negyP)) < 0) return ret;
  if ((ret = Fp2_init(miller12_aux->negyP_Fp2)) < 0) return ret;
  if ((ret=Fp2_wp_init (miller12_aux->R)) < 0) return ret;
  if ((ret=Fp2_wp_init (miller12_aux->negQ)) < 0) return ret;
  if ((ret=Fp2_wp_init (miller12_aux->Qp)) < 0) return ret;
  if ((ret=Fp2_wp_init (miller12_aux->Qp2)) < 0) return ret;
  if ((ret=Fp12_init (miller12_aux->f)) < 0) return ret;

  Fp12_set_zero (miller12_aux->f);

  if (PAIR_CURVE == BN12) {
    if (Fp2_init(bn_config.bt3) < 0) printf ("bn_config memory error.\n");
    Fp2_add (bn_config.bt3, bt, bt);
    Fp2_add (bn_config.bt3, bn_config.bt3, bt);

    bn_config.lenu = 63;
    bn_config.nafu = NULL;
    bn_config.nafu = (int *) malloc (bn_config.lenu  * sizeof (int));
    if (bn_config.nafu == NULL) return ERR_OUT_OF_MEMORY;

    for (i = 0; i < bn_config.lenu; i++) bn_config.nafu[i] = 0;
    bn_config.nafu[0] = 1;
    bn_config.nafu[55] = 1;
    bn_config.nafu[62] = 1;
  
    bn_config.signu = -1;

    bn_config.len6u2 = 66;
    bn_config.naf6u2 = NULL;
    bn_config.naf6u2 = (int *) malloc (bn_config.len6u2  * sizeof (int));
    if (bn_config.naf6u2 == NULL) return ERR_OUT_OF_MEMORY;
  
    for (i = 0; i < bn_config.len6u2; i++) bn_config.naf6u2[i] = 0;
    bn_config.naf6u2[2] = 1;
    bn_config.naf6u2[56] = -1;
    bn_config.naf6u2[58] = 1;
    bn_config.naf6u2[63] = -1;
    bn_config.naf6u2[65] = 1;

	bn_config.nafwt6u2 = 5;

	miller12_aux->premem = (uint64_t *)malloc(4*Fp_config.m->n * (bn_config.len6u2 + bn_config.nafwt6u2) * sizeof (uint64_t));
	if (miller12_aux->premem == NULL) return ERR_OUT_OF_MEMORY;
  } else if (PAIR_CURVE == BN12CP) {
    if (Fp2_init(bn_config.bt3) < 0) printf ("bn_config memory error.\n");
    Fp2_add (bn_config.bt3, bt, bt);
    Fp2_add (bn_config.bt3, bn_config.bt3, bt);

    bn_config.lenu = 63;
    bn_config.nafu = NULL;
    bn_config.nafu = (int *) malloc (bn_config.lenu  * sizeof (int));
    if (bn_config.nafu == NULL) return ERR_OUT_OF_MEMORY;

    for (i = 0; i < bn_config.lenu; i++) bn_config.nafu[i] = 0;
    bn_config.nafu[0] = -1;
    bn_config.nafu[15] = 1;
    bn_config.nafu[34] = 1;
    bn_config.nafu[37] = 1;
    bn_config.nafu[62] = 1;
  
    bn_config.signu = 1;

    bn_config.len6u2 = 66;
    bn_config.naf6u2 = NULL;
    bn_config.naf6u2 = (int *) malloc (bn_config.len6u2  * sizeof (int));
    if (bn_config.naf6u2 == NULL) return ERR_OUT_OF_MEMORY;
  
    for (i = 0; i < bn_config.len6u2; i++) bn_config.naf6u2[i] = 0;
    bn_config.naf6u2[2] = -1;
    bn_config.naf6u2[16] = -1;
    bn_config.naf6u2[18] = 1;
    bn_config.naf6u2[35] = -1;
    bn_config.naf6u2[37] = -1;
    bn_config.naf6u2[40] = 1;
    bn_config.naf6u2[63] = -1;
    bn_config.naf6u2[65] = 1;
  } 
  else if (PAIR_CURVE == BN12tiny) {
	  if (Fp2_init(bn_config.bt3) < 0) printf("bn_config memory error.\n");
	  Fp2_add(bn_config.bt3, bt, bt);
	  Fp2_add(bn_config.bt3, bn_config.bt3, bt);

	  bn_config.lenu = 13;
	  bn_config.nafu = NULL;
	  bn_config.nafu = (int *)malloc(bn_config.lenu  * sizeof (int));
	  if (bn_config.nafu == NULL) return ERR_OUT_OF_MEMORY;

	  for (i = 0; i < bn_config.lenu; i++) bn_config.nafu[i] = 0;
	  bn_config.nafu[0] = 1;
	  bn_config.nafu[2] = -1;
	  bn_config.nafu[5] = 1;
	  bn_config.nafu[9] = -1;
	  bn_config.nafu[12] = 1;
	  
	  bn_config.signu = -1;

	  bn_config.len6u2 = 15;
	  bn_config.naf6u2 = NULL;
	  bn_config.naf6u2 = (int *)malloc(bn_config.len6u2  * sizeof (int));
	  if (bn_config.naf6u2 == NULL) return ERR_OUT_OF_MEMORY;

	  for (i = 0; i < bn_config.len6u2; i++) bn_config.naf6u2[i] = 0;
	  bn_config.naf6u2[2] = -1;
	  bn_config.naf6u2[4] = -1;
	  bn_config.naf6u2[6] = -1;
	  bn_config.naf6u2[8] = 1;
	  bn_config.naf6u2[10] = 1;
	  bn_config.naf6u2[12] = 1;
	  bn_config.naf6u2[14] = 1;
  } 
  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }
  return ERR_SUCCESS;
}

void bn_bls12_miller_free_config (void) {
  Fp2_free (miller12_aux->t0);
  Fp2_free (miller12_aux->t1);
  Fp2_free (miller12_aux->t2);
  Fp2_free (miller12_aux->t3);
  Fp2_free (miller12_aux->t4);
  Fp2_free (miller12_aux->t5);
  Fp2_free (miller12_aux->t6);
  Fp2_free (miller12_aux->t7);
  Fp2_free (miller12_aux->t8);
  Fp2_free (miller12_aux->l00);
  Fp2_free (miller12_aux->l01);
  Fp2_free (miller12_aux->l10);
  Fp_free (miller12_aux->negxP);
  Fp_free (miller12_aux->negyP);
  Fp2_free(miller12_aux->negyP_Fp2);
  Fp2_wp_free (miller12_aux->R);
  Fp2_wp_free (miller12_aux->negQ);
  Fp2_wp_free (miller12_aux->Qp);
  Fp2_wp_free (miller12_aux->Qp2);
  Fp12_free (miller12_aux->f);
  free(miller12_aux->premem);
  miller12_aux->premem = NULL;

  if (PAIR_CURVE == BN12) {
    free (bn_config.nafu);
    free (bn_config.naf6u2);
    Fp2_free (bn_config.bt3);
  } 
}


int cp6_miller_initialize_config(Fp_t bt) {
	int i, ret;

	if ((ret = Fp_init(miller6_aux->t0)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t1)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t2)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t3)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t4)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t5)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t6)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t7)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->t8)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->l00)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->l01)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->l10)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->negxP)) < 0) return ret;
	if ((ret = Fp_init(miller6_aux->negyP)) < 0) return ret;
	if ((ret = Fp_wp_init(miller6_aux->R)) < 0) return ret;
	if ((ret = Fp_wp_init(miller6_aux->U)) < 0) return ret;
	if ((ret = Fp_wp_init(miller6_aux->V)) < 0) return ret;
	if ((ret = Fp_wp_init(miller6_aux->negQ)) < 0) return ret;
	if ((ret = Fp_wp_init(miller6_aux->negU)) < 0) return ret;
	if ((ret = Fp_wp_init(miller6_aux->negV)) < 0) return ret;
	if ((ret = Fp6q_init(miller6_aux->f)) < 0) return ret;
	if ((ret = Fp6q_init(miller6_aux->f1)) < 0) return ret;
	if ((ret = Fp6q_init(miller6_aux->f1inv)) < 0) return ret;
	if ((ret = Fp6q_init(miller6_aux->F)) < 0) return ret;

	Fp6q_set_zero(miller6_aux->f);
	Fp6q_set_zero(miller6_aux->f1);
	Fp6q_set_zero(miller6_aux->f1inv);
	Fp6q_set_zero(miller6_aux->F);

	if (Fp_init(cp6_config.bt3) < 0) printf("cp6_config memory error.\n");
	Fp_add(cp6_config.bt3, bt, bt);
	Fp_add(cp6_config.bt3, cp6_config.bt3, bt);

	if (PAIR_CURVE == CP6tiny) {
		cp6_config.lenm = 13;
	}
	else {
		cp6_config.lenm = 63;
	}
	cp6_config.nafm = NULL;
	cp6_config.nafm = (int *)malloc(cp6_config.lenm  * sizeof (int));
	if (PAIR_CURVE == CP6tiny) {
		cp6_config.len3m = 14;
	}
	else {
		cp6_config.len3m = 65;
	}
	cp6_config.naf3m = NULL;
	cp6_config.naf3m = (int *)malloc(cp6_config.len3m  * sizeof (int));
	if (cp6_config.naf3m == NULL) return ERR_OUT_OF_MEMORY;

	if (PAIR_CURVE == CP6) {
		for (i = 0; i < cp6_config.lenm; i++) cp6_config.nafm[i] = 0;
		cp6_config.nafm[0] = -1;
		cp6_config.nafm[15] = 1;
		cp6_config.nafm[34] = 1;
		cp6_config.nafm[37] = 1;
		cp6_config.nafm[62] = 1;
	}
	else if (PAIR_CURVE == CP6b) {
		for (i = 0; i < cp6_config.lenm; i++) cp6_config.nafm[i] = 0;
		cp6_config.nafm[0] = 1;
		cp6_config.nafm[55] = 1;
		cp6_config.nafm[62] = 1;
		cp6_config.nafwtm = 3;
		for (i = 0; i < cp6_config.len3m; i++) cp6_config.naf3m[i] = 0;
		cp6_config.naf3m[0] = -1;
		cp6_config.naf3m[2] = 1;
		cp6_config.naf3m[55] = -1;
		cp6_config.naf3m[57] = 1;
		cp6_config.naf3m[62] = -1;
		cp6_config.naf3m[64] = 1;
		cp6_config.nafwt3m = 6;
	}
	else if (PAIR_CURVE == CP6tiny) {
		for (i = 0; i < cp6_config.lenm; i++) cp6_config.nafm[i] = 0;
		cp6_config.nafm[0] = 1;
		cp6_config.nafm[2] = -1;
		cp6_config.nafm[5] = 1;
		cp6_config.nafm[9] = -1;
		cp6_config.nafm[12] = 1;
		for (i = 0; i < cp6_config.len3m; i++) cp6_config.naf3m[i] = 0;
		cp6_config.naf3m[0] = -1;
		cp6_config.naf3m[3] = -1;
		cp6_config.naf3m[5] = -1;
		cp6_config.naf3m[7] = 1;
		cp6_config.naf3m[9] = 1;
		cp6_config.naf3m[11] = 1;
		cp6_config.naf3m[13] = 1;
	}
	else {
		printf("ERROR: PAIR_CURVE not valid!");
	}
	return ERR_SUCCESS;
}

void cp6_miller_free_config(void) {
	Fp_free(miller6_aux->t0);
	Fp_free(miller6_aux->t1);
	Fp_free(miller6_aux->t2);
	Fp_free(miller6_aux->t3);
	Fp_free(miller6_aux->t4);
	Fp_free(miller6_aux->t5);
	Fp_free(miller6_aux->t6);
	Fp_free(miller6_aux->t7);
	Fp_free(miller6_aux->t8);
	Fp_free(miller6_aux->l00);
	Fp_free(miller6_aux->l01);
	Fp_free(miller6_aux->l10);
	Fp_free(miller6_aux->negxP);
	Fp_free(miller6_aux->negyP);
	Fp_wp_free(miller6_aux->R);
	Fp_wp_free(miller6_aux->U);
	Fp_wp_free(miller6_aux->V);
	Fp_wp_free(miller6_aux->negQ);
	Fp_wp_free(miller6_aux->negU);
	Fp_wp_free(miller6_aux->negV);
	Fp6q_free(miller6_aux->f);
	Fp6q_free(miller6_aux->f1);
	Fp6q_free(miller6_aux->f1inv);
	Fp6q_free(miller6_aux->F);

	free(cp6_config.nafm);
  free(cp6_config.naf3m);
	Fp_free(cp6_config.bt3);
}


int cp3_miller_initialize_config(Fp_t bt) {
	int i, ret;

	uint64_t c1[16] = { 0x89386F3B36718579, 0x5738A7002012617C, 0xFD1C20225D1BE84A, 0xFF0F4822A516FC21,
		0x07AB0721353DCE86, 0xE3E1CA50E454A0B8, 0xF1B1375F12CFA0C3, 0x6F71BCA78681D410,
		0xF10DC78679A4D5FD, 0x0DF370FC1504A859, 0x56D5BFBD0C409D3D, 0x7027471F4A4D6696,
		0xAD5A99D7A97100AF, 0x35D73AC0C1EEE801, 0xAB32BB084F2828CE, 0x085E112BF69C514A };
	uint64_t s11[16] = { 0x053EC93EE38148A6, 0x53EF981551816F15, 0x7A39C8665C4E8DDF, 0xFC2CE5068DE788CD,
		0x9B59619139DC1942, 0x73587A2BCBD2727D, 0x0CDB57661D362D34, 0xF793FD7B158DECB9,
		0xA5C98738D135AD38, 0x5AA5EFB7A9C14939, 0xDC4A7EE56ACB2CF8, 0xA16544F8086F6228,
		0xFEE2F49D23A9C745, 0xA3A6007A195B2EE6, 0x2C3A1EBDD66DE4BD, 0x135A18EBF53E09AD };
	uint64_t T[8] = { 0x3500000000002EF3, 0x931130000000609A, 0xD78537E800005ADF, 0x05754EF2480032D7,
		0xB98B526D59541277, 0x1813463651A21073, 0x8F48C426BFBAB862, 0xB1CB3E3B7C040CF7 };

	if ((ret = Fp_init(miller3_aux->t0)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t1)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t2)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t3)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t4)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t5)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t6)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t7)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->t8)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->l00)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->l01)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->l10)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->l20)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->xP2)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->negxP)) < 0) return ret;
	if ((ret = Fp_init(miller3_aux->negyP)) < 0) return ret;
	if ((ret = Fp_wp_init(miller3_aux->R)) < 0) return ret;
	if ((ret = Fp_wp_init(miller3_aux->Ptmp)) < 0) return ret;
	if ((ret = Fp_wp_init(miller3_aux->negQ)) < 0) return ret;
	if ((ret = Fp3_init(miller3_aux->f)) < 0) return ret;
	if ((ret = Fp3_init(miller3_aux->f1)) < 0) return ret;
	if ((ret = Fp3_init(miller3_aux->f1inv)) < 0) return ret;
	if ((ret = Fp3_init(miller3_aux->F)) < 0) return ret;

	Fp3_set_zero(miller3_aux->f);
	Fp3_set_zero(miller3_aux->f1);
	Fp3_set_zero(miller3_aux->f1inv);
	Fp3_set_zero(miller3_aux->F);

	if (Fp_init(cp3_config.bt) < 0) printf("cp3_config memory error.\n");
	Fp_copy(cp3_config.bt, bt);

	if (Fp_init(cp3_config.c1) < 0) printf("cp3_config memory error.\n");
	if (Fp_init(cp3_config.s11) < 0) printf("cp3_config memory error.\n");
	Fp_set(cp3_config.c1, c1);
	Fp_set(cp3_config.s11, s11);

	cp3_config.lenT = 512;
	cp3_config.T = NULL;
	cp3_config.T = (uint64_t *)malloc(8 * sizeof (uint64_t));
	for (i = 0; i < 8; i++) cp3_config.T[i] = T[i];

	return ERR_SUCCESS;
}

void cp3_miller_free_config(void) {
	Fp_free(miller3_aux->t0);
	Fp_free(miller3_aux->t1);
	Fp_free(miller3_aux->t2);
	Fp_free(miller3_aux->t3);
	Fp_free(miller3_aux->t4);
	Fp_free(miller3_aux->t5);
	Fp_free(miller3_aux->t6);
	Fp_free(miller3_aux->t7);
	Fp_free(miller3_aux->t8);
	Fp_free(miller3_aux->l00);
	Fp_free(miller3_aux->l01);
	Fp_free(miller3_aux->l10);
  Fp_free(miller3_aux->l20);
	Fp_free(miller3_aux->xP2);
	Fp_free(miller3_aux->negxP);
	Fp_free(miller3_aux->negyP);
	Fp_wp_free(miller3_aux->R);
	Fp_wp_free(miller3_aux->Ptmp);
	Fp_wp_free(miller3_aux->negQ);

	Fp3_free(miller3_aux->f);
	Fp3_free(miller3_aux->f1);
	Fp3_free(miller3_aux->f1inv);
	Fp3_free(miller3_aux->F);

  Fp_free(cp3_config.c1);
  Fp_free(cp3_config.s11);

	free(cp3_config.nafT);
	Fp_free(cp3_config.bt);

  free(cp3_config.T);
}

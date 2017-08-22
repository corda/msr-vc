/* pairing/finalexpo.c
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

typedef struct {
	Fp3_t t0;
	Fp3_t t1;
	Fp3_t t2;
	Fp3_t t3;
	Fp3_t t4;
	Fp3_t t5;
	Fp3_t t6;
	Fp3_t t7;
	Fp3_t t8;
	Fp3_t t9;
} aux_finalexpo3_t;

typedef aux_finalexpo3_t Aux_finalexpo3_t[1];

static Aux_finalexpo3_t finalexpo3_aux;


typedef struct {
  Fp6q_t t0;
  Fp6q_t t1;
  Fp6q_t t2;
  Fp6q_t t3;
  Fp6q_t t4;
  Fp6q_t t5;
  Fp6q_t t6;
  Fp6q_t t7;
  Fp6q_t t8;
  Fp6q_t t9;
  Fp6q_t tmp0;
  Fp6q_t tmp1;
} aux_finalexpo6_t;

typedef aux_finalexpo6_t Aux_finalexpo6_t[1];

static Aux_finalexpo6_t finalexpo6_aux;

typedef struct {
  Fp12_t t0;
  Fp12_t t1;
  Fp12_t t2;
  Fp12_t t3;
  Fp12_t t4;
  Fp12_t t5;
  Fp12_t t6;
  Fp12_t t7;
  Fp12_t t8;
  Fp12_t t9;
} aux_finalexpo12_t;

typedef aux_finalexpo12_t Aux_finalexpo12_t[1];

static Aux_finalexpo12_t finalexpo12_aux;


typedef struct {
  Fp3_t tt0;
  Fp3_t tt1;
  Fp3_t tt2;
  Fp_t tmp0;
  Fp_t tmp1;
} aux_cyclotomic6_t;

typedef aux_cyclotomic6_t Aux_cyclotomic6_t[1];

static Aux_cyclotomic6_t cyclotomic6_aux;

typedef struct {
  Fp6_t tt0;
  Fp6_t tt1;
  Fp6_t tt2;
} aux_cyclotomic12_t;

typedef aux_cyclotomic12_t Aux_cyclotomic12_t[1];

static Aux_cyclotomic12_t cyclotomic12_aux;




static void Fp12_squ_cyclotomic_mem (Fp12_t c, Fp12_t a, Aux_cyclotomic12_t mem) {
  Fp4_squ_sep (mem->tt0->a0, mem->tt1->a1, a->a0->a0, a->a1->a1);
  Fp4_squ_sep (mem->tt0->a1, mem->tt1->a2, a->a1->a0, a->a0->a2);
  Fp4_squ_sep (mem->tt0->a2, mem->tt1->a0, a->a0->a1, a->a1->a2);
  Fp2_mulxi (mem->tt1->a0, mem->tt1->a0);
  Fp6_add (mem->tt2, mem->tt0, mem->tt0);
  Fp6_add (mem->tt0, mem->tt2, mem->tt0);
  Fp6_add (mem->tt2, mem->tt1, mem->tt1);
  Fp6_add (mem->tt1, mem->tt2, mem->tt1);
  Fp6_add (c->a0, a->a0, a->a0);
  Fp6_neg (c->a0, c->a0);
  Fp6_add (c->a1, a->a1, a->a1);
  Fp6_add (c->a0, c->a0, mem->tt0);
  Fp6_add (c->a1, c->a1, mem->tt1);
}

void Fp12_squ_cyclotomic (Fp12_t c, Fp12_t a) {
  Fp12_squ_cyclotomic_mem (c, a, cyclotomic12_aux);
}

static void bn_exp_by_u (Fp12_t c, Fp12_t a, Fp12_t tmp0, Fp12_t tmp1) {
  int i;

  Fp12_copy (c, a);
  Fp12_copy (tmp0, a);
  Fp12_p6pow (tmp1, a);
  for (i = bn_config.lenu - 2; i>=0; i--) {
    Fp12_squ_cyclotomic(c, c);
    if (bn_config.nafu[i] == 1) Fp12_mul_lazy (c, c, tmp0);
    if (bn_config.nafu[i] == -1) Fp12_mul_lazy (c, c, tmp1);    
  }
}

static void bn_finalexpo_neg_mem (Fp12_t c, Fp12_t a, Aux_finalexpo12_t mem) {
  Fp12_inv (mem->t0, a);
  Fp12_p6pow (mem->t1, a);
  Fp12_mul_lazy (mem->t0, mem->t0, mem->t1); // t0 = a^(p^6-1)
  Fp12_p2pow (mem->t1, mem->t0);
  Fp12_mul_lazy (mem->t0, mem->t1, mem->t0); // t0 = (a^(p^6-1))^(p^2 + 1) = f

  bn_exp_by_u (mem->t1, mem->t0, mem->t5, mem->t6);
  Fp12_p6pow (mem->t1, mem->t1);      // fx
  Fp12_squ_cyclotomic (mem->t2, mem->t1);   // f2x
  Fp12_squ_cyclotomic (mem->t1, mem->t2);   // f4x
  Fp12_mul_lazy (mem->t1, mem->t2, mem->t1);// f6x
  bn_exp_by_u (mem->t3, mem->t1, mem->t5, mem->t6); 
  Fp12_p6pow (mem->t3, mem->t3); // f6x2
  Fp12_squ_cyclotomic (mem->t4, mem->t3); // f12x2
  bn_exp_by_u (mem->t4, mem->t4, mem->t5, mem->t6); 
  Fp12_p6pow (mem->t4, mem->t4);  // f12x3
  Fp12_mul_lazy (mem->t4, mem->t4, mem->t3); // a
  Fp12_mul_lazy (mem->t4, mem->t4, mem->t1); // a
  Fp12_p6pow (mem->t2, mem->t2); // b
  Fp12_mul_lazy (mem->t2, mem->t2, mem->t4); // b
  Fp12_mul_lazy (c, mem->t4, mem->t3); // fd
  Fp12_mul_lazy (c, c, mem->t0); // fd
  Fp12_ppow (mem->t1, mem->t2); // y0
  Fp12_p2pow (mem->t3, mem->t4); // y1
  Fp12_p6pow (mem->t0, mem->t0); // y2
  Fp12_mul_lazy (mem->t0, mem->t0, mem->t2); // y2
  Fp12_p3pow (mem->t0, mem->t0); // y2
  Fp12_mul_lazy (c, c, mem->t1);
  Fp12_mul_lazy (c, c, mem->t3);
  Fp12_mul_lazy (c, c, mem->t0);
}

void bn_finalexpo_neg (Fp12_t c, Fp12_t a) {
  bn_finalexpo_neg_mem (c, a, finalexpo12_aux);
}

static void bn_finalexpo_pos_mem (Fp12_t c, Fp12_t a, Aux_finalexpo12_t mem) {
  Fp12_inv (mem->t0, a);
  Fp12_p6pow (mem->t1, a);
  Fp12_mul_lazy (mem->t0, mem->t0, mem->t1); // t0 = a^(p^6-1)
  Fp12_p2pow (mem->t1, mem->t0);
  Fp12_mul_lazy (mem->t0, mem->t1, mem->t0); // t0 = (a^(p^6-1))^(p^2 + 1) = f

  bn_exp_by_u (mem->t1, mem->t0, mem->t5, mem->t6); // fx
  Fp12_squ_cyclotomic (mem->t2, mem->t1);   // f2x
  Fp12_squ_cyclotomic (mem->t1, mem->t2);   // f4x
  Fp12_mul_lazy (mem->t1, mem->t2, mem->t1);// f6x
  bn_exp_by_u (mem->t3, mem->t1, mem->t5, mem->t6); // f6x2
  Fp12_squ_cyclotomic (mem->t4, mem->t3); // f12x2
  bn_exp_by_u (mem->t4, mem->t4, mem->t5, mem->t6); // f12x3
  Fp12_mul_lazy (mem->t4, mem->t4, mem->t3); // a
  Fp12_mul_lazy (mem->t4, mem->t4, mem->t1); // a
  Fp12_p6pow (mem->t2, mem->t2); // b
  Fp12_mul_lazy (mem->t2, mem->t2, mem->t4); // b
  Fp12_mul_lazy (c, mem->t4, mem->t3); // fd
  Fp12_mul_lazy (c, c, mem->t0); // fd
  Fp12_ppow (mem->t1, mem->t2); // y0
  Fp12_p2pow (mem->t3, mem->t4); // y1
  Fp12_p6pow (mem->t0, mem->t0); // y2
  Fp12_mul_lazy (mem->t0, mem->t0, mem->t2); // y2
  Fp12_p3pow (mem->t0, mem->t0); // y2
  Fp12_mul_lazy (c, c, mem->t1);
  Fp12_mul_lazy (c, c, mem->t3);
  Fp12_mul_lazy (c, c, mem->t0);
}

void bn_finalexpo_pos (Fp12_t c, Fp12_t a) {
  bn_finalexpo_pos_mem (c, a, finalexpo12_aux);
}


/* Square an Fp2 element given by coefficients separately, return coefficients separately. */
static void cp6_Fp2_squ_sep (Fp_t c0, Fp_t c1, Fp_t a0, Fp_t a1, Fp_t t0, Fp_t t1) {	
	Fp_mul (t0, a0, a0);         // t0 = a0^2
  Fp_mul (t1, a1, a1);         // t1 = a1^2
  Fp_add (c1, a0, a1);         // c1 = a0 + a1
  Fp_add (c0, t1, t1);
  Fp_add (c0, c0, t1);         // c0 = 3*a1^2
  Fp_add (c0, t0, c0);         // c0 = a0^2 + 3*a1^2
  Fp_mul (c1, c1, c1);         // (a0 + a1)^2 = a0^2 + 2*a0*a1 + a1^2
  Fp_sub (c1, c1, t0);         // 2*a0*a1 + a1^2
  Fp_sub (c1, c1, t1);         // 2*a0*a1
}

static void Fp6q_squ_cyclotomic_mem (Fp6q_t c, Fp6q_t a, Aux_cyclotomic6_t mem) {
  cp6_Fp2_squ_sep (mem->tt0->a0, mem->tt1->a1, a->a0->a0, a->a1->a1, mem->tmp0, mem->tmp1);
  cp6_Fp2_squ_sep (mem->tt0->a1, mem->tt1->a2, a->a1->a0, a->a0->a2, mem->tmp0, mem->tmp1);
  cp6_Fp2_squ_sep (mem->tt0->a2, mem->tt1->a0, a->a0->a1, a->a1->a2, mem->tmp0, mem->tmp1);
  Fp_add (mem->tmp0, mem->tt1->a0, mem->tt1->a0);
  Fp_add (mem->tt1->a0, mem->tmp0, mem->tt1->a0);
  Fp3_add (mem->tt2, mem->tt0, mem->tt0);
  Fp3_add (mem->tt0, mem->tt2, mem->tt0);
  Fp3_add (mem->tt2, mem->tt1, mem->tt1);
  Fp3_add (mem->tt1, mem->tt2, mem->tt1);
  Fp3_add (c->a0, a->a0, a->a0);
  Fp3_neg (c->a0, c->a0);
  Fp3_add (c->a1, a->a1, a->a1);
  Fp3_add (c->a0, c->a0, mem->tt0);
  Fp3_add (c->a1, c->a1, mem->tt1);
}

void Fp6q_squ_cyclotomic (Fp6q_t c, Fp6q_t a) {
  Fp6q_squ_cyclotomic_mem (c, a, cyclotomic6_aux);
}

static void cp6_exp_by_u (Fp6q_t c, Fp6q_t a, Fp6q_t tmp0, Fp6q_t tmp1) {
  int i;

  Fp6q_copy (c, a);
  Fp6q_copy (tmp0, a);
  Fp6q_p3pow (tmp1, a);
  for (i = cp6_config.lenm - 2; i>=0; i--) {
    Fp6q_squ_cyclotomic (c, c);
    if (cp6_config.nafm[i] == 1) Fp6q_mul_lazy (c, c, tmp0);
    if (cp6_config.nafm[i] == -1) Fp6q_mul_lazy (c, c, tmp1);    
  }
}

static void cp6_finalexpo_mem (Fp6q_t c, Fp6q_t a, Aux_finalexpo6_t mem) {
  Fp6q_inv (mem->t0, a);
  Fp6q_p3pow (mem->t1, a);
  Fp6q_mul_lazy  (mem->t0, mem->t0, mem->t1); // t0 = a^(p^3-1)
  Fp6q_ppow (mem->t1, mem->t0);
  Fp6q_mul_lazy  (mem->t0, mem->t1, mem->t0); // t0 = (a^(p^3-1))^(p + 1) = f

  cp6_exp_by_u (mem->t1, mem->t0, mem->tmp0, mem->tmp1); 
  Fp6q_squ_cyclotomic (mem->t2, mem->t1);
  Fp6q_squ_cyclotomic (mem->t2, mem->t2);
  Fp6q_squ_cyclotomic (mem->t2, mem->t2);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t2);
  Fp6q_mul_lazy (mem->t2, mem->t0, mem->t1);
  Fp6q_squ_cyclotomic (mem->t0, mem->t0);
  cp6_exp_by_u (mem->t1, mem->t1, mem->tmp0, mem->tmp1);
  Fp6q_mul_lazy (mem->t0, mem->t0, mem->t1);
  Fp6q_squ_cyclotomic (mem->t3, mem->t1);
  Fp6q_squ_cyclotomic (mem->t3, mem->t3);
  Fp6q_mul_lazy (mem->t3, mem->t3, mem->t1);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t3);
  cp6_exp_by_u (mem->t1, mem->t1, mem->tmp0, mem->tmp1);
  Fp6q_mul_lazy (mem->t2, mem->t1, mem->t2);
  Fp6q_squ_cyclotomic (mem->t1, mem->t1);
  Fp6q_mul_lazy (mem->t0, mem->t0, mem->t1);
  Fp6q_squ_cyclotomic (mem->t3, mem->t1);
  Fp6q_squ_cyclotomic (mem->t3, mem->t3);
  Fp6q_squ_cyclotomic (mem->t3, mem->t3);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t3);
  cp6_exp_by_u (mem->t1, mem->t1, mem->tmp0, mem->tmp1);
  Fp6q_squ_cyclotomic (mem->t3, mem->t1);
  Fp6q_mul_lazy (mem->t0, mem->t0, mem->t3);
  Fp6q_squ_cyclotomic (mem->t4, mem->t3);
  Fp6q_squ_cyclotomic (mem->t5, mem->t4);
  Fp6q_squ_cyclotomic (mem->t5, mem->t5);
  Fp6q_mul_lazy (mem->t5, mem->t5, mem->t4);
  Fp6q_mul_lazy (mem->t3, mem->t3, mem->t5);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t3);
  cp6_exp_by_u (mem->t1, mem->t1, mem->tmp0, mem->tmp1);
  Fp6q_squ_cyclotomic (mem->t4, mem->t1);
  Fp6q_squ_cyclotomic (mem->t4, mem->t4);
  Fp6q_squ_cyclotomic (mem->t4, mem->t4);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t4);
  Fp6q_squ_cyclotomic (mem->t3, mem->t1);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t1);
  Fp6q_squ_cyclotomic (mem->t4, mem->t3);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t4);
  cp6_exp_by_u (mem->t3, mem->t3, mem->tmp0, mem->tmp1);
  Fp6q_squ_cyclotomic (mem->t4, mem->t3);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t4);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t3);
  cp6_exp_by_u (mem->t4, mem->t4, mem->tmp0, mem->tmp1);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t4);
  Fp6q_p3pow (mem->t2, mem->t2);
  Fp6q_ppow (mem->t0, mem->t0);
  Fp6q_mul_lazy (c, mem->t0, mem->t2);
}

void cp6_finalexpo (Fp6q_t c, Fp6q_t a) {
  cp6_finalexpo_mem (c, a, finalexpo6_aux);
}

static void cp6b_finalexpo_mem (Fp6q_t c, Fp6q_t a, Aux_finalexpo6_t mem) {
  Fp6q_inv (mem->t0, a);
  Fp6q_p3pow (mem->t1, a);
  Fp6q_mul_lazy  (mem->t0, mem->t0, mem->t1); // t0 = a^(p^3-1)
  Fp6q_ppow (mem->t1, mem->t0);
  Fp6q_mul_lazy  (mem->t0, mem->t1, mem->t0); // t0 = (a^(p^3-1))^(p + 1) = f

  Fp6q_squ_cyclotomic (mem->t1, mem->t0);
  Fp6q_squ_cyclotomic (mem->t2, mem->t1);
  Fp6q_mul_lazy (mem->t3, mem->t2, mem->t0); 
  Fp6q_mul_lazy (mem->t0, mem->t3, mem->t2);  
  Fp6q_mul_lazy (mem->t2, mem->t3, mem->t1);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1); 
  
  Fp6q_squ_cyclotomic (mem->t4, mem->t0); 
  Fp6q_squ_cyclotomic (mem->t1, mem->t4);  
  Fp6q_squ_cyclotomic (mem->t1, mem->t1); 
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t4);
  Fp6q_mul_lazy (mem->t4, mem->t0, mem->t4);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1);

  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_squ_cyclotomic (mem->t6, mem->t5);
  Fp6q_squ_cyclotomic (mem->t6, mem->t6);
  Fp6q_mul_lazy (mem->t5, mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t6, mem->t6, mem->t5);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t5);
  Fp6q_mul_lazy (mem->t3, mem->t3, mem->t6);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t6);
  Fp6q_squ_cyclotomic (mem->t6, mem->t6);  
  Fp6q_squ_cyclotomic (mem->t6, mem->t6);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t6);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1);

  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t6, mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t6);
  Fp6q_squ_cyclotomic (mem->t5, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t5);
  Fp6q_mul_lazy (mem->t4, mem->t4, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t5);
  Fp6q_squ_cyclotomic (mem->t6, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t6);
  Fp6q_mul_lazy (mem->t5, mem->t5, mem->t6);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t5);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1);

  Fp6q_squ_cyclotomic (mem->t0, mem->t0);
  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_p3pow (mem->t5, mem->t5);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t0, mem->t0, mem->t5);
  Fp6q_squ_cyclotomic (mem->t6, mem->t0);
  Fp6q_mul_lazy (mem->t6, mem->t6, mem->t5);
  Fp6q_mul_lazy (mem->t3, mem->t3, mem->t6);
  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t0, mem->t0, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_squ_cyclotomic (mem->t5, mem->t5);
  Fp6q_squ_cyclotomic (mem->t6, mem->t5);
  Fp6q_squ_cyclotomic (mem->t5, mem->t6);
  Fp6q_mul_lazy (mem->t5, mem->t5, mem->t6);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t5);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1);

  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t0);
  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t5);
  Fp6q_squ_cyclotomic (mem->t0, mem->t0);
  Fp6q_squ_cyclotomic (mem->t0, mem->t0);
  Fp6q_squ_cyclotomic (mem->t0, mem->t0);
  Fp6q_squ_cyclotomic (mem->t0, mem->t0);
  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t5);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1);

  Fp6q_squ_cyclotomic (mem->t5, mem->t0);
  Fp6q_mul_lazy (mem->t2, mem->t2, mem->t5);
  cp6_exp_by_u (mem->t0, mem->t0, mem->tmp0, mem->tmp1);

  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t0);
  Fp6q_p3pow (mem->t4, mem->t4);
  Fp6q_mul_lazy (mem->t3, mem->t3, mem->t4);
  Fp6q_ppow (mem->t3, mem->t3);
  Fp6q_p3pow (mem->t2, mem->t2);
  Fp6q_mul_lazy (mem->t1, mem->t1, mem->t2);
  Fp6q_mul_lazy (c, mem->t3, mem->t1);
}

void cp6b_finalexpo (Fp6q_t c, Fp6q_t a) {
  cp6b_finalexpo_mem (c, a, finalexpo6_aux);
}

static void cp3_finalexpo_mem(Fp3_t c, Fp3_t a, Aux_finalexpo3_t mem) {
	Fp3_inv(mem->t0, a);
	Fp3_ppow(mem->t1, a);
	Fp3_mul(mem->t0, mem->t0, mem->t1); // t0 = f0:=a^(p-1)
	Fp3_ppow(mem->t1, mem->t0);              // t1 = f1:=f0^p
	Fp3_ppow(mem->t2, mem->t1);              // t2 = f0^(p^2)
	Fp3_mul(mem->t2, mem->t2, mem->t1); // t2 = f0^(p^2 + p) = f0^(-1)
	Fp3_squ(mem->t3, mem->t0);
	if (cp3_config.nafa2[cp3_config.naflena2 - 2] == 1) Fp3_mul(mem->t3, mem->t3, mem->t0);
	if (cp3_config.nafa2[cp3_config.naflena2 - 2] == -1) Fp3_mul(mem->t3, mem->t3, mem->t2);
	for (int i = cp3_config.naflena2 - 3; i >= 0; i--) {
		Fp3_squ(mem->t3, mem->t3);
		if (cp3_config.nafa2[i] == 1) Fp3_mul(mem->t3, mem->t3, mem->t0);
		if (cp3_config.nafa2[i] == -1) Fp3_mul(mem->t3, mem->t3, mem->t2);
	}
	Fp3_mul(mem->t1, mem->t1, mem->t3); // t1 = f1:=f0^(a2+p) = a^((p-1)*(a2+p))
	Fp3_mul(mem->t2, mem->t0, mem->t1); // t2 = ff:=f0*f1
	// Double exponentiation below needs to compute f0^a0*f1^a1.
	Fp3_copy(c, mem->t1);
	for (int i = cp3_config.lena1 - 2; i > cp3_config.lena0 - 1; i--) {
		Fp3_squ(c, c);
		if (BIT(cp3_config.a1, i) == 1) Fp3_mul(c, c, mem->t1);
	}
	for (int i = cp3_config.lena0 - 1; i >= 0; i--) {
		Fp3_squ(c, c);
		if (BIT(cp3_config.a1, i) == 0 && BIT(cp3_config.a0, i) == 1) Fp3_mul(c, c, mem->t0);
		if (BIT(cp3_config.a1, i) == 1 && BIT(cp3_config.a0, i) == 0) Fp3_mul(c, c, mem->t1);
		if (BIT(cp3_config.a1, i) == 1 && BIT(cp3_config.a0, i) == 1) Fp3_mul(c, c, mem->t2);
	}
}

void cp3_finalexpo(Fp3_t c, Fp3_t a) {
	cp3_finalexpo_mem(c, a, finalexpo3_aux);
}



int bn_bls12_finalexpo_initialize_config (void) {
  int ret=1;

  if ((ret=Fp12_init (finalexpo12_aux->t0)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t1)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t2)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t3)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t4)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t5)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t6)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t7)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t8)) < 0) return ret;
  if ((ret=Fp12_init (finalexpo12_aux->t9)) < 0) return ret;
  if ((ret=Fp6_init (cyclotomic12_aux->tt0)) < 0) return ret;
  if ((ret=Fp6_init (cyclotomic12_aux->tt1)) < 0) return ret;
  if ((ret=Fp6_init (cyclotomic12_aux->tt2)) < 0) return ret;
  return ret;
}

void bn_bls12_finalexpo_free_config (void) {
  Fp12_free (finalexpo12_aux->t0);
  Fp12_free (finalexpo12_aux->t1);
  Fp12_free (finalexpo12_aux->t2);
  Fp12_free (finalexpo12_aux->t3);
  Fp12_free (finalexpo12_aux->t4);
  Fp12_free (finalexpo12_aux->t5);
  Fp12_free (finalexpo12_aux->t6);
  Fp12_free (finalexpo12_aux->t7);
  Fp12_free (finalexpo12_aux->t8);
  Fp12_free (finalexpo12_aux->t9);
  Fp6_free (cyclotomic12_aux->tt0);
  Fp6_free (cyclotomic12_aux->tt1);
  Fp6_free (cyclotomic12_aux->tt2);
}


int cp6_finalexpo_initialize_config (void) {
  int ret = 1;

  if ((ret=Fp6q_init (finalexpo6_aux->t0)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t1)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t2)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t3)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t4)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t5)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t6)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t7)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t8)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->t9)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->tmp0)) < 0) return ret;
  if ((ret=Fp6q_init (finalexpo6_aux->tmp1)) < 0) return ret;
  if ((ret=Fp3_init (cyclotomic6_aux->tt0)) < 0) return ret;
  if ((ret=Fp3_init (cyclotomic6_aux->tt1)) < 0) return ret;
  if ((ret=Fp3_init (cyclotomic6_aux->tt2)) < 0) return ret;
  if ((ret=Fp_init (cyclotomic6_aux->tmp0)) < 0) return ret;
  if ((ret=Fp_init (cyclotomic6_aux->tmp1)) < 0) return ret;
  return ret;
}

void cp6_finalexpo_free_config (void) {
  Fp6q_free (finalexpo6_aux->t0);
  Fp6q_free (finalexpo6_aux->t1);
  Fp6q_free (finalexpo6_aux->t2);
  Fp6q_free (finalexpo6_aux->t3);
  Fp6q_free (finalexpo6_aux->t4);
  Fp6q_free (finalexpo6_aux->t5);
  Fp6q_free (finalexpo6_aux->t6);
  Fp6q_free (finalexpo6_aux->t7);
  Fp6q_free (finalexpo6_aux->t8);
  Fp6q_free (finalexpo6_aux->t9);
  Fp6q_free (finalexpo6_aux->tmp0);
  Fp6q_free (finalexpo6_aux->tmp1);
  Fp3_free (cyclotomic6_aux->tt0);
  Fp3_free (cyclotomic6_aux->tt1);
  Fp3_free (cyclotomic6_aux->tt2);
  Fp_free  (cyclotomic6_aux->tmp0);
  Fp_free  (cyclotomic6_aux->tmp1);
}

int cp3_finalexpo_initialize_config (void) {
	int i, ret = 1;

	uint64_t a0[8] = { 0xA80000000000002B, 0xBF31F00000000054, 0xE7E8D23800000051, 0x9767A965DB000030,
		0x56B8A2E871520013, 0xF63776CA7CB7BF05, 0x5EBACA27219E58A8, 0x0158CFFC9264E111 };
	uint64_t a1[9] = { 0x040000000001040C, 0x930DA000000216A7, 0xC27F18800001F720, 0x9C0344DB430119CB,
		0x482CA96A81166681, 0x3D4FB452307C83C4, 0x1BA9804B76774979, 0xE6120613BDFF7756, 0x0000000000000003 };
	uint64_t a2[9] = { 0x6D000000000000FC, 0x3E9B3000000001FA, 0x6FC263E8000001E2, 0xCA65630848000117,
		0x6650AD919D54006A, 0x97D404BD57BC8C1B, 0x859C4645C71F1FB4, 0x05633FF249938445 };

	if ((ret=Fp3_init (finalexpo3_aux->t0)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t1)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t2)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t3)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t4)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t5)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t6)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t7)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t8)) < 0) return ret;
	if ((ret=Fp3_init (finalexpo3_aux->t9)) < 0) return ret;


	cp3_config.lena0 = 505;
	cp3_config.a0 = NULL;
	cp3_config.a0 = (uint64_t *)malloc(8 * sizeof(uint64_t));
	for (i = 0; i < 8; i++) cp3_config.a0[i] = a0[i];

	cp3_config.lena1 = 514;
	cp3_config.a1 = NULL;
	cp3_config.a1 = (uint64_t *)malloc(9 * sizeof(uint64_t));
	for (i = 0; i < 9; i++) cp3_config.a1[i] = a1[i];

	cp3_config.lena2 = 507;
	cp3_config.a2 = NULL;
	cp3_config.a2 = (uint64_t *)malloc(8 * sizeof(uint64_t));
	for (i = 0; i < 8; i++) cp3_config.a2[i] = a2[i];

	cp3_config.naflena2 = 508;
	cp3_config.nafa2 = NULL;
	cp3_config.nafa2 = (int *)malloc(cp3_config.naflena2  * sizeof (int));

	for (i = 0; i < cp3_config.naflena2; i++) cp3_config.nafa2[i] = 0;
	cp3_config.nafa2[2] = -1;
	cp3_config.nafa2[8] = 1;
	cp3_config.nafa2[56] = 1;
	cp3_config.nafa2[58] = -1;
	cp3_config.nafa2[60] = -1;
	cp3_config.nafa2[63] = 1;
	cp3_config.nafa2[65] = 1;
	cp3_config.nafa2[67] = -1;
	cp3_config.nafa2[73] = 1;
	cp3_config.nafa2[108] = -1;
	cp3_config.nafa2[110] = 1;
	cp3_config.nafa2[112] = -1;
	cp3_config.nafa2[114] = -1;
	cp3_config.nafa2[117] = 1;
	cp3_config.nafa2[119] = 1;
	cp3_config.nafa2[121] = -1;
	cp3_config.nafa2[126] = 1;
	cp3_config.nafa2[129] = 1;
	cp3_config.nafa2[133] = -1;
	cp3_config.nafa2[137] = 1;
	cp3_config.nafa2[163] = 1;
	cp3_config.nafa2[165] = -1;
	cp3_config.nafa2[170] = 1;
	cp3_config.nafa2[173] = -1;
	cp3_config.nafa2[175] = 1;
	cp3_config.nafa2[177] = 1;
	cp3_config.nafa2[182] = -1;
	cp3_config.nafa2[188] = -1;
	cp3_config.nafa2[191] = -1;
	cp3_config.nafa2[195] = -1;
	cp3_config.nafa2[197] = 1;
	cp3_config.nafa2[200] = 1;
	cp3_config.nafa2[219] = 1;
	cp3_config.nafa2[222] = 1;
	cp3_config.nafa2[227] = 1;
	cp3_config.nafa2[232] = -1;
	cp3_config.nafa2[234] = 1;
	cp3_config.nafa2[237] = -1;
	cp3_config.nafa2[239] = -1;
	cp3_config.nafa2[241] = -1;
	cp3_config.nafa2[243] = 1;
	cp3_config.nafa2[245] = -1;
	cp3_config.nafa2[247] = 1;
	cp3_config.nafa2[249] = 1;
	cp3_config.nafa2[251] = 1;
	cp3_config.nafa2[254] = -1;
	cp3_config.nafa2[256] = -1;
	cp3_config.nafa2[258] = -1;
	cp3_config.nafa2[260] = -1;
	cp3_config.nafa2[263] = 1;
	cp3_config.nafa2[274] = 1;
	cp3_config.nafa2[276] = 1;
	cp3_config.nafa2[278] = 1;
	cp3_config.nafa2[280] = 1;
	cp3_config.nafa2[282] = -1;
	cp3_config.nafa2[285] = 1;
	cp3_config.nafa2[287] = -1;
	cp3_config.nafa2[289] = 1;
	cp3_config.nafa2[292] = 1;
	cp3_config.nafa2[295] = -1;
	cp3_config.nafa2[297] = -1;
	cp3_config.nafa2[300] = -1;
	cp3_config.nafa2[302] = -1;
	cp3_config.nafa2[304] = 1;
	cp3_config.nafa2[308] = 1;
	cp3_config.nafa2[310] = 1;
	cp3_config.nafa2[313] = -1;
	cp3_config.nafa2[315] = 1;
	cp3_config.nafa2[317] = -1;
	cp3_config.nafa2[319] = -1;
	cp3_config.nafa2[322] = -1;
	cp3_config.nafa2[325] = 1;
	cp3_config.nafa2[330] = -1;
	cp3_config.nafa2[332] = 1;
	cp3_config.nafa2[335] = 1;
	cp3_config.nafa2[338] = -1;
	cp3_config.nafa2[342] = -1;
	cp3_config.nafa2[347] = -1;
	cp3_config.nafa2[349] = -1;
	cp3_config.nafa2[351] = -1;
	cp3_config.nafa2[353] = -1;
	cp3_config.nafa2[358] = -1;
	cp3_config.nafa2[360] = 1;
	cp3_config.nafa2[362] = 1;
	cp3_config.nafa2[370] = 1;
	cp3_config.nafa2[372] = 1;
	cp3_config.nafa2[374] = -1;
	cp3_config.nafa2[379] = -1;
	cp3_config.nafa2[381] = 1;
	cp3_config.nafa2[383] = 1;
	cp3_config.nafa2[386] = 1;
	cp3_config.nafa2[388] = -1;
	cp3_config.nafa2[390] = -1;
	cp3_config.nafa2[397] = 1;
	cp3_config.nafa2[400] = -1;
	cp3_config.nafa2[405] = 1;
	cp3_config.nafa2[408] = -1;
	cp3_config.nafa2[411] = 1;
	cp3_config.nafa2[414] = -1;
	cp3_config.nafa2[417] = -1;
	cp3_config.nafa2[419] = 1;
	cp3_config.nafa2[422] = 1;
	cp3_config.nafa2[425] = -1;
	cp3_config.nafa2[427] = 1;
	cp3_config.nafa2[430] = 1;
	cp3_config.nafa2[434] = -1;
	cp3_config.nafa2[437] = 1;
	cp3_config.nafa2[439] = -1;
	cp3_config.nafa2[441] = -1;
	cp3_config.nafa2[443] = 1;
	cp3_config.nafa2[447] = -1;
	cp3_config.nafa2[449] = -1;
	cp3_config.nafa2[451] = 1;
	cp3_config.nafa2[454] = 1;
	cp3_config.nafa2[458] = 1;
	cp3_config.nafa2[463] = -1;
	cp3_config.nafa2[466] = 1;
	cp3_config.nafa2[468] = 1;
	cp3_config.nafa2[471] = -1;
	cp3_config.nafa2[473] = 1;
	cp3_config.nafa2[475] = 1;
	cp3_config.nafa2[478] = 1;
	cp3_config.nafa2[481] = 1;
	cp3_config.nafa2[484] = -1;
	cp3_config.nafa2[494] = 1;
	cp3_config.nafa2[496] = -1;
	cp3_config.nafa2[498] = 1;
	cp3_config.nafa2[501] = -1;
	cp3_config.nafa2[503] = -1;
	cp3_config.nafa2[505] = -1;
	cp3_config.nafa2[507] = 1;

	return ret;
}

void cp3_finalexpo_free_config(void) {
	Fp3_free(finalexpo3_aux->t0);
	Fp3_free(finalexpo3_aux->t1);
	Fp3_free(finalexpo3_aux->t2);
	Fp3_free(finalexpo3_aux->t3);
	Fp3_free(finalexpo3_aux->t4);
	Fp3_free(finalexpo3_aux->t5);
	Fp3_free(finalexpo3_aux->t6);
	Fp3_free(finalexpo3_aux->t7);
	Fp3_free(finalexpo3_aux->t8);
	Fp3_free(finalexpo3_aux->t9);

  free(cp3_config.a0);
  free(cp3_config.a1);
  free(cp3_config.a2);
	free(cp3_config.nafa2);
}



#include "Qappairing.h"
#include "QapFp4.h"

typedef struct {
  QapFp12_t t0;
  QapFp12_t t1;
  QapFp12_t t2;
  QapFp12_t t3;
  QapFp12_t t4;
  QapFp12_t t5;
  QapFp12_t t6;
  QapFp12_t t7;
  QapFp12_t t8;
  QapFp12_t t9;
} aux_finalexpo12_t;

typedef aux_finalexpo12_t Aux_finalexpo12_t[1];

static Aux_finalexpo12_t finalexpo12_aux;

typedef struct {
  QapFp6_t tt0;
  QapFp6_t tt1;
  QapFp6_t tt2;
} aux_cyclotomic12_t;

typedef aux_cyclotomic12_t Aux_cyclotomic12_t[1];

static Aux_cyclotomic12_t cyclotomic12_aux;

static void Fp12_squ_cyclotomic_mem (QapFp12_t c, QapFp12_t a, Aux_cyclotomic12_t mem) {
  QapFp4_squ_sep (mem->tt0->a0, mem->tt1->a1, a->a0->a0, a->a1->a1);
  QapFp4_squ_sep (mem->tt0->a1, mem->tt1->a2, a->a1->a0, a->a0->a2);
  QapFp4_squ_sep (mem->tt0->a2, mem->tt1->a0, a->a0->a1, a->a1->a2);
  QapFp2_mulxi (mem->tt1->a0, mem->tt1->a0);
  QapFp6_add (mem->tt2, mem->tt0, mem->tt0);
  QapFp6_add (mem->tt0, mem->tt2, mem->tt0);
  QapFp6_add (mem->tt2, mem->tt1, mem->tt1);
  QapFp6_add (mem->tt1, mem->tt2, mem->tt1);
  QapFp6_add (c->a0, a->a0, a->a0);
  QapFp6_neg (c->a0, c->a0);
  QapFp6_add (c->a1, a->a1, a->a1);
  QapFp6_add (c->a0, c->a0, mem->tt0);
  QapFp6_add (c->a1, c->a1, mem->tt1);
}

void QapFp12_squ_cyclotomic (QapFp12_t c, QapFp12_t a) {
  Fp12_squ_cyclotomic_mem (c, a, cyclotomic12_aux);
}

static void bn_exp_by_u (QapFp12_t c, QapFp12_t a, QapFp12_t tmp0, QapFp12_t tmp1) {
  int i;

  QapFp12_copy (c, a);
  QapFp12_copy (tmp0, a);
  QapFp12_p6pow (tmp1, a);
  for (i = bn_config.lenu - 2; i>=0; i--) {
    QapFp12_squ_cyclotomic(c, c);
    if (bn_config.nafu[i] == 1) QapFp12_mul_lazy (c, c, tmp0);
    if (bn_config.nafu[i] == -1) QapFp12_mul_lazy (c, c, tmp1);    
  }
}

static void bn_finalexpo_neg_mem (QapFp12_t c, QapFp12_t a, Aux_finalexpo12_t mem) {
  QapFp12_inv (mem->t0, a);
  QapFp12_p6pow (mem->t1, a);
  QapFp12_mul_lazy (mem->t0, mem->t0, mem->t1); // t0 = a^(p^6-1)
  QapFp12_p2pow (mem->t1, mem->t0);
  QapFp12_mul_lazy (mem->t0, mem->t1, mem->t0); // t0 = (a^(p^6-1))^(p^2 + 1) = f

  bn_exp_by_u (mem->t1, mem->t0, mem->t5, mem->t6);
  QapFp12_p6pow (mem->t1, mem->t1);      // fx
  QapFp12_squ_cyclotomic (mem->t2, mem->t1);   // f2x
  QapFp12_squ_cyclotomic (mem->t1, mem->t2);   // f4x
  QapFp12_mul_lazy (mem->t1, mem->t2, mem->t1);// f6x
  bn_exp_by_u (mem->t3, mem->t1, mem->t5, mem->t6); 
  QapFp12_p6pow (mem->t3, mem->t3); // f6x2
  QapFp12_squ_cyclotomic (mem->t4, mem->t3); // f12x2
  bn_exp_by_u (mem->t4, mem->t4, mem->t5, mem->t6); 
  QapFp12_p6pow (mem->t4, mem->t4);  // f12x3
  QapFp12_mul_lazy (mem->t4, mem->t4, mem->t3); // a
  QapFp12_mul_lazy (mem->t4, mem->t4, mem->t1); // a
  QapFp12_p6pow (mem->t2, mem->t2); // b
  QapFp12_mul_lazy (mem->t2, mem->t2, mem->t4); // b
  QapFp12_mul_lazy (c, mem->t4, mem->t3); // fd
  QapFp12_mul_lazy (c, c, mem->t0); // fd
  QapFp12_ppow (mem->t1, mem->t2); // y0
  QapFp12_p2pow (mem->t3, mem->t4); // y1
  QapFp12_p6pow (mem->t0, mem->t0); // y2
  QapFp12_mul_lazy (mem->t0, mem->t0, mem->t2); // y2
  QapFp12_p3pow (mem->t0, mem->t0); // y2
  QapFp12_mul_lazy (c, c, mem->t1);
  QapFp12_mul_lazy (c, c, mem->t3);
  QapFp12_mul_lazy (c, c, mem->t0);
}

void Qapbn_finalexpo_neg (QapFp12_t c, QapFp12_t a) {
  bn_finalexpo_neg_mem (c, a, finalexpo12_aux);
}

static void bn_finalexpo_pos_mem (QapFp12_t c, QapFp12_t a, Aux_finalexpo12_t mem) {
  QapFp12_inv (mem->t0, a);
  QapFp12_p6pow (mem->t1, a);
  QapFp12_mul_lazy (mem->t0, mem->t0, mem->t1); // t0 = a^(p^6-1)
  QapFp12_p2pow (mem->t1, mem->t0);
  QapFp12_mul_lazy (mem->t0, mem->t1, mem->t0); // t0 = (a^(p^6-1))^(p^2 + 1) = f

  bn_exp_by_u (mem->t1, mem->t0, mem->t5, mem->t6); // fx
  QapFp12_squ_cyclotomic (mem->t2, mem->t1);   // f2x
  QapFp12_squ_cyclotomic (mem->t1, mem->t2);   // f4x
  QapFp12_mul_lazy (mem->t1, mem->t2, mem->t1);// f6x
  bn_exp_by_u (mem->t3, mem->t1, mem->t5, mem->t6); // f6x2
  QapFp12_squ_cyclotomic (mem->t4, mem->t3); // f12x2
  bn_exp_by_u (mem->t4, mem->t4, mem->t5, mem->t6); // f12x3
  QapFp12_mul_lazy (mem->t4, mem->t4, mem->t3); // a
  QapFp12_mul_lazy (mem->t4, mem->t4, mem->t1); // a
  QapFp12_p6pow (mem->t2, mem->t2); // b
  QapFp12_mul_lazy (mem->t2, mem->t2, mem->t4); // b
  QapFp12_mul_lazy (c, mem->t4, mem->t3); // fd
  QapFp12_mul_lazy (c, c, mem->t0); // fd
  QapFp12_ppow (mem->t1, mem->t2); // y0
  QapFp12_p2pow (mem->t3, mem->t4); // y1
  QapFp12_p6pow (mem->t0, mem->t0); // y2
  QapFp12_mul_lazy (mem->t0, mem->t0, mem->t2); // y2
  QapFp12_p3pow (mem->t0, mem->t0); // y2
  QapFp12_mul_lazy (c, c, mem->t1);
  QapFp12_mul_lazy (c, c, mem->t3);
  QapFp12_mul_lazy (c, c, mem->t0);
}

void Qapbn_finalexpo_pos (QapFp12_t c, QapFp12_t a) {
  bn_finalexpo_pos_mem (c, a, finalexpo12_aux);
}

int Qapbn_bls12_finalexpo_initialize_config (void) {
  int ret=1;

  if ((ret=QapFp12_init (finalexpo12_aux->t0)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t1)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t2)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t3)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t4)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t5)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t6)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t7)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t8)) < 0) return ret;
  if ((ret=QapFp12_init (finalexpo12_aux->t9)) < 0) return ret;

  if ((ret=QapFp6_init (cyclotomic12_aux->tt0)) < 0) return ret;
  if ((ret=QapFp6_init (cyclotomic12_aux->tt1)) < 0) return ret;
  if ((ret=QapFp6_init (cyclotomic12_aux->tt2)) < 0) return ret;
  return ret;
}

void Qapbn_bls12_finalexpo_free_config (void) {
  QapFp12_free (finalexpo12_aux->t0);
  QapFp12_free (finalexpo12_aux->t1);
  QapFp12_free (finalexpo12_aux->t2);
  QapFp12_free (finalexpo12_aux->t3);
  QapFp12_free (finalexpo12_aux->t4);
  QapFp12_free (finalexpo12_aux->t5);
  QapFp12_free (finalexpo12_aux->t6);
  QapFp12_free (finalexpo12_aux->t7);
  QapFp12_free (finalexpo12_aux->t8);
  QapFp12_free (finalexpo12_aux->t9);

  QapFp6_free (cyclotomic12_aux->tt0);
  QapFp6_free (cyclotomic12_aux->tt1);
  QapFp6_free (cyclotomic12_aux->tt2);
}

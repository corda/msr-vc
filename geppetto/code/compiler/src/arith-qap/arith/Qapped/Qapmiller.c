#include "Qappairing.h"
#include "QapFp2.h"
#include "QapFp4.h"
#include "QapFp6.h"
#include "QapFp12.h"

bn_config_t bn_config;

#define PREMEM_SIZE	388
	// measured in Fp2s.  2*(3+(bn_config.len6u2-3)*3+2)

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
  QapFp2_t l00;
  QapFp2_t l01;
  QapFp2_t l10;
  QapFp_t negxP;
  QapFp_t negyP;
  QapFp2_t negyP_Fp2;
  QapPoint_wp2_t R;
  QapPoint_wp2_t negQ;
  QapPoint_wp2_t Qp;
  QapPoint_wp2_t Qp2;
  QapFp12_t f;
  QapFp12_t premem[PREMEM_SIZE];
} aux_miller12_t;

typedef aux_miller12_t Aux_miller12_t[1];

static Aux_miller12_t miller12_aux;


/* Line function for a doubling step on BN curves using affine coordinates.
 * The coefficients returned by this function are the coefficients of the line
 * l: = -yQ*v + (lambda*xQ*v*w + nu*v ^ 2 * w), note that this assumes that yQ:=-yQ has been precomputed.
 */
static void Qapbn_bls12_line_affdbl_mem(QapFp2_t lambda, QapFp2_t nu, QapPoint_wp2_t R, Aux_miller12_t mem) {
	QapFp2_squ(lambda, R->X); //lambda: = x1 ^ 2;
	QapFp2_add(mem->t1, lambda, lambda); //c: = lambda + lambda;
	QapFp2_add(lambda, lambda, mem->t1); //lambda: = c + lambda;
	QapFp2_add(mem->t1, R->Y, R->Y); //c: = 2 * y1;
	// printf("\nQapFp2 Inv: "); QapFp2_print(mem->t1);
	QapFp2_inv(mem->t1, mem->t1); //c: = 1 / c;
	QapFp2_mul(lambda, lambda, mem->t1); //lambda: = lambda*c;
	QapFp2_mul(mem->t1, lambda, R->X); //c: = lambda*x1;
	QapFp2_sub(nu, R->Y, mem->t1);//nu: = y1 - c;
	QapFp2_squ(mem->t0, lambda); //sqr: = lambda ^ 2;
	QapFp2_add(mem->t1, R->X, R->X); //c: = 2 * x1;
	QapFp2_sub(R->X, mem->t0, mem->t1);//x3: = sqr - c;
	QapFp2_mul(R->Y, R->X, lambda); //y3: = lambda*x3;
	QapFp2_add(R->Y, R->Y, nu); //y3: = y3 + nu;
	QapFp2_neg(R->Y, R->Y); //y3: = -y3;
}

void Qapbn_bls12_line_affdbl(QapFp2_t lambda, QapFp2_t nu, QapPoint_wp2_t R) {
	Qapbn_bls12_line_affdbl_mem(lambda, nu, R, miller12_aux);
}

/* Line function for an addition step on BN curves using affine coordinates.
* The coefficients returned by this function are the coefficients of the line
* l: = -yQ*v + (lambda*xQ*v*w + nu*v ^ 2 * w), not that this assumes that xQ:=-xQ has been precomputed.
*/
static void Qapbn_bls12_line_affadd_mem(QapFp2_t lambda, QapFp2_t nu, QapPoint_wp2_t R, QapPoint_wp2_t P, Aux_miller12_t mem) {
	QapFp2_sub(lambda, P->Y, R->Y); //lambda: = y2 - y1;
	QapFp2_sub(mem->t0, P->X, R->X); //c: = x2 - x1;
	QapFp2_inv(mem->t0, mem->t0); //c: = 1 / c;
	QapFp2_mul(lambda, lambda, mem->t0); //lambda: = lambda*c;
	QapFp2_squ(mem->t1, lambda); //sqr: = lambda ^ 2;
	QapFp2_mul(nu, R->X, lambda); //nu: = x1*lambda;
	QapFp2_sub(R->X, mem->t1, R->X); //x3: = sqr - x1;
	QapFp2_sub(R->X, R->X, P->X); //x3: = x3 - x2;
	QapFp2_sub(nu, R->Y, nu); //nu: = y1 - nu;
	QapFp2_mul(R->Y, lambda, R->X); //y3: = lambda*x3;
	QapFp2_add(R->Y, R->Y, nu); //y3: = y3 + nu;
	QapFp2_neg(R->Y, R->Y); //y3: = -y3;
}

void Qapbn_bls12_line_affadd(QapFp2_t lambda, QapFp2_t nu, QapPoint_wp2_t R, QapPoint_wp2_t P) {
	Qapbn_bls12_line_affadd_mem(lambda, nu, R, P, miller12_aux);
}

void Qapbn_Fp2_wp_ppow_untwist(QapPoint_wp2_t Qp, QapPoint_wp2_t Q) {
  QapFp2_ppow(Qp->X, Q->X);
  QapFp2_ppow(Qp->Y, Q->Y);
  QapFp2_mul(Qp->Y, Qp->Y, Fp12_config.wppow3);
  QapFp2_muli(Qp->X, Qp->X);
  QapFp2_mulFp(Qp->X, Qp->X, Fp6_config.vppow);
}

static void Qapbn_Fp2_wpaff_ppow_untwist(QapPoint_wp2_t Qp, QapPoint_wp2_t Q) {
	QapFp2_ppow(Qp->X, Q->X);
	QapFp2_ppow(Qp->Y, Q->Y);
	QapFp2_mul(Qp->Y, Qp->Y, Fp12_config.wppow3);
	QapFp2_muli(Qp->X, Qp->X);
	QapFp2_mulFp(Qp->X, Qp->X, Fp6_config.vppow);
}

void Qapbn_Fp2_wp_p2pow_untwist(QapPoint_wp2_t Qp, QapPoint_wp2_t Q) {
  QapFp2_mulFp(Qp->X, Q->X, Fp6_config.vp2pow);
  QapFp2_mulFp(Qp->Y, Q->Y, Fp12_config.wp2pow3);
}


//#define TICKER 0
#if TICKER
#define TICK log0 = log; log = nRoot(); assert(100000 + log - log0)
#else 
#define TICK {}
#endif

static void Qapbn_optimal_ate_miller_aff_mem(QapFp12_t f, QapPoint_wp2_t Q, QapPoint_wp_t P, Aux_miller12_t mem) {
	int i;

#if TICKER
	int log, log0;
	log = 0;
#endif

	QapFp12_set_zero(f);
	QapFp_neg(mem->negyP, P->Y);
	QapFp_set_ui(mem->negyP_Fp2->a1, 0);
	QapFp_copy(mem->negyP_Fp2->a0, mem->negyP);
	QapFp2_wp_copy(mem->R, Q);
	QapFp2_wp_neg(mem->negQ, Q);

	Qapbn_bls12_line_affdbl(mem->l01, mem->l10, mem->R);
	QapFp_set_ui(f->a0->a1->a1, 0);
	QapFp_copy (f->a0->a1->a0, mem->negyP);
	QapFp2_mulFp(f->a1->a1, mem->l01, P->X);
	QapFp2_copy(f->a1->a2, mem->l10);
	TICK;

	if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
		Qapbn_bls12_line_affadd(mem->l01, mem->l10, mem->R, Q);
		QapFp2_mulFp(mem->l01, mem->l01, P->X);
		QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
	}
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
		Qapbn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->negQ);
		QapFp2_mulFp(mem->l01, mem->l01, P->X);
		QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
	}
	//printf ("1: f = "); QapFp12_print (f);
	for (i = bn_config.len6u2 - 3; i >= 0; i--) { // 36* + 24* or 56* (Left/Both) for each i (Left)
		QapFp12_squ_lazy(f, f);
		TICK; // 36* (left)

		//if (i>62) {printf ("1squ.%d: f = ", i); QapFp12_print (f);}
		Qapbn_bls12_line_affdbl(mem->l01, mem->l10, mem->R);
		QapFp2_mulFp(mem->l01, mem->l01, P->X);
		QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
		//if (i>62) {printf ("1squstep.%d: f = ", i); QapFp12_print (f);}
		if (bn_config.naf6u2[i] == 1) {
			Qapbn_bls12_line_affadd(mem->l01, mem->l10, mem->R, Q);
			QapFp2_mulFp(mem->l01, mem->l01, P->X);
			QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
		}
		if (bn_config.naf6u2[i] == -1) {
			Qapbn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->negQ);
			QapFp2_mulFp(mem->l01, mem->l01, P->X);
			QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
		}
		TICK; // 24* or 48* (Left) or 56* (Both)
		//if (i>62) {printf ("i = %d: f = ", i); QapFp12_print (f);}  
	}
	//printf ("2: f = "); QapFp12_print (f);
	Qapbn_Fp2_wpaff_ppow_untwist(mem->Qp, Q);
	Qapbn_Fp2_wp_p2pow_untwist(mem->Qp2, Q);
	QapFp2_wp_neg(mem->Qp2, mem->Qp2);

	if (bn_config.signu == -1) {
		QapFp2_wp_neg(mem->R, mem->R);
		QapFp12_p6pow(f, f);
	}

	Qapbn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->Qp);
	QapFp2_mulFp(mem->l01, mem->l01, P->X);
	QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);

	Qapbn_bls12_line_affadd(mem->l01, mem->l10, mem->R, mem->Qp2);
	QapFp2_mulFp(mem->l01, mem->l01, P->X);
	QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, mem->l10);
	//printf ("3: f = "); QapFp12_print (f);
}

void Qapbn_optimal_ate_miller_aff(QapFp12_t f, QapPoint_wp2_t Q, QapPoint_wp_t P) {
	Qapbn_optimal_ate_miller_aff_mem(f, Q, P, miller12_aux);
}

#define UPDATE_POINTERS(j)	{ pl01 = &premem[j*2+0]; pl10 = &premem[j*2+1]; }

static void Qapbn_optimal_ate_miller_precompute_aff_mem(QapPoint_wp2_t Q, QapFp2_t premem[PREMEM_SIZE], Aux_miller12_t mem) {
	int i, j;
	QapFp2_t *pl01, *pl10;

	QapFp2_wp_copy(mem->R, Q);
	QapFp2_wp_neg(mem->negQ, Q);

	j = 0;
	UPDATE_POINTERS(j);
	Qapbn_bls12_line_affdbl(*pl01, *pl10, mem->R);
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
		j++;
		UPDATE_POINTERS(j);
		Qapbn_bls12_line_affadd(*pl01, *pl10, mem->R, Q);
	}
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
		j++;
		UPDATE_POINTERS(j);
		Qapbn_bls12_line_affadd(*pl01, *pl10, mem->R, mem->negQ);
	}
	for (i = bn_config.len6u2 - 3; i >= 0; i--) {
		j++;
		UPDATE_POINTERS(j);
		Qapbn_bls12_line_affdbl(*pl01, *pl10, mem->R);
		if (bn_config.naf6u2[i] == 1) {
			j++;
			UPDATE_POINTERS(j);
			Qapbn_bls12_line_affadd(*pl01, *pl10, mem->R, Q);
		}
		if (bn_config.naf6u2[i] == -1) {
			j++;
			UPDATE_POINTERS(j);
			Qapbn_bls12_line_affadd(*pl01, *pl10, mem->R, mem->negQ);
		}
	}
	Qapbn_Fp2_wpaff_ppow_untwist(mem->Qp, Q);
	Qapbn_Fp2_wp_p2pow_untwist(mem->Qp2, Q);
	QapFp2_wp_neg(mem->Qp2, mem->Qp2);

	if (bn_config.signu == -1) {
		QapFp2_wp_neg(mem->R, mem->R);
	}

	j++;
	UPDATE_POINTERS(j);
	Qapbn_bls12_line_affadd(*pl01, *pl10, mem->R, mem->Qp);

	j++;
	UPDATE_POINTERS(j);
	Qapbn_bls12_line_affadd(*pl01, *pl10, mem->R, mem->Qp2);
}

void Qapbn_optimal_ate_miller_precompute_aff(QapPoint_wp2_t Q, QapFp2_t premem[PREMEM_SIZE]) {
	Qapbn_optimal_ate_miller_precompute_aff_mem(Q, premem, miller12_aux);
}

static void Qapbn_optimal_ate_miller_useprecomputed_aff_mem(QapFp12_t f, QapPoint_wp_t P, QapFp2_t premem[PREMEM_SIZE], Aux_miller12_t mem) {
	int i, j;
	QapFp2_t *pl01, *pl10;

	QapFp12_set_zero(f);
	QapFp_neg(mem->negyP, P->Y);
	QapFp_set_ui(mem->negyP_Fp2->a1, 0);
	QapFp_copy(mem->negyP_Fp2->a0, mem->negyP);

	j = 0;
	UPDATE_POINTERS(j);

	QapFp_set_ui(f->a0->a1->a1, 0);
	QapFp_copy(f->a0->a1->a0, mem->negyP);
	QapFp2_mulFp(f->a1->a1, *pl01, P->X);
	QapFp2_copy(f->a1->a2, *pl10);
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == 1) {
		j++;
		UPDATE_POINTERS(j);
		QapFp2_mulFp(*pl01, *pl01, P->X);
		QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, *pl01, *pl10);
	}
	if (bn_config.naf6u2[bn_config.len6u2 - 2] == -1) {
		j++;
		UPDATE_POINTERS(j);
		QapFp2_mulFp(mem->l01, *pl01, P->X);
		QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, *pl10);
	}
	for (i = bn_config.len6u2 - 3; i >= 0; i--) {
		QapFp12_squ_lazy(f, f);
		j++;
		UPDATE_POINTERS(j);
		QapFp2_mulFp(mem->l01, *pl01, P->X);
		QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, *pl10);
		if (bn_config.naf6u2[i] == 1) {
			j++;
			UPDATE_POINTERS(j);
			QapFp2_mulFp(mem->l01, *pl01, P->X);
			QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, *pl10);
		}
		if (bn_config.naf6u2[i] == -1) {
			j++;
			UPDATE_POINTERS(j);
			QapFp2_mulFp(mem->l01, *pl01, P->X);
			QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, *pl10);
		}
	}

	if (bn_config.signu == -1) {
		QapFp12_p6pow(f, f);
	}
	
	j++;
	UPDATE_POINTERS(j);
	QapFp2_mulFp(mem->l01, *pl01, P->X);
	QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, *pl10);

	j++;
	UPDATE_POINTERS(j);
	QapFp2_mulFp(mem->l01, *pl01, P->X);
	QapFp12_mul_sparse_untwist_lazy(f, f, mem->negyP_Fp2, mem->l01, *pl10);
}

void Qapbn_optimal_ate_miller_useprecomputed_aff(QapFp12_t f, QapPoint_wp_t P, QapFp2_t premem[PREMEM_SIZE]) {
	Qapbn_optimal_ate_miller_useprecomputed_aff_mem(f, P, premem, miller12_aux);
}


/********************************************************************************************/

int Qapbn_bls12_miller_initialize_config (QapFp2_t bt) {
  int i, ret;

  if ((ret=QapFp2_init (miller12_aux->t0)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t1)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t2)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t3)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t4)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t5)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t6)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t7)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->t8)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->l00)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->l01)) < 0) return ret;
  if ((ret=QapFp2_init (miller12_aux->l10)) < 0) return ret;
  if ((ret=QapFp_init (miller12_aux->negxP)) < 0) return ret;
  if ((ret=QapFp_init (miller12_aux->negyP)) < 0) return ret;
  if ((ret=QapFp2_init(miller12_aux->negyP_Fp2)) < 0) return ret;
  if ((ret=QapFp2_wp_init (miller12_aux->R)) < 0) return ret;
  if ((ret=QapFp2_wp_init (miller12_aux->negQ)) < 0) return ret;
  if ((ret=QapFp2_wp_init (miller12_aux->Qp)) < 0) return ret;
  if ((ret=QapFp2_wp_init (miller12_aux->Qp2)) < 0) return ret;
  if ((ret=QapFp12_init (miller12_aux->f)) < 0) return ret;

  QapFp12_set_zero (miller12_aux->f);

  if (PAIR_CURVE == BN12) {
    if (QapFp2_init(bn_config.bt3) < 0) printf ("bn_config memory error.\n");
    QapFp2_add (bn_config.bt3, bt, bt);
    QapFp2_add (bn_config.bt3, bn_config.bt3, bt);

    bn_config.lenu = 63;
#if 0
    bn_config.nafu = NULL;
    bn_config.nafu = (int *) malloc (bn_config.lenu  * sizeof (int));
    if (bn_config.nafu == NULL) return ERR_OUT_OF_MEMORY;
#endif
	QapDbgAssert(sizeof(bn_config.nafu)==bn_config.lenu*sizeof(bn_config.nafu[0]));

    for (i = 0; i < bn_config.lenu; i++) bn_config.nafu[i] = 0;
    bn_config.nafu[0] = 1;
    bn_config.nafu[55] = 1;
    bn_config.nafu[62] = 1;
  
    bn_config.signu = -1;

    bn_config.len6u2 = 66;
#if 0
    bn_config.naf6u2 = NULL;
    bn_config.naf6u2 = (int *) malloc (bn_config.len6u2  * sizeof (int));
    if (bn_config.naf6u2 == NULL) return ERR_OUT_OF_MEMORY;
#endif
	QapDbgAssert(sizeof(bn_config.naf6u2)==bn_config.len6u2*sizeof(bn_config.naf6u2[0]));
  
    for (i = 0; i < bn_config.len6u2; i++) bn_config.naf6u2[i] = 0;
    bn_config.naf6u2[2] = 1;
    bn_config.naf6u2[56] = -1;
    bn_config.naf6u2[58] = 1;
    bn_config.naf6u2[63] = -1;
    bn_config.naf6u2[65] = 1;
  } else if (PAIR_CURVE == BN12tiny) {
    if (QapFp2_init(bn_config.bt3) < 0) printf("bn_config memory error.\n");
    QapFp2_add(bn_config.bt3, bt, bt);
    QapFp2_add(bn_config.bt3, bn_config.bt3, bt);

    bn_config.lenu = 13;
#if 0
    bn_config.nafu = NULL;
    bn_config.nafu = (int *)malloc(bn_config.lenu  * sizeof (int));
    if (bn_config.nafu == NULL) return ERR_OUT_OF_MEMORY;
#endif
    QapDbgAssert(sizeof(bn_config.nafu) >= bn_config.lenu*sizeof(bn_config.nafu[0]));

    for (i = 0; i < bn_config.lenu; i++) bn_config.nafu[i] = 0;
    bn_config.nafu[0] = 1;
    bn_config.nafu[2] = -1;
    bn_config.nafu[5] = 1;
    bn_config.nafu[9] = -1;
    bn_config.nafu[12] = 1;

    bn_config.signu = -1;

    bn_config.len6u2 = 15;
#if 0
    bn_config.naf6u2 = NULL;
    bn_config.naf6u2 = (int *)malloc(bn_config.len6u2  * sizeof (int));
    if (bn_config.naf6u2 == NULL) return ERR_OUT_OF_MEMORY;
#endif
    QapDbgAssert(sizeof(bn_config.naf6u2) >= bn_config.len6u2*sizeof(bn_config.naf6u2[0]));

    for (i = 0; i < bn_config.len6u2; i++) bn_config.naf6u2[i] = 0;
    bn_config.naf6u2[2] = -1;
    bn_config.naf6u2[4] = -1;
    bn_config.naf6u2[6] = -1;
    bn_config.naf6u2[8] = 1;
    bn_config.naf6u2[10] = 1;
    bn_config.naf6u2[12] = 1;
    bn_config.naf6u2[14] = 1;
  }  else {
    printf ("ERROR: PAIR_CURVE not valid!");
  }
  return ERR_SUCCESS;
}

void Qapbn_bls12_miller_free_config (void) {
  QapFp2_free (miller12_aux->t0);
  QapFp2_free (miller12_aux->t1);
  QapFp2_free (miller12_aux->t2);
  QapFp2_free (miller12_aux->t3);
  QapFp2_free (miller12_aux->t4);
  QapFp2_free (miller12_aux->t5);
  QapFp2_free (miller12_aux->t6);
  QapFp2_free (miller12_aux->t7);
  QapFp2_free (miller12_aux->t8);
  QapFp2_free (miller12_aux->l00);
  QapFp2_free (miller12_aux->l01);
  QapFp2_free (miller12_aux->l10);
  QapFp_free (miller12_aux->negxP);
  QapFp_free (miller12_aux->negyP);
  QapFp2_free(miller12_aux->negyP_Fp2);
  QapFp2_wp_free (miller12_aux->R);
  QapFp2_wp_free (miller12_aux->negQ);
  QapFp2_wp_free (miller12_aux->Qp);
  QapFp2_wp_free (miller12_aux->Qp2);
  QapFp12_free (miller12_aux->f);

  if (PAIR_CURVE == BN12) {
#if 0
    free (bn_config.nafu);
    free (bn_config.naf6u2);
#endif
    QapFp2_free (bn_config.bt3);
  }
}


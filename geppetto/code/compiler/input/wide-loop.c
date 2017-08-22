#define GEPPETTO_NAME     "wide-loop"  // Must be defined before including geppetto.h
#define GEPPETTO_VERIFY_H "build/wide-loop.verify.h"  // Must be defined before including geppetto.h
#define GEPPETTO_VERIFY_C "build/wide-loop.verify.c"  // Must be defined before including geppetto.h
#define NUM_BANKS 3
#include "geppetto.h"

// -----------------------

// adapted from loop.c, then v-loop.c to explore trade-offs on the size of the internal variables.
// we try to bootstrap this one at a large scale

#include "pinocchio.h" 
#ifndef MQAP
#include <stdio.h>
#include <string.h>
#include "pinocchio.c"
#endif

// Set to 0 for wide-loop Classic; 1 for adjacency-matrix mode
#define ADJACENCY_MATRIX_MODE 0

// ---------------- roughly adapted from condition.c 

// Computation settings
#define SIZE 5    // with 10 we pay 400 K* just for the bootstrapped IO. 
#define INNER 10  // #iteration within each loop; avoid q-overflow; if (32 + LOG2 SIZE) > 253 we need intermediate truncations.
#define ITER 4    // #quapped iteration; must be even.

typedef struct Matrix { int x[SIZE][SIZE]; } matrix;

// another parameter to tweak; for now we have a reasonable mix of 1, 16, and 32 bit arithmetic

#if ADJACENCY_MATRIX_MODE 
#define MAXELT 1
#else
#define MAXELT 0xFFFF
#endif

int wltrunc(int n) { return (arith(n & MAXELT)); };
void wlbound(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		R->x[i][j] = bound(R->x[i][j], 0, MAXELT);
}

void safe(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		R->x[i][j] += 1;
}
void zero(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		R->x[i][j] = 0;
}
void prn(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		printf(" x[%3d][%3d] = %d\n", i, j, R->x[i][j]);;
}
void nonzero(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
	if (R->x[i][j] == 0) {
		printf("setting [%d][%d]\n", i, j);
		R->x[i][j] = 1;
	}
}
void one(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		R->x[i][j] = (i == j);
}
void mataddone(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
		R->x[i][i] = R->x[i][i] + 1;
}
void rndm(int n, struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		R->x[i][j] = wltrunc(i * 11 + j * 13 + n);
	nonzero(R);
}
void trc(struct Matrix *R) {
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++)
		R->x[i][j] = wltrunc(R->x[i][j]);
}
void matmul(struct Matrix *X, struct Matrix *Y, struct Matrix *R) {
	// using a temp matrix for simplicity
	struct Matrix t;
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++) {
		t.x[i][j] = 0;
		for (int k = 0; k < SIZE; k++)
			t.x[i][j] += X->x[i][k] * Y->x[k][j];
	}
	trc(&t);
	*R = t;
}
void matadd(struct Matrix *X, struct Matrix *Y, struct Matrix *R) {
	// using a temp matrix for simplicity
	//struct Matrix t;
	for (int i = 0; i < SIZE; i++)
	for (int j = 0; j < SIZE; j++) {
		R->x[i][j] = X->x[i][j] + Y->x[i][j];
	}
}

// -------------------------

struct Outer {
	matrix A;
};

struct Loop {
	matrix X;
}; 

void body(struct Outer *o, struct Loop *l) {
	for (int i = 0; i < INNER; i++) {
		matmul(&(l->X), &(o->A), &(l->X));
#if !ADJACENCY_MATRIX_MODE 
		mataddone(&(l->X));
#endif
	}
	// prn(&(l->X);
	printf("X[3,1] = %d \n", l->X.x[3][1]);
};


void spec(int n, struct Outer* o, struct Loop* init, struct Loop* result) {
	int i;
	struct Loop l = *init;
	for (i = 0; i < n; i++) body(o, &l);
	*result = l;
};

//-------------------------------------------------------------

BANK(Odd,struct Loop) 

BANK(Even,struct Loop)

BANK(Const,struct Outer)

#ifdef VERIFY
#include GEPPETTO_VERIFY_C
#endif

Odd even(Const c, Even l0) {
	struct Outer o; load_Const(c, &o);
	struct Loop l; load_Even(l0, &l);
	body(&o, &l);
	return(save_Odd(&l));
}
Even odd(Const c, Odd l1) {
	struct Outer o; load_Const(c, &o);
	struct Loop l; load_Odd(l1, &l);
	body(&o, &l);
	return(save_Even(&l));
}

void outer(int n, struct Outer *constants, struct Loop* initials, struct Loop* finals) {
	Const c = save_Const(constants);
	// printf("a = %d\n", constants->a);
	Even x0 = save_Even(initials);
	Odd  x1 = 0;
	for (int i = 0; i < n; i++) {
		if (i % 2) x0 = OUTSOURCE( odd, c, x1);
		else       x1 = OUTSOURCE(even, c, x0);
	}
	if (n % 2) load_Odd (x1, finals);
	else       load_Even(x0, finals);
}

// a bit simpler assuming n is even
void outer2(int n, struct Outer *o, struct Loop* x, struct Loop* r) {
	int i;
	Const c = save_Const(o);
	Even x0 = save_Even(x);
	Odd x1 = 0;
	for (i = 0; i < n; i += 2) {
		x1 = OUTSOURCE(even, c, x0);
		x0 = OUTSOURCE(odd, c, x1);
	}
	load_Even(x0, r);
}


//----------------------------------------------------------------

#define BOOTSTRAP 0 // meaning that the program is bootstrapped, not necessarily compiling in bootstrap mode
#ifdef BOOTSTRAP
// a mostly generic template for bootstrapping single-boot-level proofs 
#define GEPPETTO_NAME_BOOT "wide-loop-verify"
#define GEPPETTO_VERIFY_BOOT_H "build/wide-loop-verify.verify.h"
#define GEPPETTO_VERIFY_BOOT_C "build/wide-loop-verify.verify.c"
#define NUM_BANKS_BOOT 3
#include "geppetto-boot.h"

BANK_BOOT(0,Odd0, struct Loop)

BANK_BOOT(1,Even0, struct Loop)

BANK_BOOT(2,Const0, struct Outer)

Odd0 outer_proxy(Const0 constants, Even0 initials) {
	init(); // initialization must occur here to be compiled by bootstrapping
	struct Outer o;
	struct Loop vars;
	load_Const0(constants, &o);
	load_Even0(initials, &vars);
	outer(ITER,&o, &vars, &vars);
	safe(&(vars.X));
	return save_Odd0(&vars);
}

// this assume a FIXED, ODD number of interations (ITER)
void outer_stub(struct Outer *constants, struct Loop* initials, struct Loop* finals) {
	wlbound(&(constants->A));
	wlbound(&(initials->X));
	Const0 cst = save_Const0(constants);
	Even0 e0 = save_Even0(initials);
	Odd0 o0 = OUTSOURCE_BOOT(0, outer_proxy, cst, e0);
	load_Odd0(o0, finals);
}
#endif

//----------------------------------------------------------------

int main() {
#ifdef BOOTSTRAP
	init_BOOT();
#else
	init();
#endif
	struct Outer o;
	struct Loop v;
	printf("Sampling A\n");
	rndm(1, &(o.A));
	printf("Sampling X\n");
	rndm(3, &(v.X));
	//prn(&(o.A));
	//prn(&(v.X));
	printf("Running spec:\n");
	spec(ITER, &o, &v, &v);
	//prn(&(vars.X));
	printf("A[3,1] = %d spec size=%d inner=%d iter=%d\n", v.X.x[3][1], SIZE, INNER, ITER);

	printf("Running impl:\n");
	rndm(3, &(v.X));
#ifdef BOOTSTRAP
	timer_start();
	outer_stub(&o, &v, &v);
	timer_end();
	timer_print();
#else
	timer_start();
	outer(ITER, &o, &v, &v);
	timer_end();
	timer_print();
#endif
	//prn(&(vars.X));
	printf("A[3,1] = %d impl size=%d inner=%d iter=%d\n", v.X.x[3][1], SIZE, INNER, ITER);
	return 0;
}

#include "one-matrix-ifc.h"
#include "pinocchio.h" 
#ifndef MQAP
#include <stdio.h>
#endif

void initM(struct M* A)
{
	int i, j;
	for (i = 0; i < SIZE; i += 1)
		for (j = 0; j < SIZE; j += 1)
			A->x[i][j] = 0;
}

void addM(struct M* A, struct M* B, struct M* R)
{
	int i, j;
	for (i = 0; i < SIZE; i += 1)
	{
		for (j = 0; j < SIZE; j += 1)
		{
			R->x[i][j] = A->x[i][j] + B->x[i][j];
		}
	}
}

void mulM(struct M* A, struct M* B, struct M* R)
{
	int i, j, k;
	int t;
	for (i = 0; i < SIZE; i += 1)
	{
		for (j = 0; j < SIZE; j += 1)
		{
			t = 0;
			for (k = 0; k < SIZE; k += 1)
			{
				t += A->x[i][k] * B->x[k][j];
			}
#if TRUNCATE==1
			R->x[i][j] = t;
#else
			R->x[i][j] = bound(t, 0, 1 << 30);
#endif
		}
	}
}

void printM(struct M* A)
{
	int i, j;
	for (i = 0; i < SIZE; i += 1)
		for (j = 0; j < SIZE; j += 1)
			printf("M[%4d,%4d] = %10d\n", i, j, A->x[i][j]);
}

void outsource(struct bank_A *a, struct M *B, struct bank_R *r)
{
	struct M* A = &(a->A);
	struct M* R = &(r->R);
	struct M tmp;

	mulM(A, B, R);
	// A*A*A*A*A*B,  ~32 splittable equations.
#if 0
	mulM(A, R, &tmp);
	mulM(A, &tmp, R);
	mulM(A, R, &tmp);
	mulM(A, &tmp, R);
#endif
}

#ifndef MQAP
int main() {
	struct bank_A A = { { 5, 6, 7, 8 } };
	struct M B = { 1, 2, 3, 4 };
	struct bank_R C;
	outsource(&A, &B, &C);
	printM(&C.R);
	return 0;
}
#endif

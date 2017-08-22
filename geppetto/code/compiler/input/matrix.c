#define GEPPETTO_NAME     "matrix"  // Must be defined before including geppetto.h
#define GEPPETTO_VERIFY_H "build/"GEPPETTO_NAME".h"        
#define GEPPETTO_VERIFY_C "build/"GEPPETTO_NAME".verify.c" 
#define NUM_BANKS 3
#include "geppetto.h"

#include "pinocchio.h"
#ifndef MQAP
#include <stdio.h>
#include <string.h>
#include "pinocchio.c"
#endif

// size of the matrix (plain multiplication takes SIZE^3 multiplications)
//#define SIZE 3
#define SIZE 40
//#define SIZE 80

// number of bits in the input matrix (controls the need for truncation)
#define WIDTH 10

#define PRINT (SIZE < 11)
#define BASE 0

typedef struct { unsigned int a[SIZE][SIZE]; } matrix;

void fill(matrix *M, int seed) {
	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++) {
			M->a[i][j] = ((i * 13 + j * 17 + SIZE * SIZE * 7) * seed) % (1 << WIDTH);
		}
}

// limits the need & costs of truncation on the result
void boundm(matrix *M) {
	if (WIDTH < 31)
		for (int i = 0; i < SIZE; i++)
			for (int j = 0; j < SIZE; j++)
				M->a[i][j] = bound(M->a[i][j], 0, (1 << WIDTH) - 1);
}

void printm(matrix *M) {
#if PRINT
	for (int j = 0; j < SIZE; j++) printf(" ----------");
	printf("\n");
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++)
			printf("| %8d ", M->a[i][j]);
		printf("|\n");
	}
	for (int j = 0; j < SIZE; j++) printf(" ----------");
	printf("\n\n");
#endif
}

void mul(matrix *M, matrix *N, matrix *R) {
	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++) {
			R->a[i][j] = 0;
			for (int k = 0; k < SIZE; k++)
				R->a[i][j] += M->a[i][k] * N->a[k][j];
		}
}

BANK(bankM, matrix);
BANK(bankN, matrix);
BANK(bankR, matrix);

bankR job(bankM m, bankN n) {
	matrix M, N, R;
	load_bankM(m, &M); 
	load_bankN(n, &N); 
	boundm(&M);
	boundm(&N);
	mul(&M, &N, &R);
	return save_bankR(&R);
}

bankR halfjob(bankN n) {
	matrix M, N, R;
	fill(&M, 5425);
	load_bankN(n, &N);
	boundm(&M);
	boundm(&N);
	mul(&M, &N, &R);
	return save_bankR(&R);
}

int main() {
	init();

	matrix M, N, R;
	printf("matrix multiplication takes %d *\n\n", SIZE*SIZE*SIZE);
	fill(&M, 5425);
	fill(&N, 7653);
	printm(&M);
	printm(&N);

	bankM m = save_bankM(&M);
	bankN n = save_bankN(&N);
	bankR r;
	r = OUTSOURCE(job, m, n);
	load_bankR(r, &R);
	printm(&R);

	r = OUTSOURCE(halfjob, n);
	load_bankR(r, &R);
	printm(&R);
}

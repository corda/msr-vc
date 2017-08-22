#include "pinocchio.h" 
#ifndef MQAP
#include <stdio.h>
#include <string.h>
#include "pinocchio.c"
#endif

#define GEPPETTO_NAME "unit-test"
#define GEPPETTO_VERIFY_H "build/unit-test.verify.h"
#define GEPPETTO_VERIFY_C "build/unit-test.verify.c"
#define NUM_BANKS 2
#include "geppetto.h"

// testing nested structures, globals, and comparisons
// (consider adding more c-semantics unit tests)

typedef struct {
	int x;
	int y;
	int z;
} triple;

typedef triple Base[3];


BANK(IN,triple)
BANK(OUT, int)

// --------------------------------------------- 
void prod(triple* r, triple* p, triple* q) {
	r->x = p->x*q->x;
	r->y = p->x*q->y;
	r->z = p->x*q->z;
}

int inner(triple* p) {
	return p->x*p->x + p->y*p->y + p->z * p->z;
}

static Base GBase;

void fill() {
	GBase[0].x = 32;
	GBase[0].y = 12;
	GBase[0].z = 43;
	GBase[1].x = 21;
	GBase[1].y = 13;
	GBase[1].z = 54;
	GBase[2].x = 23;
	GBase[2].y = 87;
	GBase[2].z = 39;
}

OUT test_struct(IN in) {
	triple t[1], r[1];
	load_IN(in, t);

	int s = 0;
	fill();
	prod(r, GBase, t);
	s += inner(r);
	prod(r, GBase+1, t);
	s += inner(r);
	prod(r, GBase+2, t);
	s += inner(r);

	return(save_OUT(&s));
}

void outer_struct() {
	triple t;
	t.x = 2;
	t.y = 3;
	t.z = 4;
	IN in = save_IN(&t);
	OUT out = OUTSOURCE(test_struct, in);
	int r;
	load_OUT(out, &r);
	printf("struct test: s = %d.\n", r);
}

// --------------------------------------------- 
OUT test_cmps(IN in) {
	triple t; load_IN(in, &t);
	int x = bound(t.x, -200, 200);
	int y = bound(t.y, -200, 200);
	printf("%3d == %3d is %d\n", x, y, (x == y));
	printf("%3d != %3d is %d\n", x, y, (x != y));
	printf("%3d <= %3d is %d\n", x, y, (x <= y));
	printf("%3d <  %3d is %d\n", x, y, (x < y));
	int s = t.z;
	return(save_OUT(&s));
}

void cmps(int x, int y) {
	triple t;
	t.x = x;
	t.y = y;
	t.z = 0;
	IN in = save_IN(&t);
	OUT out = OUTSOURCE(test_cmps, in);
}

void outer_cmps() {

	cmps(10, 1);
	cmps(10, 10);
	cmps(10, 100);
	cmps(-5, -7);
	cmps(-5, -5);
	cmps(-5, 6);
}

// --------------------------------------------- 
int main() {
	init();
	outer_struct(); 
	outer_cmps();
	return 0;
}

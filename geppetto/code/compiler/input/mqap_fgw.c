#define GEPPETTO_NAME     "mqap_fgw"  // must be defined before including geppetto.h
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

// a tiny multi-qap example, maybe used as a running example

struct input {
	int x;
	int y;
};

BANK(INPUT, struct input)
BANK(INTER, int)
BANK(OUTPUT, int)

const int a = 4;
#define NOTRUNC(T) (bound((T), 0, 1 << 30))

INTER f(INPUT xy)
{
	struct input in[1]; load_INPUT(xy, in);
	int t = NOTRUNC(in->x * in->x + in->y + a);
	return(save_INTER(&t));

}

OUTPUT g(INTER z) 
{
	int r; load_INTER(z, &r);
	r = NOTRUNC(r * r - 1);
	return(save_OUTPUT(&r));
}

int main()
{
	init();

	struct input sample;
	sample.x = 2;
	sample.y = 3;

	INPUT xy = save_INPUT(&sample);
	printf("calling f({%d,%d}).\n", sample.x, sample.y);
	INTER z = OUTSOURCE(f, xy);

#ifdef OPAQUE 
	printf("calling g on intermediate result.\n");

#else
	int s; load_INTER(z, &s);
	printf("calling g(%d).\n", s);
#endif

	OUTPUT r = OUTSOURCE(g, z);
	int t; load_OUTPUT(r, &t);
	printf("the result is %d.\n", t);
}

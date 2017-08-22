#include "pinocchio.h"
// a tiny multi-qap example

struct bank_input {
	int x;
	int y;
};

struct bank_inter {
	int z;
};

struct bank_output {
	int r;
};

#define MAX 100

void f(int* a, struct bank_input* in, struct bank_inter* inter)
{
	int t = in->x * in->x;
	inter->z = bound(t + in->y + *a, 0, 1 << 30);

}

void g(struct bank_output* out, struct bank_inter* inter)
{
	out->r = bound(inter->z * inter->z - 1, 0, 1 << 30);
}

#if NATIVE
//#include <stdio.h>
//#include <string.h>


int main()
{
	struct bank_input in = { 2, 3 };
	struct bank_inter inter = { 0 };
	struct bank_output out = { 0 };
	int a = 4;
	printf("in = {%3d,%3d }, inter = {%3d }, out = {%3d }\n", in.x, in.y, inter.z, out.r);
	printf("calling f(%d,&in,&inter)\n", a);
	f(&a, &in, &inter);
	printf("in = {%3d,%3d }, inter = {%3d }, out = {%3d }\n", in.x, in.y, inter.z, out.r);
	printf("calling g(&inter,&out)\n");
	g(&out, &inter);
	printf("in = {%3d,%3d }, inter = {%3d }, out = {%3d }\n", in.x, in.y, inter.z, out.r);
}
#endif


//#include <stdio.h>
//#include <string.h>
#include "pinocchio.h"
#include "trivial-ifc.h"
#ifndef PARAM
#define PARAM 0
#endif
// a bunch of simple test cases for debugging our restricted semantics for C.
// use the config parameter to select one of them.

int ten(int a)
{
	int i;
	int r;
	r = 0;
	for (i = 0; i < 10; i++) { r = r + a; }
	return r;
}

int hundred(int* a, unsigned int foo)
{
	int i;
	int r;
	r = 0;
	for (i = 0; i < 100; i++) { r = foo * r + *a; }
	return r;
}

void dummy(int a)
{
	int r;
	r = a * a;
}

void sample(int a, int b, struct bank_Out* z)
{
	z->r = a + b;
}
// arguments at +1 +3 +5 

void sample2(int a, int b, int * z)
{
	int r = a + *z;
	*z = b + r;
}
// arguments at +1 +3 +5 

void idem(int *config, struct bank_In *input, struct bank_Out *output)
{
	int z = *config;
	output->r = input->x * input->y;
}

int shaloop(int b, int c, int d)
{
	// TODO: check optimality
	// expected: 3*32x for splitting the inputs
	return ((b & c) | ((~b) & d))          // expected at most 3*32x, can be coded in 2*32 as b*c + (1-b)*d
		+ (b ^ c ^ d)                     // expected 2*32
		+ ((b & c) | (b & d) | (c & d));  // expected at most 5*32x; we get 4*32.
}

void outsource(struct Config *config, struct bank_In *input, struct bank_Out *output)
{
	int x = input->x;
	int y = input->y;
	int z = input->z;

	if (PARAM == 0) { // testing comparisons (many inlined by clang)
		int z;
		z = 0;
		z += (3 < 4); z *= 2;
		z += (4 < 4); z *= 2;
		z += (5 < 4); z *= 2;
		z += (3 == 4); z *= 2;
		z += (4 == 4); z *= 2;
		z += (5 == 4); z *= 2;
		z += (3 <= 4); z *= 2;
		z += (4 <= 4); z *= 2;
		z += (5 <= 4); z *= 2;
		int maxint = 2147483647;
		int minint = (0 - maxint) - 1;
		z += (minint < maxint);
		output->r = z;
		z = 0;
		z += (x < y); z *= 2;
		z += (y < z); z *= 2;
		z += (y < x); z *= 2;
		z += (x >= y); z *= 2;
		z += (x >= y); z *= 2;
		z += (y >= x); z *= 2;
		printf("comparing %d and %d\n", x, y);
		z += (x == y); z *= 2;
		z += (x == z); z *= 2;
		z += (x != y); z *= 2;
		output->r *= z;
	}
	if (PARAM == 1) { // testing arithmetic and bitwise operations  
		output->r = 7 * x * y + 11 + z; // 1x
		output->r += x & y;              // 32*(2 for decomposing + 3 for computing) = 160x
		output->r += x ^ y;
		output->r += x | y;              // 67-bit final truncation
	}
	if (PARAM == 2) {
		output->r = shaloop(x, y, z);
	}
	if (PARAM == 3) { // this test depends on val being unsigned! 
		unsigned int val = x ^ 3456234563;
		int amount = 5;
		int i;
		for (i = 0; i < 5 * 32; i++) val = (val << amount) | (val >> (32 - amount));
		output->r = (val ^ 3456234563) - x + y + z;
	}

	if (PARAM == 4) { // testing loops & function calls
		output->r = 1 + ten(y) + hundred(&x, 1) + z;
	}

	if (PARAM == 5) { // testing truncations
		output->r = x * y + z + config->a;
	}

	if (PARAM == 6) { // adapted from a regular example for partitioning; 
		int a[10], b[10];
		int i, n;
		int t = y + z + 3;
		output->r = 0;
		for (i = 0; i < 10; i++) a[i] = (x & (1 << i)) >> i;
		for (i = 0; i < 10; i++) b[i] = (t & (1 << i)) >> i;
		// for (n = 0; n < 500; n++)
			for (i = 0; i < 10; i++) {
				// printf("a[%d]=%d, b[%d]=%d\n", i, a[i], i, b[i]);
				output->r += a[i] * b[i];
			}

	}
	if (PARAM == 7) {
		// testing rotations and shifts a la SHA256
		unsigned int val = x + y + z;
		int amount;
		output->r = 0;
		for (amount = 1; amount < 16; amount++) {
			unsigned int r = (val >> amount) | (val << (16 - amount));
			unsigned int s = (val >> amount);
			output->r += r + s;
		}
	}
	if (PARAM == 8) { // basic example for partitioning
		int a = x * y;
		int b = x * z;
		int c = y* x;
		int d = a* c;
		int e = b * b;
		int f = d * e;
		int g = c * d;
		output->r = bound(f + g, 0, 1 << 30);
	}
	printf("test %d: %d.\n", PARAM, output->r);
}

#if _WIN32
int main() {
	struct Config cfg;
	struct bank_In input;
	struct bank_Out output;
	cfg.a = 2;
	input.x = 3;
	input.y = 5;
	input.z = 7;
	outsource(&cfg, &input, &output);
	return 0;
}
#endif



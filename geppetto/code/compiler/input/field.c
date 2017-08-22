//#include "pinocchio.h"
#include "../src/arith-qap/arith/Qapped/PrimitiveIfc.h"

struct bank_Input {
	long int a;
	long int b;
};

struct bank_Output {
	long int r;
};

void elem_width(Elem out, Elem in, int b);

void outsource(struct bank_Input *in)
{
	Elem x; elem_init(x);
	Elem y; elem_init(y);
	Elem z; elem_init(z);
	elem_set_ui(x, in->a);
	elem_set_ui(y, in->b); // elem_set_u32(y, 0,0,0,0,0,0,0,2);

	// z = a * b
	// x = b * a
	// r = (x == z)
	elem_mul(z, x, y);
	
	Elem t; elem_init(t);
	elem_width(t, z, 9*8 + 2);
	// elem_mul(x,y,x);
	printf("\n x = ");  elem_dbg_print(x);
	printf("\n y = "); elem_dbg_print(y);
	printf("\n z = "); elem_dbg_print(z);
	// printf("\n z = ...%d%d%d%d\n", elem_get_bit(z, 3), elem_get_bit(z, 2), elem_get_bit(z, 1), elem_get_bit(z, 0));
	printf("\n t = "); elem_dbg_print(t);
	printf("\n t = ...%d%d%d%d\n", elem_get_bit(t, 3), elem_get_bit(t, 2), elem_get_bit(t, 1), elem_get_bit(t, 0));

	// currently not supported on secrets:
    // out->r = elem_cmp(x,z);
}

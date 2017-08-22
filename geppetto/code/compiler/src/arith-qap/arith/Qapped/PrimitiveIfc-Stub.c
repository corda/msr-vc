#include "PrimitiveIfc.h"

struct s_elem
{
	int e;
};

void elem_static_init()
{
}

Elem* elem_init()
{
	return (Elem*) 0;
}

void elem_free(Elem* e)
{
}

void elem_modinv(Elem* out, Elem *a)
{
}

void elem_mul(Elem* out, Elem* a, Elem* b)
{
}

void elem_add(Elem* out, Elem* a, Elem* b)
{
}

void elem_sub(Elem* out, Elem* a, Elem* b)
{
}

void elem_div2(Elem* out, Elem* a)
{
}

int elem_cmp(Elem* a, Elem* b)
{
	return 0;
}

void elem_mul3(Elem* out, Elem* a)
{
}

void elem_set_ui(Elem* out, uint64_t x)
{
}

void elem_set_u32 (Elem* out,
	uint32_t a7, uint32_t a6, uint32_t a5, uint32_t a4,
	uint32_t a3, uint32_t a2, uint32_t a1, uint32_t a0)
{
}

void elem_copy(Elem* out, Elem* a)
{
}

void elem_neg(Elem* out, Elem* a)
{
}

void elem_rand(Elem* out)
{
}

void elem_dbg_print(Elem* a)
{
}

void elem_assert(int condition)
{
}

int elem_get_bit(Elem* a, int bit)
{
}

int elem_eq_zero(Elem* a)
{
}
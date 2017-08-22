// An implementation of PrimitiveArith on top of Arith Fp.
// This uses FpGen, so that we can instantiate it with the P prime
// in the encoding math; the sa-field.c module uses the same
// FpGen instatiated with Q, the size of the base curve group.

//#if NOASSERT
//#else
#include <assert.h>
//#endif

#include "PrimitiveIfc.h"
#include "Fp/FpGen.h"

FpGen_config_t g_cfg, *cfg = &g_cfg;

void elem_static_init()
{
	uint64_t *m;
	int n;
#ifndef TINY
	n = 4;
	m = (uint64_t *)malloc(n * sizeof (uint64_t));
	m[3] = 0x2523648240000001;
	m[2] = 0xBA344D8000000008;
	m[1] = 0x6121000000000013;
	m[0] = 0xA700000000000013; // 2523648240000001BA344D80000000086121000000000013A700000000000013
	/* This is the 254-bit BN prime p = 3 (mod 4). */
	/* A point on the BN curve is composed of two values mod this prime */
#else
	n=1;
	m = (uint64_t *) malloc (n * sizeof (uint64_t));
	m[0] = 0x15C9B0796B71DB;
	/* This is the 53-bit BN prime p = 3 (mod 4). */
#endif // TINY
	FpGen_initialize_config(cfg, m, n);
}

void elem_init(Elem e)
{
	FpGen_init(cfg, e->base_fp);
}

void elem_free(Elem e)
{
	FpGen_free(cfg, e->base_fp);
}

void elem_modinv(Elem out, Elem a)
{
	FpGen_modinv(cfg, out->base_fp, a->base_fp);
}

void elem_mul(Elem out, Elem a, Elem b)
{
	FpGen_mul(cfg, out->base_fp, a->base_fp, b->base_fp);
}

void elem_add(Elem out, Elem a, Elem b)
{
	FpGen_add(cfg, out->base_fp, a->base_fp, b->base_fp);
}

void elem_sub(Elem out, Elem a, Elem b)
{
	FpGen_sub(cfg, out->base_fp, a->base_fp, b->base_fp);
}

void elem_div2(Elem out, Elem a)
{
	FpGen_div2(cfg, out->base_fp, a->base_fp);
}

int elem_cmp(Elem a, Elem b)
{
	return FpGen_cmp(cfg, a->base_fp, b->base_fp);
}

int elem_eq_zero(Elem a)
{
	// This is kind of a ridiculous implementation.
  Elem z = { NULL };
  elem_init(z);
	elem_set_ui(z, 0);
	int rc = FpGen_cmp(cfg, a->base_fp, z->base_fp);
	elem_free(z);
	return rc==0;
}

void elem_mul3(Elem out, Elem a)
{
	FpGen_mul3(cfg, out->base_fp, a->base_fp);
}

void elem_set_ui(Elem out, uint64_t x)
{
	FpGen_set_ui(cfg, out->base_fp, x);
}

void elem_set_u32 (Elem out,
	uint32_t a7, uint32_t a6, uint32_t a5, uint32_t a4,
	uint32_t a3, uint32_t a2, uint32_t a1, uint32_t a0)
{
#ifndef TINY
	uint64_t storage[4];
	storage[3] = ((uint64_t)a7)<<32 | a6;
	storage[2] = ((uint64_t)a5)<<32 | a4;
	storage[1] = ((uint64_t)a3)<<32 | a2;
	storage[0] = ((uint64_t)a1)<<32 | a0;
	FpGen_set(cfg, out->base_fp, storage);

#else 
  uint64_t storage[1];
  storage[0] = ((uint64_t)a1) << 32 | a0;
  FpGen_set(cfg, out->base_fp, storage);
#endif // TINY
}

void elem_copy(Elem out, Elem a)
{
	FpGen_copy(cfg, out->base_fp, a->base_fp);
}

void elem_neg(Elem out, Elem a)
{
	FpGen_neg(cfg, out->base_fp, a->base_fp);
}

void elem_rand(Elem out)
{
	FpGen_rand(cfg, out->base_fp);
}

void elem_dbg_print(Elem a)
{
	FpGen_print(cfg, a->base_fp);
}

void elem_assert(int condition)
{
	assert(condition);
}

int elem_get_bit(Elem a, int bit)
{
#ifndef TINY
	uint64_t ta[4];
	FpGen_get(cfg, ta, a->base_fp);
	int rc = (((ta[bit >> 6] >> (bit & 0x3f))) & 0x01) != 0;
	return rc;
#else 
	assert(bit < 54);
	uint64_t ta[1];
	FpGen_get(cfg, ta, a->base_fp);
	int rc = ((ta[0] >> bit) & 0x01) != 0;
	return rc;
#endif // TINY
}

void elem_to_mont(Elem a)
{
  FpGen_set(cfg, a->base_fp, a->base_fp->limbs);
}
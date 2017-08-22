#pragma once

#include <stdint.h>

//struct s_elem;
//typedef struct s_elem Elem;

#ifdef MQAP
typedef int elem;
// typedef struct { int e; } elem;

#else
#include "Fp/FpGen.h"
typedef struct 
{
  FpGen_t base_fp;
} elem;
#endif
typedef elem Elem[1];

void elem_static_init(void);

void elem_init(Elem e);
void elem_free(Elem e);

void elem_modinv(Elem out, Elem a);
void elem_div2(Elem out, Elem a);
void elem_mul3(Elem out, Elem a);
void elem_copy(Elem out, Elem a);
void elem_neg(Elem out, Elem a);

void elem_add(Elem out, Elem a, Elem b);
void elem_sub(Elem out, Elem a, Elem b);
void elem_mul(Elem out, Elem a, Elem b);

int elem_cmp(Elem a, Elem b);	// we shouldn't need this, which is good because it's hard to implement inexpensively.
int elem_eq_zero(Elem a);	// this should have a cheap constraint-based implementation. Returns 1 if a==0
void elem_set_ui(Elem out, uint64_t x);
void elem_set_u32 (Elem out,
	uint32_t a7, uint32_t a6, uint32_t a5, uint32_t a4,
	uint32_t a3, uint32_t a2, uint32_t a1, uint32_t a0);
void elem_rand(Elem out);

void elem_dbg_print(Elem a);
void elem_assert(int condition);
int elem_get_bit(Elem a, int bit);

#ifndef MQAP
void elem_to_mont(Elem a);
#endif








                                   
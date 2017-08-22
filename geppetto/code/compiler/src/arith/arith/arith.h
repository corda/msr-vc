/* arith.h
 * by Joppe W. Bos and Michael Naehrig (mnaehrig@microsoft.com), 
 * Cryptography Research Group, MSR Redmond
 * 
 * This file is part of the ARITH library version 0.01
 *
 * DISCLAIMER:
 * This is RESEARCH code.
 * Parts of the code might be incomplete.
 * Please use at your own risk.
 * DO NOT USE IN PRODUCTION CODE!!!
 * This code is still under active development.
 */

#define _CRT_SECURE_NO_WARNINGS

// Hunting for memory leaks
//#define MEM_LEAK 1
#ifdef MEM_LEAK
#define CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
// Print file names for the locations of memory leaks
#define DEBUG_NEW new(_NORMAL_BLOCK ,__FILE__, __LINE__)
#define new DEBUG_NEW
// See https://msdn.microsoft.com/en-us/library/x98tx3cf.aspx for useful tips
#endif
// Done with code for hunting memory leaks

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* Enable the Windows intrinsics so we get asm support on x86_64 */
#include <intrin.h>

#include "arith_error.h"

void binary_gcd (uint64_t *S, uint64_t *U, uint64_t *V, uint64_t *mem, int n);
void divrem (uint64_t *r, uint64_t *a, uint64_t *b, int n);


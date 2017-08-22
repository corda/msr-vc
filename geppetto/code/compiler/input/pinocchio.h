/* primitive C library for Pinocchio programs (with runtime support) */

#pragma once

void zeroAssert(int x); 
// when x <> 0, the prover fails.
// convenient for proving/verifying conditions on intermediate values.

void print(int x); 
// prints x on the prover; does nothing for the verifier

int runtime(int x);
// returns 1 if x is available only at runtime, 0 otherwise.
// currently not natively available.

int bound(int x, long int min, long int max);
// UNSAFE, unless followed by a binary decomposition
// returns x unchanged, and assumes it is within min..max
// convenient for introducing pre-conditions on inputs & for debugging

int arith(int x); 
// forces x's bits to be merged  (to ensure sharing)

int endorse(int x); 
// take a compile time constant and turns it into an input. 
// safe, but inefficient; convenient for writing tests. 

int extract(int x);
// treat a runtime value as an untracked compile-time constant 
// (use within endorse; many limitations)

int select(int x, int a, int b);
// implements x ? a : b (use to prevent clang optimizations on source selects)

int nRoot(void); 
// returns the number of equations generated so far; used for performance debugging and adaptive QAP filling
// *not* correctly implemented natively.

#if MQAP == 1
int printf(const char*, ...);
#endif

#ifdef NATIVE
void timer_start(); 
void timer_end(); 
void timer_print(); 
#include <stdio.h>

#else // !NATIVE
void timer_start() { ; }
void timer_end()   { ; }
void timer_print() { ; } 
#endif // NATIVE

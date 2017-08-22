#pragma once

/* primitive C library for Pinocchio programs (with runtime support) */

#include <stdlib.h> //MK needed to be able to use pinocchio.c in Quapped
#include "pinocchio.h"

void zeroAssert(int x) {
    if (x != 0) {
        printf("\n\n\n zeroAssert failed.\n\n\n");
        //assert(0);
        exit(-1);
    }
}

static int np = 0;

void print(int x) { printf("%6d: %10d  %10u\n",np++,x,x); }

int bound(int x, long int min, long int max) { 
	// if (x < min | x > max) exit(-2);
	return x;
}

int arith(int x) { return x; }
// forces x's bits to be merged  (to ensure sharing)

int endorse(int x) { return x; } 

int extract(int x) { return x; }

int select(int x, int a, int b) { return (x ? a : b);  }

int nRoot() { return 0;  }

#ifdef NATIVE
unsigned long long _timer_start;
unsigned long long _timer_end;

unsigned long long NativeRDTSC() {
  unsigned long long rv;
  __asm__ (
      "push           %%ebx;"
      "cpuid;"
      "pop            %%ebx;"
      :::"%eax","%ecx","%edx");
  __asm__ __volatile__ ("rdtsc" : "=A" (rv));
  return (rv);
}
 
void timer_start() { _timer_start  = NativeRDTSC(); }
void timer_end()   { _timer_end    = NativeRDTSC(); }
void timer_print()  { 
  float seconds = (_timer_end - _timer_start) / (float)(3.6*1000*1000*1000);
  float ms = seconds * 1000;
  printf("Native execution time: %f ms\n", ms); 
}

#endif // NATIVE

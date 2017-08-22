/* Fp/main.c
 * by Joppe W. Bos and Michael Naehrig (mnaehrig@microsoft.com), 
 * Cryptography Research Group, MSR Redmond, 2014
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

#include "Fp.h"
#include "../Fp2/Fp2.h"
#include "../Fp4/Fp4.h"
#include "../Fp6/Fp6.h"
#include "../Fp12/Fp12.h"
#include "../Fp12/Fp12overFp4.h"
#include "../uint.h"
#include "../arith.h"
#include "../pairing/pairing.h"
#include <math.h>


static void show_cpu_info () {
  int CPUInfo[4] = {-1};
  char CPUBrandString[0x40];
  unsigned int nExIds, i;

  /* Calling __cpuid with 0x80000000 as the InfoType argument
   * gets the number of valid extended IDs. */
  __cpuid (CPUInfo, 0x80000000);
  nExIds = CPUInfo[0];
  memset (CPUBrandString, 0, sizeof(CPUBrandString));

  /* Get the information associated with each extended ID. */
  for (i=0x80000000; i<=nExIds; ++i) {
    __cpuid (CPUInfo, i);
        
    /* Interpret CPU brand string and cache information. */
    if  (i == 0x80000002)
      memcpy (CPUBrandString, CPUInfo, sizeof(CPUInfo));
    else if  (i == 0x80000003)
      memcpy (CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
    else if  (i == 0x80000004)
      memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
  }

  if  (nExIds >= 0x80000004) 
    printf ("CPU information:\n\t%s\n", CPUBrandString);
}

__declspec(dllexport) test() {
  show_cpu_info();
}

void test_bingcd () {
  /* Format (c,a,p) where c = a^-1 mod p. */
  uint64_t limb2[2][3][2] = {{ {0x06974CABE4B37700, 0x60BB517F57385EC7}, {0xFFF0000000000010, 0x7FFFFFFFFFFFFFFF}, {0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF} },
                             { {0xEA67F1935A3A5EAB, 0x62781A36C4FD6C24}, {0x4E9B15698C86AEF9, 0x13701878EAE3A561}, {0x58D2754B39BCA883, 0x7647F3C382902D2F} }};
  uint64_t limb10[2][3][10] = {{ {0xD7CEF0F7CBB0FFDD, 0xF07F534B40493221, 0xF5F760F39B4D535B, 0x728BC5E59B9674C4, 0x1A450DAE41A4924D, \
                                  0x6049E5A6DD6B194F, 0xD1BD8AB684C80904, 0xCF242A145FB06554, 0x48BC5F2BAF6AECF0, 0xA7C7132DF67A451A},
                                 {0xFFFFFFFFFFFFFFFF, 0x0000000000FFFFFF, 0xFFFFFFFC00000000, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, \
                                  0x000000000000003F, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFFFFFFFFFFFFFFF0},
                                 {0x0000000000000295, 0xFFFFFFFFFFFFFFFF, 0x000000000000000F, 0x0000000000000000, 0x0000000000000000, \
                                  0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFFF0000000000000, 0xFFFFFFFFFFFFFFFF} },
                               { {0xEB057D230FF19333, 0x76FFD7BD774C98A4, 0x58700DF7B18AC6B6, 0x197294E7A4213F67, 0xD9CD720A011738B2, \
                                  0x6140E72C45045C88, 0xAC3558C5440F8284, 0xB0EB0ECDCE16172A, 0x9D289D1D1F44EA31, 0x43AE549C3FA9AE98},
                                 {0xA72D8AB486435573, 0x89B80C3C7D6FD2D0, 0xEB483F1E97F44237, 0x9A7F2C8227B578C4, 0xE1008176408D879C, \
                                  0xC3F2603055C23E74, 0x0A7E892A85CC8DC9, 0xE70DFB0197B9D149, 0x71081ED0996B318A, 0x65F0AC2B5774C1AA},
                                 {0x58D2754B39BCAA8D, 0x7647F3C382902D2F, 0x14B9C0E1680BBDC8, 0x6580D17DD84A873B, 0x1EFF7E89BF767863, \
                                  0x3C0D9FCFAA3DC18B, 0xF58174D57A337236, 0x18F204FE6846AEB6, 0x8EF7612F6694CE75, 0x9A0F53D4A88B3E55} }};  
  uint64_t res[10], *mem;
  int i, j, passed=0;

  mem = NULL;
  mem = (uint64_t *) malloc (9 * 11 * sizeof (uint64_t));
  if (mem == NULL) {
    printf ("Memory alloc error\n");
  }

  for (j=0; j < 2; j++) {
    binary_gcd (res, limb2[j][1], limb2[j][2], mem, 2);
  
    for (i=0; i < 2; i++) {
    if (res[i] != limb2[j][0][i]) {
      printf ("Failure in modular inverse2 in run %d in limb %d out of 2\n", j, i);
      printf ("%08X%08X\n", (uint32_t) (res[i]>>(uint64_t)32), (uint32_t)res[i]);
      printf ("%08X%08X\n", (uint32_t) (limb2[j][0][i]>>(uint64_t)32), (uint32_t)limb2[j][0][i]);
      return;
      }
    }
    passed++;
  }

  for (j=0; j < 2; j++) {
    binary_gcd (res, limb10[j][1], limb10[j][2], mem, 10);
  
    for (i=0; i < 2; i++) {
    if (res[i] != limb10[j][0][i]) {
      printf ("Failure in modular inverse10 in run %d in limb %d out of 10\n", j, i);
      printf ("%08X%08X\n", (uint32_t) (res[i]>>(uint64_t)32), (uint32_t)res[i]);
      printf ("%08X%08X\n", (uint32_t) (limb10[j][0][i]>>(uint64_t)32), (uint32_t)limb10[j][0][i]);
      return;
      }
    }
    passed++;
  }
  printf ("All %d modular inverse tests passed\n", passed);
  free (mem);
}

void karatsuba_threshold () {
  uint64_t *c, *a, *b, *tmp;
  int size, i, j;
  uint64_t start, end;
  uint32_t schoolbook_c, karatsuba_c, karatsuba_c_min, threshold;

  for (size = 10; size < 100; size++) {
    karatsuba_c_min = UINT32_MAX;

    a   = (uint64_t *) malloc (  size * sizeof (uint64_t));
    b   = (uint64_t *) malloc (  size * sizeof (uint64_t));
    tmp = (uint64_t *) malloc ( (2*(size-4+3*((uint64_t)(log((double)size-3.0)/log(2.0))))) * sizeof (uint64_t));
    c   = (uint64_t *) malloc (2*size * sizeof (uint64_t));
  
    for (i=0; i < size; i++) {
      a[i] = rand ();
      b[i] = rand ();
    }

    start = __rdtsc();
    for (i=0; i < 1000; i++) {
      mul_schoolbook (c,a,b,size);  
    }
    end = __rdtsc();
    end = end - start;

    schoolbook_c = (uint32_t) ((double) end / i);

    for (j=8; j < size; j++) {
      SCHOOLBOOK_THRESHOLD = j;
      start = __rdtsc();
      for (i=0; i < 1000; i++) {
        //mul_karatsuba (c,a,b,size);  
        uint_karatsuba (c,a,b,tmp,size,size);
      }
      end = __rdtsc();
      end = end - start;

      karatsuba_c = (uint32_t) ((double) end / i);
      if (karatsuba_c < karatsuba_c_min) {
        karatsuba_c_min = karatsuba_c;
        threshold = j;
      }
    }

    if (karatsuba_c_min < schoolbook_c) {
      printf ("Average number of cycles for %d-limb schoolbook: %d\n", size, (uint32_t) schoolbook_c);
      printf ("Average number of cycles for %d-limb Karatsuba with %d-limb threshold: %d\n", size, threshold, karatsuba_c_min);
      printf ("\n");
    } else {
      printf ("For %d-limb numbers schoolbook is faster (%d vs %d for threshold %d).\n", size, schoolbook_c, karatsuba_c_min, threshold);
    }

    free (a);
    free (b);
    free (c);
    free (tmp);
  }
}

void test_multiply () {
  uint64_t *c1, *c2, *a, *b, *tmp;
  int size, i, passed = 0, k;

  for (size = 3; size < 100; size++) {
    size = 2*size;
    a   = (uint64_t *) malloc (  size    * sizeof (uint64_t));
    b   = (uint64_t *) malloc (  size    * sizeof (uint64_t));
    tmp = (uint64_t *) malloc ( (2*(size-4+3*((uint64_t)(log((double)size-3.0)/log(2.0))))) * sizeof (uint64_t));
    c1  = (uint64_t *) malloc (2*size    * sizeof (uint64_t));
    c2  = (uint64_t *) malloc (2*size    * sizeof (uint64_t));
  
    size = size/2;

    for (k=0; k < 1000; k++) {
      for (i=0; i < size; i++) {
        a[i] = random_uint64_t();
        b[i] = random_uint64_t();
      }
    
      mul_schoolbook (c1,a,b,size);  
      //mul_karatsuba (c2,a,b,size);  
      for (i=0; i < 2*size; i++) c2[i] = 0;
      uint_karatsuba (c2,a,b,tmp,size,size);  
    
      if (uint_cmp (c1,c2,2*size) != 0) {
        printf ("Error in multiplication test for size %d.\n", size);
        Print (a,size);
        Print (b,size);
        Print (c1,2*size);
        Print (c2,2*size);
        return;
      }
      passed++;
    }
    free (a);
    free (b);
    free (c1);
    free (c2);
    free (tmp);
  }
  printf ("Multiplication phase passed %d tests successful\n", passed);
}

int Fp_stress_test () {
  Fp_t a, b, c, d;
  Fpl_t C;
  uint64_t *t;
  int i,j,k;

  for (i=2; i < 50; i++) {
    t = NULL;
    t = (uint64_t *) malloc (i * sizeof (uint64_t));
    if (t == NULL) return ERR_OUT_OF_MEMORY;

    for (j=0; j < i; j++) {
      t[j] = random_uint64_t();
    }
    t[0] |= 1; /* We want an odd modulus. */
    if (Fp_initialize_config (t,i) < 0) printf ("Stress test init error\n");    

    Fp_init (a);
    Fp_init (b);
    Fp_init (c);
    Fp_init (d);
    Fpl_init (C);

    for (k=0; k < 1000; k++) {
      for (j=0; j < i; j++) {
        t[j] =  random_uint64_t();
      }
      Fp_set (a, t);
        
      for (j=0; j < i; j++) {
        t[j] =  random_uint64_t();
      }
      Fp_set (b, t);

      a->limbs[i-1] %= (Fp_config.m->p[i-1]-1); /* Make sure the values are reduced mod p */
      b->limbs[i-1] %= (Fp_config.m->p[i-1]-1);

      Fp_add (c, a, b);
      Fp_sub (d, a, b);
      Fp_mul (c, c, d);

      Fp_mul (a, a, a);
      Fp_mul (b, b, b);
      Fp_sub (a, a, b);

      if (Fp_cmp (a, c) != 0) {
        printf ("Fp stress-test error 1\n");
        printf ("i = %d, k = %d\n", i, k);

        Print (a->limbs, i);
        Print (c->limbs, i);
        getchar ();
      } 

      Fp_mul (c, a, a);
      
      Fp_mul_no_red (C, a, a);
      Fpl_red (b, C);

      if (Fp_cmp (c, b) != 0) {
        printf ("Fp Mul stress-test error 2\n");
        printf ("i = %d, k = %d\n", i, k);

        Print (a->limbs, i);       
        getchar ();
      } 

      Fp_rand (a);
      Fp_set (b, a->limbs);
      Fp_get (c->limbs, b);

      if (Fp_cmp (a, c) != 0) {
        printf ("Fp stress test error: conversion\n");
        Print (Fp_config.m->p, i);
        Print (Fp_config.m->R2, i);
        getchar ();
      }

#if 0
      Fp_modinv (b, a);
      Fp_mul (c, a, b); /* Montmul(a*R, a^-1 * R) = R */
      Fp_get (d->limbs, c); /* d = 1 */
      m = (uint64_t) 0;
      for (j=1; j < i; j++) m |= d->limbs[j];

      if (m != 0 || d->limbs[0] != (uint64_t) 1) {
        printf ("Fp stress test: inversion error\n");
        Print (Fp_config.m->p, i);
        Print (Fp_config.m->R2, i);
        Print (Fp_config.m->R3, i);

        Print (a->limbs, i);
        Print (b->limbs, i);
        Print (d->limbs, i);

        getchar ();
      }
#endif

    }
    Fp_free (a);
    Fp_free (b);
    Fp_free (c);
    Fp_free (d);
    Fpl_free (C);
    Fp_free_config ();
    free (t);
    //free (tmp);
  }
  return 1;
}


int Fp2_stress_test (uint64_t *m, int n) {
  Fp2_t a, b, c, d;
  Fp2l_t A, B, C, D;
  int k;
  
  if (Fp_initialize_config (m,n) < 0) printf ("Stress test init error\n");
  if (Fp2_initialize_config () < 0) printf ("Stress test init error\n");

  if (Fp2_init (a) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2_init (b) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2_init (c) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2_init (d) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2l_init (A) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2l_init (B) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2l_init (C) < 0) printf ("Fp2 stress test memory error.\n");
  if (Fp2l_init (D) < 0) printf ("Fp2 stress test memory error.\n");

  for (k=0; k < 1000; k++) {
    Fp2_rand (a);
    Fp2_rand (b);
    
    Fp2_add (c, a, b);
    Fp2_sub (d, a, b);
    Fp2_mul (d, c, d);

    Fp2_squ (a, a);
    Fp2_squ (b, b);
    Fp2_sub (c, a, b);

    if (Fp2_cmpeq (c, d) != 1) {
      printf ("Fp2 stress-test error 1\n");
      getchar ();
    }

    Fp2_rand (a);
    Fp2_rand (b);

    Fp2_neg (c, a);
    Fp2_add(c, c, b);
    Fp2_sub(d, b, a);

    if (Fp2_cmpeq (c, d) != 1) {
      printf ("Fp2 stress-test error 2\n");
      getchar ();
    }

    Fp2_inv (b, a);
    Fp2_inv (b, b);

    if (Fp2_cmpeq (a, b) != 1) {
      printf ("Fp2 stress-test inversion error\n");
      
      Print (a->a0->limbs, 4);
      Print (b->a0->limbs, 4);

      Print (a->a1->limbs, 4);    
      Print (b->a1->limbs, 4);

      getchar ();
    }

    Fp2_rand (a);
    Fp2_rand (b);
    
    Fp2_add (c, a, b);
    Fp2_sub (d, a, b);
    Fp2_mul_no_red_o2 (C, c, d);
    Fp2_squ_no_red (A, a);
    Fp2_squ_no_red (B, b);
    Fp2l_sub_o2 (D, A, B);
    
    Fp2l_red(c, C);
    Fp2l_red(d, D);

    if (Fp2_cmpeq (c, d) != 1) {
      printf ("Fp2l stress-test error 3\n");
      getchar ();
    }
        
    Fp2_rand (a);
    Fp2_rand (b);

    Fp2_add (c, a, b);
    Fp2_sub (d, a, b);
    Fp2_mul_no_red_o1 (C, c, d, 1);

    Fp2_squ_no_red (A, a);
    Fp2_squ_no_red (B, b);
    Fp2l_sub_o1 (D, A, B, 1);

    Fp2l_red(c, C);
    Fp2l_red(d, D);

    if (Fp2_cmpeq (c, d) != 1) {
      printf ("Fp2l stress-test error 4\n");
      getchar ();
    }
    
    Fp2_rand (a);
    Fp2_rand (b);

    /*t = NULL;
    t = (uint64_t *) malloc (5 * sizeof (uint64_t));
    t[4] = 0x0000000000000000; 
    t[3] = 0x0000000000000000;
    t[2] = 0x0000000000000000;
    t[1] = 0x0000000000000000;
    t[0] = 0x0000000000000001;

    Fp_set(a->a0, t+1);
    Fp_set(a->a1, t);
    Fp_set(b->a0, t+1);
    Fp_set(b->a1, t+1);
    */

    Fp2_add (c, a, b);
    Fp2_sub (d, a, b);
    Fp2_mul_no_red_o1 (C, c, d, 2);

    Fp2_squ_no_red (A, a);
    Fp2_squ_no_red (B, b);
    Fp2l_sub_o1 (D, A, B, 1);

    Fp2l_red(c, C);
    Fp2l_red(d, D);

    if (Fp2_cmpeq (c, d) != 1) {
      printf ("Fp2l stress-test error 5\n");
      getchar ();
    }
  }

  Fp2_free (a);
  Fp2_free (b);
  Fp2_free (c);
  Fp2_free (d);
  Fp2_free_config ();
  Fp_free_config ();
  return 0;
}




#if 0
int main (int argc, char *argv) {
  Fp_t a, b, c;
  uint64_t *m;
  int n, ret, i;
  int pair_curve;

  printf ("Welcome to the sanity checker for the Finite Field Library (FFL)\n");

  printf ("\n");
  show_cpu_info ();
  printf ("\n");
  

  
  /*
    if (Fp_stress_test () < 0) {
      printf ("Fp stress test failed.\n");
    } else {
      printf ("Fp stress test successful.\n");
    }
  */

  m = NULL;
  if (PAIR_CURVE == BN12) {
    n = 4;
    m = (uint64_t *) malloc (4 * sizeof (uint64_t));
    m[3] = 0x2523648240000001;
    m[2] = 0xBA344D8000000008;
    m[1] = 0x6121000000000013;
    m[0] = 0xA700000000000013; 
    // 2523648240000001BA344D80000000086121000000000013A700000000000013
    /* This is the 254-bit BN prime p = 3 (mod 4). */

    if (Fp2_stress_test (m,n) < 0) {
      printf ("Fp2 stress test failed.\n");
    } else {
      printf ("Fp2 stress test successful.\n");
    }

    if (Fp6_stress_test (m,n) < 0) {
      printf ("Fp6 stress test failed.\n");
    } else {
      printf ("Fp6 stress test successful.\n");
    }
 
    if (Fp12_stress_test (m,n) < 0) {
      printf ("Fp12 stress test failed.\n");
    } else {
      printf ("Fp12 stress test successful.\n");
    }
  }

/*
  p[1] = 0x7FFFFFFFFFFFFFFF;
  p[0] = 0xFFFFFFFFFFFFFFFF;

  ret = Fp_initialize_config (p, 2);
  printf ("Return code: %d\n", ret);
*/

  //printf ("Testing the modular inversion code...\n");
  //test_bingcd ();
  //getchar ();

  //printf ("Testing the various multiplication routines...\n");
  //test_multiply ();
  //getchar ();


 // karatsuba_threshold ();
  getchar ();

  return 0;
}
#endif

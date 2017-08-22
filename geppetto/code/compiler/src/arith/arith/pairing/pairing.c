/* pairing/pairing.c
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

#include "../Fp/Fp.h"
#include "../Fp2/Fp2.h"
#include "../Fp4/Fp4.h"
#include "../Fp6/Fp6.h"
#include "../Fp6/Fp6overFp3.h"
#include "../Fp12/Fp12.h"
#include "../Fp12/Fp12overFp4.h"
#include "../Curve/Fp_weierstrass.h"
#include "../Curve/Fp2_weierstrass.h"
#include "pairing.h"
#include "../uint.h"
#include "../arith.h"
#include <math.h>

#define BIT(a,i) (((a)[i/64] >> (uint64_t) (i%64))&1)

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

#if 1
int main () {
  uint64_t *m, *b0, *b1;
  Fp_t b;
  Point_wp_t P0, Q0;
  Point_wp_t P0aff, Q0aff;
  int i, n;
  uint64_t start_mil, end_mil, start_mil_aff, end_mil_aff, start_mil_pc, end_mil_pc, start_mil_usepc, end_mil_usepc, start_fe, end_fe;
 
  if (PAIR_CURVE == BN12) {
    Fp2_t B, B3, l00, l01, l10;
    Point_wp2_t P, Q, R, Paff, Qaff, Raff;
    Fp12_t f1, f2, f3;
    Fp_t k1;
	Fp_t alph, bet, gam;
	uint64_t *premem;
   
    printf ("Running BN12 pairing...\n");
    n = 4;
    m = (uint64_t *) malloc (n * sizeof (uint64_t));
    b0 = (uint64_t *) malloc (n * sizeof (uint64_t));
    b1 = (uint64_t *) malloc (n * sizeof (uint64_t));
    m[3] = 0x2523648240000001;
    m[2] = 0xBA344D8000000008;
    m[1] = 0x6121000000000013;
    m[0] = 0xA700000000000013; // 2523648240000001BA344D80000000086121000000000013A700000000000013
    /* This is the 254-bit BN prime p = 3 (mod 4). */

    Fp_initialize_config (m, n);
    Fp_init (b);

	Fp_init(alph);
	Fp_init(bet);
	Fp_init(gam);
    Fp_wp_init (P0);
    Fp_wp_init (Q0);
	Fp_wp_init(P0aff);
	Fp_wp_init(Q0aff);
    
    Fp2_initialize_config (); 
    Fp4_initialize_config ();
    Fp6_initialize_config ();
    Fp12_initialize_config ();

    Fp2_init (B);
    Fp2_init (B3);
    Fp2_init (l00);
    Fp2_init (l01);
    Fp2_init (l10);
    Fp12_init (f1);
    Fp12_init (f2);
    Fp12_init (f3);
    Fp12_set_zero (f1);
    Fp12_set_zero (f2);
    Fp12_set_zero (f3);

    Fp_set_ui (b, 2);

    //m[3] = 0x10E8BA3F60E7FB98;
    //m[2] = 0x6D6589B8879FBF64;
    //m[1] = 0xA772FFFCFAC39939;
    //m[0] = 0xF64B7BE9639BD57D;
	m[3] = 0x2523648240000001;
	m[2] = 0xBA344D8000000008;
	m[1] = 0x6121000000000013;
	m[0] = 0xA700000000000012;
	Fp_set (P0aff->X, m);

    //m[3] = 0x0C5D16EEAFF21A84;
    //m[2] = 0x5929E8BE8D0117F8;
    //m[1] = 0xB951CFA243F6850B;
    //m[0] = 0x8FD5BC0A844AFFA5;
	m[3] = 0x0000000000000000;
	m[2] = 0x0000000000000000;
	m[1] = 0x0000000000000000;
	m[0] = 0x0000000000000001;
	Fp_set(P0aff->Y, m);

	Fp_weierstrass_aff2proj(P0, P0aff);

    b1[3] = 0x2523648240000001;
    b1[2] = 0xBA344D8000000008;
    b1[1] = 0x6121000000000013;
    b1[0] = 0xA700000000000012;
  
    b0[3] = 0x0000000000000000;
    b0[2] = 0x0000000000000000;
    b0[1] = 0x0000000000000000;
    b0[0] = 0x0000000000000001;
    Fp2_set (B, b0, b1);

    Fp_initialize_weierstrass (&b, 1);
    Fp2_initialize_weierstrass (B);
    Fp2_wp_init (P);
    Fp2_wp_init (Q);
    Fp2_wp_init (R);
	Fp2_wp_init(Paff);
	Fp2_wp_init(Qaff);
	Fp2_wp_init(Raff);
    bn_bls12_miller_initialize_config (B);
    bn_bls12_finalexpo_initialize_config ();
  
    Fp2_copy (B3, bn_config.bt3);

       /*b1[3] = 0x1C9EF62920C85936;
    b1[2] = 0x371AED1A1EB30C83;
    b1[1] = 0x006A310CBA6D1E26;
    b1[0] = 0xE429E3CBDB314928;*/
	b1[3] = 0x0516AAF9BA737833;
	b1[2] = 0x310AA78C5982AA5B;
	b1[1] = 0x1F4D746BAE3784B7;
	b1[0] = 0x0D8C34C1E7D54CF3;
	
    /*b0[3] = 0x1ABE36D06DD07358;
    b0[2] = 0xB7CBB532CF6961DF;
    b0[1] = 0x743BD07961A287A3;
    b0[0] = 0x59B0D4C307B5C914;*/
	b0[3] = 0x061A10BB519EB62F;
	b0[2] = 0xEB8D8C7E8C61EDB6;
	b0[1] = 0xA4648BBB4898BF0D;
	b0[0] = 0x91EE4224C803FB2B;
	Fp2_set (Paff->X, b0, b1);
  
   /* b1[3] = 0x085B1A1EF635AD65;
    b1[2] = 0x29F489BAF930310E;
    b1[1] = 0x70B25864B6178492;
    b1[0] = 0x5E66EEAB664F1C4B;*/
	b1[3] = 0x0EBB2B0E7C8B1526;
	b1[2] = 0x8F6D4456F5F38D37;
	b1[1] = 0xB09006FFD739C957;
	b1[0] = 0x8A2D1AEC6B3ACE9B;

    /*b0[3] = 0x079A4686E85D1054;
    b0[2] = 0xFEEAE258FB2103C8;
    b0[1] = 0xD6A34D4C2FDD968F;
    b0[0] = 0x6526DE28C4CC4EF3;*/
	b0[3] = 0x021897A06BAF9343;
	b0[2] = 0x9A90E096698C8223;
	b0[1] = 0x29BD0AE6BDBE09BD;
	b0[0] = 0x19F0E07891CD2B9A;
	Fp2_set (Paff->Y, b0, b1);

	Fp2_weierstrass_aff2proj(P, Paff);

    printf ("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff (P0));
    printf ("Affine point P on curve?: %d\n", Fp2_weierstrass_oncurve_aff (P));
    Fp2_weierstrass_aff2proj (P, P);
    printf ("Projective point P on curve?: %d\n", Fp2_weierstrass_oncurve_jacproj (P));
    Fp2_copy (R->X, P->X);
    Fp2_copy (R->Y, P->Y);
    Fp2_copy (R->Z, P->Z);
    
    printf("\nNAF(6u+2) = ");
    for (i=0; i<bn_config.len6u2; i++) {
      printf ("%d", bn_config.naf6u2[i]);
    }
    printf("\n");

	printf("\n(p+1)/4 = ");
	for (i = 0; i<64*Fp_config.m->n; i++) {
		printf("%d", BIT(Fp_config.sqrt_exp_p3mod4,i));
	}
	printf("\n");
	printf("%d bits\n", Fp_config.sqrt_exp_len_p3mod4);

	for (i = 0; i < 1000; i++) {
		Fp_rand(alph);
		//Fp_set_ui(alph, (uint64_t) i);
		//Fp_neg(alph, alph);
		//printf("\nalpha = "); Fp_print(alph);
		Fp_sqrt(bet, alph);
		//printf("\nbeta = "); Fp_print(bet);
		Fp_mul(gam, bet, bet);
		//printf("\nbeta^2 = "); Fp_print(gam);
		if (Fp_cmp(alph, gam) != 0) {
			Fp_neg(gam, gam);
			if (Fp_cmp(alph, gam) != 0) {
				printf("Sqrt error!");
			}
		}
	}
	

    Fp_weierstrass_aff2proj (P0, P0);
    Fp_weierstrass_dbl (Q0, P0);
    Fp_weierstrass_add (Q0, Q0, P0);
    Fp_weierstrass_dbl (Q0, Q0);
    Fp_weierstrass_dbl (Q0, Q0);
    Fp_weierstrass_add (Q0, Q0, P0);
    Fp_weierstrass_dbl (Q0, Q0);
	Fp_weierstrass_add (Q0, Q0, P0);
	Fp_weierstrass_dbl (Q0, Q0);
    Fp_weierstrass_add (Q0, Q0, P0);
    if (Fp_weierstrass_oncurve_jacproj (Q0) != 1) {
      printf ("G1 error\n");
      getchar ();
    }
	
    
    Fp2_weierstrass_dbl (Q, P);
    Fp2_weierstrass_add (Q, Q, P);
    Fp2_weierstrass_dbl (Q, Q);
    Fp2_weierstrass_dbl (Q, Q);
    Fp2_weierstrass_add (Q, Q, P);
    Fp2_weierstrass_dbl (Q, Q);
	Fp2_weierstrass_add (Q, Q, P);
	Fp2_weierstrass_dbl (Q, Q);
    Fp2_weierstrass_add (Q, Q, P);
    Fp2_wp_neg (Q, Q);
    if (Fp2_weierstrass_oncurve_jacproj (Q) != 1) {
	  printf ("G2 error\n");
      getchar ();
    }
	Fp_weierstrass_jacproj2aff(Q0, Q0);
    Fp2_weierstrass_jacproj2aff (Q, Q);

    bn_optimal_ate_miller_aff (f1, Q, P0);
    bn_optimal_ate_miller_aff (f2, P, Q0);
	bn_finalexpo_neg (f1, f1);
    bn_finalexpo_neg (f2, f2);
    Fp12_inv(f2, f2);
    printf ("projective pairing: e([k]Q, P) = e(Q, [k]P)? %d\n", Fp12_cmpeq (f1, f2));
    if (Fp12_cmpeq (f1, f2) != 1) {
      printf ("BN pairing error 1\n");
      getchar ();
    }

	Fp_weierstrass_affdbl(Q0aff, P0aff);
	Fp_weierstrass_affdbl(Q0aff, P0aff);
	Fp_weierstrass_affdbl(Q0aff, P0aff);
	Fp_weierstrass_affdbl(Q0aff, P0aff);
	Fp_weierstrass_affadd(Q0aff, Q0aff, P0aff);
	Fp_weierstrass_affdbl(Q0aff, Q0aff);
	Fp_weierstrass_affdbl(Q0aff, Q0aff);
	Fp_weierstrass_affadd(Q0aff, Q0aff, P0aff);
	Fp_weierstrass_affdbl(Q0aff, Q0aff);
	Fp_weierstrass_affadd(Q0aff, Q0aff, P0aff);
	Fp_weierstrass_affdbl(Q0aff, Q0aff);
	if (Fp_weierstrass_oncurve_aff(Q0aff) != 1) {
		printf("G1 affine error\n");
		getchar();
	}

	Fp2_weierstrass_affdbl(Qaff, Paff);
	Fp2_weierstrass_affdbl(Qaff, Paff);
	Fp2_weierstrass_affdbl(Qaff, Paff);
	Fp2_weierstrass_affdbl(Qaff, Paff);
	Fp2_weierstrass_affadd(Qaff, Qaff, Paff);
	Fp2_weierstrass_affdbl(Qaff, Qaff);
	Fp2_weierstrass_affdbl(Qaff, Qaff);
	Fp2_weierstrass_affadd(Qaff, Qaff, Paff);
	Fp2_weierstrass_affdbl(Qaff, Qaff);
	Fp2_weierstrass_affadd(Qaff, Qaff, Paff);
	Fp2_weierstrass_affdbl(Qaff, Qaff);
	Fp2_wp_neg(Qaff, Qaff);
	if (Fp2_weierstrass_oncurve_aff(Qaff) != 1) {
		printf("G2 affine error\n");
		getchar();
	}
	
	bn_optimal_ate_miller_aff(f1, Qaff, P0aff);
	bn_optimal_ate_miller_aff(f2, Paff, Q0aff);
	bn_finalexpo_neg(f1, f1);
	bn_finalexpo_neg(f2, f2);
	Fp12_inv(f2, f2);
	printf("affine pairing: e([k]Q, P) = e(Q, [k]P)? %d\n", Fp12_cmpeq(f1, f2));
	if (Fp12_cmpeq(f1, f2) != 1) {
		printf("BN affine pairing error 1\n");
		getchar();
	}

	premem = (uint64_t *)malloc(4 * Fp_config.m->n * (bn_config.len6u2 + bn_config.nafwt6u2) * sizeof (uint64_t));
	if (premem == NULL) return ERR_OUT_OF_MEMORY;

	bn_optimal_ate_miller_precompute_aff(Qaff, premem);
	bn_optimal_ate_miller_useprecomputed_aff(f1, P0aff, premem);

	bn_optimal_ate_miller_precompute_aff(Paff, premem);
	bn_optimal_ate_miller_useprecomputed_aff(f2, Q0aff, premem);

	bn_finalexpo_neg(f1, f1);
	bn_finalexpo_neg(f2, f2);
	Fp12_inv(f2, f2);
	printf("affine pairing: e([k]Q, P) = e(Q, [k]P)? %d\n", Fp12_cmpeq(f1, f2));
	if (Fp12_cmpeq(f1, f2) != 1) {
		printf("BN affine precomputed pairing error 1\n");
		getchar();
	}


    start_mil = __rdtsc();
    for (i=0; i < 1000; i++) {
      bn_optimal_ate_miller (f1, Q, P0);
    }
    end_mil = __rdtsc();
    start_fe = __rdtsc();
    for (i=0; i < 1000; i++) {
      bn_finalexpo_neg (f3, f1);
    }
    end_fe = __rdtsc();
	start_mil_aff = __rdtsc();
	for (i = 0; i < 1000; i++) {
		bn_optimal_ate_miller_aff(f1, Q, P0);
	}
	end_mil_aff = __rdtsc();
	start_mil_pc = __rdtsc();
	for (i = 0; i < 1000; i++) {
		bn_optimal_ate_miller_precompute_aff(Q, premem);
	}
	end_mil_pc = __rdtsc();
	start_mil_usepc = __rdtsc();
	for (i = 0; i < 1000; i++) {
		bn_optimal_ate_miller_useprecomputed_aff(f1, P0, premem);
	}
	end_mil_usepc = __rdtsc();
    
    bn_optimal_ate_miller (f2, P, Q0);
    bn_finalexpo_neg (f2, f2);
	//printf("f2 = "); Fp12_print(f2);
	//printf("f3 = "); Fp12_print(f3);
	/*if (Fp12_cmpeq (f3, f2) != 1) {
      printf ("BN pairing error 2\n");
      getchar ();
    }//*/

    printf ("Time BN12 OPTIMAL_ATE:  %fM cycles\n", ((uint32_t) (((double) end_mil-start_mil) / 1000.0)) / 1000000.);
    printf ("Time BN12 FINAL_EXPO:  %fM cycles\n", ((uint32_t) (((double) end_fe-start_fe) / 1000.0)) / 1000000.);
    printf ("Time BN12 PAIRING:  %fM cycles\n", ((uint32_t) (((double) end_fe-start_mil) / 1000.0)) / 1000000.);
	printf ("Time BN12 OPTIMAL_ATE (affine):  %fM cycles\n", ((uint32_t)(((double)end_mil_aff - start_mil_aff) / 1000.0)) / 1000000.);
	printf ("Time BN12 OPTIMAL_ATE (precompute):  %fM cycles\n", ((uint32_t)(((double)end_mil_pc - start_mil_pc) / 1000.0)) / 1000000.);
	printf ("Time BN12 OPTIMAL_ATE (use precomputation):  %fM cycles\n", ((uint32_t)(((double)end_mil_usepc - start_mil_usepc) / 1000.0)) / 1000000.);


	Fp_free (b);
    Fp_free (k1);
    Fp2_free (B);
    Fp2_free (B3);
    Fp2_free (l00);
    Fp2_free (l01);
    Fp2_free (l10);
    Fp12_free (f1);
    Fp12_free (f2);
    Fp12_free (f3);
    Fp2_wp_free (P);
    Fp2_wp_free (Q);
	Fp2_wp_free(R);
	Fp2_wp_free(Qaff);
	Fp2_wp_free(Paff);
	Fp2_wp_free(Raff);
	Fp_wp_free(P0aff);
	Fp_wp_free(Q0aff);
    Fp_wp_free (P0);
    Fp_wp_free (Q0);
    bn_bls12_finalexpo_free_config ();
    bn_bls12_miller_free_config ();
    Fp12_free_config (); 
    Fp6_free_config ();
    Fp4_free_config ();
    Fp2_free_weierstrass ();
    Fp2_free_config (); 
    Fp_free_config ();
  } 
  else if (PAIR_CURVE == BN12CP) {
    Fp2_t B, B3, l00, l01, l10;
    Point_wp2_t P, Q, R;
    Fp12_t f1, f2, f3;
   
    printf ("Running BN12CP pairing...\n");
    n = 4;
    m = (uint64_t *) malloc (n * sizeof (uint64_t));
    b0 = (uint64_t *) malloc (n * sizeof (uint64_t));
    b1 = (uint64_t *) malloc (n * sizeof (uint64_t));
    m[3] = 0x2400005100016456;
    m[2] = 0x51FF9E2B74EEEDF7;
    m[1] = 0x7654B7908E8B90CA;
    m[0] = 0xD2827528FFD90013; // 0x240000510001645651FF9E2B74EEEDF77654B7908E8B90CAD2827528FFD90013
    /* This is the 254-bit BN prime p = 3 (mod 4) such that the corresponding BN group order
       divides the Cocks-Pinch curve CP6. */

    Fp_initialize_config (m, n);
    Fp_init (b);
    Fp_wp_init (P0);
    Fp_wp_init (Q0);
    
    Fp2_initialize_config (); 
    Fp4_initialize_config ();
    Fp6_initialize_config ();
    Fp12_initialize_config ();

    Fp2_init (B);
    Fp2_init (B3);
    Fp2_init (l00);
    Fp2_init (l01);
    Fp2_init (l10);
    Fp12_init (f1);
    Fp12_init (f2);
    Fp12_init (f3);
    Fp12_set_zero (f1);
    Fp12_set_zero (f2);
    Fp12_set_zero (f3);
    
    Fp_set_ui (b, 2);

    m[3] = 0x190E728DCF438BE6;
    m[2] = 0x940064DBA531727C;
    m[1] = 0x9E6997443354EAA9;
    m[0] = 0x61B8915866546C1B;
    Fp_set (P0->X, m);

    m[3] = 0x12C0780D35162927;
    m[2] = 0x04F93C81D2583F24;
    m[1] = 0x0FDBABF86B664960;
    m[0] = 0x4488E11F5290673F;
    Fp_set (P0->Y, m);

    b1[3] = 0x2400005100016456;
    b1[2] = 0x51FF9E2B74EEEDF7;
    b1[1] = 0x7654B7908E8B90CA;
    b1[0] = 0xD2827528FFD90012;
  
    b0[3] = 0x0000000000000000;
    b0[2] = 0x0000000000000000;
    b0[1] = 0x0000000000000000;
    b0[0] = 0x0000000000000001;
    Fp2_set (B, b0, b1);

    Fp_initialize_weierstrass (&b, 1);
    Fp2_initialize_weierstrass (B);
    Fp2_wp_init (P);
    Fp2_wp_init (Q);
    Fp2_wp_init (R);
    bn_bls12_miller_initialize_config (B);
    bn_bls12_finalexpo_initialize_config ();
  
    Fp2_copy (B3, bn_config.bt3);

    b1[3] = 0x22878CF192B8B65E;
    b1[2] = 0x217BE0F00A60EBF2;
    b1[1] = 0x4536F25DA42FFA2B;
    b1[0] = 0xED2C8A3664B14F22;
  
    b0[3] = 0x1EFE5990A27712E3;
    b0[2] = 0xCA2E2F1F0BBA3F18;
    b0[1] = 0x5CB9B80B5402B2B3;
    b0[0] = 0x8E07986D0325775A;
    Fp2_set (P->X, b0, b1);
  
    b1[3] = 0x17AA1849A47EB48F;
    b1[2] = 0x3A6A0BEEE7810263;
    b1[1] = 0x01727155287C9BCB;
    b1[0] = 0x8A0FA76176CB3646;
  
    b0[3] = 0x175657516621CD68;
    b0[2] = 0x867C0F5ED3C36D94;
    b0[1] = 0x5E2956368412A4C8;
    b0[0] = 0x02E2C7671BAD1C73;
    Fp2_set (P->Y, b0, b1);

    printf ("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff (P0));
    printf ("Affine point P on curve?: %d\n", Fp2_weierstrass_oncurve_aff (P));
    Fp2_weierstrass_aff2proj (P, P);
    printf ("Projective point P on curve?: %d\n", Fp2_weierstrass_oncurve_jacproj (P));
    Fp2_copy (R->X, P->X);
    Fp2_copy (R->Y, P->Y);
    Fp2_copy (R->Z, P->Z);

    printf("\nNAF(6u+2) = ");
    for (i=0; i<bn_config.len6u2; i++) {
      printf ("%d", bn_config.naf6u2[i]);
    }
    printf("\n");
    
    Fp_weierstrass_aff2proj (P0, P0);
    Fp_weierstrass_dbl (Q0, P0);
    Fp_weierstrass_add (Q0, Q0, P0);
    Fp_weierstrass_dbl (Q0, Q0);
    Fp_weierstrass_dbl (Q0, Q0);
    Fp_weierstrass_add (Q0, Q0, P0);
    Fp_weierstrass_dbl (Q0, Q0);
    Fp_weierstrass_add (Q0, Q0, P0);
    if (Fp_weierstrass_oncurve_jacproj (Q0) != 1) {
      printf ("G1 error\n");
      getchar ();
    }
    Fp_weierstrass_jacproj2aff (Q0, Q0);
    Fp2_weierstrass_dbl (Q, P);
    Fp2_weierstrass_add (Q, Q, P);
    Fp2_weierstrass_dbl (Q, Q);
    Fp2_weierstrass_dbl (Q, Q);
    Fp2_weierstrass_add (Q, Q, P);
    Fp2_weierstrass_dbl (Q, Q);
    Fp2_weierstrass_add (Q, Q, P);
    Fp2_wp_neg (Q, Q);
    if (Fp2_weierstrass_oncurve_jacproj (Q) != 1) {
      printf ("G2 error\n");
      getchar ();
    }
    Fp2_weierstrass_jacproj2aff (Q, Q);

    bn_optimal_ate_miller (f1, Q, P0);
    bn_optimal_ate_miller (f2, P, Q0);
    bn_finalexpo_pos (f1, f1);
    bn_finalexpo_pos (f2, f2);
    Fp12_inv(f2, f2);
    printf ("e([k]Q, P) = e(Q, [k]P)? %d\n", Fp12_cmpeq (f1, f2));
    if (Fp12_cmpeq (f1, f2) != 1) {
      printf ("BN pairing error 1\n");
      getchar ();
    }
    
    start_mil = __rdtsc();
    for (i=0; i < 1000; i++) {
      bn_optimal_ate_miller (f1, Q, P0);
    }
    end_mil = __rdtsc();
    start_fe = __rdtsc();
    for (i=0; i < 1000; i++) {
      bn_finalexpo_pos (f3, f1);
    }
    end_fe = __rdtsc();
    
    
    bn_optimal_ate_miller (f2, P, Q0);
    bn_finalexpo_pos (f2, f2);
    /*if (Fp12_cmpeq (f3, f2) != 1) {
      printf ("BN pairing error 2\n");
      getchar ();
    }//*/

    printf ("Time BN12CP OPTIMAL_ATE:  %fM cycles\n", ((uint32_t) (((double) end_mil-start_mil) / 1000.0)) / 1000000.);
    printf ("Time BN12CP FINAL_EXPO:  %fM cycles\n", ((uint32_t) (((double) end_fe-start_fe) / 1000.0)) / 1000000.);
    printf ("Time BN12CP PAIRING:  %fM cycles\n", ((uint32_t) (((double) end_fe-start_mil) / 1000.0)) / 1000000.);

    Fp_free (b);
    Fp2_free (B);
    Fp2_free (B3);
    Fp2_free (l00);
    Fp2_free (l01);
    Fp2_free (l10);
    Fp12_free (f1);
    Fp12_free (f2);
    Fp12_free (f3);
    Fp2_wp_free (P);
    Fp2_wp_free (Q);
    Fp_wp_free (P0);
    Fp_wp_free (Q0);
    bn_bls12_finalexpo_free_config ();
    bn_bls12_miller_free_config ();
    Fp12_free_config (); 
    Fp6_free_config ();
    Fp4_free_config ();
    Fp2_free_weierstrass ();
    Fp2_free_config (); 
    Fp_free_config ();
  } 
  else if (PAIR_CURVE == BN12tiny) {
	  Fp2_t B, B3, l00, l01, l10;
	  Point_wp2_t P, Q, R;
	  Fp12_t f1, f2, f3;
	  Fp_t alph, bet, gam;

	  printf("Running BN12tiny pairing...\n");
	  n = 1;
	  m = (uint64_t *)malloc(n * sizeof (uint64_t));
	  b0 = (uint64_t *)malloc(n * sizeof (uint64_t));
	  b1 = (uint64_t *)malloc(n * sizeof (uint64_t));
	  m[0] = 0x15C9B0796B71DB;
	  /* This is the 53-bit BN prime p = 3 (mod 4). */

	  Fp_initialize_config(m, n);
	  Fp_init(b);

	  Fp_init(alph);
	  Fp_init(bet);
	  Fp_init(gam);
	  Fp_wp_init(P0);
	  Fp_wp_init(Q0);

	  Fp2_initialize_config();
	  Fp4_initialize_config();
	  Fp6_initialize_config();
	  Fp12_initialize_config();

	  Fp2_init(B);
	  Fp2_init(B3);
	  Fp2_init(l00);
	  Fp2_init(l01);
	  Fp2_init(l10);
	  Fp12_init(f1);
	  Fp12_init(f2);
	  Fp12_init(f3);
	  Fp12_set_zero(f1);
	  Fp12_set_zero(f2);
	  Fp12_set_zero(f3);

	  Fp_set_ui(b, 2);

	  m[0] = 0x321347FC84CEE;
	  Fp_set(P0->X, m);

	  m[0] = 0xC5D3F93EEC09B;
	  Fp_set(P0->Y, m);

	  b1[0] = 0x15C9B0796B71DA;

	  b0[0] = 0x0000000000000001;
	  Fp2_set(B, b0, b1);

	  Fp_initialize_weierstrass(&b, 1);
	  Fp2_initialize_weierstrass(B);
	  Fp2_wp_init(P);
	  Fp2_wp_init(Q);
	  Fp2_wp_init(R);
	  bn_bls12_miller_initialize_config(B);
	  bn_bls12_finalexpo_initialize_config();

	  Fp2_copy(B3, bn_config.bt3);

	  b1[0] = 0x11B4D9B9078C9C;

	  b0[0] = 0x1183E42BE63029;
	  Fp2_set(P->X, b0, b1);

	  b1[0] = 0x44923E174AB60;

	  b0[0] = 0x10B2AC7252D760;
	  Fp2_set(P->Y, b0, b1);

	  printf("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff(P0));
	  printf("Affine point P on curve?: %d\n", Fp2_weierstrass_oncurve_aff(P));
	  Fp2_weierstrass_aff2proj(P, P);
	  printf("Projective point P on curve?: %d\n", Fp2_weierstrass_oncurve_jacproj(P));
	  Fp2_copy(R->X, P->X);
	  Fp2_copy(R->Y, P->Y);
	  Fp2_copy(R->Z, P->Z);

	  printf("\nNAF(6u+2) = ");
	  for (i = 0; i<bn_config.len6u2; i++) {
		  printf("%d", bn_config.naf6u2[i]);
	  }
	  printf("\n");

	  printf("\n(p+1)/4 = ");
	  for (i = 0; i<64 * Fp_config.m->n; i++) {
		  printf("%d", BIT(Fp_config.sqrt_exp_p3mod4, i));
	  }
	  printf("\n");
	  printf("%d bits\n", Fp_config.sqrt_exp_len_p3mod4);

	  for (i = 0; i < 1000; i++) {
		  Fp_rand(alph);
		  //Fp_set_ui(alph, (uint64_t) i);
		  //Fp_neg(alph, alph);
		  //printf("\nalpha = "); Fp_print(alph);
		  Fp_sqrt(bet, alph);
		  //printf("\nbeta = "); Fp_print(bet);
		  Fp_mul(gam, bet, bet);
		  //printf("\nbeta^2 = "); Fp_print(gam);
		  if (Fp_cmp(alph, gam) != 0) {
			  Fp_neg(gam, gam);
			  if (Fp_cmp(alph, gam) != 0) {
				  printf("Sqrt error!");
			  }
		  }
	  }


	  Fp_weierstrass_aff2proj(P0, P0);
	  Fp_weierstrass_dbl(Q0, P0);
	  Fp_weierstrass_add(Q0, Q0, P0);
	  Fp_weierstrass_dbl(Q0, Q0);
	  Fp_weierstrass_dbl(Q0, Q0);
	  Fp_weierstrass_add(Q0, Q0, P0);
	  Fp_weierstrass_dbl(Q0, Q0);
	  Fp_weierstrass_add(Q0, Q0, P0);
	  if (Fp_weierstrass_oncurve_jacproj(Q0) != 1) {
		  printf("G1 error\n");
		  getchar();
	  }
	  Fp_weierstrass_jacproj2aff(Q0, Q0);
	  Fp2_weierstrass_dbl(Q, P);
	  Fp2_weierstrass_add(Q, Q, P);
	  Fp2_weierstrass_dbl(Q, Q);
	  Fp2_weierstrass_dbl(Q, Q);
	  Fp2_weierstrass_add(Q, Q, P);
	  Fp2_weierstrass_dbl(Q, Q);
	  Fp2_weierstrass_add(Q, Q, P);
	  Fp2_wp_neg(Q, Q);
	  if (Fp2_weierstrass_oncurve_jacproj(Q) != 1) {
		  printf("G2 error\n");
		  getchar();
	  }
	  Fp2_weierstrass_jacproj2aff(Q, Q);

	  bn_optimal_ate_miller(f1, Q, P0);
	  bn_optimal_ate_miller(f2, P, Q0);
	  bn_finalexpo_neg(f1, f1);
	  bn_finalexpo_neg(f2, f2);
	  Fp12_inv(f2, f2);
	  printf("e([k]Q, P) = e(Q, [k]P)? %d\n", Fp12_cmpeq(f1, f2));
	  if (Fp12_cmpeq(f1, f2) != 1) {
		  printf("BN pairing error 1\n");
		  getchar();
	  }

	  start_mil = __rdtsc();
	  for (i = 0; i < 1000; i++) {
		  bn_optimal_ate_miller(f1, Q, P0);
	  }
	  end_mil = __rdtsc();
	  start_fe = __rdtsc();
	  for (i = 0; i < 1000; i++) {
		  bn_finalexpo_neg(f3, f1);
	  }
	  end_fe = __rdtsc();

	  bn_optimal_ate_miller(f2, P, Q0);
	  bn_finalexpo_neg(f2, f2);
	  /*if (Fp12_cmpeq (f3, f2) != 1) {
	  printf ("BN pairing error 2\n");
	  getchar ();
	  }//*/

	  printf("Time BN12tiny OPTIMAL_ATE:  %fM cycles\n", ((uint32_t)(((double)end_mil - start_mil) / 1000.0)) / 1000000.);
	  printf("Time BN12tiny FINAL_EXPO:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_fe) / 1000.0)) / 1000000.);
	  printf("Time BN12tiny PAIRING:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_mil) / 1000.0)) / 1000000.);

	  Fp_free(b);
	  Fp2_free(B);
	  Fp2_free(B3);
	  Fp2_free(l00);
	  Fp2_free(l01);
	  Fp2_free(l10);
	  Fp12_free(f1);
	  Fp12_free(f2);
	  Fp12_free(f3);
	  Fp2_wp_free(P);
	  Fp2_wp_free(Q);
	  Fp_wp_free(P0);
	  Fp_wp_free(Q0);
	  bn_bls12_finalexpo_free_config();
	  bn_bls12_miller_free_config();
	  Fp12_free_config();
	  Fp6_free_config();
	  Fp4_free_config();
	  Fp2_free_weierstrass();
	  Fp2_free_config();
	  Fp_free_config();
  }
  else if (PAIR_CURVE == CP6b) {
    Fp_t B[2], B3;
    Point_wp_t P, Q, R;
    Fp6q_t f1, f2, f3;
	Fp_t alph, bet, gam;
	uint64_t *premem;
    
    uint64_t pCP6b[8] = { 0xF9000000000005BF,0xAA8EC00000000BD3,
                          0x0CF85A8000000B1F,0x8761FD7D40000638,
                          0x0A67549B77800241,0xB007E82F1F3CB08B,
                          0x41358FBC1F137315,0x158CFFC9264E1116 };
    /* This is the 509-bit CP6b prime p = 3 (mod 4). */

    uint64_t Px[8] = { 0x549400D7139BD8C3,0xF63D76CC9E5074C4,
                       0x9ECA293D846AD7F4,0xBCFDE6DF8A0528AC,
                       0x7412962AED6CB7E3,0x8CA78CFBE18D6C03,
                       0x9D41ADD3F736DF40,0x103E84948482F3AC };
    uint64_t Py[8] = { 0xA49BC75E21A206DA,0x7EC88B25BB3E39DF,
                       0xC54F80610A16417D,0xF7CB2BB06F683525,
                       0xB76C9F7CC52E9739,0xBAB17A19BBB5F211,
                       0x176CE782E1557D4E,0x0521E32E8D0B8F03 };
    /* This is the point over Fp in G1. */
   
    uint64_t Qx[8] = { 0x4179C530BAF016F7,0x7D29FB703C56B725,
                       0x8D8AA7AEAEBB06D3,0xA6CA4B97C79D31C1,
                       0xA49308E0F3BB3294,0xB1ECCE42F22C7800,
                       0x9B245837BEC9D3B5,0x113E15D42AB7CF02 };
    uint64_t Qy[8] = { 0x0E83D7FB80304008,0x1A0B1C4E85A3FC34,
                       0x731C345CFA535AB3,0x507571381258BBEF,
                       0x3B427423987C4B28,0x2535C8DC6F89FBF2,
                       0x0E34A6FB987994A0,0x00FCE2C1051EE0D1 };

    n = 8;
    printf ("Running CP6b pairing...\n");

    Fp_initialize_config (pCP6b, n);
    Fp3_initialize_config ();
    Fp6q_initialize_config ();

    Fp_init (B[0]);
    Fp_init (B[1]);
    Fp_init (B3);
	Fp_init(alph);
	Fp_init(bet);
	Fp_init(gam);
    Fp_wp_init (P0);
    Fp_wp_init (Q0);

    Fp_set_ui (B[0], 21);
    Fp_set_ui (B[1], 63);

    Fp_initialize_weierstrass (B, 2);
    Fp_wp_init (P);
    Fp_wp_init (Q);
    Fp_wp_init (R);

    Fp_set (P0->X, Px);
    Fp_set (P0->Y, Py);

    Fp_set (P->X, Qx);
    Fp_set (P->Y, Qy);

    printf ("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b (P0, 0)); 
    printf ("Affine point P on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b (P, 1));
    Fp_weierstrass_aff2proj (P, P);
    printf ("Projective point P on curve?: %d\n", Fp_weierstrass_oncurve_jacproj_select_b (P, 1));
    Fp_copy (R->X, P->X);
    Fp_copy (R->Y, P->Y);
    Fp_copy (R->Z, P->Z);

    
    cp6_miller_initialize_config (B[1]);
    cp6_finalexpo_initialize_config ();

    printf("\nNAF(m) = ");
    for (i=0; i<cp6_config.lenm; i++) {
      printf ("%d", cp6_config.nafm[i]);
    }
    printf("\n");
    printf("\nNAF(3m) = ");
    for (i=0; i<cp6_config.len3m; i++) {
      printf ("%d", cp6_config.naf3m[i]);
    }
    printf("\n");

	printf("\n(p+1)/4 = ");
	for (i = 0; i<64 * Fp_config.m->n; i++) {
		printf("%d", BIT(Fp_config.sqrt_exp_p3mod4, i));
	}
	printf("\n");
	printf("%d bits\n", Fp_config.sqrt_exp_len_p3mod4);

	for (i = 0; i < 1000; i++) {
		Fp_rand(alph);
		//Fp_set_ui(alph, (uint64_t) i);
		//Fp_neg(alph, alph);
		//printf("\nalpha = "); Fp_print(alph);
		Fp_sqrt(bet, alph);
		//printf("\nbeta = "); Fp_print(bet);
		Fp_mul(gam, bet, bet);
		//printf("\nbeta^2 = "); Fp_print(gam);
		if (Fp_cmp(alph, gam) != 0) {
			Fp_neg(gam, gam);
			if (Fp_cmp(alph, gam) != 0) {
				printf("Sqrt error!");
			}
		}
	}

    for (i=0; i<1; i++) {
      Fp_weierstrass_aff2proj (P0, P0);
      Fp_weierstrass_dbl (Q0, P0);
      Fp_weierstrass_add (Q0, Q0, P0);
      Fp_weierstrass_dbl (Q0, Q0);
      Fp_weierstrass_dbl (Q0, Q0);
      Fp_weierstrass_add (Q0, Q0, P0);
      if (Fp_weierstrass_oncurve_jacproj_select_b (Q0, 0) != 1) {
        printf ("G1 error\n");
        getchar ();
      }
      Fp_weierstrass_jacproj2aff (Q0, Q0);
      Fp_weierstrass_dbl (Q, P);
      Fp_weierstrass_add (Q, Q, P);
      Fp_weierstrass_dbl (Q, Q);
      Fp_weierstrass_dbl (Q, Q);
      Fp_weierstrass_add (Q, Q, P);
      //Fp_wp_neg (Q, Q);
      if (Fp_weierstrass_oncurve_jacproj_select_b (Q, 1) != 1) {
        printf ("G2 error\n");
        getchar ();
      }
      Fp_weierstrass_jacproj2aff (Q, Q);
    }

    Fp6q_init (f1);
    Fp6q_init (f2);
    Fp6q_init (f3);
    Fp6q_set_zero (f1);
    Fp6q_set_zero (f2);
    Fp6q_set_zero (f3);

    //printf ("f1 ="); Fp6q_print (f1);
    //printf ("f2 ="); Fp6q_print (f2);
    
    cp6b_optimal_ate_miller (f1, Q, P0);
    cp6b_optimal_ate_miller (f2, P, Q0);
    //printf ("f1 ="); Fp6q_print (f1);
    //printf ("f2 ="); Fp6q_print (f2);

    cp6b_finalexpo (f1, f1);
    cp6b_finalexpo (f2, f2);
    //Fp6q_inv (f1,f1);

    printf ("e([k]Q, P) = e(Q, [k]P)? %d\n", Fp6q_cmpeq (f1, f2));
    if (Fp6q_cmpeq (f1, f2) != 1) {
      printf ("CP6 pairing error\n");
      getchar ();
    }
	//printf ("f1 ="); Fp6q_print (f1);
	//printf ("f2 ="); Fp6q_print (f2);

	cp6b_optimal_ate_miller_aff(f1, Q, P0);
	cp6b_optimal_ate_miller_aff(f2, P, Q0);
	//printf ("f1 ="); Fp6q_print (f1);
	//printf ("f2 ="); Fp6q_print (f2);

	cp6b_finalexpo(f1, f1);
	cp6b_finalexpo(f2, f2);
	//Fp6q_inv (f1,f1);

	printf("affine e([k]Q, P) = e(Q, [k]P)? %d\n", Fp6q_cmpeq(f1, f2));
	if (Fp6q_cmpeq(f1, f2) != 1) {
		printf("CP6b affine pairing error\n");
		getchar();
	}

	premem = (uint64_t *)malloc(3 * Fp_config.m->n * (cp6_config.lenm + cp6_config.len3m + cp6_config.nafwtm + cp6_config.nafwt3m) * sizeof (uint64_t));
	if (premem == NULL) return ERR_OUT_OF_MEMORY;

	cp6b_optimal_ate_miller_precompute(Q, premem);
	cp6b_optimal_ate_miller_useprecomputed(f1, P0, premem);

	cp6b_optimal_ate_miller_precompute(P, premem);
	cp6b_optimal_ate_miller_useprecomputed(f2, Q0, premem);
	//printf ("f1 ="); Fp6q_print (f1);
	//printf ("f2 ="); Fp6q_print (f2);

	cp6b_finalexpo(f1, f1);
	cp6b_finalexpo(f2, f2);
	//Fp6q_inv (f1,f1);

	printf("with precomputation: e([k]Q, P) = e(Q, [k]P)? %d\n", Fp6q_cmpeq(f1, f2));
	if (Fp6q_cmpeq(f1, f2) != 1) {
		printf("CP6b pairing error (with precomputation)\n");
		getchar();
	}

	start_mil = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6b_optimal_ate_miller(f1, Q, P0);
	}
	end_mil = __rdtsc();
	start_fe = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6b_finalexpo(f3, f1);
	}
	end_fe = __rdtsc();

	start_mil_pc = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6b_optimal_ate_miller_precompute(Q, premem);
	}
	end_mil_pc = __rdtsc();
	start_mil_usepc = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6b_optimal_ate_miller_useprecomputed(f1, P0, premem);
	}
	end_mil_usepc = __rdtsc();
	start_mil_aff = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6b_optimal_ate_miller_aff(f1, Q, P0);
	}
	end_mil_aff = __rdtsc();

	printf("Time CP6b OPTIMAL_ATE:  %fM cycles\n", ((uint32_t)(((double)end_mil - start_mil) / 1000.0)) / 1000000.);
	printf("Time CP6b FINAL_EXPO:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_fe) / 1000.0)) / 1000000.);
	printf("Time CP6b PAIRING:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_mil) / 1000.0)) / 1000000.);
	printf("Time CP6b OPTIMAL_ATE (affine):  %fM cycles\n", ((uint32_t)(((double)end_mil_aff - start_mil_aff) / 1000.0)) / 1000000.);
	printf("Time CP6b OPTIMAL_ATE (precompute):  %fM cycles\n", ((uint32_t)(((double)end_mil_pc - start_mil_pc) / 1000.0)) / 1000000.);
	printf("Time CP6b OPTIMAL_ATE (use precomputation):  %fM cycles\n", ((uint32_t)(((double)end_mil_usepc - start_mil_usepc) / 1000.0)) / 1000000.);

    //*/
    Fp_wp_free (P);
    Fp_wp_free (Q);
    Fp_wp_free (R);
    Fp_free_weierstrass ();
    Fp_wp_free (P0);
    Fp_wp_free (Q0);
    Fp_free (B[0]);
    Fp_free (B[1]);
    Fp_free (B3);

    cp6_finalexpo_free_config ();
    cp6_miller_free_config ();
    Fp6q_free_config ();
    Fp3_free_config ();
    Fp_free_config ();
  } else if (PAIR_CURVE == CP6) {
    Fp_t B[2], B3;
    Point_wp_t P, Q, R;
    Fp6q_t f1, f2, f3;
	Fp_t alph, bet, gam;
    
    uint64_t pCP6[8] = { 0x38B36CB83785822F,0xB45319A45547C003,
                         0xF7B9CAC8E6B9CD9C,0x43F0CC2E901010B1,
                         0x44BC84D632FB1F1B,0xE987D67EBDB08639,
                         0x3DD162F62EA4AA91,0x05100016C8007DD9 };
    /* This is the 507-bit CP6 prime p = 3 (mod 4). */

    uint64_t Px[8] = { 0x6EB1875827D6E976,0x6C6D1458F347DA36,
                       0x1F73DAF6AE18E938,0x3D8E78F285B83E43,
                       0x0D6A82277E4E58F4,0x713F078E69DEE165,
                       0x042E2F616DAEC635,0x03B308587093811C };
    uint64_t Py[8] = { 0xE38520961CAB7143,0x9A7DB2E37D898A88,
                       0xBD8940B6B7E05CE4,0x7AFDB9DFA3C4E33D,
                       0xEC06DD4BE9A06AA4,0xBCCE96F4B2BF3F1C,
                       0xDC310C00C6C19693,0x001E1A098CA4FE17 };
    /* This is the point over Fp in G1. */
   
    uint64_t Qx[8] = { 0x33981B1322FCECD5,0x8AB9E8C90804F0BC,
                       0xBFA6DE6369789A37,0xF7A7EE2E83C513B1,
                       0xE066A9431C0BB993,0xBD4AD1162A7B9029,
                       0x2F24B81F4986E22D,0x046CCB7CED12D7C7 };
    uint64_t Qy[8] = { 0xB86AAE30B8905FB8,0xECACCCA41E0574E9,
                       0x818BB90A7AFC04A3,0xCD7F8AD9081C8035,
                       0x3F74302CB786A01F,0x9A5CAD58CADC625E,
                       0xC5E06A7D5CCC612D,0x02EFC9B152AE226B };

    n = 8;
    printf ("Running CP6 pairing...\n");

    Fp_initialize_config (pCP6, n);
    Fp3_initialize_config ();
    Fp6q_initialize_config ();

    Fp_init (B[0]);
    Fp_init (B[1]);
    Fp_init (B3);
	Fp_init(alph);
	Fp_init(bet);
	Fp_init(gam);
    Fp_wp_init (P0);
    Fp_wp_init (Q0);

    Fp_set_ui (B[0], 2);
    Fp_set_ui (B[1], 6);

    Fp_initialize_weierstrass (B, 2);
    Fp_wp_init (P);
    Fp_wp_init (Q);
    Fp_wp_init (R);

    Fp_set (P0->X, Px);
    Fp_set (P0->Y, Py);

    Fp_set (P->X, Qx);
    Fp_set (P->Y, Qy);

    printf ("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b (P0, 0)); 
    printf ("Affine point P on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b (P, 1));
    Fp_weierstrass_aff2proj (P, P);
    printf ("Projective point P on curve?: %d\n", Fp_weierstrass_oncurve_jacproj_select_b (P, 1));
    Fp_copy (R->X, P->X);
    Fp_copy (R->Y, P->Y);
    Fp_copy (R->Z, P->Z);

    cp6_miller_initialize_config (B[1]);
    cp6_finalexpo_initialize_config ();

    printf("\nNAF(m) = ");
    for (i=0; i<cp6_config.lenm; i++) {
      printf ("%d", cp6_config.nafm[i]);
    }
    printf("\n");

	printf("\n(p+1)/4 = ");
	for (i = 0; i<64 * Fp_config.m->n; i++) {
		printf("%d", BIT(Fp_config.sqrt_exp_p3mod4, i));
	}
	printf("\n");
	printf("%d bits\n", Fp_config.sqrt_exp_len_p3mod4);

	for (i = 0; i < 1000; i++) {
		Fp_rand(alph);
		//Fp_set_ui(alph, (uint64_t) i);
		//Fp_neg(alph, alph);
		//printf("\nalpha = "); Fp_print(alph);
		Fp_sqrt(bet, alph);
		//printf("\nbeta = "); Fp_print(bet);
		Fp_mul(gam, bet, bet);
		//printf("\nbeta^2 = "); Fp_print(gam);
		if (Fp_cmp(alph, gam) != 0) {
			Fp_neg(gam, gam);
			if (Fp_cmp(alph, gam) != 0) {
				printf("Sqrt error!");
			}
		}
	}

    for (i=0; i<1; i++) {
      Fp_weierstrass_aff2proj (P0, P0);
      Fp_weierstrass_dbl (Q0, P0);
      Fp_weierstrass_add (Q0, Q0, P0);
      Fp_weierstrass_dbl (Q0, Q0);
      Fp_weierstrass_dbl (Q0, Q0);
      Fp_weierstrass_add (Q0, Q0, P0);
      if (Fp_weierstrass_oncurve_jacproj_select_b (Q0, 0) != 1) {
        printf ("G1 error\n");
        getchar ();
      }
      Fp_weierstrass_jacproj2aff (Q0, Q0);
      Fp_weierstrass_dbl (Q, P);
      Fp_weierstrass_add (Q, Q, P);
      Fp_weierstrass_dbl (Q, Q);
      Fp_weierstrass_dbl (Q, Q);
      Fp_weierstrass_add (Q, Q, P);
      //Fp_wp_neg (Q, Q);
      if (Fp_weierstrass_oncurve_jacproj_select_b (Q, 1) != 1) {
        printf ("G2 error\n");
        getchar ();
      }
      Fp_weierstrass_jacproj2aff (Q, Q);
    }

    Fp6q_init (f1);
    Fp6q_init (f2);
    Fp6q_init (f3);
    Fp6q_set_zero (f1);
    Fp6q_set_zero (f2);
    Fp6q_set_zero (f3);

    //printf ("f1 ="); Fp6q_print (f1);
    //printf ("f2 ="); Fp6q_print (f2);

    cp6_optimal_ate_miller (f1, Q, P0);
    cp6_optimal_ate_miller (f2, P, Q0);
    //printf ("f1 ="); Fp6q_print (f1);
    //printf ("f2 ="); Fp6q_print (f2);

    cp6_finalexpo (f1, f1);
    cp6_finalexpo (f2, f2);
    //Fp6q_inv (f1,f1);

    printf ("e([k]Q, P) = e(Q, [k]P)? %d\n", Fp6q_cmpeq (f1, f2));
    if (Fp6q_cmpeq (f1, f2) != 1) {
      printf ("CP6 pairing error\n");
      getchar ();
    }

    start_mil = __rdtsc();
    for (i=0; i < 1000; i++) {
      cp6_optimal_ate_miller (f1, Q, P0);
    }
    end_mil = __rdtsc();
    start_fe = __rdtsc();
    for (i=0; i < 1000; i++) {
      cp6_finalexpo (f3, f1);
    }
    end_fe = __rdtsc();
    
    printf ("Time CP6 OPTIMAL_ATE:  %fM cycles\n", ((uint32_t) (((double) end_mil-start_mil) / 1000.0)) / 1000000.);
    printf ("Time CP6 FINAL_EXPO:  %fM cycles\n", ((uint32_t) (((double) end_fe-start_fe) / 1000.0)) / 1000000.);
    printf ("Time CP6 PAIRING:  %fM cycles\n", ((uint32_t) (((double) end_fe-start_mil) / 1000.0)) / 1000000.);


    Fp_wp_free (P);
    Fp_wp_free (Q);
    Fp_wp_free (R);
    Fp_free_weierstrass ();
    Fp_wp_free (P0);
    Fp_wp_free (Q0);
    Fp_free (B[0]);
    Fp_free (B[1]);
    Fp_free (B3);

    cp6_finalexpo_free_config ();
    cp6_miller_free_config ();
    Fp6q_free_config ();
    Fp3_free_config ();
    Fp_free_config ();
  }
  else if (PAIR_CURVE == CP6tiny) {
	Fp_t B[2], B3;
	Point_wp_t P, Q, R;
	Fp6q_t f1, f2, f3;
	Fp_t alph, bet, gam;

	uint64_t pCP6[2] = { 0x2108F75F57168FE7, 0x76AD76A42A2 };
	/* This is the 107-bit CP6 prime p = 3 (mod 4). */

	uint64_t Px[2] = { 0xF25743F1668655EF, 0x555807C161F };
	uint64_t Py[2] = { 0x8D2D33A7D83B9E10, 0x44BCCF08BBC };
	/* This is the point over Fp in G1. */

	uint64_t Qx[2] = { 0x182624C621C68307, 0x3665C10BD94 };
	uint64_t Qy[2] = { 0x2F8E180B89361A5D, 0x32BC9CF3EDE };

	n = 2;
	printf("Running CP6tiny pairing...\n");

	Fp_initialize_config(pCP6, n);
	Fp3_initialize_config();
	Fp6q_initialize_config();

	Fp_init(B[0]);
	Fp_init(B[1]);
	Fp_init(B3);
	Fp_init(alph);
	Fp_init(bet);
	Fp_init(gam);
	Fp_wp_init(P0);
	Fp_wp_init(Q0);

	Fp_set_ui(B[0], 13);
	Fp_set_ui(B[1], 39);

	Fp_initialize_weierstrass(B, 2);
	Fp_wp_init(P);
	Fp_wp_init(Q);
	Fp_wp_init(R);

	Fp_set(P0->X, Px);
	Fp_set(P0->Y, Py);

	Fp_set(P->X, Qx);
	Fp_set(P->Y, Qy);

	printf("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b(P0, 0));
	printf("Affine point P on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b(P, 1));
	Fp_weierstrass_aff2proj(P, P);
	printf("Projective point P on curve?: %d\n", Fp_weierstrass_oncurve_jacproj_select_b(P, 1));
	Fp_copy(R->X, P->X);
	Fp_copy(R->Y, P->Y);
	Fp_copy(R->Z, P->Z);

	cp6_miller_initialize_config(B[1]);
	cp6_finalexpo_initialize_config();

	printf("\nNAF(m) = ");
	for (i = 0; i<cp6_config.lenm; i++) {
		printf("%d", cp6_config.nafm[i]);
	}
	printf("\n");

	printf("\n(p+1)/4 = ");
	for (i = 0; i<64 * Fp_config.m->n; i++) {
		printf("%d", BIT(Fp_config.sqrt_exp_p3mod4, i));
	}
	printf("\n");
	printf("%d bits\n", Fp_config.sqrt_exp_len_p3mod4);

	for (i = 0; i < 1000; i++) {
		Fp_rand(alph);
		//Fp_set_ui(alph, (uint64_t) i);
		//Fp_neg(alph, alph);
		//printf("\nalpha = "); Fp_print(alph);
		Fp_sqrt(bet, alph);
		//printf("\nbeta = "); Fp_print(bet);
		Fp_mul(gam, bet, bet);
		//printf("\nbeta^2 = "); Fp_print(gam);
		if (Fp_cmp(alph, gam) != 0) {
			Fp_neg(gam, gam);
			if (Fp_cmp(alph, gam) != 0) {
				printf("Sqrt error!");
			}
		}
	}

	for (i = 0; i<1; i++) {
		Fp_weierstrass_aff2proj(P0, P0);
		Fp_weierstrass_dbl(Q0, P0);
		Fp_weierstrass_add(Q0, Q0, P0);
		Fp_weierstrass_dbl(Q0, Q0);
		Fp_weierstrass_dbl(Q0, Q0);
		Fp_weierstrass_add(Q0, Q0, P0);
		if (Fp_weierstrass_oncurve_jacproj_select_b(Q0, 0) != 1) {
			printf("G1 error\n");
			getchar();
		}
		Fp_weierstrass_jacproj2aff(Q0, Q0);
		Fp_weierstrass_dbl(Q, P);
		Fp_weierstrass_add(Q, Q, P);
		Fp_weierstrass_dbl(Q, Q);
		Fp_weierstrass_dbl(Q, Q);
		Fp_weierstrass_add(Q, Q, P);
		//Fp_wp_neg (Q, Q);
		if (Fp_weierstrass_oncurve_jacproj_select_b(Q, 1) != 1) {
			printf("G2 error\n");
			getchar();
		}
		Fp_weierstrass_jacproj2aff(Q, Q);
	}

	Fp6q_init(f1);
	Fp6q_init(f2);
	Fp6q_init(f3);
	Fp6q_set_zero(f1);
	Fp6q_set_zero(f2);
	Fp6q_set_zero(f3);

	//printf ("f1 ="); Fp6q_print (f1);
	//printf ("f2 ="); Fp6q_print (f2);

	cp6b_optimal_ate_miller(f1, Q, P0);
	cp6b_optimal_ate_miller(f2, P, Q0);
	//printf ("f1 ="); Fp6q_print (f1);
	//printf ("f2 ="); Fp6q_print (f2);

	cp6b_finalexpo(f1, f1);
	cp6b_finalexpo(f2, f2);
	//Fp6q_inv (f1,f1);

	printf("e([k]Q, P) = e(Q, [k]P)? %d\n", Fp6q_cmpeq(f1, f2));
	if (Fp6q_cmpeq(f1, f2) != 1) {
		printf("CP6 pairing error\n");
		getchar();
	}

	start_mil = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6_optimal_ate_miller(f1, Q, P0);
	}
	end_mil = __rdtsc();
	start_fe = __rdtsc();
	for (i = 0; i < 1000; i++) {
		cp6_finalexpo(f3, f1);
	}
	end_fe = __rdtsc();

	printf("Time CP6tiny OPTIMAL_ATE:  %fM cycles\n", ((uint32_t)(((double)end_mil - start_mil) / 1000.0)) / 1000000.);
	printf("Time CP6tiny FINAL_EXPO:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_fe) / 1000.0)) / 1000000.);
	printf("Time CP6tiny PAIRING:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_mil) / 1000.0)) / 1000000.);

	getchar();

	Fp_wp_free(P);
	Fp_wp_free(Q);
	Fp_wp_free(R);
	Fp_free_weierstrass();
	Fp_wp_free(P0);
	Fp_wp_free(Q0);
	Fp_free(B[0]);
	Fp_free(B[1]);
	Fp_free(B3);

	cp6_finalexpo_free_config();
	cp6_miller_free_config();
	Fp6q_free_config();
	Fp3_free_config();
	Fp_free_config();
  }
 else if (PAIR_CURVE == CP3) {
	 Fp_t B[2];
	 Point_wp_t P, Q, R;
	 Fp3_t f1, f2, f3;
	 //uint64_t *premem;

	 uint64_t pCP3[16] = { 0x1500000005D641EF, 0xF27A90001804027E,
		 0xDF0B05082F4B95A0, 0x0ACAC924F81EC1DA,
		 0xA8EF36C595CF5E4F, 0xEB5F7B2090C91C8C,
		 0xFDAE1FBE2E9344F8, 0x6E7B70C7C1E1ED72,
		 0xA39DF221BDDF6026, 0x2BF6589BEC58AF41,
		 0xA90B68306AB9785C, 0x7301607FD616A3FD,
		 0xCEB7B890F7C56B98, 0x2DEF7BDF8234EAE9,
		 0xB7368AE71D33BC04, 0x540531971EED5543 };

	 /* This is the 1023-bit CP3 prime p = 3 (mod 4). */

	 uint64_t Px[16] = { 0x758CE21AD1C2BE3E, 0x0C820B5BC4666991,
		 0x32E11228028216A5, 0x37940D6905381AF1,
		 0x6D44464C1E9B33A9, 0xCD0A244DF694893F,
		 0x6071A8DFDEE42983, 0x54B5E323C7DED5F4,
		 0x5FA7BC828317C453, 0x0C56050AF0713A5E,
		 0x9A9F35C0FF516FAA, 0x670D3319658D4763,
		 0x0AD3A71077DACE33, 0x86CEF75311F2B675,
		 0xBCC91FAABDAEC15E, 0x1CE3057DA88BBCBA };
	 uint64_t Py[16] = { 0x4BBBA36B15E3351E, 0xE2323BE5C035DFAF,
		 0x39E74D8B8E656050, 0x998F123AE159D451,
		 0x3F9285C69B85C8AE, 0x6C8BBED349A2BBE3,
		 0x1B44946ADA9E7A7A, 0x3660AC7781F16A5E,
		 0x8D103674B090D9C6, 0xBBA04490A82FFEC0,
		 0x3DA002C85F8C83C7, 0x77B05A4EB3A02B7F,
		 0x9A9CF058FAE676C4, 0xABD0D34D251CE709,
		 0x1794BED9C920ED53, 0x29852D209F1373A5 };

	 uint64_t Qx[16] = { 0xC5929CFB23259A13, 0xBECC2D7FD8E166C3,
		 0x6AD6BEA4E8EC2385, 0x6266AE6B695D7533,
		 0xEEABA40CC791ED22, 0x975FE9EE29D69B73,
		 0x6F6AAB22E490C9B0, 0x34DF16327FF1B1B5,
		 0xA37E0C4D542D99D1, 0xBC39D34D3C64089B,
		 0x20AB141639629A0C, 0x539C4C7EF16062CD,
		 0xE3692158F74F7A9C, 0x047BCE81FD8CACDD,
		 0x8629D4E07C56595E, 0x0BADCF277AB45C44 };
	 uint64_t Qy[16] = { 0x4E41A9F542CD76B7, 0xD3AE6839C73230FB,
		 0xC89BA875B3579302, 0x460F55F967EB2EC0,
		 0x8BC3B3AC34E61A53, 0xB619994E08BABD9D,
		 0x9708A88BB55A3F37, 0x7755F188258AFC24,
		 0xAD9E0E9C3E01AE76, 0xE80CA9292A5770EB,
		 0x5324028AA7AB6170, 0xF39FAC955522EDCF,
		 0x1A8EB2D3A0A9F7ED, 0xD5F0368FB00D9F08,
		 0xB8C9CDB6E13D0407, 0x21E83F4C17A96BEC };

	 n = 16;
	 printf("Running CP3 pairing...\n");

	 Fp_initialize_config(pCP3, n);
	 Fp3_initialize_config();

	 Fp_init(B[0]);
	 Fp_init(B[1]);

	 Fp_set_ui(B[0], 1);
	 Fp_set_ui(B[1], 11);

	 Fp_initialize_weierstrass(B, 2);
	 Fp_wp_init(P0);
	 Fp_wp_init(Q0);
	 Fp_wp_init(P);
	 Fp_wp_init(Q);
	 Fp_wp_init(R);

	 Fp_set(P0->X, Px);
	 Fp_set(P0->Y, Py);
	 Fp_set_ui(P0->Z, 1);

	 Fp_set(P->X, Qx);
	 Fp_set(P->Y, Qy);
	 Fp_set_ui(P->Z, 1);

	 printf("Affine point P0 on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b(P0, 0));
	 printf("Affine point P on curve?: %d\n", Fp_weierstrass_oncurve_aff_select_b(P, 1));
	 Fp_weierstrass_aff2proj(P, P);
	 printf("Projective point P on curve?: %d\n", Fp_weierstrass_oncurve_jacproj_select_b(P, 1));
	 Fp_copy(R->X, P->X);
	 Fp_copy(R->Y, P->Y);
	 Fp_copy(R->Z, P->Z);

	 for (i = 0; i<1; i++) {
		 Fp_weierstrass_aff2proj(P0, P0);
		 Fp_weierstrass_dbl(Q0, P0);
		 Fp_weierstrass_add(Q0, Q0, P0);
		 Fp_weierstrass_dbl(Q0, Q0);
		 Fp_weierstrass_dbl(Q0, Q0);
		 Fp_weierstrass_add(Q0, Q0, P0);
		 Fp_weierstrass_dbl(Q0, P0);
		 Fp_weierstrass_add(Q0, Q0, P0);
		 Fp_weierstrass_dbl(Q0, Q0);
		 Fp_weierstrass_dbl(Q0, Q0);
		 Fp_weierstrass_add(Q0, Q0, P0);
		 if (Fp_weierstrass_oncurve_jacproj_select_b(Q0, 0) != 1) {
			 printf("G1 error\n");
			 getchar();
		 }
		 Fp_weierstrass_jacproj2aff(Q0, Q0);

		 Fp_weierstrass_dbl(Q, P);
		 Fp_weierstrass_add(Q, Q, P);
		 Fp_weierstrass_dbl(Q, Q);
		 Fp_weierstrass_dbl(Q, Q);
		 Fp_weierstrass_add(Q, Q, P);
		 Fp_weierstrass_dbl(Q, P);
		 Fp_weierstrass_add(Q, Q, P);
		 Fp_weierstrass_dbl(Q, Q);
		 Fp_weierstrass_dbl(Q, Q);
		 Fp_weierstrass_add(Q, Q, P);
		 //Fp_wp_neg (Q, Q);
		 if (Fp_weierstrass_oncurve_jacproj_select_b(Q, 1) != 1) {
			 printf("G2 error\n");
			 getchar();
		 }
		 Fp_weierstrass_jacproj2aff(Q, Q);
	 }

	 cp3_miller_initialize_config(B[1]);
	 cp3_finalexpo_initialize_config();

	 Fp3_init(f1);
	 Fp3_init(f2);
	 Fp3_init(f3);
	 /*
	 printf("T = ");
	 Print(cp3_config.T, 8); printf("\n\n");

	 printf("a0 = ");
	 Print(cp3_config.a0, 8); printf("\n");
	 printf("a1 = ");
	 Print(cp3_config.a1, 9); printf("\n");
	 printf("a2 = ");
	 Print(cp3_config.a2, 8); printf("\n");
	 printf("c1 = ");
	 Fp_print(cp3_config.c1); printf("\n");
	 printf("s11 = ");
	 Fp_print(cp3_config.s11); printf("\n");
	 */
	 printf("NAF(a2) = ");
	 for (i = 0; i<cp3_config.naflena2; i++) {
		 printf("%d", cp3_config.nafa2[i]);
	 }
	 printf("\n\n\n");

	 //printf("\nPoint P0: "); Fp_wp_print(P0);
	 //printf("\nPoint P: "); Fp_wp_print(P);
	 //printf("\n\n");

	 cp3_ate_miller(f1, Q, P0);
	 cp3_ate_miller(f2, P, Q0);
	 //printf ("f1 ="); Fp3_print (f1);
	 //printf ("f2 ="); Fp3_print (f2);

	 cp3_finalexpo(f1, f1);
	 cp3_finalexpo(f2, f2);

	 //printf("f1 ="); Fp3_print(f1);
	 //printf("f2 ="); Fp3_print(f2);

	 printf("e([k]Q, P) = e(Q, [k]P)? %d\n", Fp3_cmpeq(f1, f2));
	 if (Fp3_cmpeq(f1, f2) != 1) {
		 printf("CP3 pairing error\n");
		 getchar();
	 }

	 start_mil = __rdtsc();
	 for (i = 0; i < 1000; i++) {
		 cp3_ate_miller(f1, Q, P0);
	 }
	 end_mil = __rdtsc();
	 start_fe = __rdtsc();
	 for (i = 0; i < 1000; i++) {
		 cp3_finalexpo(f3, f1);
	 }
	 end_fe = __rdtsc();

	 printf("Time CP3 ATE:  %fM cycles\n", ((uint32_t)(((double)end_mil - start_mil) / 1000.0)) / 1000000.);
	 printf("Time CP3 FINAL_EXPO:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_fe) / 1000.0)) / 1000000.);
	 printf("Time CP3 PAIRING:  %fM cycles\n", ((uint32_t)(((double)end_fe - start_mil) / 1000.0)) / 1000000.);

	 cp3_finalexpo_free_config();
	 cp3_miller_free_config();

	 Fp_wp_free(P);
	 Fp_wp_free(Q);
	 Fp_wp_free(R);
	 Fp_free_weierstrass();

	 Fp_wp_free(P0);
	 Fp_wp_free(Q0);
	 Fp_free(B[0]);
	 Fp_free(B[1]);

	 Fp3_free_config();
	 Fp_free_config();
 } else {
      printf ("ERROR: PAIR_CURVE not valid!\n");
  }
  
  getchar ();
}
#endif
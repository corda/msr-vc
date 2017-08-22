/* curves/Fp2_weierstrass.h
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
 
#ifndef __FP2_WEIERSTRASS
#define __FP2_WEIERSTRASS

#include "../uint.h"
#include "../Fp2/Fp2.h"

typedef struct {
  Fp2_t X;
  Fp2_t Y;
  Fp2_t Z;
} point_wp2_t;

typedef point_wp2_t Point_wp2_t[1];

typedef struct {
	Fp2_t X;
	Fp2_t Y;
} point_wp2aff_t;

typedef point_wp2aff_t Point_wp2aff_t[1];


int Fp2_initialize_weierstrass (Fp2_t b);
void Fp2_free_weierstrass (void);

void Fp2_weierstrass_dbl (Point_wp2_t c, Point_wp2_t a);
void Fp2_weierstrass_add (Point_wp2_t c, Point_wp2_t a, Point_wp2_t b);
void Fp2_weierstrass_mix (Point_wp2_t c, Point_wp2_t a, Point_wp2_t b);
void Fp2_weierstrass_aff (Point_wp2_t c, Point_wp2_t a, Point_wp2_t b);
void Fp2_weierstrass_affdbl(Point_wp2_t c, Point_wp2_t a);
void Fp2_weierstrass_affadd(Point_wp2_t c, Point_wp2_t a, Point_wp2_t b);
int Fp2_weierstrass_oncurve_jacproj (Point_wp2_t x);
int Fp2_weierstrass_oncurve_homproj (Point_wp2_t x);
int Fp2_weierstrass_oncurve_aff (Point_wp2_t x);
int Fp2_weierstrass_equal_aff (Point_wp2_t x, Point_wp2_t y);
int Fp2_weierstrass_equal_jacproj(Point_wp2_t a, Point_wp2_t b);
void Fp2_weierstrass_aff2proj (Point_wp2_t c, Point_wp2_t a);  
void Fp2_weierstrass_jacproj2aff (Point_wp2_t c, Point_wp2_t a);
void Fp2_weierstrass_homproj2aff (Point_wp2_t c, Point_wp2_t a);

int Fp2_wp_init (Point_wp2_t p);
void Fp2_wp_free (Point_wp2_t p);
void Fp2_wp_copy (Point_wp2_t c, Point_wp2_t a);
void Fp2_wp_neg (Point_wp2_t c, Point_wp2_t a);

int Fp2_wpaff_init(Point_wp2aff_t p);
void Fp2_wpaff_free(Point_wp2aff_t p);
void Fp2_wpaff_copy(Point_wp2aff_t c, Point_wp2aff_t a);
void Fp2_wpaff_neg(Point_wp2aff_t c, Point_wp2aff_t a);

#endif /* __FP2_WEIERSTRASS */
/* curves/Fp_twistededwards.h
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
 
#ifndef __FP_WEIERSTRASS
#define __FP_WEIERSTRASS

#include "../uint.h"
#include "../Fp/Fp.h"

typedef struct {
  Fp_t X;
  Fp_t Y;
  Fp_t Z;
} point_wp_t;

typedef point_wp_t Point_wp_t[1];

typedef struct {
	Fp_t X;
	Fp_t Y;
} point_wpaff_t;

typedef point_wpaff_t Point_wpaff_t[1];

int Fp_initialize_weierstrass (Fp_t *b, int num);
void Fp_free_weierstrass (void);

void Fp_weierstrass_dbl (Point_wp_t c, Point_wp_t a);
void Fp_weierstrass_add (Point_wp_t c, Point_wp_t a, Point_wp_t b);
void Fp_weierstrass_mix (Point_wp_t c, Point_wp_t a, Point_wp_t b); // b is affine
void Fp_weierstrass_aff (Point_wp_t c, Point_wp_t a, Point_wp_t b);
void Fp_weierstrass_affdbl(Point_wp_t c, Point_wp_t a);
void Fp_weierstrass_affadd(Point_wp_t c, Point_wp_t a, Point_wp_t b);
int Fp_weierstrass_oncurve_jacproj (Point_wp_t x);
int Fp_weierstrass_oncurve_aff (Point_wp_t x);
int Fp_weierstrass_oncurve_jacproj_select_b (Point_wp_t x, int num);
int Fp_weierstrass_oncurve_aff_select_b (Point_wp_t x, int num);
void Fp_weierstrass_aff2proj (Point_wp_t c, Point_wp_t a);  
void Fp_weierstrass_jacproj2aff (Point_wp_t c, Point_wp_t a);
int Fp_weierstrass_equal_aff(Point_wp_t a, Point_wp_t b);
int Fp_weierstrass_equal_jacproj(Point_wp_t a, Point_wp_t b);

int Fp_wp_init (Point_wp_t p);
void Fp_wp_free (Point_wp_t p);
void Fp_wp_copy (Point_wp_t c, Point_wp_t a);
void Fp_wp_neg (Point_wp_t c, Point_wp_t a);

int Fp_wpaff_init(Point_wpaff_t p);
void Fp_wpaff_free(Point_wpaff_t p);
void Fp_wpaff_copy(Point_wpaff_t c, Point_wpaff_t a);
void Fp_wpaff_neg(Point_wpaff_t c, Point_wpaff_t a);


#endif /* __FP_WEIERSTRASS */
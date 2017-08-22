#ifndef __FP_WEIERSTRASS
#define __FP_WEIERSTRASS

#include "Qapuint.h"
#include "QapFp.h"

typedef struct {
  QapFp_t X;
  QapFp_t Y;
} point_wp_t;

typedef point_wp_t QapPoint_wp_t[1];

int QapFp_initialize_weierstrass (QapFp_t *b, int num);
void QapFp_free_weierstrass (void);

void QapFp_weierstrass_mix (QapPoint_wp_t c, QapPoint_wp_t a, QapPoint_wp_t b); // b is affine
void QapFp_weierstrass_affdbl(QapPoint_wp_t c, QapPoint_wp_t a);
void QapFp_weierstrass_affadd(QapPoint_wp_t c, QapPoint_wp_t a, QapPoint_wp_t b);
int QapFp_weierstrass_oncurve_aff (QapPoint_wp_t x);
int QapFp_weierstrass_equal_aff(QapPoint_wp_t x, QapPoint_wp_t y);

int QapFp_wp_init (QapPoint_wp_t p);
void QapFp_wp_free (QapPoint_wp_t p);
void QapFp_wp_copy (QapPoint_wp_t c, QapPoint_wp_t a);
void QapFp_wp_neg (QapPoint_wp_t c, QapPoint_wp_t a);

#endif /* __FP_WEIERSTRASS */
#ifndef __FP2_WEIERSTRASS
#define __FP2_WEIERSTRASS

#include "Qapuint.h"
#include "QapFp2.h"

typedef struct {
  QapFp2_t X;
  QapFp2_t Y;
} point_wp2_t;

typedef point_wp2_t QapPoint_wp2_t[1];

int QapFp2_initialize_weierstrass (QapFp2_t b);
void QapFp2_free_weierstrass (void);

void QapFp2_weierstrass_affdbl(QapPoint_wp2_t c, QapPoint_wp2_t a);
void QapFp2_weierstrass_affadd(QapPoint_wp2_t c, QapPoint_wp2_t a, QapPoint_wp2_t b);
int QapFp2_weierstrass_oncurve_aff (QapPoint_wp2_t x);
int QapFp2_weierstrass_equal_aff (QapPoint_wp2_t x, QapPoint_wp2_t y);

int QapFp2_wp_init (QapPoint_wp2_t p);
void QapFp2_wp_free (QapPoint_wp2_t p);
void QapFp2_wp_copy (QapPoint_wp2_t c, QapPoint_wp2_t a);
void QapFp2_wp_neg (QapPoint_wp2_t c, QapPoint_wp2_t a);

#endif /* __FP2_WEIERSTRASS */
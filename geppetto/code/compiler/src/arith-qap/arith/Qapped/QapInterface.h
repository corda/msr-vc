#pragma once

#include "QapFp.h"
#include "QapFp2.h"
#include "QapFp_weierstrass.h"
#include "QapFp2_weierstrass.h"

typedef struct {
	QapFp_t m;
	QapFp_t b;
	int n;
	QapFp2_t B, B3, l00, l01, l10;

	// Projective versions of the Base and Twist generators (also a version of g^1)
	QapPoint_wp_t Lg;
	QapPoint_wp2_t Rg;

	// Frequently encoded elts
	QapPoint_wp_t Lzero;
	QapPoint_wp2_t Rzero;
} QapInterface;

void QapInterface_init(QapInterface* interface);

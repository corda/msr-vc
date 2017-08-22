#pragma once
#include "Poly.h"
#include <memory>

/* Performs polynomial interpolation. Given a set of degree+1 points (x_i,y_i),
     finds a polynomial of degree 'degree' that passes through them all. The 
     set of x_i values are constrained to be in a geometric progression.
   In typical use cases, x values are likely to be fixed ahead of time, with 
     many different y sequences being used over the same x sequence. So we 
     precompute as much as possible using just the x values, which is then 
     shared over many interpolations.
  In case x values are not reused, a convenience wrapper is provided in
     Poly::interpolateGeometric in Poly.cpp

  Excuse the long name, but I don't expect it to be typed out very frequently
*/

class PolyGeometricInterpolator 
{
	int degree, size;
	Field* field;
	FieldElt* a;
	std::unique_ptr<FieldEltArray> r, q;
	std::unique_ptr<FieldEltArray> ui;
	Poly tInv, wiPoly;
	std::unique_ptr<FieldEltArray> uc, zc;
	mutable Poly wcPoly;

	FieldElt *one, *zero;
public:
	// Constructors also perform precomputation based on the x-values
	// Note: xr^i values must be unique for i in [0,degree]
	PolyGeometricInterpolator(Field* field, FieldElt* x0, FieldElt* xr, int degree);
	PolyGeometricInterpolator(const PolyGeometricInterpolator&) = delete;
  ~PolyGeometricInterpolator();

	// The y values for the actual interpolation. y->length() must be degree+1
	Poly* operator()(FieldEltArray* y) const;
};
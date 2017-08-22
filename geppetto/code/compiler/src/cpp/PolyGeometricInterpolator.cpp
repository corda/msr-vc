#include "PolyGeometricInterpolator.h"


// Slow O(n^2) polynomial multiplication used for debugging
void slowMul(Field* field, Poly& a, Poly& b, Poly& c)
{
	FieldElt* t = field->newElt();
	FieldEltArray *ac = a.coefficients, *bc = b.coefficients,
		          *cc = c.coefficients;

	// Poly::zero(c, a.getDegree() + b.getDegree());
	for (int i = 0; i <= a.getDegree() + b.getDegree(); ++i)
		field->zero(cc->elt(i));

	for (int i = 0; i <= a.getDegree(); ++i)
		for (int j = 0; j <= b.getDegree(); ++j)
		{
			field->mul(ac->elt(i), bc->elt(j), t);
			field->add(t, cc->elt(i + j), cc->elt(i + j));
		}

  field->delElt(t);
}

//#define SLOW_MUL
#ifndef SLOW_MUL
#define polyMul(a,b,c) Poly::mul((a),(b),(c))
#else
// Ugly macro, but I needed to capture "field"
#define polyMul(a,b,c) slowMul(field,(a),(b),(c))
#endif

// Implementation is as explained in the following paper:
//   Alin Bostan and Eric Schost. "Polynomial evaluation and interpolation 
//   on special sets of points". Journal of Complexity.

PolyGeometricInterpolator::PolyGeometricInterpolator(
	Field* field, FieldElt* x0, FieldElt* xr, int degree)
	: degree(degree), size(degree+1), field(field)
	, a(field->newElt())                  // x0
	, r (field->newEltArray(size, true))  // Powers of xr
	, q (field->newEltArray(size, true))  // q[i] = xr^(i*(i-1)/2) = q[i-1]*r[i-1]

	, ui(field->newEltArray(size, true))  // Used in interpolation 
	, wiPoly(field,degree)                //   wrt Newtonian basis
	, tInv(field, degree)                 //

	, uc(field->newEltArray(size, true))  // Used for converting to
	, zc(field->newEltArray(size, true))  //   the monomial basis
	, wcPoly(field,degree)
  , one(field->newElt()), zero(field->newElt())
{
	FieldEltArray* tPolyCoeff = tInv.coefficients;

	// Grr, Field api uses FieldElt*
  FieldElt* t  = field->newElt();
  FieldElt* t2 = field->newElt();
	field->one(one);
	field->zero(zero);
	field->copy(x0, a);

	field->one(r->elt(0));
	field->one(q->elt(0));
	field->one(ui->elt(0));
	field->one(uc->elt(0));
	field->one(zc->elt(0));
	field->one(tPolyCoeff->elt(0));

	for (int i = 1; i < size; ++i)
	{
		// r[i] = r[i-1]*xr; q[i] = q[i-1]*r[i-1]
		field->mul(r->elt(i - 1), xr, r->elt(i));
		field->mul(q->elt(i - 1), r->elt(i - 1), q->elt(i));

		// ui[i] = ui[i-1] * (r[i]-1)
		field->sub(r->elt(i), one, t); 
		field->mul(ui->elt(i - 1), t, ui->elt(i));

		// tPolyCoeff[i] = (-1)^i * q[i] / u[i]
		FieldElt* tpi = tPolyCoeff->elt(i);
		field->div(q->elt(i), ui->elt(i), tpi);
		if (i % 2) field->sub(zero, tpi, tpi);

		// uc[i] = uc[i-1]r[i] / (1-r[i])
		field->mul(uc->elt(i - 1), r->elt(i), t);
		field->sub(one, r->elt(i), t2);
		field->div(t, t2, uc->elt(i));

		// zc[i] = (-1)^i * u[i]/q[i]
		FieldElt* zci = zc->elt(i);
		field->div(uc->elt(i), q->elt(i), zci);
		if (i % 2) field->sub(zero, zci, zci);
	}

	/*
	// Debugging block:
	for (int i = 0; i < size; ++i)
	{
		printf("q:");
		field->print(q->elt(i)); printf(" u:");
		field->print(ui->elt(i)); printf(" t:");
		field->print(tInv.coefficients->elt(i)); printf("\n");
	}
	*/
  field->delElt(t);
  field->delElt(t2);
}


PolyGeometricInterpolator::~PolyGeometricInterpolator() {
  field->delElt(a);
  field->delElt(one);
  field->delElt(zero);
}


Poly* PolyGeometricInterpolator::operator()(FieldEltArray* y) const
{
	assert(degree+1 == y->length());

  FieldElt* telt = field->newElt();  // temp var
	Poly gPoly(field, 2 * degree);
	Poly ucRev(field, degree), gcPoly(field, 2 * degree);
	unique_ptr<Poly> result(new Poly(field, degree));

	FieldEltArray* wi = wiPoly.coefficients;
	for (int i = 0; i < size; ++i)
		field->div(y->elt(i), ui->elt(i), wi->elt(i));

	polyMul((Poly&)wiPoly, (Poly&)tInv, gPoly);    // <-- Expensive step 1

	/*
	// Debugging block
	printf(" ---- y-values ----\n");
	for (int i = 0; i < size; ++i)
	{
		field->print(y->elt(i)); 
		printf(" ");
	}
	printf("\n");

	printf(" ---- Newtonian coefficients ---\n");
	for (int i = 0; i < size; ++i)
	{
		field->div(gPoly.coefficients->elt(i), q->elt(i), t);
		field->print(t); printf(" ");
	}
	printf("\n");
	*/

	FieldEltArray *vc = gPoly.coefficients, *wc = wcPoly.coefficients;
	for (int i = 0; i < size; ++i)
	{
		assert(i < vc->length());
		if (i % 2) field->sub(zero, vc->elt(i), vc->elt(i));
		field->div(vc->elt(i), uc->elt(i), wc->elt(i));
		// reversed because mul^t
		field->copy(uc->elt(i), ucRev.coefficients->elt(degree - i));
	}
	wcPoly.setCoefficients(wc, degree);
	polyMul(ucRev, wcPoly, gcPoly);                // <-- Expensive step 2

	FieldElt* ap = field->newElt();
	int ilim = gcPoly.getDegree()-degree;
	field->one(ap);
	for (int i = 0; i < size; ++i)
	{
		// Poly::mul(a,b,c) does not ensure c.getDegree()=a.getDegree()+b.getDegree(), 
		// specially if a or b has zero as leading coefficients
		assert(i>ilim || i<gcPoly.coefficients->length());
		field->mul(i>ilim?zero:gcPoly.coefficients->elt(i + degree), zc->elt(i), telt);
    field->div(telt, ap, result->coefficients->elt(i));
		field->mul(a, ap, ap);
	}
  field->delElt(telt);
  field->delElt(ap);
	return result.release();
}
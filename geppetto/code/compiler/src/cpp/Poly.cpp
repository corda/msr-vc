#include "Poly.h"
#include "PolyGeometricInterpolator.h"
#include "Util.h"
#include "EnigmaWrappers.h"
#include <stdio.h>
#include <assert.h>

// Microsoft's bignum library expects NDEBUG to be defined (or badness occurs!)
#ifndef NDEBUG
  #define UNDEFINENDEBUG
  #define NDEBUG
#endif
#include <msbignum.h>
#include <field.h>
#include "pfc.h"
#include "modfft.h"
#ifdef UNDEFINENDEBUG
  #undef NDEBUG
  #undef UNDEFINENDEBUG
#endif

void Poly::init(Field* field, int degree) {
	this->degree = degree;
	this->size = degree + 1;
	this->coefficients = field->newEltArray(this->size, true); 
	this->field = field;
}

Poly::Poly(Field* field, int degree) {
	init(field, degree);
}

Poly::Poly(Field* field, int degree, int* coeffs) {
	init(field, degree);
	
	for (int i = 0; i < this->size; i++) {
 		field->set(this->coefficients->elt(i), coeffs[i]);
	}
}

Poly::~Poly() {
	delete coefficients;
}

void Poly::setDegree(int d) {
	reserveSpace(d);
	this->degree = d;
}

// Set the coefficients to the array provided
void Poly::setCoefficients(FieldEltArray* coeff, int degree) {
	if (coeff != coefficients) {
		delete coefficients;
		coefficients = coeff;
		size = degree + 1;		
	}

	this->degree = degree;
}

// Get a pointer to an array large enough to hold coeff for a degree d poly
FieldEltArray* Poly::getCoefficientArray(int degree) {
	if (degree + 1 > size) {
		return field->newEltArray(degree+1, true);
	} else {
		return coefficients;
	}
}


// Reserve enough space to ensure we can represent a poly of degree d
void Poly::reserveSpace(int d) {
	if (d + 1 > size) {	// Need more room!
		delete coefficients;
		coefficients = field->newEltArray(d + 1, true);
		size = d + 1;
	} 
}

void Poly::print() {
	printf("Degree %d, (", degree);
	for (int i = 0; i < degree + 1; i++) {
		field->print(coefficients->elt(i));
		printf("\t");
	}
	printf(")");
}
	
bool Poly::equals(Poly& p) {
	int minLen = min(degree, p.getDegree()) + 1;

	for (int i = 0; i < minLen; i++) {
		if (!field->equal(coefficients->elt(i), p.coefficients->elt(i))) {
			return false;
		}
	}

	// Check that the rest of the coefficients are 0 for the poly with larger degree
	FieldElt* zero = field->newElt();
	field->zero(zero);
	for (int i = minLen; i < max(degree, p.getDegree()) + 1; i++) {
		if (degree < p.getDegree()) {
			if (!field->equal(p.coefficients->elt(i), zero)) {
				return false;
			}
		} else {
			if (!field->equal(coefficients->elt(i), zero)) {
				return false;
			}
		}
	}
	field->delElt(zero);

	return true;	
}


void Poly::firstDerivative(Poly& poly) {
	firstDerivative(poly, poly);
}

void Poly::firstDerivative(Poly& poly, Poly& result) {
	FieldEltArray* resultCoeff = result.getCoefficientArray(poly.getDegree()-1);

	FieldElt* degree = poly.field->newElt();
	for (int d = 1; d < poly.getDegree() + 1; d++) {
		// result.coefficients->elt(d-1) = poly.coefficients->elt(d) * d		
		poly.field->set(degree, d);
		poly.field->mul(poly.coefficients->elt(d), degree, resultCoeff->elt(d-1));
	}
	poly.field->delElt(degree);
	//zero(resultCoeff[poly.getDegree()]);
	result.setCoefficients(resultCoeff, poly.getDegree()-1);
}

// Implements high-school-style multiplication as a sanity check on the FFT version
void Poly::mulSlow(Poly& a, Poly& b, Poly& result) {
	assert(a.field == b.field && b.field == result.field);
	int resultDegree = a.getDegree() + b.getDegree();
	FieldEltArray* resultCoeff = result.getCoefficientArray(resultDegree);

	a.field->zero(resultCoeff, resultDegree + 1);

	FieldElt* temp = a.field->newElt();
	for (int da = 0; da < a.getDegree() + 1; da++) {
		for (int db = 0; db < b.getDegree() + 1; db++) {
			a.field->mul(a.coefficients->elt(da), b.coefficients->elt(db), temp);
			a.field->add(resultCoeff->elt(da+db), temp, resultCoeff->elt(da+db));
		}
	}
	a.field->delElt(temp);

	result.setCoefficients(resultCoeff, resultDegree);
}

// Multiply a polynomial p by a constant c
void Poly::constMul(Poly& p, FieldElt* c, Poly& result) {
	assert(p.field == result.field);
	FieldEltArray* resultCoeff = result.getCoefficientArray(p.getDegree());

	for (int d = 0; d < p.getDegree() + 1; d++) {
		p.field->mul(p.coefficients->elt(d), c, resultCoeff->elt(d));
	}

	result.setCoefficients(resultCoeff, p.getDegree());
}

// Makes a polynomial monic.  Sets c to be coeff of x^maxdegree
// TODO: Optimize if p is already monic
void Poly::makeMonic(Poly& p, Poly& result, FieldElt* leadingCoeff) {
	assert(p.field == result.field);
	FieldElt* inverse, *uno, *nada;
	inverse = p.field->newElt();
	uno = p.field->newElt();
	nada = p.field->newElt();

	p.field->one(uno);
	p.field->zero(nada);

	// Trim any leading 0s
	int resultDegree = p.getDegree();
	while (p.field->equal(p.coefficients->elt(resultDegree), nada) && resultDegree > 0) {
		resultDegree--;
	}

	result.setDegree(resultDegree);

	for (int d = 0; d < resultDegree + 1; d++) {
		p.field->copy(p.coefficients->elt(d), result.coefficients->elt(d));
	}
	
	
	if (p.field->equal(result.coefficients->elt(result.getDegree()), uno)) {
		// Already monic
		p.field->one(leadingCoeff);
	} else {
		// Monicize by dividing all coefficients by the leading coefficient
		p.field->copy(result.coefficients->elt(result.getDegree()), leadingCoeff);
		p.field->div(uno, leadingCoeff, inverse);
		
		constMul(result, inverse, result);
	}

	p.field->delElt(uno);
	p.field->delElt(inverse);
	p.field->delElt(nada);
}

bool Poly::isZero(Poly& p) {
	FieldElt* nada = p.field->newElt();
	p.field->zero(nada);

	for (int d = 0; d < p.getDegree() + 1; d++) {
		if (!p.field->equal(p.coefficients->elt(d), nada)) {
			p.field->delElt(nada);
			return false;
		}
	}
	p.field->delElt(nada);
	return true;
}

// Computes: result <- a * b
void Poly::mul(Poly& a, Poly& b, Poly& result) {
	//mulSlow(a, b, result);
	assert(a.field == b.field && b.field == result.field);

	FieldElt* one = a.field->newElt();
	a.field->one(one);

	bool aMonic = a.field->equal(a.coefficients->elt(a.getDegree()), one);
	bool bMonic = a.field->equal(b.coefficients->elt(b.getDegree()), one);
	bool monic = aMonic && bMonic;

	a.field->delElt(one);

	// Special case the easiest cases
	if (a.getDegree() == 0 && b.getDegree() == 0) {
		result.setDegree(0);
		a.field->mul(a.coefficients->elt(0), b.coefficients->elt(0), result.coefficients->elt(0));
	} else if (isZero(a) || isZero(b)) {
		result.setDegree(a.getDegree() + b.getDegree());
		zero(result, result.getDegree());
	} else {
		// Make both inputs monic
		Poly monicA(a.field), monicB(b.field);
		FieldElt* monicCorrection, *aCorrection, *bCorrection;
		monicCorrection = a.field->newElt();
		aCorrection = a.field->newElt();
		bCorrection = b.field->newElt();

		makeMonic(a, monicA, aCorrection);			
		makeMonic(b, monicB, bCorrection);
		a.field->mul(aCorrection, bCorrection, monicCorrection);

		result.setDegree(monicA.getDegree() + monicB.getDegree());

		// Special case the small, easy values
		if (monicA.getDegree() == 1 && monicB.getDegree() == 1) {  // (x+a)(x+b)
			a.field->one(result.coefficients->elt(2));                                                   // x^2
			a.field->add(monicA.coefficients->elt(0), monicB.coefficients->elt(0), result.coefficients->elt(1));	// (a+b)x^1
			a.field->mul(monicA.coefficients->elt(0), monicB.coefficients->elt(0), result.coefficients->elt(0));	// (a*b)x^0		
		} else if (monicA.getDegree() == 2 && monicB.getDegree() == 1) {  // (x^2+ax+b)(x+c) = x^3 + ax^2 + cx^2 + bx + acx + bc
			a.field->one(result.coefficients->elt(3));	                                            	// x^3
			a.field->add(monicA.coefficients->elt(1), monicB.coefficients->elt(0), result.coefficients->elt(2));	// (a+c)x^2
			a.field->mul(monicA.coefficients->elt(1), monicB.coefficients->elt(0), result.coefficients->elt(1));
			a.field->add(result.coefficients->elt(1), monicA.coefficients->elt(0), result.coefficients->elt(1));	// (ac+b)x^1
			a.field->mul(monicA.coefficients->elt(0), monicB.coefficients->elt(0), result.coefficients->elt(0));	// (b*c)x^0		
		} else if (monicA.getDegree() == 2 && monicB.getDegree() == 2) { // (x^2+ax+b)(x^2+cx+d) = x^4 + ax^3 + bx^2 + cx^3 + acx^2 + bcx + dx^2 + adx + bd
			FieldElt* temp;
			temp = a.field->newElt();
			a.field->one(result.coefficients->elt(4));	                                            	// x^4
			a.field->add(monicA.coefficients->elt(1), monicB.coefficients->elt(1), result.coefficients->elt(3));	// (a+c)x^3
			a.field->mul(monicA.coefficients->elt(1), monicB.coefficients->elt(1), result.coefficients->elt(2));
			a.field->add(result.coefficients->elt(2), monicA.coefficients->elt(0), result.coefficients->elt(2));
			a.field->add(result.coefficients->elt(2), monicB.coefficients->elt(0), result.coefficients->elt(2));	// (ac + b + d)x^2
			a.field->mul(monicA.coefficients->elt(0), monicB.coefficients->elt(1), temp);
			a.field->mul(monicA.coefficients->elt(1), monicB.coefficients->elt(0), result.coefficients->elt(1));
			a.field->add(result.coefficients->elt(1), temp, result.coefficients->elt(1));	// (bc+ad)x^1
			a.field->mul(monicA.coefficients->elt(0), monicB.coefficients->elt(0), result.coefficients->elt(0)); // (bd)x^0
			a.field->delElt(temp);
		} else {	// Use an FFT	
			bigctx_t BignumCtx;
			bigctx_t* pbigctx = &BignumCtx;
			memset(pbigctx, 0, sizeof(bigctx_t));

			int lenA = monicA.getDegree() + 1;
			int lenB = monicB.getDegree() + 1;

			const int lg_convol_lng = (int)significant_bit_count((digit_t)(lenA + lenB - 1));
			int convol_lng = (int)1 << lg_convol_lng;
			convol_lng *= (int)a.field->get_fft_info()->info.nprime;

			// TODO: Cache this memory to avoid allocating every time
			FieldEltArray* scratch1 = a.field->newEltArray(convol_lng, true);
			FieldEltArray* scratch2 = a.field->newEltArray(convol_lng, true);

			lenA--;
			lenB--;

			if (lenA + lenB - 1 > result.getDegree()) {
				assert(0);
				// We need more space to hold the intermediate result
				// TODO: Find a more elegant solution
				FieldEltArray* temp = a.field->newEltArray(lenA + lenB, true);
				modpoly_monic_product(monicA.coefficients->get_digits(), lenA, monicB.coefficients->get_digits(), lenB, temp->get_digits(), scratch1->get_digits(), scratch2->get_digits(), &a.field->get_fft_info()->info, PBIGCTX_PASS);
				// Take the result mod x^(result.degree+1)
				for (int i = 0; i < result.getDegree() + 1; i++) {
					a.field->copy(temp->elt(i), result.coefficients->elt(i));
				}	
				delete temp;
			} else {
				modpoly_monic_product(monicA.coefficients->get_digits(), lenA, monicB.coefficients->get_digits(), lenB, result.coefficients->get_digits(), scratch1->get_digits(), scratch2->get_digits(), &a.field->get_fft_info()->info, PBIGCTX_PASS);
				a.field->one(result.coefficients->elt(result.getDegree()));	// modpoly_monic_product doesn't bother with the x.getDegree() coefficient, since it assumes it's monic
			}		

			delete scratch1;
			delete scratch2;
		}

		if (!monic) { // Need to reverse the corrections we made earlier to convert a and b to monic polys
			constMul(result, monicCorrection, result);
		}

		a.field->delElt(aCorrection);
		a.field->delElt(bCorrection);
		a.field->delElt(monicCorrection);		
	}

}


// result <- a + b
void Poly::add(Poly& a, Poly& b, Poly& result) {
	assert(a.field == b.field && b.field == result.field);

	int resultDegree = max(a.getDegree(), b.getDegree());
	FieldEltArray* resultCoeff = result.getCoefficientArray(resultDegree);

	for (int i = 0; i < min(a.getDegree(), b.getDegree()) + 1; i++) {
		a.field->add(a.coefficients->elt(i), b.coefficients->elt(i), resultCoeff->elt(i));
	}

	if (a.getDegree() < b.getDegree()) {
		for (int i = a.getDegree()+1; i < b.getDegree() + 1; i++) {
			a.field->copy(b.coefficients->elt(i), resultCoeff->elt(i));
		}
	} else if (b.getDegree() < a.getDegree() + 1) {
		for (int i = b.getDegree()+1; i < a.getDegree()+1; i++) {
			a.field->copy(a.coefficients->elt(i), resultCoeff->elt(i));
		}
	}

	result.setCoefficients(resultCoeff, resultDegree);
}

// result <- a - b
void Poly::sub(Poly& a, Poly& b, Poly& result) {
	assert(a.field == b.field && b.field == result.field);

	int resultDegree = max(a.getDegree(), b.getDegree());
	FieldEltArray* resultCoeff = result.getCoefficientArray(resultDegree);

	for (int i = 0; i < min(a.getDegree(), b.getDegree()) + 1; i++) {
		a.field->sub(a.coefficients->elt(i), b.coefficients->elt(i), resultCoeff->elt(i));
	}

	if (a.getDegree() < b.getDegree()) {
		FieldElt* nada;
		nada = a.field->newElt();
		a.field->zero(nada);

		for (int i = a.getDegree()+1; i < b.getDegree() + 1; i++) {
			a.field->sub(nada, b.coefficients->elt(i), resultCoeff->elt(i));
		}
		a.field->delElt(nada);
	} else if (b.getDegree() < a.getDegree()) {
		for (int i = b.getDegree()+1; i < a.getDegree() + 1; i++) {
			a.field->copy(a.coefficients->elt(i), resultCoeff->elt(i));
		}
	}

	result.setCoefficients(resultCoeff, resultDegree);
}

// Reverse into a fresh polynomial
void Poly::reverse(Poly& f, Poly& fRev) { 
	assert(f.field == fRev.field);
	fRev.setDegree(f.getDegree());
	
	for (int i = 0; i < f.getDegree()+1; i++) {
		f.field->copy(f.coefficients->elt(i), fRev.coefficients->elt((f.getDegree()+1) - i - 1));
	}
}

// Reverse in place
void Poly::reverse(Poly& f) { 
	FieldElt* tmp = f.field->newElt();
	for (int i = 0; i < (f.getDegree()+1) / 2; i++) {
		f.field->copy(f.coefficients->elt(i), tmp);
		f.field->copy(f.coefficients->elt((f.getDegree()+1) - i - 1), f.coefficients->elt(i));
		f.field->copy(tmp, f.coefficients->elt((f.getDegree()+1) - i - 1));
	}
	f.field->delElt(tmp);
}

// Compute f mod x^degree
void Poly::reduce(Poly& f, int degree) {
	assert(degree > 0);

	if (f.getDegree() >= degree) {
		f.setDegree(degree - 1);
	}
	// Else, no effect
}

// Create a polynomial that has all 0 coefficients
void Poly::zero(Poly& f, int degree) {
	f.setDegree(degree);	
	f.field->zero(f.coefficients, degree+1);
}

// Compute: result <- f / x^degree
void Poly::div(Poly& f, int degree, Poly& result) {
	assert(f.field == result.field);

	if (f.getDegree() - degree < 0) {	// Should only happen if we're dividing the 0 poly
		FieldElt* nada = f.field->newElt();
		f.field->zero(nada);
		for (int i = 0; i < f.getDegree(); i++) {
			assert(f.field->equal(f.coefficients->elt(i), nada));
		}
		result.setDegree(0);
		f.field->zero(result.coefficients->elt(0));
		f.field->delElt(nada);
	} else {
		result.setDegree(f.getDegree() - degree);	

		for (int i = 0; i < result.getDegree()+1; i++) {
			f.field->copy(f.coefficients->elt(i+degree), result.coefficients->elt(i));
		}
	}
}

// Compute result <- f^-1 mod x.getDegree()
// Algorithm based on http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf
void Poly::invert(Poly& f, int degree, Poly& result) {
	assert(f.field == result.field);
#ifdef _DEBUG
	FieldElt* uno = f.field->newElt(); f.field->one(uno);
	assert(f.field->equal(f.coefficients->elt(0), uno));
	f.field->delElt(uno);
#endif // _DEBUG

	if (degree == 0) {
		assert(0); // TODO
	} else if (degree == 1) {	// Because we know the constant coeff is 1, Inv(f) mod x = 1
		result.setDegree(0);
		f.field->one(result.coefficients->elt(0));
	} else {
		int splitDegree;

		if (degree % 2 == 1) {
			splitDegree = degree - 1;
		} else {
			splitDegree = degree / 2;
		}

		Poly inverse0(f.field), inverse0sqr(f.field), inverse1(f.field);
		invert(f, splitDegree, inverse0);	// inverse0 * f = 1 mod x^{d/2}
		mul(inverse0, inverse0, inverse0sqr);

		// Let f = f0 + f1 * x^{d/2}
		Poly f0(f.field), f1(f.field);
		f0.setDegree(splitDegree - 1);
		f1.setDegree(degree - splitDegree);		
		for (int i = 0; i < f0.getDegree() + 1; i++) {
			if (i < f.getDegree() + 1) {
				f.field->copy(f.coefficients->elt(i), f0.coefficients->elt(i));
			} else { // f is 0 here
				f.field->zero(f0.coefficients->elt(i));
			}
		}
		for (int i = 0; i < f1.getDegree() + 1; i++) {
			int fIndex = i + splitDegree;
			if (fIndex < f.getDegree() + 1) {
				f.field->copy(f.coefficients->elt(fIndex), f1.coefficients->elt(i));
			} else { // f is 0 here
				f.field->zero(f1.coefficients->elt(i));
			}			
		}		

		// inverse1 = -(inverse0^2 * f1 + inverse0 * b), where b = x^{-d/2} (inverse0 * f0 - 1)
		// => inverse1 = -(inverse0^2 * f1 + x^{-d/2} (inverse0^2 * f0 - inverse0))		
		Poly temp(f.field), nada(f.field);
		zero(nada, splitDegree);
		mul(inverse0sqr, f0, temp);		// temp = inverse0^2 * f0
		sub(temp, inverse0, temp);		// temp = (inverse0^2 * f0 - inverse0)
		div(temp, splitDegree, temp);	// temp = x^{-d/2} (inverse0^2 * f0 - inverse0)
		mul(inverse0sqr, f1, inverse1);	// inverse1 = inverse0^2 * f1
		add(inverse1, temp, inverse1);	// inverse1 = (inverse0^2 * f1 + x^{-d/2} (inverse0^2 * f0 - inverse0))
		zero(nada, inverse1.getDegree());
		sub(nada, inverse1, inverse1);	// inverse1 *= -1
		
		// Finally, result = inverse0 + x^{splitDegree} inverse1
		result.setDegree(min(degree - 1, inverse1.getDegree() + splitDegree));
		for (int i = 0; i < splitDegree; i++) {
			int inverse0index = i;
			if (inverse0index < inverse0.getDegree() + 1) {
				f.field->copy(inverse0.coefficients->elt(i), result.coefficients->elt(i));
			} else {
				f.field->zero(result.coefficients->elt(i));
			}
		}
		for (int i = splitDegree; i < result.getDegree()+1; i++) {
			int inverse1index = i - splitDegree;
			if (inverse1index < inverse1.getDegree() + 1) {
				f.field->copy(inverse1.coefficients->elt(inverse1index), result.coefficients->elt(i));
			} else {
				f.field->zero(result.coefficients->elt(i));
			}
		}
	}
}

// Computes r, q such that f = q*g + r
// Algorithm based on http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf
// Note: Rev(q) = Rev(g)^-1 * Rev(f) mod (x^(deg(f)+deg(g)))
// Requires that g is monic
void Poly::mod(Poly& f, Poly& g, Poly& q, Poly& r) {
	assert(f.field == g.field && g.field == q.field && q.field == r.field);

	Poly revF(f.field), revG(f.field);

#ifdef _DEBUG
	FieldElt* uno = f.field->newElt(); f.field->one(uno);
	assert(f.field->equal(g.coefficients->elt(g.getDegree()), uno));
	f.field->delElt(uno);
#endif // _DEBUG
	assert(g.getDegree() >= 1);	// Can't take mods of a constant

	if (f.getDegree() < g.getDegree()) {  // No work needs to be done
		FieldEltArray* rCoeff = r.getCoefficientArray(f.getDegree());		
		for (int i = 0; i < f.getDegree() + 1; i++) {
			f.field->copy(f.coefficients->elt(i), rCoeff->elt(i));
		}
		zero(q, g.getDegree());
		r.setCoefficients(rCoeff, f.getDegree());	
	} else {
		reverse(f, revF);
		reverse(g, revG);

		invert(revG, f.getDegree()-g.getDegree()+1, revG);	
		mul(revF, revG, q);
		reduce(q, f.getDegree()-g.getDegree()+1);

		reverse(q);

		Poly temp(f.field);
		mul(g, q, temp);	
		sub(f, temp, r);	// r <- f - q*g
		reduce(r, g.getDegree());
	}
}

// Generates a random polynomial of specified degree, or random degree if degree=-1
void Poly::polyRand(Field* field, Poly& r, int max_degree, int degree) {
	if (degree == -1) {
		r.setDegree((int)((rand() / (double) RAND_MAX) * max_degree));
	} else {
		r.setDegree(degree);
	}

	for (int i = 0; i < r.getDegree() + 1; i++) {
		field->assignRandomElt(r.coefficients->elt(i));
	}

	if (rand() % 2 == 0) {
		// Make it monic, since there are some odd corner cases there
		field->one(r.coefficients->elt(r.getDegree()));
	}
}

void Poly::multiEval(Poly& f, FieldEltArray* x, int len, FieldEltArray* result) {
	PolyTree* tree = new PolyTree(f.field, x, len);
	multiEval(f, x, len, result, tree, 0, 0, tree->height - 1);
	delete tree;
}


// Evaluate f at x, which contains len points, and place the result in the result array, starting at index resultIndex.
// Uses a precomputed tree of polynomials, which we traverse top down
void Poly::multiEval(Poly& f, FieldEltArray* x, int len, FieldEltArray* result, PolyTree* tree, int resultIndex, int polyIndex, int depth) {
	if (len == 0) {
		assert(0);
	} else if (depth == 0) {
		f.field->copy(f.coefficients->elt(0), result->elt(resultIndex));
	} else {
		int newLen = (int) pow(2.0, (int)(my_log2(len)-1));

		if (len < pow(2.0, depth) + pow(2.0, depth-1)) {	// Divide in two
			Poly r0(f.field), r1(f.field), q(f.field);
			mod(f, *tree->polys[depth-1][polyIndex], q, r0);
			mod(f, *tree->polys[depth-1][polyIndex+1], q, r1);

			multiEval(r0, x, newLen, result, tree, resultIndex, 2*polyIndex, depth-1);
			multiEval(r1, &x[newLen], len - newLen, result, tree, resultIndex+newLen, 2*(polyIndex+1), depth-1);
		} else {	// Divide in three
			Poly r0(f.field), r1(f.field), r2(f.field), q(f.field);
			mod(f, *tree->polys[depth-1][polyIndex], q, r0);
			mod(f, *tree->polys[depth-1][polyIndex+1], q, r1);
			mod(f, *tree->polys[depth-1][polyIndex+2], q, r2);

			multiEval(r0, x, newLen, result, tree, resultIndex, 2*polyIndex, depth-1);
			multiEval(r1, &x[newLen], newLen, result, tree, resultIndex+newLen, 2*(polyIndex+1), depth-1);
			multiEval(r2, &x[2*newLen], len - 2*newLen, result, tree, resultIndex+2*newLen, 2*(polyIndex+2), depth-1);
		}
	}
}

// Interpolate the len points (x_i,y'_i), using the auxilliary info in the polynomial tree, placing the result in Poly& result
// Interpolate with respect to the polynomials at level depth in the tree, and index polyIndex within that particular level
void Poly::interpolateRecursive(Field* field, FieldEltArray* y, int y_index, PolyTree* tree, Poly& result, int len, int depth, int polyIndex) {
	if (len == 0) {
		assert(0);
	} else if (depth == 0) {
		result.setDegree(0);
		field->copy(y->elt(y_index), result.coefficients->elt(0));
	} else {
		int newLen = (int) pow(2.0, (int)(my_log2(len)-1));

		if (len < pow(2.0, depth) + pow(2.0, depth-1)) {	// Divide in two
			Poly r0(field), r1(field);
		
			interpolateRecursive(field, y, y_index, tree, r0, newLen, depth-1, 2*polyIndex);
			interpolateRecursive(field, y, y_index + newLen, tree, r1, len - newLen, depth-1, 2*(polyIndex+1));

			Poly result0(field), result1(field);
			mul(r0, *tree->polys[depth-1][polyIndex+1], result0);
			mul(r1, *tree->polys[depth-1][polyIndex], result1);
			add(result0, result1, result);
		} else {	// Divide in three
			Poly r0(field), r1(field), r2(field);
		
			interpolateRecursive(field, y, y_index + 0		  , tree, r0, newLen, depth-1, 2*polyIndex);
			interpolateRecursive(field, y, y_index + newLen   , tree, r1, newLen, depth-1, 2*(polyIndex+1));
			interpolateRecursive(field, y, y_index + 2*newLen , tree, r2, len - 2*newLen, depth-1, 2*(polyIndex+2));

			Poly result0(field), result1(field), result2(field);
			//mul(r0, *tree->polys[depth-1][levelindex+1], result0);
			//mul(r1, *tree->polys[depth-1][levelindex], result1);
			//mul(r2, *tree->polys[depth-1][levelindex+2], result2);
			//add(result0, result1, result);
			//add(result, result2, result);

			// Proposed fix.  TODO: Should precompute these extra multiplications
			mul(*tree->polys[depth-1][polyIndex+1], *tree->polys[depth-1][polyIndex+2], result0);
			mul(r0, result0, result0);

			mul(*tree->polys[depth-1][polyIndex], *tree->polys[depth-1][polyIndex+2], result1);
			mul(r1, result1, result1);

			mul(*tree->polys[depth-1][polyIndex], *tree->polys[depth-1][polyIndex+1], result2);
			mul(r2, result2, result2);

			add(result0, result1, result);
			add(result, result2, result);
		}
	}
}

//void Poly::interpolateRecursive(Field* field, FieldEltArray* y, PolyTree* tree, Poly& result, int len, int depth, int polyIndex) {
//	if (len == 0) {
//		assert(0);
//	} else if (depth == 0) {
//		result.setDegree(0);
//		field->copy(y->elt(0), result.coefficients->elt(0));
//	} else {
//		int newLen = (int) pow(2.0, (int)(my_log2(len)-1));
//
//		if (len < pow(2.0, depth) + pow(2.0, depth-1)) {	// Divide in two
//			Poly r0(field), r1(field);
//		
//			interpolateRecursive(field, y->elt(0), tree, r0, newLen, depth-1, 2*polyIndex);
//			interpolateRecursive(field, y->elt(newLen), tree, r1, len - newLen, depth-1, 2*(polyIndex+1));
//
//			Poly result0(field), result1(field);
//			mul(r0, *tree->polys[depth-1][polyIndex+1], result0);
//			mul(r1, *tree->polys[depth-1][polyIndex], result1);
//			add(result0, result1, result);
//		} else {	// Divide in three
//			Poly r0(field), r1(field), r2(field);
//		
//			interpolateRecursive(field, y->elt(0), tree, r0, newLen, depth-1, 2*polyIndex);
//			interpolateRecursive(field, y->elt(newLen), tree, r1, newLen, depth-1, 2*(polyIndex+1));
//			interpolateRecursive(field, y->elt(2*newLen), tree, r2, len - 2*newLen, depth-1, 2*(polyIndex+2));
//
//			Poly result0(field), result1(field), result2(field);
//			//mul(r0, *tree->polys[depth-1][levelindex+1], result0);
//			//mul(r1, *tree->polys[depth-1][levelindex], result1);
//			//mul(r2, *tree->polys[depth-1][levelindex+2], result2);
//			//add(result0, result1, result);
//			//add(result, result2, result);
//
//			// Proposed fix.  TODO: Should precompute these extra multiplications
//			mul(*tree->polys[depth-1][polyIndex+1], *tree->polys[depth-1][polyIndex+2], result0);
//			mul(r0, result0, result0);
//
//			mul(*tree->polys[depth-1][polyIndex], *tree->polys[depth-1][polyIndex+2], result1);
//			mul(r1, result1, result1);
//
//			mul(*tree->polys[depth-1][polyIndex], *tree->polys[depth-1][polyIndex+1], result2);
//			mul(r2, result2, result2);
//
//			add(result0, result1, result);
//			add(result, result2, result);
//		}
//	}
//}


// Given a set of len points (x_i,y_i), find the degree n+1 polynomial that they correspond to
// using the precalculated polynomial tree and Lagrange denominators ( denominator_i = \prod_{j!=i} (x_i-x_j) )
// Note that: f(u) = \sum (y_i / denominator_i) * \prod_{j!=i} (u - x_j)
Poly* Poly::interpolate(Field* field, FieldEltArray* y, int len, PolyTree* polyTree, FieldEltArray* denominators) {
	// Divide to get y'_i, where f(u) = \sum y'_i * \prod_{j!=i} (u - x_j)
	FieldEltArray* yPrime = field->newEltArray(len, true); 
	for (int i = 0; i < len; i++) {
		field->div(y->elt(i), denominators->elt(i), yPrime->elt(i));
		//mul(divisors[i], y[i], yPrime[i]);
	}

	Poly* result = new Poly(field);	
	interpolateRecursive(field, yPrime, 0, polyTree, *result, len, polyTree->height - 1, 0);

	delete yPrime;
	return result;
}

FieldEltArray* Poly::genLagrangeDenominators(Field* field, Poly& func, PolyTree* tree, FieldEltArray* x, int len) {
	// Differentiate P_{k0}
	Poly fPrime(field);
	firstDerivative(func, fPrime);

	// Multipoint eval P'_{k0} at x_i to obtain the Lagrange denominators, where denominator_i = \prod_{j!=i} (x_i-x_j)
	FieldEltArray* denominators = field->newEltArray(len, true); 
	multiEval(fPrime, x, len, denominators, tree, 0, 0, tree->height - 1);

	/*	// Don't bother with this.  Be lazy.  Let the interpolate step do the division
	// Invert to get divisor_i = \prod_{j!=i} [1 / (x_i-x_j)]
	FieldElt uno;
	one(uno);
	for (int i = 0; i < len; i++) {
		div(uno, divisors[i], divisors[i]);
	}
	*/

	return denominators;
}

// Given a set of len points (x_i,y_i), find the degree len+1 polynomial that they correspond to
// Implements Algorithm 2 from http://www.mpi-inf.mpg.de/~csaha/lectures/lec6.pdf
Poly* Poly::interpolate(Field* field, FieldEltArray* x, FieldEltArray* y, int len) {
	PolyTree* tree = new PolyTree(field, x, len);

	FieldEltArray* denominators = genLagrangeDenominators(field, *tree->polys[tree->height-1][0], tree, x, len);

	Poly* result = interpolate(field, y, len, tree, denominators);
	
	// Cleanup
	delete tree;
	delete denominators;

	return result;
}

// Computes the coefficient representation of the Lagrange basis polynomial indicated by index
// Observe that: L_index(x) = \prod_{i != index} 1/(s_index - s_i) * \prod_{i != index} (x - s_i)
// The left term has no dependence on x, so we can compute that separately and then apply to the
// coefficients that come from the right term
FieldEltArray* Poly::genLagrangeCoeffs(Field* field, FieldEltArray* evalPts, int len, int index) {
	FieldElt* zero, *tmp;
	zero = field->newElt();
	tmp  = field->newElt();
	field->zero(zero);
	FieldEltArray* coefficients = field->newEltArray(len, true);
	field->zero(coefficients, len);
	
	field->one(coefficients->elt(0));	
	// Compute (x-evalPt[0])(x-evalPt[1]) ... (x-evalPt[d]), one term at a time
	int rootCount = 0;	// Separate from i, since we skip one of the roots
	for (int i = 0; i < len; i++) { // For each root
		if (i == index) {
			continue;	
		} 
		rootCount++;
		// Current coefficient = prevValShiftedOne - root*prevVal
		for (int j = rootCount; j >= 0; j--) {
			if (j == 0) {
				field->sub(zero, evalPts->elt(i), tmp);	// Multiply by -1*r_i
				field->mul(tmp, coefficients->elt(0), coefficients->elt(0));
			} else {
				field->mul(evalPts->elt(i), coefficients->elt(j), coefficients->elt(j));
				field->sub(coefficients->elt(j-1), coefficients->elt(j), coefficients->elt(j));
			}
		}
	}

	// Compute the divisor Prod_{i <> j} (x_i - x_j)
	FieldElt* prod;
	prod = field->newElt();
	field->one(prod);
	for (int i = 0; i < len; i++) { // For each root
		if (i == index) {
			continue;	
		} 
		field->sub(evalPts->elt(index), evalPts->elt(i), tmp);
		field->mul(tmp, prod, prod);		
	}

	// Apply the divisor to all of the coefficients
	for (int i = 0; i < len; i++) { // For each root
		field->div(coefficients->elt(i), prod, coefficients->elt(i));	
	}

	field->delElt(zero);
	field->delElt(tmp);
	field->delElt(prod);	
	return coefficients;
}


//void Poly::Initialize(Field* field) {
//	pbigctx = &BignumCtx;
//	memset(pbigctx, 0, sizeof(bigctx_t));
//
//	modfft_init(field->get_modulus(), &fftinfo, PBIGCTX_PASS);
//}
//
//void Poly::Cleanup() {
//	modfft_uninit(&fftinfo, PBIGCTX_PASS);
//}


Poly* Poly::interpolateGeometric(Field* field, FieldElt* x0, FieldElt* xr, FieldEltArray* y, int len)
{
	PolyGeometricInterpolator pgi(field, x0, xr, len-1);
	return pgi(y);
}

Poly* Poly::interpolateGeometric(PolyGeometricInterpolator& pgi, FieldEltArray* y)
{
	return pgi(y);
}
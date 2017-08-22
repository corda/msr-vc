#include "Field.h"
#include "Encoding.h"
#include "Poly.h"
#include <iostream>
#include <cassert>
#include <time.h>
#include "timer.h"
#include "FileArchiver.h"
#include "Archive.h"

#include <boost/program_options.hpp>
using namespace std;
using namespace boost::program_options;

void randomFieldAddSubTest(Field* field, int num_tests) {
	FieldElt* f1 = field->newElt();
	FieldElt* f2 = field->newElt();
	FieldElt* f3 = field->newElt();
	FieldElt* result = field->newElt();

	for (int trial = 0; trial < num_tests; trial++) {
		field->assignRandomElt(f1);
		field->assignRandomElt(f2);
		field->assignRandomElt(f3);
		
		// Test reversability: f1 - f2 - f3 + f2 + f3 == f1
		field->sub(f1, f2, result);
		field->sub(result, f3, result);
		field->add(result, f2, result);
		field->add(result, f3, result);

		assert(field->equal(f1, result));
	}

	field->delElt(f1);
	field->delElt(f2);
	field->delElt(f3);
	field->delElt(result);
}

void randomFieldAssocTest(Field* field, int num_tests) {
	FieldElt* f1 = field->newElt();
	FieldElt* f2 = field->newElt();
	FieldElt* f3 = field->newElt();
	FieldElt* tmp = field->newElt();	
	FieldElt* result1 = field->newElt();
	FieldElt* result2 = field->newElt();

	for (int trial = 0; trial < num_tests; trial++) {
		field->assignRandomElt(f1);
		field->assignRandomElt(f2);
		field->assignRandomElt(f3);
		
		// Test associativity of addition: f1 * (f2 + f3) == f1 * f2 + f1 * f3
		field->add(f2, f3, result1);
		field->mul(f1, result1, result1);

		field->mul(f1, f2, tmp);
		field->mul(f1, f3, result2);
		field->add(tmp, result2, result2);

		assert(field->equal(result1, result2));

		// Test associativity of subtraction: f1 * (f2 - f3) == f1 * f2 - f1 * f3
		field->sub(f2, f3, result1);
		field->mul(f1, result1, result1);

		field->mul(f1, f2, tmp);
		field->mul(f1, f3, result2);
		field->sub(tmp, result2, result2);

		assert(field->equal(result1, result2));
	}

	field->delElt(f1);
	field->delElt(f2);
	field->delElt(f3);
	field->delElt(tmp);
	field->delElt(result1);
	field->delElt(result2);
}

void randomFieldArrayTest(Field* field, int num_tests) {
	FieldEltArray* arr = field->newEltArray(num_tests, true);
	FieldElt* add_result1 = field->newElt();
	FieldElt* add_result2 = field->newElt();
	FieldElt* mul_result1 = field->newElt();
	FieldElt* mul_result2 = field->newElt();

	for (int i = 0; i < num_tests; i++) {
		// Randomize
		for (int j = 0; j < num_tests; j++) {
			field->assignRandomElt(arr->elt(j));
		}

		// Make sure we get the same result no matter how we access the array
		field->zero(add_result1);
		field->zero(add_result2);
		field->zero(mul_result1);
		field->zero(mul_result2);

		for (int j = 0; j < num_tests; j++) {
			field->add(add_result1, arr->elt(j),				 add_result1);
			field->add(add_result2, arr->elt(num_tests - j - 1), add_result2);

			field->mul(mul_result1, arr->elt(j),				 mul_result1);
			field->mul(mul_result2, arr->elt(num_tests - j - 1), mul_result2);
		}

		assert(field->equal(add_result1, add_result2));
		assert(field->equal(mul_result1, mul_result2));
	}

	delete arr;
	field->delElt(add_result1);
	field->delElt(add_result2);
	field->delElt(mul_result1);
	field->delElt(mul_result2);	
}

void randomPolyMulTest(Field* field, int num_tests, int max_degree) {
	Poly p1(field), p2(field), resultFast(field), resultSlow(field);

	for (int trial = 0; trial < num_tests; trial++) {
		Poly::polyRand(field, p1, max_degree);
		Poly::polyRand(field, p2, max_degree);		

		Poly::mul(p1, p2, resultFast);
		Poly::mulSlow(p1, p2, resultSlow);

		assert(resultFast.equals(resultSlow));
	}
}

void randomPolyAddSubTest(Field* field, int num_tests, int max_degree) {
	Poly p1(field), p2(field), p3(field), result(field);

	for (int trial = 0; trial < num_tests; trial++) {
		Poly::polyRand(field, p1, max_degree);
		Poly::polyRand(field, p2, max_degree);		
		Poly::polyRand(field, p3, max_degree);

		// Test reversability: p1 - p2 - p3 + p2 + p3 == p1
		Poly::sub(p1, p2, result);
		Poly::sub(result, p3, result);
		Poly::add(result, p2, result);
		Poly::add(result, p3, result);

		assert(result.equals(p1));
	}
}

void randomPolyAddMulTest(Field* field, int num_tests, int max_degree) {
	Poly p1(field), p2(field), p3(field), result1(field), result2(field), temp(field);

	for (int trial = 0; trial < num_tests; trial++) {
		Poly::polyRand(field, p1, max_degree);
		Poly::polyRand(field, p2, max_degree);		
		Poly::polyRand(field, p3, max_degree);

		// Test associativity: p1 * (p2 + p3) == p1 * p2 + p1 * p3
		Poly::add(p2, p3, result1);
		Poly::mul(p1, result1, result1);
		
		Poly::mul(p1, p2, temp);
		Poly::mul(p1, p3, result2);
		Poly::add(result2, temp, result2);

		assert(result1.equals(result2));
	}
}

void randomPolySubMulTest(Field* field, int num_tests, int max_degree) {
	Poly p1(field), p2(field), p3(field), result1(field), result2(field), temp(field);

	for (int trial = 0; trial < num_tests; trial++) {
		Poly::polyRand(field, p1, max_degree);
		Poly::polyRand(field, p2, max_degree);		
		Poly::polyRand(field, p3, max_degree);

		// Test associativity: p1 * (p2 - p3) == p1 * p2 - p1 * p3
		Poly::sub(p2, p3, result1);
		Poly::mul(p1, result1, result1);
		
		Poly::mul(p1, p2, temp);
		Poly::mul(p1, p3, result2);
		Poly::sub(temp, result2, result2);

		assert(result1.equals(result2));
	}
}

void randomPolyModTest(Field* field, int num_tests, int max_degree) {
	Poly f(field), m(field), q(field), r(field), result(field);

	for (int trial = 0; trial < num_tests; trial++) {
		int modDegree = (int)((rand() / (double) RAND_MAX) * max_degree) + 1;	// +1, b/c we need mod to be non-constant
		assert(modDegree >= 1);

		Poly::polyRand(field, f, max_degree);
		Poly::polyRand(field, m, max_degree, modDegree);		

		// Ensure m is monic
		field->one(m.coefficients->elt(m.getDegree()));

		Poly::mod(f, m, q, r);
		
		Poly::mul(q, m, result);
		Poly::add(result, r, result);

		assert(result.equals(f));
	}
}

// Setting geom to true fixes x values to powers of xr
void randomPolyEvalInterpTest(Field* field, int max_degree, int degree = -1, int numTests=100, bool geom=false) {

  printf("Interpolating polynomials with degree ");
  if (degree == -1) {
    printf("chosen at random\n");
  } else {
    printf("%d\n", degree);
  }

	Poly f(field);
	for (int trial = 0; trial < numTests; trial++) {
		// Pick a random function
		Poly::polyRand(field, f, max_degree, degree);

		// Pick some random points to evaluate at
		int numPts = f.getDegree() + 1;
		FieldEltArray* x = field->newEltArray(numPts, true);
		FieldEltArray* y = field->newEltArray(numPts, true);
		FieldElt* xr = field->newElt();
		field->set(xr, 3);

		for (int i = 0; i < numPts; i++) {
			if (!geom) field->assignRandomElt(x->elt(i));
			else if (i == 0) field->set(x->elt(i), 1);
			else field->mul(x->elt(i - 1), xr, x->elt(i));

			// Make sure this value is unique, since we'll be using it for interpolation
			bool unique;
			do {
				unique = true;
				for (int j = 0; j < i; j++) {
					if (field->equal(x->elt(i), x->elt(j))) {
						unique = false;
						field->assignRandomElt(x->elt(i));
						break;
					}
				}
			} while (!unique);
		}

		Poly::multiEval(f, x, numPts, y);
		Poly* fprime;
		if (!geom) fprime = Poly::interpolate(field, x, y, numPts);
		else fprime = Poly::interpolateGeometric(field, x->elt(0), xr, y, numPts);
		assert(f.equals(*fprime));
		    
		delete x;
		delete y;
    delete xr;
		delete fprime;
	}
}

void randomEncAddSubTest(Encoding* encoding, elt_t t, int num_tests) {
	Field* field = encoding->getSrcField();
	FieldElt* f1 = field->newElt();
	FieldElt* f2 = field->newElt();
	FieldElt* f3 = field->newElt();

	EncodedElt* e1 = encoding->new_elt(t);
	EncodedElt* e2 = encoding->new_elt(t);
	EncodedElt* e3 = encoding->new_elt(t);
  EncodedElt* t1 = encoding->new_elt(t);
  EncodedElt* t2 = encoding->new_elt(t);
  EncodedElt* tmp1 = encoding->new_elt(t);
  EncodedElt* tmp2 = encoding->new_elt(t);
	EncodedElt* result = encoding->new_elt(t);

	for (int trial = 0; trial < num_tests; trial++) {
		field->assignRandomElt(f1);
		field->assignRandomElt(f2);
		field->assignRandomElt(f3);
		
		encoding->encode(t, f1, e1);
		encoding->encode(t, f2, e2);
		encoding->encode(t, f3, e3);

		// Basic equality checks
		assert(encoding->equal(t, e1, e1));
		assert(encoding->equal(t, e2, e2));
		assert(encoding->equal(t, e3, e3));

		assert(!encoding->equal(t, e1, e2));
		assert(!encoding->equal(t, e2, e3));

		// Check that ordering doesn't affect result    
		encoding->add(t, e1, e2, result);
		encoding->add(t, result, e3, tmp1);

		encoding->add(t, e1, e3, result);
		encoding->add(t, result, e2, tmp2);
		assert(encoding->equal(t, tmp1, tmp2));

		encoding->sub(t, e1, e2, result);
		encoding->sub(t, result, e3, tmp1);

		encoding->sub(t, e1, e3, result);
		encoding->sub(t, result, e2, tmp2);
		assert(encoding->equal(t, tmp1, tmp2));

		// Test reversability: E(f1) - E(f2) - E(f3) + E(f2) + E(f3) == E(f1)
		encoding->sub(t, e1, e2, result);
		encoding->sub(t, result, e3, result);
		encoding->add(t, result, e2, result);
		encoding->add(t, result, e3, result);

		assert(encoding->equal(t, e1, result));
	}

	field->delElt(f1);
	field->delElt(f2);
	field->delElt(f3);

	encoding->del_elt(t, e1);
	encoding->del_elt(t, e2);
	encoding->del_elt(t, e3);
  encoding->del_elt(t, t1);
  encoding->del_elt(t, t2);
  encoding->del_elt(t, tmp1);
  encoding->del_elt(t, tmp2);
	encoding->del_elt(t, result);
}

void randomEncAssocTest(Encoding* encoding, elt_t t, int num_tests) {
	Field* field = encoding->getSrcField();
	FieldElt* f1 = field->newElt();
	FieldElt* f2 = field->newElt();
	FieldElt* f3 = field->newElt();
	FieldElt* tmp = field->newElt();	
	FieldElt* result1 = field->newElt();

	EncodedElt* e1 = encoding->new_elt(t);
	EncodedElt* e2 = encoding->new_elt(t);
	EncodedElt* e3 = encoding->new_elt(t);
	EncodedElt* resultA = encoding->new_elt(t);
	EncodedElt* resultB = encoding->new_elt(t);

	for (int trial = 0; trial < num_tests; trial++) {
		field->assignRandomElt(f1);
		field->assignRandomElt(f2);
		field->assignRandomElt(f3);

		encoding->encode(t, f1, e1);
		encoding->encode(t, f2, e2);
		encoding->encode(t, f3, e3);
		
		// Test associativity of addition: f1 * (E(f2) + E(f3)) == f1 * E(f2) + f1 * E(f3) == E(f1*(f2+f3))
		encoding->add(t, e2, e3, resultA);
		encoding->mul(t, resultA, f1, resultA);

		encoding->mul(t, e2, f1, e2);
		encoding->mul(t, e3, f1, e3);
		encoding->add(t, e2, e3, resultB);

		assert(encoding->equal(t, resultA, resultB));

		field->add(f2, f3, result1);
		field->mul(f1, result1, result1);
		encoding->encode(t, result1, resultB);
		assert(encoding->equal(t, resultA, resultB));

		// Test associativity of subtraction: f1 * (E(f2) - E(f3)) == f1 * E(f2) - f1 * E(f3) == E(f1*(f2-f3))
		encoding->encode(t, f1, e1);
		encoding->encode(t, f2, e2);
		encoding->encode(t, f3, e3);

		encoding->sub(t, e2, e3, resultA);
		encoding->mul(t, resultA, f1, resultA);

		encoding->mul(t, e2, f1, e2);
		encoding->mul(t, e3, f1, e3);
		encoding->sub(t, e2, e3, resultB);

		assert(encoding->equal(t, resultA, resultB));

		field->sub(f2, f3, result1);
		field->mul(f1, result1, result1);
		encoding->encode(t, result1, resultB);

		assert(encoding->equal(t, resultA, resultB));
	}

	field->delElt(f1);
	field->delElt(f2);
	field->delElt(f3);
	field->delElt(tmp);
	field->delElt(result1);

	encoding->del_elt(t, e1);
	encoding->del_elt(t, e2);
	encoding->del_elt(t, e3);
	encoding->del_elt(t, resultA);
	encoding->del_elt(t, resultB);
}


void randomEncMulTest(Encoding* encoding, int num_tests) {
	Field* field = encoding->getSrcField();
	FieldElt* f1 = field->newElt();
	FieldElt* f2 = field->newElt();
	FieldElt* f3 = field->newElt();
	FieldElt* f4 = field->newElt();
	FieldElt* f5 = field->newElt();

	EncodedElt* L1 = encoding->new_elt(L);
	EncodedElt* R1 = encoding->new_elt(R);
	EncodedElt* L2 = encoding->new_elt(L);
	EncodedElt* R2 = encoding->new_elt(R);
	EncodedElt* Lsum = encoding->new_elt(L);
	EncodedElt* Rsum = encoding->new_elt(R);

	EncodedProduct* result1 = encoding->new_prod();
	EncodedProduct* result2 = encoding->new_prod();
	EncodedProduct* result3 = encoding->new_prod();

	for (int trial = 0; trial < num_tests; trial++) {
		field->assignRandomElt(f1);
		field->assignRandomElt(f2);
		field->assignRandomElt(f3);
		field->assignRandomElt(f4);

		// Very basic equality test for target group	
		encoding->encode(L, f1, L1);  // L1 = g^f1
		encoding->encode(R, f2, R1);  // R1 = h^f2

		encoding->mul(L1, R1, result1);
		encoding->mul(L1, R1, result2);
		assert(encoding->equals(result1, result2));

		// Even simpler pairing check
		field->one(f5);
		encoding->encode(L, f1, L1);  // L1 = g^f1
		encoding->encode(R, f5, R1);  // R1 = h

		encoding->encode(L, f5, L2);  // L2 = g
		encoding->encode(R, f1, R2);  // R2 = h^f1
			
		encoding->mul(L1, R1, result1);	// result1 = e(g^f1, h)
		encoding->mul(L2, R2, result2); // result2 = e(g, h^f1)	

		assert(encoding->equals(result1, result2));

		// Simple pairing check
		encoding->encode(L, f1, L1);  // L1 = g^f1
		encoding->encode(R, f2, R1);  // R1 = h^f2

		encoding->encode(L, f2, L2);  // L2 = g^f2
		encoding->encode(R, f1, R2);  // R2 = h^f1

		encoding->mul(L1, R1, result1);	// result1 = e(g^f1, h^f2)
		encoding->mul(L2, R2, result2); // result2 = e(g^f2, h^f1)
		assert(encoding->equals(result1, result2));

		// More complex pairing check		
		encoding->encode(L, f1, L1);  // L1 = g^f1
		encoding->encode(R, f2, R1);  // R1 = h^f2

		encoding->mul(L, L1, f3, L2); // L2 = g^f1f3
		encoding->mul(R, R1, f4, R2); // R2 = h^f2f4

		encoding->mul(L2, R2, result1);	// e(g,h)^f1f2f3f4

		field->mul(f3, f4, f5);		
		encoding->mul(L, L1, f5, L2);	  // L2 = g^f1f3f4
		encoding->mul(L2, R1, result2); // e(g,h)^f1f2f3f4

		encoding->mul(R, R1, f5, R2);	  // R2 = h^f2f3f4
		encoding->mul(L1, R2, result3); // e(g,h)^f1f2f3f4

		assert(encoding->equals(result1, result2));
		assert(encoding->equals(result2, result3));

		// One more test of operations in the target group
		field->assignRandomElt(f3);
		field->assignRandomElt(f4);

		encoding->encode(L, f1, L1);  // L1 = g^f1
		encoding->encode(R, f2, R1);  // R1 = h^f2

		encoding->encode(L, f3, L2);  // L2 = g^f3'
		encoding->encode(R, f4, R2);  // R2 = h^f4'

		encoding->add(L, L1, L2, Lsum);	// g^f1+f3'
		encoding->add(R, R1, R2, Rsum); // h^f2+f4'
		encoding->mul(Lsum, Rsum, result1);	// e(g,h)^{(f1+f3')(f2+f4')}

		encoding->mul(L1, R1, result2);	// e(g,h)^{f1f2}
		encoding->mul(L1, R2, result3);	// e(g,h)^{f1f4'}
		encoding->add(result2, result3, result3);
		encoding->mul(L2, R1, result2); // e(g,h)^{f3'f2}
		encoding->add(result2, result3, result3);
		encoding->mul(L2, R2, result2);	// e(g,h)^{f3'f4'}
		encoding->add(result2, result3, result3);		// e(g,h)^{(f1f2 + f1f4' + f3'f2 + f3'f4')

		assert(encoding->equals(result1, result3));
	}

	field->delElt(f1);
	field->delElt(f2);
	field->delElt(f3);
	field->delElt(f4);
	field->delElt(f5);

	encoding->del_elt(L, L1);
	encoding->del_elt(L, L2);
	encoding->del_elt(R, R1);
	encoding->del_elt(R, R2);
	encoding->del_elt(L, Lsum);
	encoding->del_elt(R, Rsum);
	delete result1;
	delete result2;
	delete result3;		
}

void randomEncInnerProductTest(Encoding* encoding, elt_t t, int num_tests) {
	Field* field = encoding->getSrcField();

	FieldEltArray* exp = field->newEltArray(num_tests, true);
	FieldEltArray* randf = field->newEltArray(num_tests, true);	
	EncodedEltArray* bases = encoding->new_elt_array(t, num_tests, true);

	EncodedElt* tmp = encoding->new_elt(t);
	EncodedElt* resultSlow = encoding->new_elt(t);
	EncodedElt* resultFast = encoding->new_elt(t);

	for (int test = 0; test < num_tests; test++) {
		// Choose all new random values
		for (int i = 0; i < num_tests; i++) {
			field->assignRandomElt(exp->elt(i));			
			//field->one(exp->elt(i));
			field->assignRandomElt(randf->elt(i));
			encoding->encode(t, randf->elt(i), bases->elt(i));	// bases[i] = g^{r_i}
		}

		// Fast inner product
    int field_elt_max_bitsize = field->eltSize()*sizeof(digit_t) * 8;
    ip_handle_t h = encoding->prepareForInnerProduct(bases, num_tests, 5, field_elt_max_bitsize, 4, false);
		encoding->innerProduct(h, exp, num_tests, resultFast);
		encoding->doneWithInnerProduct(h);

		// Slow inner product
		encoding->zero(t, resultSlow);
		for (int i = 0; i < num_tests; i++) {
			encoding->mul(t, bases->elt(i), exp->elt(i), tmp);
			encoding->add(t, tmp, resultSlow, resultSlow);
		}

		assert(encoding->equal(t, resultSlow, resultFast));	
	}

	delete exp;
	delete randf;
	delete bases;
	encoding->del_elt(t, tmp);
	encoding->del_elt(t, resultSlow);
	encoding->del_elt(t, resultFast);

}

void runFieldCorrectnessTest(Field* field, int num_tests) {	
	randomFieldAddSubTest(field, num_tests);
	randomFieldAssocTest(field, num_tests);
	randomFieldArrayTest(field, num_tests);
}

void runPolyCorrectnessTest(Field* field, int num_tests, int max_degree) {	
	randomPolyAddSubTest(field, num_tests, max_degree);
	randomPolyMulTest(field, num_tests, max_degree);
	randomPolyAddMulTest(field, num_tests, max_degree);
	randomPolySubMulTest(field, num_tests, max_degree);
	randomPolyModTest(field, num_tests, max_degree);
	randomPolyEvalInterpTest(field, max_degree, -1, num_tests);
	randomPolyEvalInterpTest(field, max_degree, -1, num_tests, true);
}

void runEncCorrectnessTestSided(Encoding* encoding, elt_t t, int num_tests) {	
	randomEncAddSubTest(encoding, t, num_tests);
	randomEncAssocTest(encoding, t, num_tests);
	randomEncInnerProductTest(encoding, t, num_tests);
}

void runEncCorrectnessTest(Encoding* encoding, int num_tests) {	
	runEncCorrectnessTestSided(encoding, L, num_tests);
	runEncCorrectnessTestSided(encoding, R, num_tests);
	randomEncMulTest(encoding, num_tests);
}

void runPrintTest(Encoding* encoding) {
	EncodedElt* eL = encoding->new_elt(L);
	EncodedElt* eR = encoding->new_elt(R);
	EncodedProduct* eP = encoding->new_prod();

	Field* field = encoding->getSrcField();
	FieldElt* f = field->newRandomElt();

	encoding->encode(L, f, eL);
	encoding->encode(R, f, eR);
	encoding->mul(eL, eR, eP);

	printf("\nField elt: ");
	field->print(f);
	printf("\nEncoded L elt: ");
	encoding->print(L, eL);
	printf("\nEncoded R elt: ");
	encoding->print(R, eR);
	printf("\nEncoded Product: ");
	encoding->print(eP);
	printf("\n");

	encoding->del_elt(L, eL);
	encoding->del_elt(R, eR);
	delete eP;
	field->delElt(f);
}



void Write(Archive* arc, Archiver* writer, Encoding* encoding, EncodedElt* eL, EncodedElt* eR, EncodedProduct* eP) {
  encoding->write(writer, L, eL, false);
  encoding->write(writer, R, eR, false);
  writer->write((uint8_t*)eP, encoding->size_of_prod());
}

void WriteTest(Encoding* encoding, int num_trials) {
  FileArchiver file_archiver("test_file.dat", FileArchiver::Write);
  Archive arc(&file_archiver, encoding);

  Field* field = encoding->getSrcField();
  FieldElt* f1 = field->newRandomElt();
  FieldElt* f2 = field->newRandomElt();
  EncodedElt* eL = encoding->new_elt(L);
  EncodedElt* eR = encoding->new_elt(R);
  EncodedProduct* eP = encoding->new_prod();

  // Write out some known values
  encoding->zero(L, eL);
  encoding->zero(R, eR);
  encoding->mul(eL, eR, eP);
  Write(&arc, &file_archiver, encoding, eL, eR, eP);
  encoding->one(L, eL);
  encoding->one(R, eR);
  encoding->mul(eL, eR, eP);
  Write(&arc, &file_archiver, encoding, eL, eR, eP);

  // Write out some random values
  for (int i = 0; i < num_trials; i++) {
    field->assignRandomElt(f1);
    field->assignRandomElt(f2);
    encoding->encode(L, f1, eL);
    encoding->encode(R, f2, eR);
    encoding->mul(eL, eR, eP);
    Write(&arc, &file_archiver, encoding, eL, eR, eP);
  }

  field->delElt(f1);
  field->delElt(f2);
  encoding->del_elt(L, eL);
  encoding->del_elt(R, eR);
  delete eP;
}


void Read(Archive* arc, Archiver* reader, Encoding* encoding, EncodedElt* eL, EncodedElt* eR, EncodedProduct* eP) {
  encoding->read(reader, L, eL);
  encoding->read(reader, R, eR);
  reader->read((uint8_t*)eP, encoding->size_of_prod());
}

void ReadTest(Encoding* encoding, int num_trials) {
  FileArchiver file_archiver("test_file.dat", FileArchiver::Read);
  Archive arc(&file_archiver, encoding);

  EncodedElt* eL = encoding->new_elt(L);
  EncodedElt* eR = encoding->new_elt(R);
  EncodedProduct* eP = encoding->new_prod();

  EncodedElt* tL = encoding->new_elt(L);
  EncodedElt* tR = encoding->new_elt(R);
  EncodedProduct* tP = encoding->new_prod();

  // Read in some known values
  Read(&arc, &file_archiver, encoding, eL, eR, eP);
  encoding->zero(L, tL);
  encoding->zero(R, tR);
  encoding->mul(tL, tR, tP);
  assert(encoding->equal(L, tL, eL));
  assert(encoding->equal(R, tR, eR));
  assert(encoding->equals(tP, eP));

  Read(&arc, &file_archiver, encoding, eL, eR, eP);
  encoding->one(L, tL);
  encoding->one(R, tR);
  encoding->mul(tL, tR, tP);
  assert(encoding->equal(L, tL, eL));
  assert(encoding->equal(R, tR, eR));
  assert(encoding->equals(tP, eP));

  // Read in some random values.  Can't check actual value,
  // but we can check that multiplication still works
  for (int i = 0; i < num_trials; i++) {
    Read(&arc, &file_archiver, encoding, eL, eR, eP);
    encoding->mul(eL, eR, tP);
    assert(encoding->equals(tP, eP));
  }

  encoding->del_elt(L, eL);
  encoding->del_elt(R, eR);
  delete eP;
  encoding->del_elt(L, tL);
  encoding->del_elt(R, tR);
  delete tP;
}

void runFileTest(Encoding* encoding, int num_trials) {
  WriteTest(encoding, num_trials);
  ReadTest(encoding, num_trials);
}

/*
void runCompressionTest(variables_map& config, Encoding* encoding, int num_trials) {
  if (!(config["encoding"].as<string>().compare("enigma-bn") == 0 || config["encoding"].as<string>().compare("arith-cp") == 0)) {
    printf("Compression not enabled for this encoding.  Skipping compression test\n");
    return;
  }

  Field* field = encoding->getSrcField();
  FieldElt* f1 = field->newElt();
  FieldElt* f2 = field->newElt(); 
  EncodedElt* eL = encoding->new_elt(L);
  EncodedElt* eR = encoding->new_elt(R);

  EncodedElt* eLc = encoding->new_elt(L);
  EncodedElt* eRc = encoding->new_elt(R);

  EncodedElt* eLd = encoding->new_elt(L);
  EncodedElt* eRd = encoding->new_elt(R);

  for (int i = 0; i < num_trials; i++) {
    field->assignRandomElt(f1);
    field->assignRandomElt(f2);

    encoding->encode(L, f1, eL);
    encoding->encode(R, f2, eR);

    encoding->copy(L, eL, eLc);
    encoding->copy(R, eR, eRc);

    encoding->compress(L, eLc);
    encoding->compress(R, eRc);

    if (config["encoding"].as<string>().compare("enigma-bn") == 0) {
      EncodedEltEnigmaBN* bnL = (EncodedEltEnigmaBN*)eLc;
      EncodedEltEnigmaBN* bnR = (EncodedEltEnigmaBN*)eRc;
      encoding->decompress(L, bnL->get_digits(), eLd);
      encoding->decompress(R, bnR->get_digits(), eRd);
    } else if (config["encoding"].as<string>().compare("arith-cp") == 0) {
      //EncodedEltArithCP* cpL = (EncodedEltArithCP*)eLc;
      //EncodedEltArithCP* cpR = (EncodedEltArithCP*)eRc;
      //encoding->decompress(L, cpL->me->X, eLd);
      //encoding->decompress(R, cpR->me->X, eRd);
      encoding->decompress(L, ((point_wp_t*) eLc)->X->limbs, eLd);
      encoding->decompress(R, ((point_wp2_t*)eRc)->X->limbs, eRd);
    }

    assert(encoding->equal(L, eL, eLd));
    assert(encoding->equal(R, eR, eRd));
  }  

}
*/

void functionality_tests(variables_map& config) {
	Encoding* encoding = NewEncoding(config["encoding"].as<string>());
	if (encoding == NULL) { return; }
	Field* field = encoding->getSrcField();
  
	int max_degree = 200;
	int num_tests = config["trials"].as<int>();

	runFieldCorrectnessTest(field, num_tests);
	runPolyCorrectnessTest(field, num_tests, max_degree);
	runEncCorrectnessTest(encoding, num_tests);
	runPrintTest(encoding);
  //runFileTest(encoding, num_tests);

	encoding->prepareForManyEncodings(max_degree*10, 4, false);
	runEncCorrectnessTest(encoding, num_tests);  // Rerun with faster encodings
	encoding->doneWithManyEncodings();

	delete encoding;
}

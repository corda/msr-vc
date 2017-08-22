#include <iostream>
#include <cassert>
#include <time.h>
#include "timer.h"
#include "Encoding.h"
#include "PolyTree.h"
#include "EncodingEnigmaBN.h"

#include <boost/program_options.hpp>
using namespace std;
using namespace boost::program_options;

void microbenchmarks(variables_map& config) {

#ifdef _DEBUG
	printf("\nRunning in DEBUG mode -- performance measurements should not be trusted!\n\n");
#endif

	string option = config["micro"].as<string>();
	int num_trials = config["trials"].as<int>();
	printf("About to run %d tests...\n", num_trials);

	Encoding* encoding = NewEncoding(config["encoding"].as<string>());
	if (encoding == NULL) { return; }
	Field* field = encoding->getSrcField();
  if (option.compare("mont-test") == 0) {
    assert(config["encoding"].as<string>().compare("enigma-bn") == 0);
    field = ((EncodingEnigmaBN*)encoding)->getMontField();

    ((EncodingEnigmaBN*)encoding)->testToFromMont(field, num_trials);
  }

	FieldEltArray* field1 = field->newEltArray(num_trials, true);
	FieldEltArray* field2 = field->newEltArray(num_trials, true);
	FieldEltArray* field32 = field->newEltArray(num_trials, true);

	EncodedEltArray* lencodings = encoding->new_elt_array(L, num_trials, true);
	EncodedEltArray* rencodings = encoding->new_elt_array(R, num_trials, true);
	EncodedEltArray* lencodings2 = encoding->new_elt_array(L, num_trials, true);
	EncodedEltArray* rencodings2 = encoding->new_elt_array(R, num_trials, true);
	EncodedProduct** prods = new EncodedProduct*[num_trials];

	FieldElt *fieldg1 = field->newElt(), *fieldg0 = field->newElt();
	// I'll assume 5^i are unique values, for i in [0,num_trials)
	field->set(fieldg1, 5);
	field->set(fieldg0, 1);
	for (int i = 0; i < num_trials; i++) {
		field->assignRandomElt(field1->elt(i));
		field->assignRandomElt(field2->elt(i));
		field->assignRandomElt(field32->elt(i));

    field->truncate(field1->elt(i), 254);
    field->truncate(field2->elt(i), 254);
		field->truncate(field32->elt(i), 32);   

		prods[i] = encoding->new_prod();
	}

	Timer* total = timers.newTimer("Total", NULL);
	total->start();

	/////////////// Field ops ///////////////////////////
  if (option.compare("all") == 0 || option.compare("field") == 0 || option.compare("mont-test") == 0) {
		TIMEREP(field->add(field1->elt(i), field2->elt(i), field1->elt(i)), num_trials, "FieldAdd", "Total");
		TIMEREP(field->sub(field1->elt(i), field2->elt(i), field1->elt(i)), num_trials, "FieldSub", "Total");
		TIMEREP(field->mul(field1->elt(i), field2->elt(i), field1->elt(i)), num_trials, "FieldMul", "Total");
		TIMEREP(field->mul(field1->elt(i), field2->elt(i), field1->elt(i)), num_trials, "FieldDiv", "Total");	

		printf("Finished field ops\n");  fflush(stdout);
	}

	/////////////// Encoding ops ///////////////////////////
	if (option.compare("all") == 0 || option.compare("encoding") == 0) {
		TIMEREP(encoding->encode(L, field1->elt(i), lencodings->elt(i)), num_trials, "EncodeSlowL", "Total");
		TIMEREP(encoding->encode(R, field1->elt(i), rencodings->elt(i)), num_trials, "EncodeSlowR", "Total");

		TIME(encoding->prepareForManyEncodings(num_trials, config["mem"].as<int>(), false), "PrepForEnc", "Total");																											 
		TIMEREP(encoding->encode(L, field1->elt(i), lencodings->elt(i)), num_trials, "EncodeFastL", "Total");
		TIMEREP(encoding->encode(R, field1->elt(i), rencodings->elt(i)), num_trials, "EncodeFastR", "Total");
		TIMEREP(encoding->encode(L, field2->elt(i), lencodings2->elt(i)), num_trials, "ignoreL", "Total");
		TIMEREP(encoding->encode(R, field2->elt(i), rencodings2->elt(i)), num_trials, "ignoreR", "Total");
		encoding->doneWithManyEncodings();

		printf("Finished encoding test values\n");  fflush(stdout);

		TIMEREP(encoding->add(L, lencodings->elt(i), lencodings2->elt(i), lencodings2->elt(i)), num_trials, "EncAddL", "Total");
		TIMEREP(encoding->add(R, rencodings->elt(i), rencodings2->elt(i), rencodings2->elt(i)), num_trials, "EncAddR", "Total");

		ip_handle_t hL, hR;
    int field_elt_max_bitsize = field->eltSize()*sizeof(digit_t) * 8;
    TIME(hL = encoding->prepareForInnerProduct(lencodings, num_trials, 1, field_elt_max_bitsize, config["mem"].as<int>(), false), "PrepInnerL", "Total");
		TIME(encoding->innerProduct(hL, field1, num_trials, lencodings2->elt(0)), "InnerProdL", "Total");
		encoding->doneWithInnerProduct(hL);

    TIME(hR = encoding->prepareForInnerProduct(rencodings, num_trials, 1, field_elt_max_bitsize, config["mem"].as<int>(), false), "PrepInnerR", "Total");
		TIME(encoding->innerProduct(hR, field1, num_trials, rencodings2->elt(0)), "InnerProdR", "Total");
		encoding->doneWithInnerProduct(hR);

		TIMEREP(encoding->mul(lencodings->elt(i), rencodings->elt(i), prods[i]), num_trials, "Pairing", "Total");

		printf("Finished large encoding tests\n");  fflush(stdout);

		// Measure the time for smaller field elements
		TIME(encoding->prepareForManyEncodings(num_trials, config["mem"].as<int>(), false), "PrepForEnc", "Total");																											 
		TIMEREP(encoding->encode(L, field32->elt(i), lencodings2->elt(i)), num_trials, "EncodeFast32L", "Total");
		TIMEREP(encoding->encode(R, field32->elt(i), rencodings2->elt(i)), num_trials, "EncodeFast32R", "Total");	
		encoding->doneWithManyEncodings();

		TIME(hL = encoding->prepareForInnerProduct(lencodings, num_trials, 1, 32, config["mem"].as<int>(), false), "PrepInnerL", "Total");
		TIME(encoding->innerProduct(hL, field32, num_trials, lencodings2->elt(0)), "InnerProd32L", "Total");
		encoding->doneWithInnerProduct(hL);

		TIME(hR = encoding->prepareForInnerProduct(rencodings, num_trials, 1, 32, config["mem"].as<int>(), false), "PrepInnerR", "Total");
		TIME(encoding->innerProduct(hR, field32, num_trials, rencodings2->elt(0)), "InnerProd32R", "Total");
		encoding->doneWithInnerProduct(hR);

		printf("Finished 32-bit encoding tests\n");  fflush(stdout);
	}
	
	/////////////// Poly Ops ///////////////////////////
  if (option.compare("all") == 0 || option.compare("poly-all") == 0 || option.compare("poly-fast") == 0 || option.compare("poly-slow") == 0 || option.compare("mont-test") == 0) {
		Timer* polyPreComp = timers.newTimer("PolyPreComp", "Total");
		polyPreComp->start();
		PolyTree* tree = new PolyTree(field, field1, num_trials);
		FieldEltArray* denominators = Poly::genLagrangeDenominators(field, *tree->polys[tree->height-1][0], tree, field1, num_trials);
		polyPreComp->stop();

		TIME(Poly::interpolate(field, field2, num_trials, tree, denominators), "PolyInterp", "Total");
		TIME(Poly::interpolateGeometric(field, fieldg0, fieldg1, field2, num_trials), "PolyInterpGeom", "Total");
		
		// Compare slow vs. optimized polynomial multiplication
		Poly p1(field), p2(field), r(field);
		Poly::polyRand(field, p1, num_trials);
		Poly::polyRand(field, p2, num_trials);
		if (option.compare("poly-all") == 0 || option.compare("poly-fast") == 0) {
			TIME(Poly::mul(p1, p2, p1), "PolyMul", "Total");
		}
		if (option.compare("poly-all") == 0 || option.compare("poly-slow") == 0) {				 
			TIME(Poly::mulSlow(p1, p2, p1), "PolyMulSlow", "Total");
		}
	}

	total->stop();

	if (config.count("raw")) {
		timers.printRaw();
	} else {
		timers.printStats();
	}
	printf("\n");

	delete field1, field2, field32;
	delete lencodings,  rencodings;
	delete lencodings2, rencodings2;
	delete [] prods;

	delete encoding;
}
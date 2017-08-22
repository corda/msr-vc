#pragma once

extern "C" {
#include "Fp/Fp.h"
#include "Fp/Fp.h"
#include "Fp2/Fp2.h"
#include "Fp4/Fp4.h"
#include "Fp6/Fp6.h"
#include "Fp6/Fp6overFp3.h"
#include "Fp12/Fp12.h"
#include "Fp12/Fp12overFp4.h"
#include "Curve/Fp_weierstrass.h"
#include "Curve/Fp2_weierstrass.h"
#include "pairing/pairing.h"
#include "uint.h"
#include "arith.h"
}

#include "Encoding.h"


class LEncodedEltArithBN : public EncodedEltArithBN { 
public: 
	static LEncodedEltArithBN* newElt();
	Point_wp_t me;
};

class REncodedEltArithBN : public EncodedEltArithBN { 
public: 
	static LEncodedEltArithBN* newElt();
	Point_wp2_t me;
};

class EncodedProductArithBN : public EncodedProduct { 
public: 
	Fp12_t me;
  ~EncodedProductArithBN();
};

class EncodingArithBNMembers {
public:
	uint64_t *m, *b0, *b1;
	Fp_t b;
	int n;
	uint64_t start_mil, end_mil, start_fe, end_fe;
	Fp2_t B, B3, l00, l01, l10;
	Encoding* enigma_encoding;

  // Projective versions of the Base and Twist generators (also a version of g^1)
  LEncodedEltArithBN Lg;
  REncodedEltArithBN Rg;

  
	// Frequently encoded elts
	LEncodedEltArithBN Lzero;
  REncodedEltArithBN Rzero;
};



// Both the base and twist curves are defined over Fp,
// so we don't need different representations or operations
class EncodedEltArithCP : public EncodedElt {
public:
  Point_wp_t me;
};

class EncodedProductArithCP : public EncodedProduct { 
public:
  bool is3;
  union Prod {
    Fp3_t  me3;  // Used for CP3
    Fp6q_t me6;  // Used for CP6 and DEBUG
  };

  Prod me;
  ~EncodedProductArithCP();
};

class EncodingArithCPMembers {
public:
  Fp_t B[2], B3;
  int n;

	Encoding* enigma_encoding;

  // Projective versions of the generator (also a version of g^1)
  EncodedEltArithCP Lg;
  EncodedEltArithCP Rg;
    
	// Frequently encoded elt
	EncodedEltArithCP Lzero;  
  EncodedEltArithCP Rzero;  
};

void write_Fp(Archiver* arc, Fp_t* fp, int num_digits);

void read_Fp(Archiver* arc, Fp_t* fp, int num_digits);

void write_Fp2(Archiver* arc, Fp2_t* fp2, int num_digits);

void read_Fp2(Archiver* arc, Fp2_t* fp2, int num_digits);
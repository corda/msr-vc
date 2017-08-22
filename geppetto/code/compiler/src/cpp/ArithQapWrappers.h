#pragma once

extern "C" {
#include "QapFp.h"
#include "QapFp.h"
#include "QapFp2.h"
#include "QapFp4.h"
#include "QapFp6.h"
//#include "QapFp6overFp3.h"
#include "QapFp12.h"
#include "QapFp_weierstrass.h"
#include "QapFp2_weierstrass.h"
#include "Qappairing.h"
#include "Qapuint.h"
#include "QapInterface.h"
}

#include "Encoding.h"


class LEncodedEltQapArith : public EncodedEltQapArith { 
public: 
	static LEncodedEltQapArith* newElt();
	QapPoint_wp_t me;
};

class REncodedEltQapArith : public EncodedEltQapArith { 
public: 
	static LEncodedEltQapArith* newElt();
	QapPoint_wp2_t me;
};

class EncodedProductQapArith : public EncodedProduct { 
public: 
	QapFp12_t me;
};

class EncodingArithQapMembers {
public:
	QapInterface qap_interface;
	Encoding* msft_encoding;

	// Projective versions of the Base and Twist generators (also a version of g^1)
	LEncodedEltQapArith Lg;
	REncodedEltQapArith Rg;

	// Frequently encoded elts
	LEncodedEltQapArith Lzero;
	REncodedEltQapArith Rzero;
};



// Both the base and twist curves are defined over Fp,
// so we don't need different representations or operations
class EncodedEltArithQap : public EncodedElt {
public:
  QapPoint_wp_t me;
};

#if 0
class EncodedProductArithQap : public EncodedProduct { 
public: 
	QapFp6q_t me;
};
#endif
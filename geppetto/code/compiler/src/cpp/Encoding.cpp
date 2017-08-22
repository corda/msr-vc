#include <assert.h>
#include <iostream>
#include "Encoding.h"
#include "EncodingDbg.h"
#include "EncodingEnigmaBN.h"
#include "EncodingArithBN.h"
#include "EncodingArithCP.h"
#include "EncodingArithQAP.h"

Encoding* Encoding::current_encoding;

const char* encoding_names[] = { "debug", "debug32", "enigma-bn", "arith-bn", "arith-cp", "arith-qap", "arith-bn-dbg", "arith-cp-dbg", "arith-cp3" };

encoding_t THE_ENCODING;  // Used for debugging purposes

Encoding* NewEncoding(encoding_t type) {
  THE_ENCODING = type;
	switch(type) {
		case ENCODING_DBG_PRIME:
			return new EncodingDbg(EncodingDbg::field_prime);
		case ENCODING_DBG_32:
			return new EncodingDbg(EncodingDbg::field32);
		case ENCODING_ENIGMA_BN:
			return new EncodingEnigmaBN();
    case ENCODING_ARITH_BN:
			return new EncodingArithBN(false);
    case ENCODING_ARITH_CP:
			return new EncodingArithCP(EncodingArithCP::CP_curve_type::CP_6);
	  case ENCODING_ARITH_QAP:
			return new EncodingArithQap();
    case ENCODING_ARITH_CP_DBG:
      return new EncodingArithCP(EncodingArithCP::CP_curve_type::DEBUG);
    case ENCODING_ARITH_BN_DBG:
      return new EncodingArithBN(true);
    case ENCODING_ARITH_CP3:
      return new EncodingArithCP(EncodingArithCP::CP_curve_type::CP_3);
		default:
			assert(0);
			return NULL;
	}  
}

Encoding* NewEncoding(string name) {
	for (int i = 0; i < sizeof(encoding_names)/sizeof(char*); i++) {
		if (name == encoding_names[i]) {
			return NewEncoding((encoding_t)(i+1));		// enum starts at 1
		}
	}
	cout << "Error: Unrecognized encoding name: " << name << endl;
	return NULL;
}

void Encoding::print_stats()
{
  printf("No stats collected\n");
}

void Encoding::collect_stats(bool)
{
}

void Encoding::reset_stats()
{
}

EncodedProduct::~EncodedProduct() { }

EncodedEltArray::EncodedEltArray() {
}

EncodedEltArray::~EncodedEltArray() {
}

#pragma once

#include "Types.h"
#include "Field.h"
#include "Archiver.h"

typedef uint32_t ip_handle_t;
enum elt_t {L = 0, R};

class EncodedElt { 
protected:
	EncodedElt() { ; }  // Ensure this ifc class can't be instantiated
};

class EncodedProduct { 
public:
	virtual ~EncodedProduct();
protected:
	EncodedProduct() { ; }  // Ensure this ifc class can't be instantiated
};

class CondensedEncElt { 
protected:
	CondensedEncElt();  // Ensure this ifc class can't be instantiated
};

class Encoding;

class EncodedEltArray { 
public:
	virtual ~EncodedEltArray();

	virtual EncodedElt* elt(uint32_t index) = 0;
	virtual void set(uint32_t index, EncodedElt* elt) = 0;
	elt_t t()    { return type; }
  int length() { return len; }
  bool pre()   { return allocated; }
  virtual void reset() = 0;
protected:
	EncodedEltArray();
	elt_t    type;	
  uint32_t len;
  bool     allocated;  // Did we allocate the elts?
};

typedef EncodedElt LEncodedElt;
typedef EncodedElt REncodedElt;
typedef EncodedEltArray LEncodedEltArray;
typedef EncodedEltArray REncodedEltArray;

class Archive;

class Encoding {
public:
	//Encoding();
	virtual ~Encoding() { ; }
	
	virtual Field* getSrcField() = 0;	// Encoding takes elements from srcField to dstField

	virtual EncodedElt* new_elt(elt_t t) = 0;
	virtual void del_elt(elt_t t, EncodedElt* elt) = 0;
	virtual EncodedProduct* new_prod() = 0;
  virtual uint32_t size_of_prod() = 0;

	virtual void print(elt_t t, EncodedElt* e) = 0;
	virtual void print(EncodedProduct* p) = 0;

	virtual EncodedEltArray* new_elt_array(elt_t t, uint32_t len, bool preallocate=true) = 0;

	virtual void encode(elt_t t, FieldElt* in, EncodedElt* out) = 0;

	virtual void copy(elt_t t, EncodedElt* src, EncodedElt* dst) = 0;

	virtual bool equal(elt_t t, EncodedElt* a, EncodedElt* b) = 0;

	// a <- 0
	virtual void zero(elt_t t, EncodedElt* a) = 0;

	// a <- 1
	virtual void one(elt_t t, EncodedElt* a) = 0;

	// Adds two encoded elements.  Result goes in r.
	virtual void add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) = 0;

	// Subtracts two encoded elements.  Result goes in r.
	virtual void sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) = 0;
	
	// b = 2 * a
	virtual void doubleIt(elt_t t, EncodedElt* a, EncodedElt* b) = 0;

	// Multiplies the value encoded in g by the constant c.  Result goes in r.
	virtual void mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r) = 0;

	// Mul two encoded values
	virtual void mul(EncodedElt* a, EncodedElt* b, EncodedProduct* r) = 0;

	// Arithmetic on the result of multiplying two encoded values
	virtual void add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) = 0;
	virtual void sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) = 0;
	virtual bool equals(EncodedProduct* a, EncodedProduct* b) = 0;

	// Checks whether the plaintexts inside the encodings obey:
	// L1*R1 - L2*R2 == L3*R3
	virtual bool mulSubEquals(EncodedElt* L1, EncodedElt* R1, EncodedElt* L2, EncodedElt* R2, EncodedElt* L3, EncodedElt* R3) = 0;

	// Checks whether the plaintexts inside the encodings obey:
	// L1*R1 == L2*R2
	virtual bool mulEquals(EncodedElt* L1, EncodedElt* R1, EncodedElt* L2, EncodedElt* R2) = 0;

	// Control caching behavior
	virtual void prepareForManyEncodings(uint64_t numEncodings, uint32_t maxMemGb, bool prepIsFree) = 0;
	virtual void doneWithManyEncodings() = 0;

  virtual ip_handle_t prepareForInnerProduct(EncodedEltArray* bases, int len, int numUses, int max_exponent_size, int maxMemGb, bool precompFree) = 0;
  virtual void doneWithInnerProduct(ip_handle_t handle) = 0;
  virtual void innerProduct(ip_handle_t handle, FieldEltArray* exp, int len, EncodedElt* r) = 0;

	virtual void write(Archiver* arc, elt_t t, EncodedElt* e, bool simple) = 0;
	virtual void read (Archiver* arc, elt_t t, EncodedElt* e) = 0;

  // Internal methods exposed only for perf testing
  virtual void compress(elt_t t, EncodedElt* elt) = 0;
  virtual void decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt) = 0;

  virtual void print_stats();
  virtual void collect_stats(bool on);
  virtual void reset_stats();

  // Used for debugging purposes only
  static Encoding* current_encoding;
};

enum encoding_t { ENCODING_DBG_PRIME = 1, ENCODING_DBG_32, ENCODING_ENIGMA_BN, ENCODING_ARITH_BN, ENCODING_ARITH_CP, ENCODING_ARITH_QAP, ENCODING_ARITH_BN_DBG, ENCODING_ARITH_CP_DBG, ENCODING_ARITH_CP3 };

Encoding* NewEncoding(encoding_t type);
Encoding* NewEncoding(string name);

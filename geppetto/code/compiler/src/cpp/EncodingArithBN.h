#pragma once

/* 
 * Implements an encoding that maps values to powers of a generator of a BN ECC,
 * built using the new ARITH crypto library
 */

#include "Types.h"
#include "Field.h"
#include "EncodingExp.h"

class EncodingArithBNMembers;

class LEncodedEltArithBN;
class REncodedEltArithBN;

class EncodedEltArithBN : public EncodedElt {

};


class EncodedProductArithBN;

class CondensedEncEltArithBN : public CondensedEncElt { 
};

// For fields with elements that consist of enigma digit_t
class EncodedEltArithBNArray : public EncodedEltArray {
public:
	// If allocate=true, we take responsibility for all memory needs
	// If not, it's assumed someone will point us at the memory needed for the FieldElts
	EncodedEltArithBNArray(elt_t t, Encoding* encoding, uint32_t len, bool allocate);
	virtual ~EncodedEltArithBNArray();

  void reset();
	inline EncodedElt* elt(uint32_t index);
	void set(uint32_t index, EncodedElt* elt);
private:
	EncodedEltArithBN** elts;  
  Encoding* encoding;
};


class EncodingArithBN : public EncodingExp {
public:
  EncodingArithBN(bool debug);
	~EncodingArithBN();	

	virtual EncodedElt* new_elt(elt_t t);
	virtual void del_elt(elt_t t, EncodedElt* elt);
	virtual EncodedProduct* new_prod();
  virtual uint32_t size_of_prod();

	virtual void print(elt_t t, EncodedElt* e);
	virtual void print(EncodedProduct* p);

	virtual EncodedEltArray* new_elt_array(elt_t t, uint32_t len, bool preallocate=true);

	//virtual void encode(elt_t t, FieldElt* in, EncodedElt* out);

	virtual void copy(elt_t t, EncodedElt* src, EncodedElt* dst);

	virtual bool equal(elt_t t, EncodedElt* a, EncodedElt* b);

	// a <- 0
	virtual void zero(elt_t t, EncodedElt* a);

	// a <- 1
	virtual void one(elt_t t, EncodedElt* a);

	// Adds two encoded elements.  Result goes in r.
  virtual void add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r);

	// Subtracts two encoded elements.  Result goes in r.
	virtual void sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r);
	
	// b = 2 * a
	virtual void doubleIt(elt_t t, EncodedElt* a, EncodedElt* b);

	// Multiplies the value encoded in g by the constant c.  Result goes in r.
	virtual void mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r);

	// Mul two encoded values
	virtual void mul(EncodedElt* a, EncodedElt* b, EncodedProduct* r);

	// Arithmetic on the result of multiplying two encoded values
	virtual void add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r);
	virtual void sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r);
	virtual bool equals(EncodedProduct* a, EncodedProduct* b);

	// Checks whether the plaintexts inside the encodings obey:
	// L1*R1 - L2*R2 == L3*R3
	virtual bool mulSubEquals(EncodedElt* L1, EncodedElt* R1, EncodedElt* L2, EncodedElt* R2, EncodedElt* L3, EncodedElt* R3);

	// Checks whether the plaintexts inside the encodings obey:
	// L1*R1 == L2*R2
	virtual bool mulEquals(EncodedElt* L1, EncodedElt* R1, EncodedElt* L2, EncodedElt* R2);
	
  virtual void write(Archiver* arc, elt_t t, EncodedElt* e, bool simple);
	virtual void read (Archiver* arc, elt_t t, EncodedElt* e);

  // Internal methods exposed only for perf testing
  virtual void compress(elt_t t, EncodedElt* elt);
  virtual void decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt);

	// ECC specific-functionality
	void testOnCurve(elt_t t, EncodedElt* a);
	bool isZero(elt_t t, EncodedElt* elt);
	void pair(LEncodedElt* L1, REncodedElt* R1, EncodedProduct* result);
	void affinize(elt_t t, EncodedElt* in, EncodedElt* out);
  void projectivize(elt_t t, EncodedElt* elt);

	virtual void encodeSlow(elt_t t, FieldElt* in, EncodedElt* out);

	virtual uint32_t bytes_per(elt_t t);
	virtual uint32_t num_exponent_bits();

	enum addsub_t { Add, Sub };
	void addsub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, addsub_t op);

	static const uint32_t num_L_digits = 3 * 4;  // Three field elements (4 digits each) to define a projective point on the base curve
	static const uint32_t num_R_digits = 3*2*4;  // Three field elements (each 2 field elements of 4 digits each) to define a projective point on the twist curve

protected:		
	EncodingArithBNMembers* members;
  bool debug;

  // Adds two encoded elements.  Result goes in r.  Set definitely distinct if you know a!=b
	virtual void add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, bool definitely_distinct);

};

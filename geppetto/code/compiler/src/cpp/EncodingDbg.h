#pragma once

/* 
 * Implements an (insecure) encoding that maps values to Z_p
 * for some small p (though it needs to be large enough to accomodate
 * degree(QAP) different roots
 */

#define DEBUG_INHERITS_FROM_EXP

#include "Types.h"
#include "Field.h"
#include "Encoding.h"
#include "EncodedDigitArray.h"
#ifdef DEBUG_INHERITS_FROM_EXP
#include "EncodingExp.h"
#endif

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

class EncodedEltDbg : public EncodedElt { 
public: 
	//EncodedEltDbg(elt_t t) { f = NULL; this->t = t; }
	EncodedEltDbg(elt_t t, FieldElt* f_in) { f = f_in; this->t = t; }
	~EncodedEltDbg() { if (f) { delete f; f = NULL; } }

	inline FieldElt* get_f() { return f; }
	inline void set_f(FieldElt* f_in) { f = f_in; }
	inline elt_t get_t() { return t; }
private:
	FieldElt* f;
	elt_t t;
};


class EncodedProductDbg : public EncodedProduct { 
public: 
	EncodedProductDbg(Field* field) { f = field->newElt(); }
	~EncodedProductDbg() { if (f) { delete f; f = NULL; } }
	FieldElt* f;
};

class CondensedEncEltDbg : public CondensedEncElt { 
};

// Simple array implementation
class EncodedEltDbgArray : public EncodedEltArray { 
public:
	EncodedEltDbgArray(elt_t type, Encoding* encoding, uint32_t len, bool preallocate=true); // Set preallocate to false to start with only pointers
	~EncodedEltDbgArray();

  void reset();
	virtual EncodedElt* elt(uint32_t index);
	virtual void set(uint32_t index, EncodedElt* elt);
	//elt_t t() { return type; }
protected:
	//elt_t type;
	EncodedEltDbg** elts;
};


#ifdef DEBUG_INHERITS_FROM_EXP
class EncodingDbg : public EncodingExp {
#else // !DEBUG_INHERITS_FROM_EXP
class EncodingDbg : public Encoding {
#endif // DEBUG_INHERITS_FROM_EXP

public:
	static enum field_type { field_prime, field32 };

	EncodingDbg(field_type type);
	~EncodingDbg();
	
	virtual Field* getSrcField();	// Encoding takes elements from srcField to dstField

	virtual EncodedElt* new_elt(elt_t t);
	virtual void del_elt(elt_t t, EncodedElt* elt);
	virtual EncodedProduct* new_prod();
  virtual uint32_t size_of_prod();

	virtual void print(elt_t t, EncodedElt* e);
	virtual void print(EncodedProduct* p);

	virtual EncodedEltArray* new_elt_array(elt_t t, uint32_t len, bool preallocate=true);

#ifdef DEBUG_INHERITS_FROM_EXP
	virtual void encodeSlow(elt_t t, FieldElt* in, EncodedElt* out);
#else // !DEBUG_INHERITS_FROM_EXP
	virtual void encode(elt_t t, FieldElt* in, EncodedElt* out);
#endif //DEBUG_INHERITS_FROM_EXP

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

	// Control caching behavior
#ifndef DEBUG_INHERITS_FROM_EXP
	virtual void prepareForManyEncodings(uint64_t numEncodings, uint32_t maxMemGb, bool prepIsFree);
	virtual void doneWithManyEncodings();

	virtual ip_handle_t prepareForInnerProduct(EncodedEltArray* bases, int len, int numUses, int maxMemGb, bool precompFree);
	virtual void doneWithInnerProduct(ip_handle_t handle);
	virtual void innerProduct(ip_handle_t handle, FieldEltArray* exp, int len, EncodedElt* r);	
#endif // DEBUG_INHERITS_FROM_EXP

  virtual void write(Archiver* arc, elt_t t, EncodedElt* e, bool simple);
	virtual void read (Archiver* arc, elt_t t, EncodedElt* e);

  // Internal methods exposed only for perf testing
  virtual void compress(elt_t t, EncodedElt* elt);
  virtual void decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt);

#ifdef DEBUG_INHERITS_FROM_EXP
	virtual uint32_t bytes_per(elt_t t);
	virtual uint32_t num_exponent_bits();
#endif

protected:		
	mp_modulus_t* p;	

#ifndef DEBUG_INHERITS_FROM_EXP	// EncodingExp provides a srcField
	Field* srcField;
#endif

	bigctx_t BignumCtx;
	bigctx_t* pbigctx;
	
	field_desc_t* msExponentField;

	ip_handle_t h_ctr;
	EncodedEltArray* bases[2];

	//int num_exponent_bits;
};





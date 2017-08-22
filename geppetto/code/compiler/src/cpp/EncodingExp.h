#pragma once

/* 
 * Implements functionality common to all encodings that take x and return g^x
 * All encoding-specific operations are left to the children
 */

#include "Types.h"
#include "Field.h"
#include "Encoding.h"

#include <unordered_map>
#include <hash_map>	
using namespace stdext;

typedef struct MultiExpPrecompTable_s {
	int numBits;  // Number of bits to take from each exponent
	int numTerms; // Number of bases to combine into a table.  Currently, code can only handle values of 1 and 2.
	int len;			// Number of bases this table covers
	int numTables;
	uint64_t sizeInMB; // Memory used	
	EncodedEltArray** precompTables;
	elt_t t;
} MultiExpPrecompTable;

typedef hash_map<uint32_t,MultiExpPrecompTable*> MultiExpTables;
typedef pair<uint32_t,MultiExpPrecompTable*> MultiExpTable;

struct fequal : std::binary_function<FieldElt*, FieldElt*, bool>
{
  static Field* field;
  bool operator()(FieldElt* const x, FieldElt* const y) const {
    return field->equal(x, y);
  }
};

struct fhash : std::unary_function<FieldElt*, std::size_t>
{
  static Field* field;
  std::size_t operator()(FieldElt* const x) const {
    return field->hash(x);
  }
};

typedef unordered_map<FieldElt*, EncodedElt*, fhash, fequal> EncodedCache;
typedef pair<FieldElt*, EncodedElt*> CachedEncoding;



class EncodingExp : public Encoding {
public:
	EncodingExp();
	virtual ~EncodingExp();
	
	virtual Field* getSrcField();	// Encoding takes an element x from srcField to g^x

	virtual EncodedElt* new_elt(elt_t t) = 0;
	virtual void del_elt(elt_t t, EncodedElt* elt) = 0;
	virtual EncodedProduct* new_prod() = 0;

	virtual void print(elt_t t, EncodedElt* e) = 0;
	virtual void print(EncodedProduct* p) = 0;

	virtual EncodedEltArray* new_elt_array(elt_t t, uint32_t len, bool preallocate=true) = 0;

	virtual void encode(elt_t t, FieldElt* in, EncodedElt* out);

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
	virtual void prepareForManyEncodings(uint64_t numEncodings, uint32_t maxMemGb, bool prepIsFree);
	virtual void doneWithManyEncodings();

  virtual ip_handle_t prepareForInnerProduct(EncodedEltArray* bases, int len, int numUses, int max_exponent_size, int maxMemGb, bool precompFree);
	virtual void doneWithInnerProduct(ip_handle_t handle);
	virtual void innerProduct(ip_handle_t handle, FieldEltArray* exp, int len, EncodedElt* r);	
	
#ifdef CONDENSE_ARCHIVES
	virtual void compress(EncodedElt* in, LCondensedEncElt& out) ;
	virtual void compress(REncodedElt* in, RCondensedEncElt& out) ;
	virtual void compressMany(EncodedEltArray<EncodedElt>* in, LCondensedEncElt* out, int count, bool overWriteOk) ;	// overWriteOk = true if okay to overwrite in
	virtual void compressMany(REncodedEltArray* in, RCondensedEncElt* out, int count, bool overWriteOk) ;	// overWriteOk = true if okay to overwrite in
	virtual void decompress(LCondensedEncElt& in, EncodedElt* out) ;
	virtual void decompress(RCondensedEncElt& in, REncodedElt* out) ;
	virtual void decompressMany(LCondensedEncElt* in, EncodedEltArray<EncodedElt>* out, int count) ;
	virtual void decompressMany(RCondensedEncElt* in, REncodedEltArray* out, int count) ;
#endif // CONDENSE_ARCHIVES

	////////  ECC specific-functionality  ////////////
	// Provided by children
	virtual uint32_t bytes_per(elt_t t) = 0;
	virtual uint32_t num_exponent_bits() = 0;
	
	// Children are expected to provide some basic encoding, so we can bootstrap ourselves
	virtual void encodeSlow(elt_t t, FieldElt* in, EncodedElt* out) = 0;
	
	// Provided by parent
	void allocatePreCalcMem(uint64_t numExps, uint32_t maxMemGb, bool precompFree);
	int computeIndex(FieldElt* in, int j, int k);

protected:		
	Field* srcField;

  EncodedCache enc_cache[2];
  int cache_hits;
  int cache_attempts;

	// Used for fast encoding calculations
	uint32_t h, v, a, b;
	EncodedEltArray** encodingTable[2];
	bool precomputedPowers;

	// Use for precomputing terms for a multi-exponentiation
	MultiExpTables PrecompTables;
	ip_handle_t h_ctr;

	virtual void prepareForManyEncodingsHelper(elt_t t);
};


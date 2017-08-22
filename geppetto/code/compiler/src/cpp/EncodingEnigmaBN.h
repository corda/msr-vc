#pragma once

/* 
 * Implements an encoding that maps values to powers of a generator of a BN ECC,
 * built using Microsoft's crypto library
 */

#include "Types.h"
#include "Field.h"
#include "EncodingExp.h"
#include "EncodedDigitArray.h"

class EncodingEnigmaBNMembers;

#pragma warning( push )
#pragma warning( disable : 4200 )	// Disable the warning about the 0-length array
class EncodedEltEnigmaBN : public EncodedElt { 
public: 
	static EncodedEltEnigmaBN* newElt(uint32_t num_digits);
	inline digit_t* get_digits() { return digits; }
private:
	digit_t digits[];

	EncodedEltEnigmaBN(EncodedEltEnigmaBN const &);            // undefined, but here to supress warnings
	EncodedEltEnigmaBN& operator=(EncodedEltEnigmaBN const &);  // undefined, but here to supress warnings
};
#pragma warning( pop )

class EncodedProductEnigmaBN : public EncodedProduct { 
public: 
	digit_t digits[48];
};

class CondensedEncEltEnigmaBN : public CondensedEncElt { 
};

typedef EncodedEltEnigmaBN LEncodedEltEnigmaBN;
typedef EncodedEltEnigmaBN REncodedEltEnigmaBN;


class EncodedEltEnigmaBNArray : public EncodedDigitArray {
public:
  EncodedEltEnigmaBNArray(elt_t t, uint32_t len, bool allocate);
};

class EncodingEnigmaBN : public EncodingExp {
public:
	EncodingEnigmaBN();
	~EncodingEnigmaBN();	

	virtual EncodedElt* new_elt(elt_t t);
	virtual void del_elt(elt_t t, EncodedElt* elt);
	virtual EncodedProduct* new_prod();
  virtual uint32_t size_of_prod();

  virtual void print(elt_t t, EncodedElt* e);
  virtual void print(elt_t t, EncodedElt* e, bool affine);
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

  // Reading and writing to files/network/etc.
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

  Field* getBaseField(bool debug=false);

	virtual void encodeSlow(elt_t t, FieldElt* in, EncodedElt* out);

	virtual uint32_t bytes_per(elt_t t);
	virtual uint32_t num_exponent_bits();

	enum addsub_t { Add, Sub };
	void addsub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, addsub_t op);

	static const uint32_t num_L_digits = 3 * 4;  // Three field elements (4 digits each) to define a projective point on the base curve
	static const uint32_t num_R_digits = 3*2*4;  // Three field elements (each 2 field elements of 4 digits each) to define a projective point on the twist curve

  static const int       aff_size_in_digits_L = 2 * 4; // 2 field elts/point @ 4 digits per field elt
  static const int condensed_size_in_digits_L = 1 * 4; // 1 field elts/point @ 4 digits per field elt

  static const int       aff_size_in_digits_R = 2 * 8; // 2 field elts/point @ 8 digits per field elt
  static const int condensed_size_in_digits_R = 1 * 8; // 1 field elts/point @ 8 digits per field elt
  
  Field* get_tiny_bn_field();

  // For performance testing purposes only
  Field* getMontField();
  void testToFromMont(Field* montField, int num_trials);

  virtual void print_stats();
  virtual void collect_stats(bool on);
  virtual void reset_stats();

  static int num_live_enc_elts[2];

protected:		
	EncodingEnigmaBNMembers* members;

	// Projective versions of the Base and Twist generators (also a version of g^1)
	EncodedEltEnigmaBN* gen[2];

	// Frequently encoded elts
	EncodedEltEnigmaBN* Zero[2];	

  FieldElt* montR;
  FieldElt* montRinv;
  Field* base_field;

  void makeBaseField(bool debug = false);

  virtual void to_mont(FieldElt* normal, FieldElt* mont);
  virtual void from_mont(FieldElt* mont, FieldElt* normal);

  virtual void compressLaff(bool elt_is_zero, digit_t* aff_in, digit_t* out);
  virtual void compressRaff(bool elt_is_zero, digit_t* aff_in, digit_t* out);
  virtual void decompressL(digit_t* in, LEncodedElt* out);
  virtual void decompressR(digit_t* in, REncodedElt* out);
  
  bool collecting_stats;
  int num_addsub[2];
  int num_double[2];
  int num_constmul[2];
  int num_pair;
};

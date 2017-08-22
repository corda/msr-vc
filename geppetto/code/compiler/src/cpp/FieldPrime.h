#pragma once
#include <assert.h>
#include <string.h>
#include "Types.h"
#include "Field.h"

typedef FieldEltDigits FieldEltPrime;

class FieldPrimeMembers;
class FieldDescription;

#ifdef __FFLIB
class __declspec(dllexport) FieldPrime : public Field {
#else
class FieldPrime : public Field {
#endif
public:  
	FieldPrime(FieldDescription* field, int bitsPerElt);
	virtual ~FieldPrime();

  static FieldPrime* get_CP6_base_field();

	// Returns a randomly selected non-zero field element
	// It's your responsibility to free it
	FieldElt* newRandomElt();

	// Sets an existing field element to a randomly selected non-zero field element
	void assignRandomElt(FieldElt* elt);

	// Sets an existing field element to a new, never returned value
	void assignFreshElt(FieldElt* elt);

	// Return a fresh field elt with no value
	// It's your responsibility to free it
	FieldElt* newElt();
	virtual void delElt(FieldElt* elt);
	FieldEltArray* newEltArray(uint32_t len, bool allocate);

	// Check if an element is zero
	bool isZero(FieldElt* elt);

	// Set an existing field element to zero
	void zero(FieldElt* elt);
	void zero(FieldEltArray* eltArray, int count) ;

	// Set an existing field element to one
	void one(FieldElt* elt);

	unsigned int bit(int index, const FieldElt* in);

	// Extracts count bits, starting at index and heading towards the MSB
	int getBits(FieldElt *elt, int index, int count);

  // Returns the most-significant position with a bit set
  virtual int highest_bit(FieldElt *elt);

	//void modulus(FieldElt *elt);

  // Number of int64s in a field element
  int eltSize();

  // Compute a non-cryptographic hash of x
  virtual size_t hash(FieldElt* x);

	// Set an existing field element to the integer value i
  void set(FieldElt* a, long long i);

  // Converts to/from integer array.  May use less than 4 ints.  w is expected to be in little endian
  virtual void getFull(const FieldElt* a, _int64 w[4]);
  virtual void setFull(FieldElt* a, const _int64 w[4]);

  // Returns the modulus the field currently uses, in little endian
  virtual void getMod(_int64 w[4]);

	// c <- a*b
	void mul(const FieldElt* a, const FieldElt* b, FieldElt* c);

	// c <- a/b
	void div(const FieldElt* a, const FieldElt* b, FieldElt* c);

	// c <- a+b
	void add(const FieldElt* a, const FieldElt* b, FieldElt* c);

	// c <- a-b
	void sub(const FieldElt* a, const FieldElt* b, FieldElt* c);

	// c <- a^exp (where ^ is exponentiation)
	void exp(const FieldElt* a, int exp, FieldElt* r);

	// Is a < b over the integers?
	//bool lt(FieldElt* a, FieldElt* b);

	// Reduce the value to be less than 2^bits, where bits <= 64
	virtual void truncate(FieldElt* a, int bits);

	// dst <- src
	inline void copy(const FieldElt* src, FieldElt* dst) { memcpy(((FieldEltDigits*)dst)->get_digits(), ((FieldEltDigits*)src)->get_digits(), num_digits_per_field_elt*sizeof(digit_t)); }

	// The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
	bool equal(const FieldElt* f1, const FieldElt* f2);
	int compare(const FieldElt* f1, const FieldElt* f2);
	
	virtual void print(FieldElt* e);	
  virtual char* to_string(FieldElt* e);

  virtual void read (Archiver* arc, FieldElt* dst);
  virtual void write(Archiver* arc, FieldElt* val);

  virtual FFT_info* get_fft_info();
	
  virtual void to_mont(FieldElt* normal, FieldElt* mont);
  virtual void from_mont(FieldElt* mont, FieldElt* normal);
	
  FieldPrimeMembers* members;

  static int num_live_field_elts;

  virtual uint32_t num_digits_per_elt() { return num_digits_per_field_elt; }

private:	
	int num_digits_per_field_elt;
	
	FieldEltDigits* ctr;
	FieldEltDigits* the_one;

	int bits_per_elt;

  FieldElt* montR;
  FieldElt* montRinv;
};

#pragma once

#include <string.h>
#include <assert.h>
#include "Types.h"
#include "Archiver.h"
#include <unordered_map>

typedef uint64_t digit_t;

class FFT_info;

class FieldEltArray;
class Field;

class FieldElt { 
protected: 
	FieldElt() {;} // Ensure this ifc class can't be instantiated
};


typedef unordered_map<unsigned long, FieldElt*> FieldElts;
typedef pair<const unsigned long, FieldElt*> FieldPair;

#pragma warning( push )
#pragma warning( disable : 4200 )	// Disable the warning about the 0-length array
class FieldEltDigits : public FieldElt { 
public: 
	static FieldEltDigits* newElt(uint32_t num_digits);
	inline digit_t* get_digits() { return digits; }
private:
	digit_t digits[];

	FieldEltDigits(FieldEltDigits const &);            // undefined, but here to supress warnings
	FieldEltDigits& operator=(FieldEltDigits const &);  // undefined, but here to supress warnings
};
#pragma warning( pop )

class FieldEltArray { 
public:
	virtual ~FieldEltArray();
	virtual FieldElt* elt(uint32_t index) = 0;
	virtual void set(uint32_t index, FieldElt* elt) = 0;
	virtual digit_t* get_digits() = 0;
	virtual void zero(Field* field) = 0;	// Zero all of the elements
  virtual int length() = 0;
  virtual bool has(uint32_t index) = 0;
  virtual bool complete() = 0;    // Does this array contain all possible elements?

  // Primitive iterator interface
  virtual void begin() = 0;
  virtual bool end() = 0;
  virtual FieldPair* next() = 0;
protected: 
	FieldEltArray(); // Ensure this ifc class can't be instantiated
};


// For fields with elements that consist of enigma digit_t
class FieldEltDigitsArray : public FieldEltArray {
public:
	// If allocate=true, we take responsibility for all memory needs
	// If not, it's assumed someone will point us at the memory needed for the FieldElts
	FieldEltDigitsArray(uint32_t size, uint32_t digits_per_elt, bool allocate);
	virtual ~FieldEltDigitsArray();

	inline FieldElt* elt(uint32_t index);
	inline digit_t* get_digits();
	void set(uint32_t index, FieldElt* elt);
	virtual void zero(Field* field);	// Zero all of the elements
  virtual int length();
  virtual bool has(uint32_t index);
  virtual bool complete();  

  virtual void begin();
  virtual bool end();
  virtual FieldPair* next();
private:
	uint32_t size;
	uint32_t digits_per_elt;
	bool allocated;
	digit_t* digits;	    // Elts we allocated
	FieldEltDigits** elts;  // Elts that someone allocated for us
  uint32_t iterator;
};


// For fields with elements that consist of enigma digit_t
class FieldEltDigitsMap : public FieldEltArray {
public:
  // If allocate=true, we take responsibility for all memory needs
  // If not, it's assumed someone will point us at the memory needed for the FieldElts
  FieldEltDigitsMap(Field* field, uint32_t digits_per_elt, bool allocate = true);
  virtual ~FieldEltDigitsMap();

  inline FieldElt* elt(uint32_t index);
  inline digit_t* get_digits();
  void set(uint32_t index, FieldElt* elt);
  virtual void zero(Field* field);	// Zero all of the elements
  virtual int length();  
  virtual bool has(uint32_t index);
  virtual bool complete();

  virtual void begin();
  virtual bool end();
  virtual FieldPair* next();
private:
  Field* field;
  uint32_t size;
  uint32_t digits_per_elt;
  bool allocate;
  FieldElts elts;
  FieldElts::iterator iterator;
};

#ifdef __FFLIB
class __declspec(dllexport) Field {
#else
class Field {
#endif
public:
	//Field();
	virtual ~Field();

	// Returns a randomly selected non-zero field element
	// It's your responsibility to free it
	virtual FieldElt* newRandomElt() = 0;

	// Sets an existing field element to a randomly selected non-zero field element
	virtual void assignRandomElt(FieldElt* elt) = 0;

	// Sets an existing field element to a new, never returned value
	virtual void assignFreshElt(FieldElt* elt) = 0;

	// Return a fresh field elt with no value
	// It's your responsibility to free it
	virtual FieldElt* newElt() = 0;
	virtual void delElt(FieldElt* elt) = 0;

	virtual FieldEltArray* newEltArray(uint32_t len, bool allocate) = 0;

	// Check if an element is zero
	virtual bool isZero(FieldElt* elt) = 0;
	// Set an existing field element to zero
	virtual void zero(FieldElt* elt) = 0;
	//static void Zero(FieldElt* elt) = 0;
	virtual void zero(FieldEltArray* eltArray, int count) = 0;

	// Set an existing field element to one
	virtual void one(FieldElt* elt) = 0;

	virtual unsigned int bit(int index, const FieldElt* in) = 0;

	// Extracts count bits, starting at index and heading towards the MSB
	virtual int getBits(FieldElt *elt, int index, int count) = 0;

  // Returns the most-significant position with a bit set
  virtual int highest_bit(FieldElt *elt) = 0;

  // Number of int64s in a field element
  virtual int eltSize() = 0;

  // Compute a non-cryptographic hash of x
  virtual size_t hash(FieldElt* x) = 0;

	// Set an existing field element to the integer value i
	virtual void set(FieldElt* a, long long i) = 0;

  // Converts to/from integer array.  May use less than 4 ints.  w is expected to be in little endian
  virtual void getFull(const FieldElt* a, _int64 w[4]) = 0;
  virtual void setFull(FieldElt* a, const _int64 w[4]) = 0;

  // Returns the modulus the field currently uses, in little endian
  virtual void getMod(_int64 w[4]) = 0;

	// c <- a*b
	virtual void mul(const FieldElt* a, const FieldElt* b, FieldElt* c) = 0;

	// c <- a/b
	virtual void div(const FieldElt* a, const FieldElt* b, FieldElt* c) = 0;

	// c <- a+b
	virtual void add(const FieldElt* a, const FieldElt* b, FieldElt* c) = 0;

	// c <- a-b
	virtual void sub(const FieldElt* a, const FieldElt* b, FieldElt* c) = 0;

	// c <- a^exp (where ^ is exponentiation)
	virtual void exp(const FieldElt* a, int exp, FieldElt* r) = 0;

	// Reduce the value to be less than 2^bits
	virtual void truncate(FieldElt* a, int bits) = 0;

	// dst <- src
	virtual void copy(const FieldElt* src, FieldElt* dst) = 0;
	virtual void copyArray(const FieldEltArray* src, FieldEltArray* dst);

	virtual bool equal(const FieldElt* f1, const FieldElt* f2) = 0;
	virtual int compare(const FieldElt* f1, const FieldElt* f2) = 0;
	
	virtual void print(FieldElt* e) = 0;
  virtual char* to_string(FieldElt* e) = 0;

  virtual void read (Archiver* arc, FieldElt* dst) = 0;
  virtual void write(Archiver* arc, FieldElt* val) = 0;

  // Used for the FFTs in Poly.*
  virtual FFT_info* get_fft_info() = 0;

  virtual void to_mont(FieldElt* normal, FieldElt* mont) = 0;
  virtual void from_mont(FieldElt* mont, FieldElt* normal) = 0;

  virtual uint32_t num_digits_per_elt() = 0;
};

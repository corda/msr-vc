#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "Field.h"
#include "FieldPrime.h"

Field::~Field() { ; }

FieldEltDigits* FieldEltDigits::newElt(uint32_t num_digits) {
	return (FieldEltDigits*) new digit_t[num_digits];
}

FieldEltArray::FieldEltArray() { ; }
FieldEltArray::~FieldEltArray() { ; }

FieldEltDigitsArray::FieldEltDigitsArray(uint32_t size, uint32_t digits_per_elt, bool allocate) {
	this->size = size;
	this->digits_per_elt = digits_per_elt;
	this->allocated = allocate;


	if (allocate) {
		// We need enough memory for all of the digits
		this->digits = new digit_t[size*digits_per_elt];
    memset(this->digits, 0, size*digits_per_elt*sizeof(digit_t));
		this->elts = NULL;
#ifdef PERF_DEBUG
    FieldPrime::num_live_field_elts += size;
#endif
	} else {
		// Just need pointers
		this->digits = NULL;
		this->elts = new FieldEltDigits*[size];
	}
}

FieldEltDigitsArray::~FieldEltDigitsArray() {
	if (allocated) {
		delete [] digits;
#ifdef PERF_DEBUG
    FieldPrime::num_live_field_elts -= this->size;
#endif
	} else {
		delete [] elts;
	}
}

digit_t* FieldEltDigitsArray::get_digits() {
	assert(allocated);	// Doesn't make sense if we're just holding pointers
	return digits; 
}

FieldElt* FieldEltDigitsArray::elt(uint32_t index) { 
	assert (index < size); 
	if (allocated) {
		return (FieldEltDigits*)&digits[index * digits_per_elt];
	} else {
		return elts[index]; 
	}
}

void FieldEltDigitsArray::set(uint32_t index, FieldElt* elt) {
	assert(index < size);
	assert(elt);
	assert(!allocated);	// If we did the allocation, should be using copy.  Set is for avoiding copies.
	FieldEltDigits* fd = (FieldEltDigits*)elt;

	elts[index] = fd;
}

void FieldEltDigitsArray::zero(Field* field) {
	if (allocated) {
		memset(digits, 0, size*digits_per_elt*sizeof(digit_t));
	} else {
		for (uint32_t i = 0; i < size; i++) {
			field->zero(elts[i]);
		}
	}
}

int FieldEltDigitsArray::length() {
  return size;
}

bool FieldEltDigitsArray::has(uint32_t index) {
  return index < size;
}

bool FieldEltDigitsArray::complete() {
  return true;
}

void FieldEltDigitsArray::begin() {
  iterator = 0;
}

bool FieldEltDigitsArray::end() {
  return iterator == size;
}
FieldPair* FieldEltDigitsArray::next() {
  assert(false);  // TODO: I think this will leak FieldPairs all over the place :(
  return new FieldPair(iterator, elt(iterator));
}

void Field::copyArray(const FieldEltArray* csrc, FieldEltArray* dst) {
	auto* src = const_cast<FieldEltArray*>(csrc);
	assert(dst->length() == src->length());
	for (int i = 0; i < src->length(); ++i) copy(src->elt(i), dst->elt(i));
}




FieldEltDigitsMap::FieldEltDigitsMap(Field* field, uint32_t digits_per_elt, bool allocate) {
  this->field = field;
  this->digits_per_elt = digits_per_elt;
  this->allocate = allocate;
}

FieldEltDigitsMap::~FieldEltDigitsMap() {
  if (allocate) { // Then we're responsible to clean up    
    for (FieldElts::iterator iter = elts.begin(); iter != elts.end(); iter++) {
      field->delElt(iter->second);
    }    
  }
}

FieldElt* FieldEltDigitsMap::elt(uint32_t index) {
  FieldElts::iterator result = elts.find(index);
  if (result == elts.end()) {
    if (!this->allocate) {
      assert(false);  // Shouldn't attempt to access something that's not here!
      return NULL;
    } else {
      FieldElt* ret = this->field->newElt();
      field->zero(ret);
      elts.insert(FieldPair(index, ret));
      return ret;
    }
  } else {
    return result->second;
  }
}

digit_t* FieldEltDigitsMap::get_digits() {
  assert(false);  // Not supported!
  return NULL;  
}

void FieldEltDigitsMap::set(uint32_t index, FieldElt* elt) {
  assert(!allocate);	// If we do the allocation, should be using copy.  Set is for avoiding copies.
  elts.insert(FieldPair(index, elt));
}

void FieldEltDigitsMap::zero(Field* field) {
  // Zero all of the elements
  for (FieldElts::iterator iter = elts.begin(); iter != elts.end(); iter++) {
    field->zero(iter->second);
  }
}

int FieldEltDigitsMap::length() {
  return (int)elts.size();
}

bool FieldEltDigitsMap::has(uint32_t index) {
  return elts.find(index) != elts.end();
}

bool FieldEltDigitsMap::complete() {
  return false;
}

void FieldEltDigitsMap::begin() {
  iterator = elts.begin();
}
bool FieldEltDigitsMap::end() {
  return iterator == elts.end();
}
FieldPair* FieldEltDigitsMap::next() {
  FieldPair* ret = &(*iterator);
  iterator++;
  return ret;
}
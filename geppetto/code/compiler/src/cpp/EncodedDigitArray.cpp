#include "EncodedDigitArray.h"
#include "EncodingEnigmaBN.h"


EncodedDigitArray::EncodedDigitArray(elt_t t, uint32_t len, bool allocate, int digits_per_elt_i) {
	this->type = t;
	this->len = len;
	this->allocated = allocate;
  this->digits_per_elt = digits_per_elt_i;

	if (allocate) {
		// We need enough memory for all of the digits
		this->digits = new digit_t[len*this->digits_per_elt];
		this->elts = NULL;

#ifdef PERF_DEBUG
    EncodingEnigmaBN::num_live_enc_elts[this->type] += this->len;
#endif 
	} else {
		// Just need pointers
		this->digits = NULL;
		this->elts = new EncodedElt*[len];
	}
}

EncodedDigitArray::~EncodedDigitArray() {
	if (allocated) {
		delete [] digits;
#ifdef PERF_DEBUG
    EncodingEnigmaBN::num_live_enc_elts[this->type] -= this->len;
#endif 
	} else {
		delete [] elts;
	}
}

digit_t* EncodedDigitArray::get_digits() {
	assert(allocated);	// Doesn't make sense if we're just holding pointers
	return digits; 
}

EncodedElt* EncodedDigitArray::elt(uint32_t index) { 
	assert (index < len); 
	if (allocated) {
		return (EncodedElt*)&digits[index * digits_per_elt];
	} else {
		return elts[index]; 
	}
}

void EncodedDigitArray::set(uint32_t index, EncodedElt* elt) {
	assert(index < len);
	assert(elt);
	assert(!allocated);	// If we did the allocation, should be using copy.  Set is for avoiding copies.	

	elts[index] = elt;
}

void EncodedDigitArray::reset() {
  assert(!allocated);
  for (uint32_t i = 0; i < len; i++) {
    elts[i] = NULL;
  }
}
#pragma once

#include "Encoding.h"


// For encodings with elements that consist of enigma digit_t's
class EncodedDigitArray : public EncodedEltArray {
public:
	// If allocate=true, we take responsibility for all memory needs
	// If not, it's assumed someone will point us at the memory needed for the Elts
	EncodedDigitArray(elt_t t, uint32_t len, bool allocate, int digits_per_elt_i);
	virtual ~EncodedDigitArray();

  void reset();
	inline EncodedElt* elt(uint32_t index);	
	void set(uint32_t index, EncodedElt* elt);

  inline digit_t* get_digits();

private:
	digit_t* digits;	    // Elts we allocated
	uint32_t digits_per_elt;
	EncodedElt** elts;  // Elts that someone allocated for us
};

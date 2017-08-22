////////////////////////////////////////////////////////
//  Code for loading encoded values in native code
////////////////////////////////////////////////////////
#pragma once

#ifndef MQAP 
#include <stdio.h>

size_t num_bytes_read = 0;

// Fix fread calls to use a buffered read
void buffered_read(void* buffer, int size, FILE* file) {
  size_t total_read = 0;
  while (total_read != size) {
    size_t read = fread(buffer, 1, size, file);
	assert(read); // unexpected EOF
    total_read += read;
  }
  num_bytes_read += total_read;
}

#ifdef TINY
#define NUM_DIGITS 1
#else
#define NUM_DIGITS 4
#endif //TINY

void load_field_elt(FILE* file, QapFp_t elt) {
	//assert(!feof(file));  
	buffered_read(elt->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e);
}

void load_int32(FILE* file, int *i) {
	//assert(!feof(file));  
	buffered_read(i, sizeof(int), file); 
}

void load_encodedL(FILE* file, struct s_EncodedEltL* elt) {
	//assert(!feof(file));  
	buffered_read(elt->e->X->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e->X->e);
	buffered_read(elt->e->Y->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e->Y->e);
	// enc_print_rawL("reading", elt);
}

void load_encodedR(FILE* file, struct s_EncodedEltR* elt) {
	//assert(!feof(file));  
	buffered_read(elt->e->X->a0->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e->X->a0->e);
	buffered_read(elt->e->X->a1->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e->X->a1->e);
	buffered_read(elt->e->Y->a0->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e->Y->a0->e);
	buffered_read(elt->e->Y->a1->e->base_fp->limbs, sizeof(uint64_t)* NUM_DIGITS, file);  elem_to_mont(elt->e->Y->a1->e);
	// enc_print_rawR("reading", elt);
}

#endif // MQAP
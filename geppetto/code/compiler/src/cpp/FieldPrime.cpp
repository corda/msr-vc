#include <stdio.h>
#include <math.h>
#include "FieldPrime.h"
#include "EnigmaWrappers.h"
#include <intrin.h>

#pragma intrinsic(_BitScanReverse64)

int FieldPrime::num_live_field_elts = 0;

FieldPrime::FieldPrime(FieldDescription* field, int bitsPerElt) {
	assert(field);
	this->bits_per_elt = bitsPerElt;
	members = new FieldPrimeMembers;

	members->msfield = field->desc;	
	num_digits_per_field_elt = (int)members->msfield->elng;
	members->pbigctx = &this->members->BignumCtx;
	memset(this->members->pbigctx, 0, sizeof(bigctx_t));

	this->ctr     = (FieldEltDigits*)newElt();
	this->the_one = (FieldEltDigits*)newElt();

	memcpy(the_one->get_digits(), field->desc->one, num_digits_per_field_elt*sizeof(digit_t));
  add(the_one, the_one, ctr); // Start at 2
  add(the_one, ctr, ctr);     // Start at 3

	modfft_init(members->msfield->modulo, &members->fft_info.info, members->PBIGCTX_PASS);

  this->montR = newElt();
  set(this->montR, 2);
  exp(this->montR, 64 * 4, this->montR);
  this->montRinv = newElt();
  FieldElt* one_f = newElt();
  one(one_f);
  div(one_f, this->montR, this->montRinv);
  delElt(one_f);
}

FieldPrime::~FieldPrime() {
	this->delElt(ctr);
	this->delElt(the_one);
  this->delElt(this->montR);
  this->delElt(this->montRinv);

	modfft_uninit(&members->fft_info.info, members->PBIGCTX_PASS);

	delete members;
}

FieldPrime* FieldPrime::get_CP6_base_field() {
  uint64_t pCP6b[8] = { 0xF9000000000005BF, 0xAA8EC00000000BD3,
    0x0CF85A8000000B1F, 0x8761FD7D40000638,
    0x0A67549B77800241, 0xB007E82F1F3CB08B,
    0x41358FBC1F137315, 0x158CFFC9264E1116 };
  /* This is the 509-bit CP6b prime p = 3 (mod 4). */

  bigctx_t BignumCtx;
  bigctx_t* pbigctx = &BignumCtx;
  mp_modulus_t* base_p = new mp_modulus_t;    // TODO: We're currently leaking this

  BOOL OK = create_modulus(pCP6b, 8, FROM_LEFT, base_p, PBIGCTX_PASS);
  assert(OK);

  FieldDescription* base_desc = new FieldDescription;   // TODO: We're currently leaking this
  base_desc->desc = new field_desc_t;
  OK = Kinitialize_prime(base_p, (field_desc_t*)base_desc->desc, PBIGCTX_PASS);
  assert(OK);

  FieldPrime* baseField = new FieldPrime(base_desc, 509);
  return baseField;
}

void FieldPrime::print(FieldElt* e) {
  char* s = to_string(e);
  printf("%s", s);
  delete[] s;
}

char* FieldPrime::to_string(FieldElt* e) {
  FieldEltDigits* ep = (FieldEltDigits*)e;
  int len = num_digits_per_field_elt * sizeof(digit_t) * 2 + 1;  // +1 for NULL
  char* ret = new char[len];  
  for (int i = num_digits_per_field_elt - 1; i >= 0; i--) {
    int str_index = num_digits_per_field_elt - 1 - i;
    sprintf_s(ret + str_index * 16, len - str_index * 16, "%0.16llX", ep->get_digits()[i]);
  }
  return ret;
}

FieldElt* FieldPrime::newElt() { 
#ifdef PERF_DEBUG
  FieldPrime::num_live_field_elts++;
#endif
	return FieldEltDigits::newElt(this->num_digits_per_field_elt);
}

void FieldPrime::delElt(FieldElt* elt) { 
#ifdef PERF_DEBUG
  FieldPrime::num_live_field_elts--;
#endif
	delete [] elt;
}

FieldEltArray* FieldPrime::newEltArray(uint32_t len, bool allocate) { 
	return (FieldEltArray*) new FieldEltDigitsArray(len, this->num_digits_per_field_elt, allocate);
}

bool FieldPrime::isZero(FieldElt* elt) {
	assert(elt);
	digit_t* d = ((FieldEltDigits*)elt)->get_digits();
	for (int i = 0; i < num_digits_per_field_elt; ++i) {
		if (d[i] != 0) return false;
	}
	return true;
}

void FieldPrime::zero(FieldElt* elt) { 
	assert(elt);
	memset(((FieldEltDigits*)elt)->get_digits(), 0, num_digits_per_field_elt*sizeof(digit_t)); 
}

void FieldPrime::zero(FieldEltArray* eltArray, int count) {
		assert(eltArray);
		eltArray->zero(this);
}

void FieldPrime::one(FieldElt* elt) { 
	assert(elt);
	copy(the_one, elt); 
}


// c <- a*b
void FieldPrime::mul(const FieldElt* a, const FieldElt* b, FieldElt* c) {
	assert(a && b && c);
	FieldEltDigits* ap = (FieldEltDigits*)a;
	FieldEltDigits* bp = (FieldEltDigits*)b;
	FieldEltDigits* cp = (FieldEltDigits*)c;

	assert(sizeof(members->mulTempSpace) >= members->msfield->ndigtemps_mul * sizeof(digit_t));
	BOOL OK = Kmul(ap->get_digits(), bp->get_digits(), cp->get_digits(), members->msfield, members->mulTempSpace, members->PBIGCTX_PASS);
	assert(OK);
}

// c <- a / b
void FieldPrime::div(const FieldElt* a, const FieldElt* b, FieldElt* c) {
	assert(a && b && c);
	FieldEltDigits* ap = (FieldEltDigits*)a;
	FieldEltDigits* bp = (FieldEltDigits*)b;
	FieldEltDigits* cp = (FieldEltDigits*)c;

	assert(sizeof(members->divTempSpace) >= members->msfield->ndigtemps_arith * sizeof(digit_t));
	BOOL OK = Kdiv(ap->get_digits(), bp->get_digits(), cp->get_digits(), members->msfield, members->divTempSpace, members->PBIGCTX_PASS);
	assert(OK);
}

// c <- a+b
void FieldPrime::add(const FieldElt* a, const FieldElt* b, FieldElt* c) {
	assert(a && b && c);
	FieldEltDigits* ap = (FieldEltDigits*)a;
	FieldEltDigits* bp = (FieldEltDigits*)b;
	FieldEltDigits* cp = (FieldEltDigits*)c;

	BOOL OK = Kadd(ap->get_digits(), bp->get_digits(), cp->get_digits(), members->msfield, members->PBIGCTX_PASS);
	assert(OK);
}

// c <- a-b
void FieldPrime::sub(const FieldElt* a, const FieldElt* b, FieldElt* c) {
	assert(a && b && c);
	FieldEltDigits* ap = (FieldEltDigits*)a;
	FieldEltDigits* bp = (FieldEltDigits*)b;
	FieldEltDigits* cp = (FieldEltDigits*)c;

	BOOL OK = Ksub(ap->get_digits(), bp->get_digits(), cp->get_digits(), members->msfield, members->PBIGCTX_PASS);
	assert(OK);
}

// c <- a^exp (where ^ is exponentiation)
void FieldPrime::exp(const FieldElt* a, int exp, FieldElt* r) {
	assert(a && r);
	FieldEltDigits* ap = (FieldEltDigits*)a;
	FieldEltDigits* rp = (FieldEltDigits*)r;

	FieldEltDigits* fexp = (FieldEltDigits*)newElt();
	set(fexp, exp);
	BOOL OK = Kexpon(ap->get_digits(), fexp->get_digits(), num_digits_per_field_elt, rp->get_digits(), members->msfield, members->PBIGCTX_PASS);
	assert(OK);
	delElt(fexp);
}

// Reduce the value to be less than 2^bits
void FieldPrime::truncate(FieldElt* a, int bits) {
	assert(a);
  assert(bits >= 0);

  int numDigitBits = sizeof(digit_t) * 8;
 
  if (bits > this->num_digits_per_field_elt * numDigitBits) {
    return;   // Nothing to truncate
  }

  int digitIndex;
  if (bits == 0) {
    digitIndex = 0;
  } else {
    digitIndex = (bits-1) / numDigitBits;
  }
  bits = bits - numDigitBits * digitIndex;
  
  // Bulk clear everything above the cutoff point
	FieldEltDigits* ap = (FieldEltDigits*)a;
	for (int i = digitIndex + 1; i < this->num_digits_per_field_elt; i++) {
		ap->get_digits()[i] = 0;
	}

  // Clear within the digit at the cutoff point
  if (bits < numDigitBits) {
    ap->get_digits()[digitIndex] &= (((uint64_t)1) << bits) - 1;
  }
}


unsigned int FieldPrime::bit(int index, const FieldElt* in) {
	assert(in);
	FieldEltDigits* inp = (FieldEltDigits*)in;
	if (index > bits_per_elt) {
		return 0;
	}

	int numDigitBits = sizeof(digit_t) * 8;
	int digitIndex = index / numDigitBits;
	int intraDigitIndex = index % numDigitBits;

	digit_t mask = ((digit_t)1) << intraDigitIndex;

	return (mask & inp->get_digits()[digitIndex]) > 0 ? 1 : 0;
}

int FieldPrime::eltSize() {
  return this->num_digits_per_field_elt;
}

// Compute a non-cryptographic hash of x
size_t FieldPrime::hash(FieldElt* in) {
  assert(in);
  FieldEltDigits* inp = (FieldEltDigits*)in;
  size_t hash = 5381;
  for (int i = 0; i < this->num_digits_per_field_elt; i++) {
    digit_t val = inp->get_digits()[i];
    hash = ((hash << 5) + hash) + val;  // hash * 33 + val = djb2's hash
  }
  return hash;
}

void FieldPrime::set(FieldElt* a, long long i) {
	assert(a);
	FieldEltDigits* ap = (FieldEltDigits*)a;
	BOOL OK = Kimmediate(i, ap->get_digits(), members->msfield, members->PBIGCTX_PASS);
	assert(OK);
}

// Converts to/from integer array.  May use less than 4 ints.  w is expected to be in little endian
void FieldPrime::getFull(const FieldElt* a, _int64 w[4]) {
  assert(a);
	FieldEltDigits* ap = (FieldEltDigits*)a;

  for (int i = 0; i < this->num_digits_per_field_elt && i < 4; i++) {
    w[i] = ap->get_digits()[i];
  }
}

void FieldPrime::setFull(FieldElt* a, const _int64 w[4]) {
  assert(a);
	FieldEltDigits* ap = (FieldEltDigits*)a;

  for (int i = 0; i < this->num_digits_per_field_elt && i < 4; i++) {
    ap->get_digits()[i] = w[i];
  }
}

void FieldPrime::getMod(_int64 w[4]) {  
  for (int i = 0; i < this->members->msfield->modulo->length && i < 4; i++) {
    w[i] = this->members->msfield->modulo->modulus[i];
  }
}

// Extracts count bits, starting at index and heading towards the MSB
int FieldPrime::getBits(FieldElt *elt, int index, int count) {	
	assert(elt);
	FieldEltDigits* eltp = (FieldEltDigits*)elt;
	const int digitSize = sizeof(digit_t)*8;
	assert(count <= digitSize);

	int interDigitIndex = index / digitSize;
	int intraDigitIndex = index % digitSize;
	int numBitsToGet = min(count, digitSize - intraDigitIndex);	// Don't try to take more bits than remain in this digit

	int mask = (1 << numBitsToGet) - 1;  // mask = 000000111, where # of 1s = numBitsToGet

	digit_t shiftedDigit = (digit_t) (eltp->get_digits()[interDigitIndex] >> intraDigitIndex);

	int ret = mask & shiftedDigit;

	if (numBitsToGet < count) {	// We were at the boundary of the prev. digit, so grab remaining bits from the next digit
		int numBitsGot = numBitsToGet;
		numBitsToGet = count - numBitsToGet;
		mask = (1 << numBitsToGet) - 1;  // mask = 000000111, where # of 1s = numBitsToGet
		ret |= (eltp->get_digits()[interDigitIndex+1] & mask) << numBitsGot;	 // Make room for the bits we did get
	}

	return ret;
}

// Returns the most-significant position with a bit set.  Returns -1 if no bits are set
int FieldPrime::highest_bit(FieldElt *elt) {
  assert(elt);
  FieldEltDigits* eltp = (FieldEltDigits*)elt;

  unsigned long ret;
  unsigned long bits_per_digit_t = sizeof(digit_t) * 8;
  long position = (num_digits_per_field_elt - 1) * bits_per_digit_t;

  for (int i = num_digits_per_field_elt - 1; i >= 0; i--) {
    unsigned char success = _BitScanReverse64(&ret, eltp->get_digits()[i]);	// Finds first 1 bit, starting from msb
    //printf("Scanning digit %llu returned status %d and ret %d\n", eltp->get_digits()[i], (int)success, ret);

    if (success != 0) {  // Found a bit set
      return ret + position;
    }

    position -= bits_per_digit_t;
  }

  printf("Failed to find any bits set\n");
  return -1;
}

// Returns a randomly selected non-zero field element
// It's your responsibility to free it
FieldElt* FieldPrime::newRandomElt() {
	FieldElt* ret = newElt();
	assignRandomElt(ret);
	return ret;
}

// Sets an existing field element to a randomly selected non-zero field element
void FieldPrime::assignRandomElt(FieldElt* elt) {
	assert(elt);
#ifdef DEBUG_RANDOM
	assignFreshElt(elt);
#else
	FieldEltDigits* eltp = (FieldEltDigits*)elt;
	BOOL OK = Krandom_nonzero(eltp->get_digits(), members->msfield, members->PBIGCTX_PASS);
	assert(OK);
#endif // DEBUG_RANDOM
}

void FieldPrime::assignFreshElt(FieldElt* elt) {
	assert(elt);
	copy(ctr, elt);
	add(ctr, the_one, ctr);	
}

int FieldPrime::compare(const FieldElt* f1, const FieldElt* f2) {
  for (int i = 0; i < num_digits_per_field_elt; i++) {
    int index = num_digits_per_field_elt - i - 1;
		digit_t a = ((FieldEltDigits*)f1)->get_digits()[index];
		digit_t b = ((FieldEltDigits*)f2)->get_digits()[index];
		if (a != b) return a < b ? -1 : 1;
	}
	return 0;
}

bool FieldPrime::equal(const FieldElt* f1, const FieldElt* f2) {
	return Kequal(((FieldEltDigits*)f1)->get_digits(), ((FieldEltDigits*)f2)->get_digits(), members->msfield, members->PBIGCTX_PASS) != 0; 
}

void FieldPrime::read(Archiver* arc, FieldElt* dst) {
  FieldEltDigits* eltp = (FieldEltDigits*)dst;
  arc->read((unsigned char*)(eltp->get_digits()), this->num_digits_per_field_elt * sizeof(digit_t));
  // Convert to the Montgomery-form that the C code expects
  //to_mont(dst, dst);
}

void FieldPrime::write(Archiver* arc, FieldElt* val) {
  FieldEltDigits* eltp = (FieldEltDigits*)val;
  arc->write((unsigned char*)(eltp->get_digits()), this->num_digits_per_field_elt * sizeof(digit_t));
}

FFT_info* FieldPrime::get_fft_info() {
  return &members->fft_info;  
}


void FieldPrime::to_mont(FieldElt* normal, FieldElt* mont) {
  mul(normal, this->montR, mont);
}

void FieldPrime::from_mont(FieldElt* mont, FieldElt* normal) {
  mul(mont, this->montRinv, normal);
}
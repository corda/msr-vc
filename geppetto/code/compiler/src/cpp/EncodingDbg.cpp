#include <assert.h>
#include "EncodingDbg.h"
#include "FieldPrime.h"
#include "EnigmaWrappers.h"

EncodingDbg::EncodingDbg(field_type type) {
	switch (type) {
	case EncodingDbg::field_prime:
		{
		pbigctx = &this->BignumCtx;
		memset(this->pbigctx, 0, sizeof(bigctx_t));

		//digit_t prime = 4294967291;
		digit_t prime = 49999;    
    //digit_t prime = 1223;
		//digit_t prime = 101;
		//digit_t prime = 11;

		// Create a small prime field
		p = new mp_modulus_t; 
		BOOL OK = create_modulus(&prime, 1, FROM_LEFT, p, PBIGCTX_PASS);
		assert(OK);

		msExponentField = new field_desc_t; 

		OK = Kinitialize_prime(p, msExponentField, PBIGCTX_PASS);
		assert(OK);

    FieldDescription* desc = new FieldDescription;
    desc->desc = msExponentField;
    srcField = new FieldPrime(desc, 16);
    delete desc;
		
		break;
		}

	}

	h_ctr = 0;
	bases[0] = NULL;
	bases[1] = NULL;
}

EncodingDbg::~EncodingDbg() {
	delete srcField;
	delete msExponentField;
	uncreate_modulus(p, PBIGCTX_PASS);
	delete p;
}

Field* EncodingDbg::getSrcField() {
	return this->srcField;
}

EncodedElt* EncodingDbg::new_elt(elt_t t) {
	return new EncodedEltDbg(t, srcField->newElt());
}

void EncodingDbg::del_elt(elt_t t, EncodedElt* elt) {
	delete elt;
}

EncodedProduct* EncodingDbg::new_prod() {
	return new EncodedProductDbg(this->srcField);
}

uint32_t EncodingDbg::size_of_prod() {
  return sizeof(digit_t);
}

void EncodingDbg::print(elt_t t, EncodedElt* e) {
	assert(e);
	EncodedEltDbg* ee = (EncodedEltDbg*)e;
	
	printf("(");
	srcField->print(ee->get_f());		
	printf(")");
}

void EncodingDbg::print(EncodedProduct* p) {
	assert(p);
	EncodedProductDbg* ep = (EncodedProductDbg*)p;

	srcField->print(ep->f);	
}


EncodedEltArray* EncodingDbg::new_elt_array(elt_t t, uint32_t len, bool preallocate) {
	return new EncodedEltDbgArray(t, this, len, preallocate);
}


#ifdef DEBUG_INHERITS_FROM_EXP
void EncodingDbg::encodeSlow(elt_t t, FieldElt* in, EncodedElt* out) {
#else // !DEBUG_INHERITS_FROM_EXP
void EncodingDbg::encode(elt_t t, FieldElt* in, EncodedElt* out) {
#endif  // DEBUG_INHERITS_FROM_EXP
	assert(in && out);
	assert(t == ((EncodedEltDbg*)out)->get_t());
	srcField->copy(in, ((EncodedEltDbg*)out)->get_f());
}

void EncodingDbg::copy(elt_t t, EncodedElt* src, EncodedElt* dst) {
	assert(src && dst);
	EncodedEltDbg* srce = (EncodedEltDbg*)src;
	EncodedEltDbg* dste = (EncodedEltDbg*)dst;
	assert(t == srce->get_t() && t == dste->get_t());
	srcField->copy(srce->get_f(), dste->get_f());
}

bool EncodingDbg::equal(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	EncodedEltDbg* be = (EncodedEltDbg*)b;

	return srcField->equal(ae->get_f(), be->get_f());
}

void EncodingDbg::zero(elt_t t, EncodedElt* a) { 
	assert(a);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	assert(t == ae->get_t());
	srcField->zero(ae->get_f());
}

void EncodingDbg::one(elt_t t, EncodedElt* a) { 
	assert(a);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	assert(t == ae->get_t());
	srcField->one(ae->get_f());
}

// Adds two encoded elements.  Result goes in r.
void EncodingDbg::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	assert(a && b && r);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	EncodedEltDbg* be = (EncodedEltDbg*)b;
	EncodedEltDbg* re = (EncodedEltDbg*)r;
	assert(t == ae->get_t() && t == be->get_t() && t == re->get_t());
	srcField->add(ae->get_f(), be->get_f(), re->get_f());
}

// Subtracts two encoded elements.  Result goes in r.
void EncodingDbg::sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	assert(a && b && r);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	EncodedEltDbg* be = (EncodedEltDbg*)b;
	EncodedEltDbg* re = (EncodedEltDbg*)r;
	assert(t == ae->get_t() && t == be->get_t() && t == re->get_t());
	srcField->sub(ae->get_f(), be->get_f(), re->get_f());
}

// b = 2 * a
void EncodingDbg::doubleIt(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	EncodedEltDbg* be = (EncodedEltDbg*)b;
	assert(t == ae->get_t() && t == be->get_t());

	FieldElt* two = srcField->newElt();
	srcField->set(two, 2);

	srcField->mul(ae->get_f(), two, be->get_f());

	srcField->delElt(two);
}

// Multiplies the value encoded in g by the constant c.  Result goes in r.
// For pairings, we compute r <- g^c
void EncodingDbg::mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r) {	
	assert(g && c && r);
	EncodedEltDbg* ge = (EncodedEltDbg*)g;
	EncodedEltDbg* re = (EncodedEltDbg*)r;
	assert(t == ge->get_t() && t == re->get_t());

	srcField->mul(ge->get_f(), c, re->get_f());
}

// Mul two encoded values
void EncodingDbg::mul(EncodedElt* a, EncodedElt* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedEltDbg* ae = (EncodedEltDbg*)a;
	EncodedEltDbg* be = (EncodedEltDbg*)b;
	assert(ae->get_t() == L && be->get_t() == R);
	EncodedProductDbg* re = (EncodedProductDbg*)r;

	srcField->mul(ae->get_f(), be->get_f(), re->f);
}

// Adds two encoded products.  Result goes in r.
void EncodingDbg::add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductDbg* ae = (EncodedProductDbg*)a;
	EncodedProductDbg* be = (EncodedProductDbg*)b;
	EncodedProductDbg* re = (EncodedProductDbg*)r;
	srcField->add(ae->f, be->f, re->f);
}

// Subtracts two encoded products.  Result goes in r.
void EncodingDbg::sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductDbg* ae = (EncodedProductDbg*)a;
	EncodedProductDbg* be = (EncodedProductDbg*)b;
	EncodedProductDbg* re = (EncodedProductDbg*)r;
	srcField->sub(ae->f, be->f, re->f);
}

bool EncodingDbg::equals(EncodedProduct* a, EncodedProduct* b) {
	assert(a && b);
	EncodedProductDbg* ae = (EncodedProductDbg*)a;
	EncodedProductDbg* be = (EncodedProductDbg*)b;
	return srcField->equal(ae->f, be->f);
}

// Checks whether the plaintexts inside the encodings obey:
// L1*R1 - L2*R2 == L3*R3
bool EncodingDbg::mulSubEquals(EncodedElt* L1, EncodedElt* R1, EncodedElt* L2, EncodedElt* R2, EncodedElt* L3, EncodedElt* R3) {
	assert(L1 && R1 && L2 && R2 && L3 && R3);
	FieldElt* tmp1 = srcField->newElt();
	FieldElt* tmp2 = srcField->newElt();
	FieldElt* tmp3 = srcField->newElt();

	EncodedEltDbg* L1e = (EncodedEltDbg*)L1;
	EncodedEltDbg* R1e = (EncodedEltDbg*)R1;
	EncodedEltDbg* L2e = (EncodedEltDbg*)L2;
	EncodedEltDbg* R2e = (EncodedEltDbg*)R2;
	EncodedEltDbg* L3e = (EncodedEltDbg*)L3;
	EncodedEltDbg* R3e = (EncodedEltDbg*)R3;
	assert(L1e->get_t() == L && L2e->get_t() == L && L3e->get_t() == L);
	assert(R1e->get_t() == R && R2e->get_t() == R && R3e->get_t() == R);
	
	srcField->mul(L1e->get_f(), R1e->get_f(), tmp1);
	srcField->mul(L2e->get_f(), R2e->get_f(), tmp2);
	srcField->mul(L3e->get_f(), R3e->get_f(), tmp3);

	srcField->sub(tmp1, tmp2, tmp1);
  bool result = srcField->equal(tmp1, tmp3);

  srcField->delElt(tmp1);
	srcField->delElt(tmp2);
	srcField->delElt(tmp3);

	return result;
}

// Checks whether the plaintexts inside the encodings obey:
// L1*R1 == L2*R2
bool EncodingDbg::mulEquals(EncodedElt* L1, EncodedElt* R1, EncodedElt* L2, EncodedElt* R2) {
	assert(L1 && R1 && L2 && R2);
	FieldElt* tmp1 = srcField->newElt();
	FieldElt* tmp2 = srcField->newElt();

	EncodedEltDbg* L1e = (EncodedEltDbg*)L1;
	EncodedEltDbg* R1e = (EncodedEltDbg*)R1;
	EncodedEltDbg* L2e = (EncodedEltDbg*)L2;
	EncodedEltDbg* R2e = (EncodedEltDbg*)R2;
	assert(L1e->get_t() == L && L2e->get_t() == L);
	assert(R1e->get_t() == R && R2e->get_t() == R);
	
	srcField->mul(L1e->get_f(), R1e->get_f(), tmp1);
	srcField->mul(L2e->get_f(), R2e->get_f(), tmp2);

	return srcField->equal(tmp1, tmp2);
}

#ifndef DEBUG_INHERITS_FROM_EXP
void EncodingDbg::prepareForManyEncodings(uint64_t numEncodings, uint32_t maxMemGb, bool prepIsFree) { ; }
void EncodingDbg::doneWithManyEncodings() { ; }

ip_handle_t EncodingDbg::prepareForInnerProduct(EncodedEltArray* bases, int len, int numUses, int maxMemGb, bool precompFree) {
	assert(bases);
	assert(0 <= h_ctr && h_ctr < 2);
	ip_handle_t loc = h_ctr;
	assert(this->bases[loc] == NULL);
	this->bases[loc] = bases;
	h_ctr++;	
	return loc;
}

void EncodingDbg::doneWithInnerProduct(ip_handle_t handle) { 
	assert(0 <= handle && handle < 2);
	assert(0 < h_ctr && h_ctr <= 2);
	this->bases[handle] = NULL;
	h_ctr--;
}

void EncodingDbg::innerProduct(ip_handle_t handle, FieldEltArray* exp, int len, EncodedElt* r) {
	assert(exp && r);
	assert(0 <= handle && handle < 2);
	EncodedEltDbg* rd = (EncodedEltDbg*)r;
	srcField->zero(rd->get_f());

	FieldElt* tmp = srcField->newElt();
	srcField->zero(tmp);

	for (int i = 0; i < len; i++) {
		srcField->mul(exp->elt(i), ((EncodedEltDbg*)bases[handle]->elt(i))->get_f(), tmp);
		srcField->add(tmp, rd->get_f(), rd->get_f());
	}

	delete tmp;
}

#endif // DEBUG_INHERITS_FROM_EXP


void EncodingDbg::write(Archiver* arc, elt_t t, EncodedElt* e, bool simple) {
	assert(arc && e);
	EncodedEltDbg* ee = (EncodedEltDbg*)e;	
	FieldEltPrime* f = (FieldEltPrime*)ee->get_f();

	arc->write((unsigned char*)f->get_digits(), this->getSrcField()->eltSize() * sizeof(digit_t));
}

void EncodingDbg::read (Archiver* arc, elt_t t, EncodedElt* e) {
	assert(arc && e);
	EncodedEltDbg* ee = (EncodedEltDbg*)e;
	FieldEltPrime* f = (FieldEltPrime*)ee->get_f();

	arc->read((unsigned char*)f->get_digits(), this->getSrcField()->eltSize() * sizeof(digit_t));
}

void EncodingDbg::compress(elt_t t, EncodedElt* e) {
  // Nothing to do
}

void EncodingDbg::decompress(elt_t t, digit_t* compressed_elt, EncodedElt* e) {
  copy(t, (EncodedElt*)compressed_elt, e);
}


#ifdef DEBUG_INHERITS_FROM_EXP
uint32_t EncodingDbg::bytes_per(elt_t t) {
	return sizeof(digit_t);
}

uint32_t EncodingDbg::num_exponent_bits() {
	return sizeof(digit_t)*8;		// A very conservative overestimate
}
#endif // DEBUG_INHERITS_FROM_EXP


EncodedEltDbgArray::EncodedEltDbgArray(elt_t type, Encoding* encoding, uint32_t len, bool preallocate) {
	this->type = type;
	this->len = len;
	elts = new EncodedEltDbg*[len];
	allocated = preallocate;
	
	for (uint32_t i = 0; i < len; i++) {
		if (preallocate) {
			elts[i] = (EncodedEltDbg*)encoding->new_elt(type);
		} else {
			elts[i] = NULL;
		}
	}
}

EncodedEltDbgArray::~EncodedEltDbgArray() {
	if (allocated) {
		for (uint32_t i = 0; i < len; i++) {
			delete elts[i];
		}
	}
	delete [] elts;
}

EncodedElt* EncodedEltDbgArray::elt(uint32_t index) {
	assert(index < len);
	assert(elts[index] != NULL);
	return elts[index];
}

void EncodedEltDbgArray::set(uint32_t index, EncodedElt* elt) {
	assert(index < len);
	assert(elts[index] == NULL);	// Otherwise we're overwriting and existing value!
	elts[index] = (EncodedEltDbg*)elt;
}


void EncodedEltDbgArray::reset() {
  assert(!allocated);
  for (uint32_t i = 0; i < len; i++) {
    elts[i] = NULL;
  }
}
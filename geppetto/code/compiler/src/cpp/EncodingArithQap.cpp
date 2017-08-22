#include <assert.h>
#include "EncodingArithQap.h"
#include <bitset>
#include <intrin.h>
#include "Poly.h"
#include "FieldPrime.h"
#include "EncodingEnigmaBN.h"
#include "ArithQapWrappers.h"

// Protect ourselves from ARITH's #define
#ifdef mul
#undef mul
#endif

EncodingArithQap::EncodingArithQap()
{
	members = new EncodingArithQapMembers;

	QapInterface_init(&members->qap_interface);

	members->msft_encoding = new EncodingEnigmaBN();
	this->srcField = members->msft_encoding->getSrcField();

	QapFp_wp_init(members->Lg.me);
	QapFp_wp_copy(members->Lg.me, members->qap_interface.Lg);
	QapFp2_wp_init(members->Rg.me);
	QapFp2_wp_copy(members->Rg.me, members->qap_interface.Rg);
	QapFp_wp_init(members->Lzero.me);
	QapFp_wp_copy(members->Lzero.me, members->qap_interface.Lzero);
	QapFp2_wp_init(members->Rzero.me);
	QapFp2_wp_copy(members->Rzero.me, members->qap_interface.Rzero);

}

EncodingArithQap::~EncodingArithQap() {
	delete members->msft_encoding;	// Delete the encoding that provided our srcField, not the srcField itself

#if 0	// TODO fill in corresponding destruction sequence
  QapFp_wp_free(members->Lzero.me);
  QapFp2_wp_free(members->Rzero.me);

  QapFp_free (members->b);
  QapFp2_free (members->B);
  QapFp2_free (members->B3);
  QapFp2_free (members->l00);
  QapFp2_free (members->l01);
  QapFp2_free (members->l10);
  QapFp2_wp_free (members->Rg.me);
  QapFp_wp_free (members->Lg.me);
  bn_bls12_finalexpo_free_config ();
  bn_bls12_miller_free_config ();
  QapFp12_free_config (); 
  QapFp6_free_config ();
  QapFp4_free_config ();
  QapFp2_free_weierstrass ();
  QapFp2_free_config (); 
  QapFp_free_config ();
#endif
}

// TODO: Adjust to count pointers as well
uint32_t EncodingArithQap::bytes_per(elt_t t) {
	if (t == L) {
		return num_L_digits*sizeof(digit_t);
	} else {
		assert(t == R);
		return num_R_digits*sizeof(digit_t);
	}
}

uint32_t EncodingArithQap::num_exponent_bits() {
	return 4*sizeof(digit_t)*8 - 1;	// BN_curve uses 4 digits, but doesn't use them all (only uses 254 bits)
}


EncodedElt* EncodingArithQap::new_elt(elt_t t) {
	if (t == L) {
    LEncodedEltQapArith* elt = new LEncodedEltQapArith;
    QapFp_wp_init(elt->me);
		return elt;
	} else {
		assert(t == R);
		REncodedEltQapArith* elt = new REncodedEltQapArith;
    QapFp2_wp_init(elt->me);
		return elt;
	}
}

EncodedEltArray* EncodingArithQap::new_elt_array(elt_t t, uint32_t len, bool preallocate) {
	return new EncodedEltQapArithArray(t, this, len, preallocate);
}


void EncodingArithQap::del_elt(elt_t t, EncodedElt* elt) {
  assert(elt);
  if (t == L) {
    QapFp_wp_free(((LEncodedEltQapArith*)elt)->me);
  } else {
    assert(t == R);
    QapFp2_wp_free(((REncodedEltQapArith*)elt)->me);
  }
	delete elt;
}

EncodedProduct* EncodingArithQap::new_prod() {
	EncodedProductQapArith* ret = new EncodedProductQapArith();
  QapFp12_init(ret->me);
  return ret;
}

uint32_t EncodingArithQap::size_of_prod() {
  return sizeof(QapFp12_t);
}

void EncodingArithQap::print_raw(elt_t t, EncodedElt* e) {
	if (t == L) {
		LEncodedEltQapArith* en = (LEncodedEltQapArith*)e;
		printf("( ");
		QapFp_print(en->me->X);
		printf(", ");
		QapFp_print(en->me->Y);
		printf(")");
	} else {
		assert(t == R);
		REncodedEltQapArith* en = (REncodedEltQapArith*)e;
		printf("( ");
		QapFp2_print(en->me->X);
		printf(", ");
    QapFp2_print(en->me->Y);
		printf(")");
	}
}

void EncodingArithQap::print(elt_t t, EncodedElt* e) {
	if (isZero(t, e))
	{
		printf("(ZERO)\n");
	}
	else
	{
		if (t == L) {
			LEncodedEltQapArith* en = (LEncodedEltQapArith*)e;
			printf("( ");
			QapFp_print(en->me->X);
			printf(", ");
			QapFp_print(en->me->Y);
			printf(")");
		}
		else {
			assert(t == R);
			REncodedEltQapArith* en = (REncodedEltQapArith*)e;
			printf("( ");
			QapFp2_print(en->me->X);
			printf(", ");
			QapFp2_print(en->me->Y);
			printf(")");
		}
	}
}

void EncodingArithQap::print(EncodedProduct* p) {
	EncodedProductQapArith* ep = (EncodedProductQapArith*)p;
	QapFp12_print(ep->me);
}

void EncodingArithQap::copy(elt_t t, EncodedElt* src, EncodedElt* dst) {
	if (t == L) {
    QapFp_wp_copy(((LEncodedEltQapArith*)dst)->me, ((LEncodedEltQapArith*)src)->me);
	} else {
		assert(t == R);
		QapFp2_wp_copy(((REncodedEltQapArith*)dst)->me, ((REncodedEltQapArith*)src)->me);
	}	
}

bool EncodingArithQap::equal(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	testOnCurve(t, a);
	testOnCurve(t, b);

  bool ret = false;

  // Short circuit 0, since affinizing 0 spins into infinity
  if (isZero(t, a)) { return isZero(t, b); }
  else if (isZero(t, b)) { return false; }

  if (t == L) {
    LEncodedEltQapArith* ae = (LEncodedEltQapArith*)a;
	  LEncodedEltQapArith* be = (LEncodedEltQapArith*)b;

    ret = 0 != QapFp_weierstrass_equal_aff(ae->me, be->me);
  } else {
    assert(t == R);

    REncodedEltQapArith* ae = (REncodedEltQapArith*)a;
	  REncodedEltQapArith* be = (REncodedEltQapArith*)b;
    
    ret = 0 != QapFp2_weierstrass_equal_aff(ae->me, be->me);
  }

  return ret;
}

void EncodingArithQap::zero(elt_t t, EncodedElt* a) { 
	if (t == L) {
	  copy(t, &members->Lzero, a);	// Use our precomputed value	
	} else {
		assert(t == R);
	  copy(t, &members->Rzero, a);	// Use our precomputed value	
	}
}

void EncodingArithQap::one(elt_t t, EncodedElt* a) {
	if (t == L) {
	  copy(t, &members->Lg, a);	// Use our precomputed value	
	} else {
		assert(t == R);
	  copy(t, &members->Rg, a);	// Use our precomputed value	
	}
}

// Adds two encoded elements.  Result goes in r.
void EncodingArithQap::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
  add(t, a, b, r, false);
}

// Adds two encoded elements.  Result goes in r.
void EncodingArithQap::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, bool definitely_distinct) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);

	if (t == L) {
    LEncodedEltQapArith* ae = (LEncodedEltQapArith*)a;
    LEncodedEltQapArith* be = (LEncodedEltQapArith*)b;
    LEncodedEltQapArith* re = (LEncodedEltQapArith*)r;
    
    if (isZero(t, a)) {
      QapFp_wp_copy(re->me, be->me);
    } else if (isZero(t, b)) {
      QapFp_wp_copy(re->me, ae->me);    
    } else if (!definitely_distinct && this->equal(t, a, b)) {
      // Library doesn't support addition of a point to itself
      QapFp_weierstrass_affdbl(re->me, ae->me);
    } else {  
      QapFp_weierstrass_affadd(re->me, ae->me, be->me);
    }
	} else {
		assert(t == R);
    REncodedEltQapArith* ae = (REncodedEltQapArith*)a;
    REncodedEltQapArith* be = (REncodedEltQapArith*)b;
    REncodedEltQapArith* re = (REncodedEltQapArith*)r;

    if (isZero(t, a)) {
      QapFp2_wp_copy(re->me, be->me);
    } else if (isZero(t, b)) {
      QapFp2_wp_copy(re->me, ae->me);    
    } else if (!definitely_distinct && this->equal(t, a, b)) {
      // Library doesn't support addition of a point to itself
      QapFp2_weierstrass_affdbl(re->me, ae->me);
    } else {
      QapFp2_weierstrass_affadd(re->me, ae->me, be->me);
    }
	}

	testOnCurve(t, r);
}

// Subtracts two encoded elements.  Result goes in r.
void EncodingArithQap::sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);

	if (t == L) {
    LEncodedEltQapArith* ae = (LEncodedEltQapArith*)a;
    LEncodedEltQapArith* be = (LEncodedEltQapArith*)b;
    LEncodedEltQapArith* re = (LEncodedEltQapArith*)r;

    LEncodedEltQapArith neg_b;
    QapFp_wp_init(neg_b.me);
    QapFp_wp_neg(neg_b.me, be->me);
//	printf("\n");
//	printf("b:"); print(t, be);
//	printf("neg_b:");  print(t, &neg_b);
    add(t, ae, &neg_b, re);
//	printf("a:"); print(t, ae);
//	printf("result:"); print(t, re);
    QapFp_wp_free(neg_b.me);
	} else {
		assert(t == R);
    REncodedEltQapArith* ae = (REncodedEltQapArith*)a;
    REncodedEltQapArith* be = (REncodedEltQapArith*)b;
    REncodedEltQapArith* re = (REncodedEltQapArith*)r;

    REncodedEltQapArith neg_b;
    QapFp2_wp_init(neg_b.me);
    QapFp2_wp_neg(neg_b.me, be->me);
    add(t, ae, &neg_b, re);
    QapFp2_wp_free(neg_b.me);
	}

	testOnCurve(t, r);
}

void EncodingArithQap::doubleIt(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	testOnCurve(t, a);

	if (t == L) {
    LEncodedEltQapArith* ae = (LEncodedEltQapArith*)a;
    LEncodedEltQapArith* be = (LEncodedEltQapArith*)b;
    if (isZero(L, ae)) {
      zero(L, be); 
    } else {
      QapFp_weierstrass_affdbl(be->me, ae->me);
    }
	} else {
		assert(t == R);
    REncodedEltQapArith* ae = (REncodedEltQapArith*)a;
    REncodedEltQapArith* be = (REncodedEltQapArith*)b;
    if (isZero(R, ae)) {
      zero(R, be);
    } else {
      QapFp2_weierstrass_affdbl(be->me, ae->me);
    }
	}

	testOnCurve(t, b);
}


// Multiplies the value encoded in g by the constant c.  Result goes in r.
// For pairings, we compute r <- g^c
void EncodingArithQap::mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r) {		
	assert(g && c && r);
	testOnCurve(t, g);

  EncodedElt* result;
  if (g == r) { // Don't want to clobber our input
    result = new_elt(t);
  } else {
    result = r;
  }
	
  // Naive square and multiply
  zero(t, result);
  for (int bit_index = num_exponent_bits(); bit_index >= 0; bit_index--) {
    doubleIt(t, result, result);
    if (srcField->getBits(c, bit_index, 1)) {
      add(t, g, result, result, true);
    }
  }

  if (g == r) { 
    // Copy over the result
    copy(t, result, r);
    del_elt(t, result);
  }

	testOnCurve(t, r);
}

void EncodingArithQap::mul(LEncodedElt* a, REncodedElt* b, EncodedProduct* r) {
	pair(a, b, r);
}

void EncodingArithQap::pair(LEncodedElt* L1, REncodedElt* R1, EncodedProduct* result) {
	assert(L1 && R1 && result);
	testOnCurve(L, L1);
	testOnCurve(R, R1);

	LEncodedEltQapArith* L1e = (LEncodedEltQapArith*)L1;
	REncodedEltQapArith* R1e = (REncodedEltQapArith*)R1;
	EncodedProductQapArith* resulte = (EncodedProductQapArith*)result;

  if (isZero(L, L1) || isZero(R, R1)) {
    // Library doesn't handle g^0 values properly, so we special-case it
    QapFp12_set_one(resulte->me);
  } else {
    Qapbn_optimal_ate_miller_aff(resulte->me, R1e->me, L1e->me);   
    Qapbn_finalexpo_neg(resulte->me, resulte->me);
  }

}

bool EncodingArithQap::isZero(elt_t t, EncodedElt* elt) {
	assert(elt);
	//testOnCurve(t, elt);

  bool ret = false;
	if (t == L) {
	  LEncodedEltQapArith* elte = (LEncodedEltQapArith*)elt;
    ret = 0 == QapFp_cmp(elte->me->Y, members->Lzero.me->Y);  // Our internal representation of 0
	} else {
		assert(t == R);
	  REncodedEltQapArith* elte = (REncodedEltQapArith*)elt;
    ret = 0 != QapFp2_cmpeq(elte->me->Y, members->Rzero.me->Y);   // Check if the Y coord is 0 ==> pt at infinity
	}

  return ret;
}


void EncodingArithQap::add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
  EncodedProductQapArith* ae = (EncodedProductQapArith*)a;
  EncodedProductQapArith* be = (EncodedProductQapArith*)b;
  EncodedProductQapArith* re = (EncodedProductQapArith*)r;

  QapFp12_mul(re->me, ae->me, be->me);
}

void EncodingArithQap::sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductQapArith* ae = (EncodedProductQapArith*)a;
	EncodedProductQapArith* be = (EncodedProductQapArith*)b;
	EncodedProductQapArith* re = (EncodedProductQapArith*)r;
  
	EncodedProductQapArith* inv = (EncodedProductQapArith*)this->new_prod();
	QapFp12_inv(inv->me, be->me);
	QapFp12_mul(re->me, ae->me, inv->me);
}

bool EncodingArithQap::equals(EncodedProduct* a, EncodedProduct* b) {
	assert(a && b);
	EncodedProductQapArith* ae = (EncodedProductQapArith*)a;
	EncodedProductQapArith* be = (EncodedProductQapArith*)b;

  return 1 == QapFp12_cmpeq(ae->me, be->me);
}

// Checks whether the plaintexts inside the encodings obey:
// L1*R1 - L2*R2 == L3*R3
bool EncodingArithQap::mulSubEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2, LEncodedElt* L3, REncodedElt* R3) {
	assert(L1 && R1 && L2 && R2 && L3 && R3);

	EncodedProductQapArith pairing1;
	EncodedProductQapArith pairing2;
	EncodedProductQapArith pairing3;

  QapFp12_init(pairing1.me);
  QapFp12_init(pairing2.me);
  QapFp12_init(pairing3.me);

	testOnCurve(L, L1); testOnCurve(R, R1);
	testOnCurve(L, L2); testOnCurve(R, R2);
	testOnCurve(L, L3); testOnCurve(R, R3);

	// Compute the three pairings
	pair(L1, R1, &pairing1);
	pair(L2, R2, &pairing2);
	pair(L3, R3, &pairing3);

	// Compute (with some notation abuse) L1*R1 - L2*R2 (in the exponents), i.e., pairing1 <- pairing1 / pairing2
	sub(&pairing1, &pairing2, &pairing1);
	bool result = equals(&pairing1, &pairing3);

  QapFp12_free(pairing1.me);
  QapFp12_free(pairing2.me);
  QapFp12_free(pairing3.me);

	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}


// Checks whether the plaintexts inside the encodings obey:
// L1*R1 == L2*R2
bool EncodingArithQap::mulEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2) {
	assert(L1 && R1 && L2 && R2);
	EncodedProductQapArith pairing1;
	EncodedProductQapArith pairing2;
  QapFp12_init(pairing1.me);
  QapFp12_init(pairing2.me);

	testOnCurve(L, L1); testOnCurve(R, R1);
	testOnCurve(L, L2); testOnCurve(R, R2);

	// Compute the two pairings
	pair(L1, R1, &pairing1);
	pair(L2, R2, &pairing2);

	BOOL result = equals(&pairing1, &pairing2); 

  QapFp12_free(pairing1.me);
  QapFp12_free(pairing2.me);
	
	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}

void EncodingArithQap::compress(elt_t t, EncodedElt* elt) {
  assert(false);		// Not yet implemented
}

void EncodingArithQap::decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt) {
  assert(false);		// Not yet implemented
}

void EncodingArithQap::write(Archiver* arc, elt_t t, EncodedElt* e, bool simple) {
	assert(false);		// Not yet implemented
}

void EncodingArithQap::read (Archiver* arc, elt_t t, EncodedElt* e) {
	assert(false);		// Not yet implemented
}

void EncodingArithQap::affinize(elt_t t, EncodedElt* in, EncodedElt* out) {
	assert(in && out);
	testOnCurve(t, in);

  copy(t, in, out);
}


void EncodingArithQap::testOnCurve(elt_t t, EncodedElt* a) {
#ifdef DEBUG_EC
	assert(a);

	if (t == L) {
    LEncodedEltQapArith* ae = (LEncodedEltQapArith*)a;
    QapFp_weierstrass_oncurve_aff(ae->me);
	} else {
		assert(t == R);
    REncodedEltQapArith* ae = (REncodedEltQapArith*)a;
    QapFp2_weierstrass_oncurve_aff(ae->me);
	}
#endif
}

void EncodingArithQap::encodeSlow(elt_t t, FieldElt* in, EncodedElt* out) {
	assert(in && out);

	if (t == L) {
    mul(t, &members->Lg, in, out);
	} else {
		assert(t == R);
    mul(t, &members->Rg, in, out);
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//			Arrays
//////////////////////////////////////////////////////////////////////////////////////////////////////


EncodedEltQapArithArray::EncodedEltQapArithArray(elt_t type, Encoding* encoding, uint32_t len, bool preallocate) {
	this->type = type;
	this->len = len;
	elts = new EncodedEltQapArith*[len];
	allocated = preallocate;
	
	for (uint32_t i = 0; i < len; i++) {
		if (preallocate) {
			elts[i] = (EncodedEltQapArith*)encoding->new_elt(type);
		} else {
			elts[i] = NULL;
		}
	}
}

EncodedEltQapArithArray::~EncodedEltQapArithArray() {
	if (allocated) {
		for (uint32_t i = 0; i < len; i++) {
			delete elts[i];
		}
	}
	delete [] elts;
}

EncodedElt* EncodedEltQapArithArray::elt(uint32_t index) {
	assert(index < len);
	assert(elts[index] != NULL);
	return elts[index];
}

void EncodedEltQapArithArray::set(uint32_t index, EncodedElt* elt) {
	assert(index < len);
	assert(elts[index] == NULL);	// Otherwise we're overwriting and existing value!
	elts[index] = (EncodedEltQapArith*)elt;
}

void EncodedEltQapArithArray::reset() {
  assert(!allocated);
  for (uint32_t i = 0; i < len; i++) {
    elts[i] = NULL;
  }
}
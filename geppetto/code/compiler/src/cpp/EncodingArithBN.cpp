#include <assert.h>
#include "EncodingArithBN.h"
#include <bitset>
#include <intrin.h>
#include "Poly.h"
#include "FieldPrime.h"
#include "EncodingEnigmaBN.h"
#include "ArithWrappers.h"

// Protect ourselves from ARITH's #define
#ifdef mul
#undef mul
#endif

EncodingArithBN::EncodingArithBN(bool debug) {
  assert(PAIR_CURVE == 0);    // Make sure no one else has set it.  
	members = new EncodingArithBNMembers;      
  this->debug = debug;

  if (!debug) { // Use the real curve
    PAIR_CURVE = BN12;
    members->n = 4;
    members->m = (uint64_t *)malloc(members->n * sizeof (uint64_t));
    members->b0 = (uint64_t *)malloc(members->n * sizeof (uint64_t));
    members->b1 = (uint64_t *)malloc(members->n * sizeof (uint64_t));
    members->m[3] = 0x2523648240000001;
    members->m[2] = 0xBA344D8000000008;
    members->m[1] = 0x6121000000000013;
    members->m[0] = 0xA700000000000013; // 2523648240000001BA344D80000000086121000000000013A700000000000013
    /* This is the 254-bit BN prime p = 3 (mod 4). */

    Fp_initialize_config(members->m, members->n);
    Fp_init(members->b);
    Fp_wp_init(members->Lg.me);
    //  Fp_wp_init (members->Q0);

    Fp2_initialize_config();
    Fp4_initialize_config();
    Fp6_initialize_config();
    Fp12_initialize_config();

    Fp2_init(members->B);
    Fp2_init(members->B3);
    Fp2_init(members->l00);
    Fp2_init(members->l01);
    Fp2_init(members->l10);

    Fp_set_ui(members->b, 2);

    // New generator that matches Enigma's
    members->m[3] = 0x2523648240000001;
    members->m[2] = 0xBA344D8000000008;
    members->m[1] = 0x6121000000000013;
    members->m[0] = 0xA700000000000012;
    Fp_set(members->Lg.me->X, members->m);

    members->m[3] = 0x0000000000000000;
    members->m[2] = 0x0000000000000000;
    members->m[1] = 0x0000000000000000;
    members->m[0] = 0x0000000000000001;
    Fp_set(members->Lg.me->Y, members->m);

    members->b1[3] = 0x2523648240000001;
    members->b1[2] = 0xBA344D8000000008;
    members->b1[1] = 0x6121000000000013;
    members->b1[0] = 0xA700000000000012;

    members->b0[3] = 0x0000000000000000;
    members->b0[2] = 0x0000000000000000;
    members->b0[1] = 0x0000000000000000;
    members->b0[0] = 0x0000000000000001;
    Fp2_set(members->B, members->b0, members->b1);

    Fp_initialize_weierstrass(&members->b, 1);
    Fp2_initialize_weierstrass(members->B);
    Fp2_wp_init(members->Rg.me);
    bn_bls12_miller_initialize_config(members->B);
    bn_bls12_finalexpo_initialize_config();

    Fp2_copy(members->B3, bn_config.bt3);

    members->b1[3] = 0x0516AAF9BA737833;
    members->b1[2] = 0x310AA78C5982AA5B;
    members->b1[1] = 0x1F4D746BAE3784B7;
    members->b1[0] = 0x0D8C34C1E7D54CF3;

    members->b0[3] = 0x061A10BB519EB62F;
    members->b0[2] = 0xEB8D8C7E8C61EDB6;
    members->b0[1] = 0xA4648BBB4898BF0D;
    members->b0[0] = 0x91EE4224C803FB2B;

    Fp2_set(members->Rg.me->X, members->b0, members->b1);

    members->b1[3] = 0x0EBB2B0E7C8B1526;
    members->b1[2] = 0x8F6D4456F5F38D37;
    members->b1[1] = 0xB09006FFD739C957;
    members->b1[0] = 0x8A2D1AEC6B3ACE9B;

    members->b0[3] = 0x021897A06BAF9343;
    members->b0[2] = 0x9A90E096698C8223;
    members->b0[1] = 0x29BD0AE6BDBE09BD;
    members->b0[0] = 0x19F0E07891CD2B9A;

    Fp2_set(members->Rg.me->Y, members->b0, members->b1);

    // For now, use Enigma's encoding to create the prime field we want
    members->enigma_encoding = new EncodingEnigmaBN();
    this->srcField = members->enigma_encoding->getSrcField();
  } else { // Use the small debugging curve
    PAIR_CURVE = BN12tiny;
    members->n = 1;
    members->m = (uint64_t *)malloc(members->n * sizeof (uint64_t));
    members->b0 = (uint64_t *)malloc(members->n * sizeof (uint64_t));
    members->b1 = (uint64_t *)malloc(members->n * sizeof (uint64_t));
    members->m[0] = 0x15C9B0796B71DB;
    /* This is the 53-bit BN prime p = 3 (mod 4). */

    Fp_initialize_config(members->m, members->n);
    Fp_init(members->b);
    Fp_wp_init(members->Lg.me);

    Fp2_initialize_config();
    Fp4_initialize_config();
    Fp6_initialize_config();
    Fp12_initialize_config();

    Fp2_init(members->B);
    Fp2_init(members->B3);
    Fp2_init(members->l00);
    Fp2_init(members->l01);
    Fp2_init(members->l10);

    Fp_set_ui(members->b, 2);

    members->m[0] = 0x321347FC84CEE;
    Fp_set(members->Lg.me->X, members->m);

    members->m[0] = 0xC5D3F93EEC09B;
    Fp_set(members->Lg.me->Y, members->m);

    members->b1[0] = 0x15C9B0796B71DA;
    members->b0[0] = 0x0000000000000001;
    Fp2_set(members->B, members->b0, members->b1);

    Fp_initialize_weierstrass(&members->b, 1);
    Fp2_initialize_weierstrass(members->B);
    Fp2_wp_init(members->Rg.me);
    bn_bls12_miller_initialize_config(members->B);
    bn_bls12_finalexpo_initialize_config();

    Fp2_copy(members->B3, bn_config.bt3);

    members->b1[0] = 0x11B4D9B9078C9C;
    members->b0[0] = 0x1183E42BE63029;
    Fp2_set(members->Rg.me->X, members->b0, members->b1);

    members->b1[0] = 0x44923E174AB60;
    members->b0[0] = 0x10B2AC7252D760;

    Fp2_set(members->Rg.me->Y, members->b0, members->b1);

    // For now, use Enigma's encoding to create the prime field we want
    members->enigma_encoding = new EncodingEnigmaBN();
    this->srcField = ((EncodingEnigmaBN*)members->enigma_encoding)->get_tiny_bn_field();
  }

  // Convert the generators to projective form
  Fp_weierstrass_aff2proj (members->Lg.me, members->Lg.me);
  Fp2_weierstrass_aff2proj (members->Rg.me, members->Rg.me);

	// Create representations of zero, since we use it a lot
  Fp_wp_init(members->Lzero.me);
  Fp_copy(members->Lzero.me->X, members->Lg.me->X);
  Fp_copy(members->Lzero.me->Y, members->Lg.me->Y);
  memset(members->Lzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

  Fp2_wp_init(members->Rzero.me);
  Fp2_copy(members->Rzero.me->X, members->Rg.me->X);
  Fp2_copy(members->Rzero.me->Y, members->Rg.me->Y);
  Fp2_set_zero(members->Rzero.me->Z);
}

EncodingArithBN::~EncodingArithBN() {
	delete members->enigma_encoding;	// Delete the encoding that provided our srcField, not the srcField itself
	
  Fp_wp_free(members->Lzero.me);
  Fp2_wp_free(members->Rzero.me);

  Fp_free (members->b);
  Fp2_free (members->B);
  Fp2_free (members->B3);
  Fp2_free (members->l00);
  Fp2_free (members->l01);
  Fp2_free (members->l10);
  Fp2_wp_free (members->Rg.me);
  Fp_wp_free (members->Lg.me);
  bn_bls12_finalexpo_free_config ();
  bn_bls12_miller_free_config ();
  Fp12_free_config (); 
  Fp6_free_config ();
  Fp4_free_config ();
  Fp_free_weierstrass();
  Fp2_free_weierstrass ();
  Fp2_free_config (); 
  Fp_free_config ();

  free(members->m );
  free(members->b0);
  free(members->b1);

  delete members;
}

// TODO: Adjust to count pointers as well
uint32_t EncodingArithBN::bytes_per(elt_t t) {
	if (t == L) {
		return num_L_digits*sizeof(digit_t);
	} else {
		assert(t == R);
		return num_R_digits*sizeof(digit_t);
	}
}

uint32_t EncodingArithBN::num_exponent_bits() {
  if (!debug) {
	  return 4*sizeof(digit_t)*8 - 1;	// BN_curve uses 4 digits, but doesn't use them all (only uses 254 bits)
  } else {
    return 54;
  }
}


EncodedElt* EncodingArithBN::new_elt(elt_t t) {
	if (t == L) {
    LEncodedEltArithBN* elt = new LEncodedEltArithBN;
    Fp_wp_init(elt->me);
		return elt;
	} else {
		assert(t == R);
		REncodedEltArithBN* elt = new REncodedEltArithBN;
    Fp2_wp_init(elt->me);
		return elt;
	}
}

EncodedEltArray* EncodingArithBN::new_elt_array(elt_t t, uint32_t len, bool preallocate) {
	return new EncodedEltArithBNArray(t, this, len, preallocate);
}


void EncodingArithBN::del_elt(elt_t t, EncodedElt* elt) {
  assert(elt);
  if (t == L) {
    Fp_wp_free(((LEncodedEltArithBN*)elt)->me);
  } else {
    assert(t == R);
    Fp2_wp_free(((REncodedEltArithBN*)elt)->me);
  }
	delete elt;
}

EncodedProduct* EncodingArithBN::new_prod() {
	EncodedProductArithBN* ret = new EncodedProductArithBN();
  Fp12_init(ret->me);
  return ret;
}

EncodedProductArithBN::~EncodedProductArithBN() {
  Fp12_free(me);
}

uint32_t EncodingArithBN::size_of_prod() {
  return sizeof(Fp_t);
}

void EncodingArithBN::print(elt_t t, EncodedElt* e) {
	if (t == L) {
		LEncodedEltArithBN* en = (LEncodedEltArithBN*)e;
		Fp_weierstrass_jacproj2aff(en->me, en->me);
		printf("( ");
		Fp_print(en->me->X);
		printf(", ");
		Fp_print(en->me->Y);
		printf(", ");
		Fp_print(en->me->Z);
		printf(")");
	} else {
		assert(t == R);
		REncodedEltArithBN* en = (REncodedEltArithBN*)e;
		Fp2_weierstrass_jacproj2aff(en->me, en->me);
		printf("( ");
		Fp2_print(en->me->X);
		printf(", ");
		Fp2_print(en->me->Y);
		printf(", ");
		Fp2_print(en->me->Z);
		printf(")");
	}
}

void EncodingArithBN::print(EncodedProduct* p) {
	EncodedProductArithBN* ep = (EncodedProductArithBN*)p;
	Fp12_print(ep->me);
}

void EncodingArithBN::copy(elt_t t, EncodedElt* src, EncodedElt* dst) {
	if (t == L) {
    Fp_wp_copy(((LEncodedEltArithBN*)dst)->me, ((LEncodedEltArithBN*)src)->me);
	} else {
		assert(t == R);
		Fp2_wp_copy(((REncodedEltArithBN*)dst)->me, ((REncodedEltArithBN*)src)->me);
	}	
}

bool EncodingArithBN::equal(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	testOnCurve(t, a);
	testOnCurve(t, b);

  bool ret = false;
  if (t == L) {
    ret = 0 != Fp_weierstrass_equal_jacproj(((LEncodedEltArithBN*)a)->me, ((LEncodedEltArithBN*)b)->me);
  } else {
    assert(t == R);
    ret = 0 != Fp2_weierstrass_equal_jacproj(((REncodedEltArithBN*)a)->me, ((REncodedEltArithBN*)b)->me);
  }
 
  return ret;
}

void EncodingArithBN::zero(elt_t t, EncodedElt* a) { 
	if (t == L) {
	  copy(t, &members->Lzero, a);	// Use our precomputed value	
	} else {
		assert(t == R);
	  copy(t, &members->Rzero, a);	// Use our precomputed value	
	}
}

void EncodingArithBN::one(elt_t t, EncodedElt* a) {
	if (t == L) {
	  copy(t, &members->Lg, a);	// Use our precomputed value	
	} else {
		assert(t == R);
	  copy(t, &members->Rg, a);	// Use our precomputed value	
	}
}

// Adds two encoded elements.  Result goes in r.
void EncodingArithBN::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
  add(t, a, b, r, false);
}

// Adds two encoded elements.  Result goes in r.
void EncodingArithBN::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, bool definitely_distinct) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);

	if (t == L) {
    LEncodedEltArithBN* ae = (LEncodedEltArithBN*)a;
    LEncodedEltArithBN* be = (LEncodedEltArithBN*)b;
    LEncodedEltArithBN* re = (LEncodedEltArithBN*)r;
    
    if (isZero(t, a)) {
      Fp_wp_copy(re->me, be->me);
    } else if (isZero(t, b)) {
      Fp_wp_copy(re->me, ae->me);    
    } else if (!definitely_distinct && this->equal(t, a, b)) {
      // Library doesn't support addition of a point to itself
      Fp_weierstrass_dbl(re->me, ae->me);
    } else {  
      Fp_weierstrass_add(re->me, ae->me, be->me);      
    }
	} else {
		assert(t == R);
    REncodedEltArithBN* ae = (REncodedEltArithBN*)a;
    REncodedEltArithBN* be = (REncodedEltArithBN*)b;
    REncodedEltArithBN* re = (REncodedEltArithBN*)r;

    if (isZero(t, a)) {
      Fp2_wp_copy(re->me, be->me);
    } else if (isZero(t, b)) {
      Fp2_wp_copy(re->me, ae->me);    
    } else if (!definitely_distinct && this->equal(t, a, b)) {
      // Library doesn't support addition of a point to itself
      Fp2_weierstrass_dbl(re->me, ae->me);
    } else {      
      Fp2_weierstrass_add(re->me, ae->me, be->me);
    }
	}

	testOnCurve(t, r);
}

// Subtracts two encoded elements.  Result goes in r.
void EncodingArithBN::sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);

	if (t == L) {
    LEncodedEltArithBN* ae = (LEncodedEltArithBN*)a;
    LEncodedEltArithBN* be = (LEncodedEltArithBN*)b;
    LEncodedEltArithBN* re = (LEncodedEltArithBN*)r;

    LEncodedEltArithBN neg_b;
    Fp_wp_init(neg_b.me);
    Fp_wp_neg(neg_b.me, be->me);
    add(t, ae, &neg_b, re);
    Fp_wp_free(neg_b.me);
	} else {
		assert(t == R);
    REncodedEltArithBN* ae = (REncodedEltArithBN*)a;
    REncodedEltArithBN* be = (REncodedEltArithBN*)b;
    REncodedEltArithBN* re = (REncodedEltArithBN*)r;

    REncodedEltArithBN neg_b;
    Fp2_wp_init(neg_b.me);
    Fp2_wp_neg(neg_b.me, be->me);
    add(t, ae, &neg_b, re);
    Fp2_wp_free(neg_b.me);
	}

	testOnCurve(t, r);
}

void EncodingArithBN::doubleIt(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	testOnCurve(t, a);

	if (t == L) {
    LEncodedEltArithBN* ae = (LEncodedEltArithBN*)a;
    LEncodedEltArithBN* be = (LEncodedEltArithBN*)b;

    Fp_weierstrass_dbl(be->me, ae->me);
	} else {
		assert(t == R);
    REncodedEltArithBN* ae = (REncodedEltArithBN*)a;
    REncodedEltArithBN* be = (REncodedEltArithBN*)b;

    Fp2_weierstrass_dbl(be->me, ae->me);
	}

	testOnCurve(t, b);
}


// Multiplies the value encoded in g by the constant c.  Result goes in r.
// For pairings, we compute r <- g^c
void EncodingArithBN::mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r) {		
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
    //if (isZero(t, r)) {
    //  printf("After index %d, r is still 0\n", bit_index);
    //}
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

void EncodingArithBN::mul(LEncodedElt* a, REncodedElt* b, EncodedProduct* r) {
	pair(a, b, r);
}

void EncodingArithBN::pair(LEncodedElt* L1, REncodedElt* R1, EncodedProduct* result) {
	assert(L1 && R1 && result);
	testOnCurve(L, L1);
	testOnCurve(R, R1);

	LEncodedEltArithBN* L1e = (LEncodedEltArithBN*)L1;
	REncodedEltArithBN* R1e = (REncodedEltArithBN*)R1;
	EncodedProductArithBN* resulte = (EncodedProductArithBN*)result;

  if (isZero(L, L1) || isZero(R, R1)) {
    // Library doesn't handle g^0 values properly, so we special-case it
    Fp12_set_one(resulte->me);
  } else {
    Point_wp_t aff_L1; 
    Point_wp2_t aff_R1;
    Fp_wp_init(aff_L1);
    Fp2_wp_init(aff_R1);

    // Convert to affine representation before doing the pairing
    Fp_weierstrass_jacproj2aff(aff_L1, L1e->me);
    Fp2_weierstrass_jacproj2aff(aff_R1, R1e->me);

    bn_optimal_ate_miller(resulte->me, aff_R1, aff_L1);   
    bn_finalexpo_neg(resulte->me, resulte->me);

    Fp_wp_free(aff_L1);
    Fp2_wp_free(aff_R1);
  }

}

bool EncodingArithBN::isZero(elt_t t, EncodedElt* elt) {
	assert(elt);

  bool ret = false;
	if (t == L) {
	  LEncodedEltArithBN* elte = (LEncodedEltArithBN*)elt;
    ret = 0 == Fp_cmp(elte->me->Z, members->Lzero.me->Z);   // Check if the Z coord is 0 ==> pt at infinity
	} else {
		assert(t == R);
	  REncodedEltArithBN* elte = (REncodedEltArithBN*)elt;
    ret = 0 != Fp2_cmpeq(elte->me->Z, members->Rzero.me->Z);   // Check if the Z coord is 0 ==> pt at infinity
	}

  return ret;
}


void EncodingArithBN::add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
  EncodedProductArithBN* ae = (EncodedProductArithBN*)a;
  EncodedProductArithBN* be = (EncodedProductArithBN*)b;
  EncodedProductArithBN* re = (EncodedProductArithBN*)r;

  Fp12_mul(re->me, ae->me, be->me);
}

void EncodingArithBN::sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductArithBN* ae = (EncodedProductArithBN*)a;
	EncodedProductArithBN* be = (EncodedProductArithBN*)b;
	EncodedProductArithBN* re = (EncodedProductArithBN*)r;
  
	EncodedProductArithBN* inv = (EncodedProductArithBN*)this->new_prod();
	Fp12_inv(inv->me, be->me);
	Fp12_mul(re->me, ae->me, inv->me);
}

bool EncodingArithBN::equals(EncodedProduct* a, EncodedProduct* b) {
	assert(a && b);
	EncodedProductArithBN* ae = (EncodedProductArithBN*)a;
	EncodedProductArithBN* be = (EncodedProductArithBN*)b;

  return 1 == Fp12_cmpeq(ae->me, be->me);
}

// Checks whether the plaintexts inside the encodings obey:
// L1*R1 - L2*R2 == L3*R3
bool EncodingArithBN::mulSubEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2, LEncodedElt* L3, REncodedElt* R3) {
	assert(L1 && R1 && L2 && R2 && L3 && R3);

	EncodedProductArithBN pairing1;
	EncodedProductArithBN pairing2;
	EncodedProductArithBN pairing3;

  Fp12_init(pairing1.me);
  Fp12_init(pairing2.me);
  Fp12_init(pairing3.me);

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

  Fp12_free(pairing1.me);
  Fp12_free(pairing2.me);
  Fp12_free(pairing3.me);

	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}


// Checks whether the plaintexts inside the encodings obey:
// L1*R1 == L2*R2
bool EncodingArithBN::mulEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2) {
	assert(L1 && R1 && L2 && R2);
	EncodedProductArithBN pairing1;
	EncodedProductArithBN pairing2;
  Fp12_init(pairing1.me);
  Fp12_init(pairing2.me);

	testOnCurve(L, L1); testOnCurve(R, R1);
	testOnCurve(L, L2); testOnCurve(R, R2);

	// Compute the two pairings
	pair(L1, R1, &pairing1);
	pair(L2, R2, &pairing2);

	BOOL result = equals(&pairing1, &pairing2); 

  Fp12_free(pairing1.me);
  Fp12_free(pairing2.me);
	
	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}

void write_zero(Archiver* arc, int num_digits) {
  size_t size = num_digits * sizeof(digit_t);
  unsigned char* zeroes = new unsigned char[size];
  memset(zeroes, 0, size);
  arc->write(zeroes, (int)size);
  delete[] zeroes;
}

void write_Fp(Archiver* arc, Fp_t* fp, int num_digits) {
  digit_t buffer[4];  
  assert(0 <= num_digits && num_digits <= 4);
  assert(arc && fp);
  Fp_get(buffer, *fp);
	arc->write((unsigned char*)buffer, num_digits * sizeof(digit_t));
}

void read_Fp(Archiver* arc, Fp_t* fp, int num_digits) {
  digit_t buffer[4];
  assert(0 <= num_digits && num_digits <= 4);
  assert(arc && fp);
	arc->read((unsigned char*)buffer, num_digits * sizeof(digit_t));
  Fp_set(*fp, buffer);
}

void write_Fp2(Archiver* arc, Fp2_t* fp2, int num_digits) {
	write_Fp(arc, &fp2[0]->a0, num_digits);
	write_Fp(arc, &fp2[0]->a1, num_digits);
}

void read_Fp2(Archiver* arc, Fp2_t* fp2, int num_digits) {
	read_Fp(arc, &fp2[0]->a0, num_digits);
	read_Fp(arc, &fp2[0]->a1, num_digits);
}

// Defined in EncodingArithCP.cpp
extern void compressAff(bool elt_is_zero, Point_wp_t aff_in, Fp_t out, Fp_t b);

// On the twist curve, we use y^2 = x^3 + (1 - i), and each point is represented as: a + i*b
void compressRaff(bool elt_is_zero, Point_wp2_t aff_in, Fp2_t out) {
  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (elt_is_zero) {
    Fp2_set_zero(out);    
    return;
  }

  // Compute x^3 + (1 - i)
  Fp2_t xCubed;
  Fp2_mul(xCubed, aff_in->X, aff_in->X);
  Fp2_mul(xCubed, xCubed, aff_in->X);

  Fp2_t oneNegOne;
  Fp_set_ui(oneNegOne->a0, 1);  // Put 1 into first slot
  Fp_set_ui(oneNegOne->a1, 0);
  Fp_sub(oneNegOne->a1, oneNegOne->a1, oneNegOne->a0); // Put -1 into the second slot

  Fp2_add(xCubed, xCubed, oneNegOne);

  // Compute a root
  Fp2_t sqrtRhs;
  assert(false);  // Need to find an implementation of sqrt
  // sqrtRhs = Fp2_sqrt(xCubed);

#ifdef DEBUG_EC
  // Check the square root is correct
  Fp2_t allegedRhs;
  Fp2_mul(allegedRhs, sqrtRhs, sqrtRhs);
  assert(Fp2_cmpeq(allegedRhs, aff_in->X)); 
#endif

  // Compute the other root
  Fp2_t negSqrtRhs;
  Fp2_neg(negSqrtRhs, sqrtRhs);

  bool equal = 0 != Fp2_cmpeq(sqrtRhs, aff_in->Y);
  
  // Since we can't actually test for negative roots, and the library gives us a random one,
  // we distinguish between them based on which last digit is smaller
  bool equalsSmaller = (equal && (sqrtRhs->a1->limbs[3] <= negSqrtRhs->a1->limbs[3])) || (!equal && (sqrtRhs->a1->limbs[3] > negSqrtRhs->a1->limbs[3]));

  // We use the top-most bit of X to store the bit indicating whether to use the smaller root
  int lastDigit = 8 / 2 - 1;	// EncodedElt contains 2 field elements
  out = aff_in->X;
  if (equalsSmaller) {
    // Set the top bit of the last digit of X to indicate that we use the smaller root
    out->a0->limbs[lastDigit] |= ((digit_t)1) << 63;
  } else {
    // Make sure the top bit is actually clear, indicating the larger root
    assert(!(out->a0->limbs[lastDigit] & (((digit_t)1) << 63)));
  }
}

// Warning: This hasn't been tested yet!
void EncodingArithBN::compress(elt_t t, EncodedElt* e) {
  assert(e);  
  bool elt_is_zero = isZero(t, e);
  affinize(t, e, e);
  if (t == L) {
    LEncodedEltArithBN* elt = (LEncodedEltArithBN*)e;
    //compressLaff(elt_is_zero, elt->me, elt->me->X);
    compressAff(elt_is_zero, elt->me, elt->me->X, members->b);
  } else {
    assert(t == R);
    REncodedEltArithBN* elt = (REncodedEltArithBN*)e;
    compressRaff(elt_is_zero, elt->me, elt->me->X);
  }
}

// Defined in EncodingArithCP.cpp
extern void decompressAff(Fp_t in, Point_wp_t out, Fp_t b, int digits_in_fp_t, Point_wp_t zero);
extern bool allZero(digit_t* value, int num_digits);

// On the twist curve, we use y^2 = x^3 + (1 - i), and each point is represented as a + i*b
void decompressAffR(Fp2_t in, Point_wp2_t out, int digits_in_fp2_t, Point_wp2_t zero) {
  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (allZero((digit_t*)in, digits_in_fp2_t)) {
    Fp2_wp_copy(out, zero);    
    return;
  }

  // Check whether the top bit of X is set.  Indicates which root to use for Y
  int lastDigit = digits_in_fp2_t / 2 - 1;	// EncodedElt contains 2 field elements
  bool useSmallerRoot = (in->a0->limbs[lastDigit] & (((digit_t)1) << 63)) > 0;

  // Either way, clear out the top bit and copy the input into the x coefficient
  in->a0->limbs[lastDigit] &= ((((digit_t)1) << 63) - 1);
  Fp2_copy(out->X, in);

  // Compute x^3 + (1 - i)
  Fp2_t xCubed;
  Fp2_mul(xCubed, out->X, out->X);
  Fp2_mul(xCubed, xCubed, out->X);

  Fp2_t oneNegOne;
  Fp_set_ui(oneNegOne->a0, 1);  // Put 1 into first slot
  Fp_set_ui(oneNegOne->a1, 0);
  Fp_sub(oneNegOne->a1, oneNegOne->a1, oneNegOne->a0); // Put -1 into the second slot

  Fp2_add(xCubed, xCubed, oneNegOne);

  // Compute a root
  assert(false);  // Need to find an implementation of sqrt
  // out->Y = Fp2_sqrt(xCubed);

#ifdef DEBUG_EC
  // Check the square root is correct
  Fp2_t allegedRhs;
  Fp2_mul(allegedRhs, out->Y, out->Y);
  assert(Fp2_cmpeq(allegedRhs, out->X));
#endif

  // Compute the other root
  Fp2_t negSqrtRhs;
  Fp2_neg(negSqrtRhs, out->Y);

  if ((useSmallerRoot && negSqrtRhs->a1->limbs[3] <= out->Y->a1->limbs[3]) ||
    (!useSmallerRoot && negSqrtRhs->a1->limbs[3] > out->Y->a1->limbs[3])) {
    Fp2_copy(out->Y, negSqrtRhs);    
  }  
}

void EncodingArithBN::decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt) {
  assert(compressed_elt);
  assert(elt);
  
  if (t == L) {
    LEncodedEltArithBN* e = (LEncodedEltArithBN*)elt;
    decompressAff(*((Fp_t*)compressed_elt), e->me, members->b, members->n, members->Lzero.me);
    Fp_weierstrass_aff2proj(e->me, e->me);    
  } else {
    assert(t == R);
    REncodedEltArithBN* e = (REncodedEltArithBN*)elt;
    decompressAffR(*((Fp2_t*)compressed_elt), e->me, members->n, members->Rzero.me);
    Fp2_weierstrass_aff2proj(e->me, e->me);
  }

  testOnCurve(t, elt);
}

void EncodingArithBN::write(Archiver* arc, elt_t t, EncodedElt* e, bool simple) {
  assert(arc && e);

  if (simple) {
    if (isZero(t, e)) {
      if (t == L) {
        write_zero(arc, members->n);  // Zero out both X and Y
        write_zero(arc, members->n);
      } else {
        assert(t == R);
        write_zero(arc, members->n);  // Zero out X 
        write_zero(arc, members->n);
        write_zero(arc, members->n);  // and Y
        write_zero(arc, members->n);
      }
      return;
    } else {
      affinize(t, e, e);   // C code expects affine coordinates
    }
  } else {
#ifdef CONDENSE_ARCHIVES
    compress(t, e);
#endif
  }

	if (t == L) {
    LEncodedEltArithBN* elt = (LEncodedEltArithBN*)e;

		write_Fp(arc, &elt->me->X, members->n);
#ifndef CONDENSE_ARCHIVES
		write_Fp(arc, &elt->me->Y, members->n);
    if (!simple) {
      write_Fp(arc, &elt->me->Z, members->n);
    }
#endif // CONDENSE_ARCHIVES
	} else {
		assert(t == R);
		REncodedEltArithBN* elt = (REncodedEltArithBN*)e;
    write_Fp2(arc, &elt->me->X, members->n);
#ifndef CONDENSE_ARCHIVES
		write_Fp2(arc, &elt->me->Y, members->n);
    if (!simple) {
      write_Fp2(arc, &elt->me->Z, members->n);
    }
#endif // CONDENSE_ARCHIVES
	}
  if (simple) {
    // Restore the point to projective form, in case we want to keep using it
    projectivize(t, e);
  }
}

void EncodingArithBN::read (Archiver* arc, elt_t t, EncodedElt* e) {
	assert(arc && e);

	if (t == L) {
    LEncodedEltArithBN* elt = (LEncodedEltArithBN*)e;
		//	Fp_init(elt->me->X);		TODO: Is this needed?
    read_Fp(arc, &elt->me->X, members->n);
#ifdef CONDENSE_ARCHIVES
    decompress(t, &elt->me->X, elt);
#else // !CONDENSE_ARCHIVES
		read_Fp(arc, &elt->me->Y, members->n);
		read_Fp(arc, &elt->me->Z, members->n);    
#endif // CONDENSE_ARCHIVES
	} else {
		assert(t == R);
		REncodedEltArithBN* elt = (REncodedEltArithBN*)e;
    read_Fp2(arc, &elt->me->X, members->n);
#ifdef CONDENSE_ARCHIVES
    decompress(t, &elt->me->X, elt);
#else // !CONDENSE_ARCHIVES
		read_Fp2(arc, &elt->me->Y, members->n);
		read_Fp2(arc, &elt->me->Z, members->n);   
#endif // CONDENSE_ARCHIVES
	}
}


void EncodingArithBN::affinize(elt_t t, EncodedElt* in, EncodedElt* out) {
	assert(in && out);
	testOnCurve(t, in);

	if (t == L) {
    LEncodedEltArithBN* ine  = (LEncodedEltArithBN*)in;
    LEncodedEltArithBN* oute = (LEncodedEltArithBN*)out;

    Fp_weierstrass_jacproj2aff(oute->me, ine->me);
	} else {
		assert(t == R);
    REncodedEltArithBN* ine  = (REncodedEltArithBN*)in;
    REncodedEltArithBN* oute = (REncodedEltArithBN*)out;

    Fp2_weierstrass_jacproj2aff(oute->me, ine->me);
	}
}

void EncodingArithBN::projectivize(elt_t t, EncodedElt* elt) {
  assert(elt);

  if (t == L) {
    LEncodedEltArithBN* elte = (LEncodedEltArithBN*)elt;
    Fp_weierstrass_aff2proj(elte->me, elte->me);
  } else {
    assert(t == R);
    REncodedEltArithBN* elte = (REncodedEltArithBN*)elt;
    Fp2_weierstrass_aff2proj(elte->me, elte->me);
  }
}


void EncodingArithBN::testOnCurve(elt_t t, EncodedElt* a) {
#ifdef DEBUG_EC
	assert(a);

	if (t == L) {
    LEncodedEltArithBN* ae = (LEncodedEltArithBN*)a;
    Fp_weierstrass_oncurve_jacproj(ae->me);
	} else {
		assert(t == R);
    REncodedEltArithBN* ae = (REncodedEltArithBN*)a;
    Fp2_weierstrass_oncurve_jacproj(ae->me);
	}
#endif
}

void EncodingArithBN::encodeSlow(elt_t t, FieldElt* in, EncodedElt* out) {
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


EncodedEltArithBNArray::EncodedEltArithBNArray(elt_t type, Encoding* encoding, uint32_t len, bool preallocate) {
	this->type = type;
  this->encoding = encoding;
	this->len = len;
	elts = new EncodedEltArithBN*[len];
	allocated = preallocate;
	
	for (uint32_t i = 0; i < len; i++) {
		if (preallocate) {
			elts[i] = (EncodedEltArithBN*)encoding->new_elt(type);
		} else {
			elts[i] = NULL;
		}
	}
}

EncodedEltArithBNArray::~EncodedEltArithBNArray() {
	if (allocated) {
		for (uint32_t i = 0; i < len; i++) {
			encoding->del_elt(this->type, elts[i]);
		}
	}
	delete [] elts;
}

EncodedElt* EncodedEltArithBNArray::elt(uint32_t index) {
	assert(index < len);
	assert(elts[index] != NULL);
	return elts[index];
}

void EncodedEltArithBNArray::set(uint32_t index, EncodedElt* elt) {
	assert(index < len);
	assert(elts[index] == NULL);	// Otherwise we're overwriting an existing value!
	elts[index] = (EncodedEltArithBN*)elt;
}

void EncodedEltArithBNArray::reset() {
  assert(!allocated);
  for (uint32_t i = 0; i < len; i++) {
    elts[i] = NULL;
  }
}
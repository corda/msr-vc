#include <assert.h>
#include "EncodingArithBN.h"
#include "EncodingArithCP.h"
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

uint32_t EncodingArithCP::num_L_digits;
uint32_t EncodingArithCP::num_R_digits;

EncodingArithCP::EncodingArithCP(CP_curve_type curve_type) {
  assert(PAIR_CURVE == 0);    // Make sure no one else has set it.
  members = new EncodingArithCPMembers;
  this->curve_type = curve_type;

  if (curve_type == CP_curve_type::CP_6) {

    PAIR_CURVE = CP6b;

    uint64_t pCP6b[8] = { 0xF9000000000005BF, 0xAA8EC00000000BD3,
      0x0CF85A8000000B1F, 0x8761FD7D40000638,
      0x0A67549B77800241, 0xB007E82F1F3CB08B,
      0x41358FBC1F137315, 0x158CFFC9264E1116 };
    /* This is the 509-bit CP6b prime p = 3 (mod 4). */

    uint64_t Px[8] = { 0x549400D7139BD8C3, 0xF63D76CC9E5074C4,
      0x9ECA293D846AD7F4, 0xBCFDE6DF8A0528AC,
      0x7412962AED6CB7E3, 0x8CA78CFBE18D6C03,
      0x9D41ADD3F736DF40, 0x103E84948482F3AC };
    uint64_t Py[8] = { 0xA49BC75E21A206DA, 0x7EC88B25BB3E39DF,
      0xC54F80610A16417D, 0xF7CB2BB06F683525,
      0xB76C9F7CC52E9739, 0xBAB17A19BBB5F211,
      0x176CE782E1557D4E, 0x0521E32E8D0B8F03 };
    /* This is the point over Fp in G1. */

    uint64_t Qx[8] = { 0x4179C530BAF016F7, 0x7D29FB703C56B725,
      0x8D8AA7AEAEBB06D3, 0xA6CA4B97C79D31C1,
      0xA49308E0F3BB3294, 0xB1ECCE42F22C7800,
      0x9B245837BEC9D3B5, 0x113E15D42AB7CF02 };
    uint64_t Qy[8] = { 0x0E83D7FB80304008, 0x1A0B1C4E85A3FC34,
      0x731C345CFA535AB3, 0x507571381258BBEF,
      0x3B427423987C4B28, 0x2535C8DC6F89FBF2,
      0x0E34A6FB987994A0, 0x00FCE2C1051EE0D1 };
    members->n = 8;
    //printf ("Running CP6b pairing...\n");

    Fp_initialize_config(pCP6b, members->n);
    Fp3_initialize_config();
    Fp6q_initialize_config();

    Fp_init(members->B[0]);
    Fp_init(members->B[1]);
    Fp_init(members->B3);

    Fp_set_ui(members->B[0], 21);
    Fp_set_ui(members->B[1], 63);

    Fp_initialize_weierstrass(members->B, 2);

    cp6_miller_initialize_config(members->B[1]);
    cp6_finalexpo_initialize_config();

    // For now, use Enigma's encoding to create the prime field we want
    members->enigma_encoding = new EncodingEnigmaBN();
    this->srcField = ((EncodingEnigmaBN*)members->enigma_encoding)->getBaseField();    

    // Convert the generators to projective form
    Fp_wp_init(members->Lg.me);
    Fp_set(members->Lg.me->X, Px);
    Fp_set(members->Lg.me->Y, Py);
    Fp_weierstrass_aff2proj(members->Lg.me, members->Lg.me);

    Fp_wp_init(members->Rg.me);
    Fp_set(members->Rg.me->X, Qx);
    Fp_set(members->Rg.me->Y, Qy);
    Fp_weierstrass_aff2proj(members->Rg.me, members->Rg.me);

    // Create representations of zero, since we use it a lot
    Fp_wp_init(members->Lzero.me);
    Fp_copy(members->Lzero.me->X, members->Lg.me->X);
    Fp_copy(members->Lzero.me->Y, members->Lg.me->Y);
    memset(members->Lzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

    Fp_wp_init(members->Rzero.me);
    Fp_copy(members->Rzero.me->X, members->Rg.me->X);
    Fp_copy(members->Rzero.me->Y, members->Rg.me->Y);
    memset(members->Rzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

    EncodingArithCP::num_L_digits = 3 * 8;  // Three field elements (8 digits each) to define a projective point on the base curve
    EncodingArithCP::num_R_digits = 3 * 8;  // Three field elements (8 digits each) to define a projective point on the twist curve
  }  else if (curve_type == CP_curve_type::CP_3) {

    PAIR_CURVE = CP3;

    uint64_t pCP3[16] = { 0x1500000005D641EF, 0xF27A90001804027E,
      0xDF0B05082F4B95A0, 0x0ACAC924F81EC1DA,
      0xA8EF36C595CF5E4F, 0xEB5F7B2090C91C8C,
      0xFDAE1FBE2E9344F8, 0x6E7B70C7C1E1ED72,
      0xA39DF221BDDF6026, 0x2BF6589BEC58AF41,
      0xA90B68306AB9785C, 0x7301607FD616A3FD,
      0xCEB7B890F7C56B98, 0x2DEF7BDF8234EAE9,
      0xB7368AE71D33BC04, 0x540531971EED5543 };

    /* This is the 1023-bit CP3 prime p = 3 (mod 4). */

    uint64_t Px[16] = { 0x758CE21AD1C2BE3E, 0x0C820B5BC4666991,
      0x32E11228028216A5, 0x37940D6905381AF1,
      0x6D44464C1E9B33A9, 0xCD0A244DF694893F,
      0x6071A8DFDEE42983, 0x54B5E323C7DED5F4,
      0x5FA7BC828317C453, 0x0C56050AF0713A5E,
      0x9A9F35C0FF516FAA, 0x670D3319658D4763,
      0x0AD3A71077DACE33, 0x86CEF75311F2B675,
      0xBCC91FAABDAEC15E, 0x1CE3057DA88BBCBA };
    uint64_t Py[16] = { 0x4BBBA36B15E3351E, 0xE2323BE5C035DFAF,
      0x39E74D8B8E656050, 0x998F123AE159D451,
      0x3F9285C69B85C8AE, 0x6C8BBED349A2BBE3,
      0x1B44946ADA9E7A7A, 0x3660AC7781F16A5E,
      0x8D103674B090D9C6, 0xBBA04490A82FFEC0,
      0x3DA002C85F8C83C7, 0x77B05A4EB3A02B7F,
      0x9A9CF058FAE676C4, 0xABD0D34D251CE709,
      0x1794BED9C920ED53, 0x29852D209F1373A5 };

    uint64_t Qx[16] = { 0xC5929CFB23259A13, 0xBECC2D7FD8E166C3,
      0x6AD6BEA4E8EC2385, 0x6266AE6B695D7533,
      0xEEABA40CC791ED22, 0x975FE9EE29D69B73,
      0x6F6AAB22E490C9B0, 0x34DF16327FF1B1B5,
      0xA37E0C4D542D99D1, 0xBC39D34D3C64089B,
      0x20AB141639629A0C, 0x539C4C7EF16062CD,
      0xE3692158F74F7A9C, 0x047BCE81FD8CACDD,
      0x8629D4E07C56595E, 0x0BADCF277AB45C44 };
    uint64_t Qy[16] = { 0x4E41A9F542CD76B7, 0xD3AE6839C73230FB,
      0xC89BA875B3579302, 0x460F55F967EB2EC0,
      0x8BC3B3AC34E61A53, 0xB619994E08BABD9D,
      0x9708A88BB55A3F37, 0x7755F188258AFC24,
      0xAD9E0E9C3E01AE76, 0xE80CA9292A5770EB,
      0x5324028AA7AB6170, 0xF39FAC955522EDCF,
      0x1A8EB2D3A0A9F7ED, 0xD5F0368FB00D9F08,
      0xB8C9CDB6E13D0407, 0x21E83F4C17A96BEC };

    members->n = 16;
    //printf("Running CP3 pairing...\n");

    Fp_initialize_config(pCP3, members->n);
    Fp3_initialize_config();

    Fp_init(members->B[0]);
    Fp_init(members->B[1]);
    Fp_init(members->B3);

    Fp_set_ui(members->B[0], 1);
    Fp_set_ui(members->B[1], 11);

    Fp_initialize_weierstrass(members->B, 2);

    cp3_miller_initialize_config(members->B[1]);
    cp3_finalexpo_initialize_config();

    //// For now, use Enigma's encoding to create the prime field we want    
    members->enigma_encoding = NULL;
    this->srcField = FieldPrime::get_CP6_base_field();

    // Convert the generators to projective form
    Fp_wp_init(members->Lg.me);
    Fp_set(members->Lg.me->X, Px);
    Fp_set(members->Lg.me->Y, Py);
    Fp_weierstrass_aff2proj(members->Lg.me, members->Lg.me);

    Fp_wp_init(members->Rg.me);
    Fp_set(members->Rg.me->X, Qx);
    Fp_set(members->Rg.me->Y, Qy);
    Fp_weierstrass_aff2proj(members->Rg.me, members->Rg.me);

    // Create representations of zero, since we use it a lot
    Fp_wp_init(members->Lzero.me);
    Fp_copy(members->Lzero.me->X, members->Lg.me->X);
    Fp_copy(members->Lzero.me->Y, members->Lg.me->Y);
    memset(members->Lzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

    Fp_wp_init(members->Rzero.me);
    Fp_copy(members->Rzero.me->X, members->Rg.me->X);
    Fp_copy(members->Rzero.me->Y, members->Rg.me->Y);
    memset(members->Rzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

    EncodingArithCP::num_L_digits = 3 * 16;
    EncodingArithCP::num_R_digits = 3 * 16;
  } else if (curve_type == CP_curve_type::DEBUG) {
    // Use the much smaller debug curve
    PAIR_CURVE = CP6tiny;

    uint64_t pCP6[2] = { 0x2108F75F57168FE7, 0x76AD76A42A2 };
    /* This is the 107-bit CP6 prime p = 3 (mod 4). */

    uint64_t Px[2] = { 0xF25743F1668655EF, 0x555807C161F };
    uint64_t Py[2] = { 0x8D2D33A7D83B9E10, 0x44BCCF08BBC };
    /* This is the point over Fp in G1. */

    uint64_t Qx[2] = { 0x182624C621C68307, 0x3665C10BD94 };
    uint64_t Qy[2] = { 0x2F8E180B89361A5D, 0x32BC9CF3EDE };

    members->n = 2;

    Fp_initialize_config(pCP6, members->n);
    Fp3_initialize_config();
    Fp6q_initialize_config();

    Fp_init(members->B[0]);
    Fp_init(members->B[1]);
    Fp_init(members->B3);

    Fp_set_ui(members->B[0], 13);
    Fp_set_ui(members->B[1], 39);

    Fp_initialize_weierstrass(members->B, 2);

    cp6_miller_initialize_config(members->B[1]);
    cp6_finalexpo_initialize_config();

    // For now, use Enigma's encoding to create the prime field we want
    members->enigma_encoding = new EncodingEnigmaBN();
    this->srcField = ((EncodingEnigmaBN*)members->enigma_encoding)->getBaseField(true);

    // Convert the generators to projective form
    Fp_wp_init(members->Lg.me);
    Fp_set(members->Lg.me->X, Px);
    Fp_set(members->Lg.me->Y, Py);
    Fp_weierstrass_aff2proj(members->Lg.me, members->Lg.me);

    Fp_wp_init(members->Rg.me);
    Fp_set(members->Rg.me->X, Qx);
    Fp_set(members->Rg.me->Y, Qy);
    Fp_weierstrass_aff2proj(members->Rg.me, members->Rg.me);

    // Create representations of zero, since we use it a lot
    Fp_wp_init(members->Lzero.me);
    Fp_copy(members->Lzero.me->X, members->Lg.me->X);
    Fp_copy(members->Lzero.me->Y, members->Lg.me->Y);
    memset(members->Lzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

    Fp_wp_init(members->Rzero.me);
    Fp_copy(members->Rzero.me->X, members->Rg.me->X);
    Fp_copy(members->Rzero.me->Y, members->Rg.me->Y);
    memset(members->Rzero.me->Z->limbs, 0, sizeof(uint64_t)*members->n);   // Should be point at infinity in projective representation

    // TODO: Why does this work?
    EncodingArithCP::num_L_digits = 3 * 8;  // Three field elements (8 digits each) to define a projective point on the base curve
    EncodingArithCP::num_R_digits = 3 * 8;  // Three field elements (8 digits each) to define a projective point on the twist curve
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
  }

}

EncodingArithCP::~EncodingArithCP() {
  Fp_wp_free(members->Lzero.me);
  Fp_wp_free(members->Rzero.me);
  Fp_wp_free(members->Lg.me);
  Fp_wp_free(members->Rg.me);

  Fp_free_weierstrass ();

  Fp_free (members->B[0]);
  Fp_free (members->B[1]);
  Fp_free (members->B3);

  cp6_finalexpo_free_config ();
  cp6_miller_free_config();
  cp3_finalexpo_free_config();
  cp3_miller_free_config();
  Fp6q_free_config ();
  Fp3_free_config ();
  Fp_free_config ();

  delete this->srcField;
  delete members->enigma_encoding;
  delete members;
}

// TODO: Adjust to count pointers as well
uint32_t EncodingArithCP::bytes_per(elt_t t) {
	if (t == L) {
		return num_L_digits*sizeof(digit_t);
	} else {
		assert(t == R);
		return num_R_digits*sizeof(digit_t);
	}
}

uint32_t EncodingArithCP::num_exponent_bits() {
  if (curve_type == CP_curve_type::CP_6) {
    return 4 * sizeof(digit_t) * 8 - 1;	// BN_curve uses 4 digits, but doesn't use them all (only uses 254 bits)
  } else if (curve_type == CP_curve_type::CP_3) {
    return 509 + 1;
  } else if (curve_type == CP_curve_type::DEBUG) {
    return 53 + 1;
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
    return 0;
  }
}


EncodedElt* EncodingArithCP::new_elt(elt_t t) {
  EncodedEltArithCP* elt = new EncodedEltArithCP;
  Fp_wp_init(elt->me);
	return elt;
}

EncodedEltArray* EncodingArithCP::new_elt_array(elt_t t, uint32_t len, bool preallocate) {
	return new EncodedEltArithCPArray(t, this, len, preallocate);
}


void EncodingArithCP::del_elt(elt_t t, EncodedElt* elt) {
  assert(elt);
  Fp_wp_free(((EncodedEltArithCP*)elt)->me);
	delete elt;
}

EncodedProduct* EncodingArithCP::new_prod() {
	EncodedProductArithCP* ret = new EncodedProductArithCP();
  ret->is3 = curve_type == CP_curve_type::CP_3;
  if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
    Fp6q_init(ret->me.me6);
  } else if (curve_type == CP_curve_type::CP_3) {
    Fp3_init(ret->me.me3);
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
  }
  return ret;
}

EncodedProductArithCP::~EncodedProductArithCP() {
  if (is3) {
    Fp3_free(me.me3);
  } else {
    Fp6q_free(me.me6);
  } 
}

uint32_t EncodingArithCP::size_of_prod() {
  return sizeof(Fp6q_t);
}

void EncodingArithCP::print(elt_t t, EncodedElt* e) {
	EncodedEltArithCP* en = (EncodedEltArithCP*)e;
	Fp_weierstrass_jacproj2aff(en->me, en->me);
	printf("( ");
	Fp_print(en->me->X);
	printf(", ");
	Fp_print(en->me->Y);
	printf(", ");
	Fp_print(en->me->Z);
	printf(")");	
}

void EncodingArithCP::print(EncodedProduct* p) {
	EncodedProductArithCP* ep = (EncodedProductArithCP*)p;	

  if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
    Fp6q_print(ep->me.me6);
  } else if (curve_type == CP_curve_type::CP_3) {
    Fp3_print(ep->me.me3);
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
  }
}

void EncodingArithCP::copy(elt_t t, EncodedElt* src, EncodedElt* dst) {
  Fp_wp_copy(((EncodedEltArithCP*)dst)->me, ((EncodedEltArithCP*)src)->me);
}

bool EncodingArithCP::equal(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	testOnCurve(t, a);
	testOnCurve(t, b);

  bool ret = 0 != Fp_weierstrass_equal_jacproj(((EncodedEltArithCP*)a)->me, ((EncodedEltArithCP*)b)->me);
  return ret;
}

void EncodingArithCP::zero(elt_t t, EncodedElt* a) { 
	if (t == L) {
	  copy(t, &members->Lzero, a);	// Use our precomputed value	
	} else {
		assert(t == R);
	  copy(t, &members->Rzero, a);	// Use our precomputed value	
	}
}

void EncodingArithCP::one(elt_t t, EncodedElt* a) {
	if (t == L) {
	  copy(t, &members->Lg, a);	// Use our precomputed value	
	} else {
		assert(t == R);
	  copy(t, &members->Rg, a);	// Use our precomputed value	
	}
}

// Adds two encoded elements.  Result goes in r.
void EncodingArithCP::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
  add(t, a, b, r, false);
}

// Adds two encoded elements.  Result goes in r.
void EncodingArithCP::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, bool definitely_distinct) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);
	assert (t == L || t == R);

  EncodedEltArithCP* ae = (EncodedEltArithCP*)a;
  EncodedEltArithCP* be = (EncodedEltArithCP*)b;
  EncodedEltArithCP* re = (EncodedEltArithCP*)r;
    
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
	
	testOnCurve(t, r);
}

// Subtracts two encoded elements.  Result goes in r.
void EncodingArithCP::sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);

  EncodedEltArithCP* ae = (EncodedEltArithCP*)a;
  EncodedEltArithCP* be = (EncodedEltArithCP*)b;
  EncodedEltArithCP* re = (EncodedEltArithCP*)r;

  EncodedEltArithCP neg_b;
  Fp_wp_init(neg_b.me);
  Fp_wp_neg(neg_b.me, be->me);
  add(t, ae, &neg_b, re);
  Fp_wp_free(neg_b.me);

	testOnCurve(t, r);
}

void EncodingArithCP::doubleIt(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	testOnCurve(t, a);

  EncodedEltArithCP* ae = (EncodedEltArithCP*)a;
  EncodedEltArithCP* be = (EncodedEltArithCP*)b;

  Fp_weierstrass_dbl(be->me, ae->me);

	testOnCurve(t, b);
}


// Multiplies the value encoded in g by the constant c.  Result goes in r.
// For pairings, we compute r <- g^c
void EncodingArithCP::mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r) {		
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

void EncodingArithCP::mul(LEncodedElt* a, REncodedElt* b, EncodedProduct* r) {
	pair(a, b, r);
}

void EncodingArithCP::pair(LEncodedElt* L1, REncodedElt* R1, EncodedProduct* result) {
	assert(L1 && R1 && result);
	testOnCurve(L, L1);
	testOnCurve(R, R1);

	EncodedEltArithCP* L1e = (EncodedEltArithCP*)L1;
	EncodedEltArithCP* R1e = (EncodedEltArithCP*)R1;
	EncodedProductArithCP* resulte = (EncodedProductArithCP*)result;

  if (isZero(L, L1) || isZero(R, R1)) {
    // Library doesn't handle g^0 values properly, so we special-case it

    if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
      Fp6q_set_one(resulte->me.me6);
    } else if (curve_type == CP_curve_type::CP_3) {
      Fp3_set_one(resulte->me.me3);
    } else {
      printf("Unexpected CP curve type\n");
      assert(false);
    }
    
  } else {
    Point_wp_t aff_L1; 
    Point_wp_t aff_R1;
    Fp_wp_init(aff_L1);
    Fp_wp_init(aff_R1);

    // Convert to affine representation before doing the pairing
    Fp_weierstrass_jacproj2aff(aff_L1, L1e->me);
    Fp_weierstrass_jacproj2aff(aff_R1, R1e->me);

    if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
      cp6b_optimal_ate_miller(resulte->me.me6, aff_R1, aff_L1);
      cp6b_finalexpo(resulte->me.me6, resulte->me.me6);
    } else if (curve_type == CP_curve_type::CP_3) {
      cp3_ate_miller(resulte->me.me3, aff_R1, aff_L1);
      cp3_finalexpo(resulte->me.me3, resulte->me.me3);
    } else {
      printf("Unexpected CP curve type\n");
      assert(false);
    }
   
    Fp_wp_free(aff_L1);
    Fp_wp_free(aff_R1);
  }

}

bool EncodingArithCP::isZero(elt_t t, EncodedElt* elt) {
	assert(elt);
	//testOnCurve(t, elt);

  bool ret = false;
  EncodedEltArithCP* elte = (EncodedEltArithCP*)elt;
	if (t == L) {	  
    ret = 0 == Fp_cmp(elte->me->Z, members->Lzero.me->Z);   // Check if the Z coord is 0 ==> pt at infinity
	} else {
		assert(t == R);
    ret = 0 == Fp_cmp(elte->me->Z, members->Rzero.me->Z);   // Check if the Z coord is 0 ==> pt at infinity
	}

  return ret;
}

void EncodingArithCP::add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
  EncodedProductArithCP* ae = (EncodedProductArithCP*)a;
  EncodedProductArithCP* be = (EncodedProductArithCP*)b;
  EncodedProductArithCP* re = (EncodedProductArithCP*)r;
  
  if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
    Fp6q_mul(re->me.me6, ae->me.me6, be->me.me6);
  } else if (curve_type == CP_curve_type::CP_3) {
    Fp3_mul(re->me.me3, ae->me.me3, be->me.me3);
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
  }
  
}

void EncodingArithCP::sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductArithCP* ae = (EncodedProductArithCP*)a;
	EncodedProductArithCP* be = (EncodedProductArithCP*)b;
	EncodedProductArithCP* re = (EncodedProductArithCP*)r;
  
  EncodedProductArithCP* inv = (EncodedProductArithCP*)this->new_prod();
  
  if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
    Fp6q_inv(inv->me.me6, be->me.me6);
    Fp6q_mul(re->me.me6, ae->me.me6, inv->me.me6);
  } else if (curve_type == CP_curve_type::CP_3) {
    Fp3_inv(inv->me.me3, be->me.me3);
    Fp3_mul(re->me.me3, ae->me.me3, inv->me.me3);
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
  }
  delete inv;
}

bool EncodingArithCP::equals(EncodedProduct* a, EncodedProduct* b) {
	assert(a && b);
	EncodedProductArithCP* ae = (EncodedProductArithCP*)a;
	EncodedProductArithCP* be = (EncodedProductArithCP*)b;
  if (curve_type == CP_curve_type::CP_6 || curve_type == CP_curve_type::DEBUG) {
    return 1 == Fp6q_cmpeq(ae->me.me6, be->me.me6);
  } else if (curve_type == CP_curve_type::CP_3) {
    return 1 == Fp3_cmpeq(ae->me.me3, be->me.me3);
  } else {
    printf("Unexpected CP curve type\n");
    assert(false);
    return false;
  }  
}

// Checks whether the plaintexts inside the encodings obey:
// L1*R1 - L2*R2 == L3*R3
bool EncodingArithCP::mulSubEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2, LEncodedElt* L3, REncodedElt* R3) {
	assert(L1 && R1 && L2 && R2 && L3 && R3);

	EncodedProductArithCP* pairing1;
	EncodedProductArithCP* pairing2;
	EncodedProductArithCP* pairing3;

  pairing1 = (EncodedProductArithCP*)this->new_prod();
  pairing2 = (EncodedProductArithCP*)this->new_prod();
  pairing3 = (EncodedProductArithCP*)this->new_prod();

	testOnCurve(L, L1); testOnCurve(R, R1);
	testOnCurve(L, L2); testOnCurve(R, R2);
	testOnCurve(L, L3); testOnCurve(R, R3);

	// Compute the three pairings
	pair(L1, R1, pairing1);
	pair(L2, R2, pairing2);
	pair(L3, R3, pairing3);

	// Compute (with some notation abuse) L1*R1 - L2*R2 (in the exponents), i.e., pairing1 <- pairing1 / pairing2
	sub(pairing1, pairing2, pairing1);
	bool result = equals(pairing1, pairing3);

  delete pairing1;
  delete pairing2;
  delete pairing3;  

	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}


// Checks whether the plaintexts inside the encodings obey:
// L1*R1 == L2*R2
bool EncodingArithCP::mulEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2) {
	assert(L1 && R1 && L2 && R2);
	EncodedProductArithCP* pairing1;
	EncodedProductArithCP* pairing2;
  pairing1 = (EncodedProductArithCP*)this->new_prod();
  pairing2 = (EncodedProductArithCP*)this->new_prod();

	testOnCurve(L, L1); testOnCurve(R, R1);
	testOnCurve(L, L2); testOnCurve(R, R2);

	// Compute the two pairings
	pair(L1, R1, pairing1);
	pair(L2, R2, pairing2);

	BOOL result = equals(pairing1, pairing2); 

  delete pairing1;
  delete pairing2;
	
	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}

// Expects out to point to a buffer with at least condensed_size_in_digits*sizeof(digit_t) space
void compressAff(bool elt_is_zero, Point_wp_t aff_in, Fp_t out, Fp_t b) {
  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (elt_is_zero) {
    Fp_set_ui(out, 0);
    return;
  }

  // Compute x^3 + b, where b is supplied
  Fp_t xCubed;
  Fp_init(xCubed);
  Fp_mul(xCubed, aff_in->X, aff_in->X);
  Fp_mul(xCubed, xCubed, aff_in->X);
  Fp_add(xCubed, xCubed, b);

  // Compute positive root
  Fp_t sqrtRhs;
  Fp_init(sqrtRhs);
  Fp_sqrt(sqrtRhs, xCubed);
  
#ifdef DEBUG_EC
  // Check the square root is correct
  Fp_t allegedRhs;
  Fp_init(allegedRhs);
  Fp_mul(allegedRhs, sqrtRhs, sqrtRhs);
  assert(0 == Fp_cmp(allegedRhs, xCubed));
#endif

  bool equal = 0 == Fp_cmp(sqrtRhs, aff_in->Y);
  out = aff_in->X;
  // We use the top-most bit of X to store the bit indicating whether to use the negative root
  int lastDigit = 8 / 2 - 1;	// EncodedElt contains 2 field elements
  if (!equal) {
    // Set the top bit of the last digit of X to indicate that we use the "negative" root		    
    out->limbs[lastDigit] |= ((digit_t)1) << 63;
  }
  else {
    // Make sure the top bit is actually clear, indicating the "positive" root		
    assert(0 == (out->limbs[lastDigit] & (((digit_t)1) << 63)));
  }
}



void EncodingArithCP::compress(elt_t t, EncodedElt* elt) {
  assert(elt);
  EncodedEltArithCP* e = (EncodedEltArithCP*)elt;
  
  affinize(t, e, e);
  if (t == L) {
    bool isZero = equal(t, elt, &members->Lzero);
    compressAff(isZero, e->me, e->me->X, members->B[0]);
  } else {
    assert(t == R);
    bool isZero = equal(t, elt, &members->Rzero);
    compressAff(isZero, e->me, e->me->X, members->B[1]);
  }
}

bool allZero(digit_t* value, int num_digits) {
  for (int i = 0; i < num_digits; i++) {
    if (value[i] != 0) {
      return false;
    }
  }
  return true;
}

void decompressAff(Fp_t in, Point_wp_t out, Fp_t b, int digits_in_fp_t, Point_wp_t zero) {
  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (allZero(in->limbs, digits_in_fp_t)) {
    Fp_wp_copy(out, zero);
    return;
  }

  // Check whether the top bit of X is set.  Indicates which root to use for Y
  int lastDigit = digits_in_fp_t - 1;
  bool useNegRoot = (in->limbs[lastDigit] & (((digit_t)1) << 63)) > 0;

  // Either way, clear out the top bit and copy the input into the x coefficient
  in->limbs[lastDigit] &= ((((digit_t)1) << 63) - 1);
  Fp_copy(out->X, in);

  // Compute x^3 + b, where b is supplied
  Fp_t xCubed;
  Fp_init(xCubed);
  Fp_mul(xCubed, out->X, out->X);
  Fp_mul(xCubed, xCubed, out->X);
  Fp_add(xCubed, xCubed, b);

  // Compute "positive" root and, if necessary, use it to compute the negative root
  Fp_sqrt(out->Y, xCubed);  

#ifdef DEBUG_EC
  // Check the square root is correct
  Fp_t allegedRhs;
  Fp_init(allegedRhs);
  Fp_mul(allegedRhs, out->Y, out->Y);
  assert(0 == Fp_cmp(allegedRhs, xCubed));
#endif

  if (useNegRoot) {
    Fp_neg(out->Y, out->Y);    
  }  
}

void EncodingArithCP::decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt) {
  assert(compressed_elt);
  assert(elt);
  EncodedEltArithCP* e = (EncodedEltArithCP*)elt;
  Fp_t* ce = (Fp_t*)compressed_elt;
  
  if (t == L) {
    decompressAff(*ce, e->me, members->B[0], members->n, members->Lzero.me);
  } else {
    assert(t == R);
    decompressAff(*ce, e->me, members->B[1], members->n, members->Rzero.me);
  }

  Fp_weierstrass_aff2proj(e->me, e->me);  
  testOnCurve(t, elt);
}

void EncodingArithCP::write(Archiver* arc, elt_t t, EncodedElt* e, bool simple) {
	assert(arc && e);
	EncodedEltArithCP* ee = (EncodedEltArithCP*)e;

#ifdef CONDENSE_ARCHIVES
  compress(t, e);
#endif // CONDENSE_ARCHIVES

	write_Fp(arc, &ee->me->X, members->n);
#ifndef CONDENSE_ARCHIVES
	write_Fp(arc, &ee->me->Y, members->n);
	write_Fp(arc, &ee->me->Z, members->n);  
#endif // CONDENSE_ARCHIVES
}

void EncodingArithCP::read(Archiver* arc, elt_t t, EncodedElt* e) {
	assert(arc && e);
	EncodedEltArithCP* ee = (EncodedEltArithCP*)e;

	read_Fp(arc, &ee->me->X, members->n);
#ifdef CONDENSE_ARCHIVES
  decompress(t, &ee->me->X, ee);
#else // !CONDENSE_ARCHIVES
	read_Fp(arc, &ee->me->Y, members->n);
	read_Fp(arc, &ee->me->Z, members->n);  
#endif // CONDENSE_ARCHIVES


}

void EncodingArithCP::affinize(elt_t t, EncodedElt* in, EncodedElt* out) {
	assert(in && out);
	testOnCurve(t, in);
	
  EncodedEltArithCP* ine  = (EncodedEltArithCP*)in;
  EncodedEltArithCP* oute = (EncodedEltArithCP*)out;

  Fp_weierstrass_jacproj2aff(oute->me, ine->me);
}



void EncodingArithCP::testOnCurve(elt_t t, EncodedElt* a) {
#ifdef DEBUG_EC
	assert(a);

  EncodedEltArithCP* ae = (EncodedEltArithCP*)a;
  Fp_weierstrass_oncurve_jacproj(ae->me);	
#endif
}

void EncodingArithCP::encodeSlow(elt_t t, FieldElt* in, EncodedElt* out) {
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


EncodedEltArithCPArray::EncodedEltArithCPArray(elt_t type, Encoding* encoding, uint32_t len, bool preallocate) {
	this->type = type;
  this->encoding = encoding;
	this->len = len;
	elts = new EncodedEltArithCP*[len];
	allocated = preallocate;
	
	for (uint32_t i = 0; i < len; i++) {
		if (preallocate) {
			elts[i] = (EncodedEltArithCP*)encoding->new_elt(type);
		} else {
			elts[i] = NULL;
		}
	}
}

EncodedEltArithCPArray::~EncodedEltArithCPArray() {
	if (allocated) {
		for (uint32_t i = 0; i < len; i++) {
      encoding->del_elt(this->type, elts[i]);
		}
	}
	delete [] elts;
}

EncodedElt* EncodedEltArithCPArray::elt(uint32_t index) {
	assert(index < len);
	assert(elts[index] != NULL);
	return elts[index];
}

void EncodedEltArithCPArray::set(uint32_t index, EncodedElt* elt) {
	assert(index < len);
	assert(elts[index] == NULL);	// Otherwise we're overwriting and existing value!
	elts[index] = (EncodedEltArithCP*)elt;
}


void EncodedEltArithCPArray::reset() {
  assert(!allocated);
  for (uint32_t i = 0; i < len; i++) {
    elts[i] = NULL;
  }
}

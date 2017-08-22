#include <assert.h>
#include "EncodingEnigmaBN.h"
#include <bitset>
#include <intrin.h>
#include "Poly.h"
#include "FieldPrime.h"
#include "EnigmaWrappers.h"


#pragma intrinsic(_BitScanReverse)

int EncodingEnigmaBN::num_live_enc_elts[2];

EncodingEnigmaBN::EncodingEnigmaBN() {
	members = new EncodingEnigmaBNMembers;

	members->pbigctx = &this->members->BignumCtx;
	memset(this->members->pbigctx, 0, sizeof(bigctx_t));

	// Setup the pairing parameters
	BOOL OK = TRUE;
	const char *bnname = named_bn[0].name;
	const char *bninfo = *(named_bn[0].bninfo);               // Pairing info 
    
	memset((void*)&members->bn, 0, sizeof(members->bn));                        
	OK = BN_build(bninfo, &members->bn, members->PBIGCTX_PASS);				  // Set up pairing struct
	assert(OK);	

	// Create a field to represent operations in the exponent
	members->p = new mp_modulus_t; 

	OK = create_modulus(members->bn.BaseCurve.gorder, members->bn.BaseCurve.lnggorder, FROM_LEFT, members->p, members->PBIGCTX_PASS);
	assert(OK);

	members->msExponentField = new field_desc_t; 

	OK = Kinitialize_prime(members->p, members->msExponentField, members->PBIGCTX_PASS);
	assert(OK);

  FieldDescription* desc = new FieldDescription;
  desc->desc = members->msExponentField;
  srcField = new FieldPrime(desc, 254);
  delete desc;

  makeBaseField();
  this->montR = base_field->newElt();
  base_field->set(this->montR, 2);
  base_field->exp(this->montR, 64 * 4, this->montR);
  this->montRinv = base_field->newElt();
  FieldElt* one_f = base_field->newElt();
  base_field->one(one_f);
  base_field->div(one_f, this->montR, this->montRinv);
  base_field->delElt(one_f);

	// Save curve pointers, so we can access them based on curve type
  members->curve[L] = &(members->bn.BaseCurve);
  members->curve[R] = &(members->bn.TwistCurve);
	
	// Build projective versions of the two members->curves' generators
	gen[L] = (EncodedEltEnigmaBN*)new_elt(L);
	gen[R] = (EncodedEltEnigmaBN*)new_elt(R);

	OK = ecaffine_projectivize(members->bn.BaseCurve.generator, gen[L]->get_digits(), &members->bn.BaseCurve, members->arithTempSpace, members->PBIGCTX_PASS);
	assert(OK);

	OK = ecaffine_projectivize(members->bn.TwistCurve.generator, gen[R]->get_digits(), &members->bn.TwistCurve, members->arithTempSpace, members->PBIGCTX_PASS);
	assert(OK);

	// Create a representation of zero, since we use it a lot
	FieldEltPrime* zero = (FieldEltPrime*)srcField->newElt();
	srcField->zero(zero);

	Zero[L] = (EncodedEltEnigmaBN*)new_elt(L);
	Zero[R] = (EncodedEltEnigmaBN*)new_elt(R);
	OK = ecproj_single_exponentiation(gen[L]->get_digits(), zero->get_digits(), members->bn.pbits, members->bn.BaseCurve.lnggorder, Zero[L]->get_digits(), &members->bn.BaseCurve, members->PBIGCTX_PASS);
	assert(OK);
	OK = ecproj_single_exponentiation(gen[R]->get_digits(), zero->get_digits(), members->bn.pbits, members->bn.TwistCurve.lnggorder, Zero[R]->get_digits(), &members->bn.TwistCurve, members->PBIGCTX_PASS);
	assert(OK);

	srcField->delElt(zero);

  this->num_addsub[0] = 0;
  this->num_double[0] = 0;
  this->num_constmul[0] = 0;
  this->num_addsub[1] = 0;
  this->num_double[1] = 0;
  this->num_constmul[1] = 0;
  this->num_pair = 0;
  this->collecting_stats = false;

  Encoding::current_encoding = this;
}

Field* EncodingEnigmaBN::getMontField() {
  // Create a field to represent operations in the exponent
  mp_modulus_t* prime = new mp_modulus_t;

  BOOL OK = create_modulus(members->bn.BaseCurve.gorder, members->bn.BaseCurve.lnggorder, FROM_RIGHT, prime, members->PBIGCTX_PASS);
  assert(OK);

  field_desc_t* msfield = new field_desc_t;

  OK = Kinitialize_prime(prime, msfield, members->PBIGCTX_PASS);
  assert(OK);

  FieldDescription* desc = new FieldDescription;
  desc->desc = msfield;
  Field* f = new FieldPrime(desc, 254);
  delete desc;

  return f;   // TODO: Leaking memory!
}

void EncodingEnigmaBN::testToFromMont(Field* montField, int num_trials) {
  Field* f = this->getSrcField();
  FieldEltArray* monts = montField->newEltArray(num_trials, true);
  FieldEltArray* normals = f->newEltArray(num_trials, true);

  for (int i = 0; i < num_trials; i++) {
    f->assignRandomElt(normals->elt(i));
  }

  TIMEREP(to_modular(((FieldEltDigits*)normals->elt(i))->get_digits(), 4, ((FieldEltDigits*)monts->elt(i))->get_digits(), ((FieldPrime*)montField)->members->msfield->modulo, members->PBIGCTX_PASS), num_trials, "ToMont", "Total");
  TIMEREP(from_modular(((FieldEltDigits*)monts->elt(i))->get_digits(), ((FieldEltDigits*)normals->elt(i))->get_digits(), ((FieldPrime*)montField)->members->msfield->modulo, members->PBIGCTX_PASS), num_trials, "FromMont", "Total");

  delete monts;
  delete normals;
}

Field* EncodingEnigmaBN::getBaseField(bool debug) {
  return this->base_field;
}

void EncodingEnigmaBN::makeBaseField(bool debug) {
  members->base_desc = new FieldDescription;
 
  // Create a field to represent operations in the base field
	members->base_p = new mp_modulus_t;   
  Field* baseField = NULL;
  if (!debug) {
    BOOL OK = create_modulus(members->bn.BaseField.modulo->modulus, 4, FROM_LEFT, members->base_p, members->PBIGCTX_PASS);
    assert(OK);

    members->base_desc->desc = new field_desc_t;
    OK = Kinitialize_prime(members->base_p, (field_desc_t*)members->base_desc->desc, members->PBIGCTX_PASS);
    assert(OK);

    baseField = new FieldPrime(members->base_desc, 254);
  } else {
    digit_t mod[1];
    mod[0] = 0x15C9B0796B71DB;
    BOOL OK = create_modulus(mod, 1, FROM_LEFT, members->base_p, members->PBIGCTX_PASS);
    assert(OK);

    members->base_desc->desc = new field_desc_t;
    OK = Kinitialize_prime(members->base_p, (field_desc_t*)members->base_desc->desc, members->PBIGCTX_PASS);
    assert(OK);

    baseField = new FieldPrime(members->base_desc, 53);
  }

  this->base_field = baseField;
}

EncodingEnigmaBN::~EncodingEnigmaBN() {
	this->del_elt(L, gen[L]);
	this->del_elt(R, gen[R]);
	this->del_elt(L, Zero[L]);
	this->del_elt(R, Zero[R]);

  this->base_field->delElt(this->montR);
  this->base_field->delElt(this->montRinv);

	delete srcField;
	delete members->msExponentField;
	uncreate_modulus(members->p, members->PBIGCTX_PASS);	
  delete members->p;

	if (members->base_p != NULL) {
		uncreate_modulus(members->base_p, members->PBIGCTX_PASS);
		delete members->base_p;
		delete members->base_desc->desc;
		delete members->base_desc;    
	}


  delete this->base_field;

	BN_free(&members->bn, members->PBIGCTX_PASS);	
  delete members;
}

uint32_t EncodingEnigmaBN::bytes_per(elt_t t) {
	if (t == L) {
		return num_L_digits*sizeof(digit_t);
	} else {
		assert(t == R);
		return num_R_digits*sizeof(digit_t);
	}
}

uint32_t EncodingEnigmaBN::num_exponent_bits() {
	return 4*sizeof(digit_t)*8 - 1;	// BN_members->curve uses 4 digits, but doesn't use them all (only uses 254 bits)
}


EncodedElt* EncodingEnigmaBN::new_elt(elt_t t) {
#ifdef PERF_DEBUG
  EncodingEnigmaBN::num_live_enc_elts[t]++;
#endif 
	if (t == L) {
		return (EncodedEltEnigmaBN*) new digit_t[num_L_digits];
	} else {
		assert(t == R);
		return (EncodedEltEnigmaBN*) new digit_t[num_R_digits];
	}
}

EncodedEltArray* EncodingEnigmaBN::new_elt_array(elt_t t, uint32_t len, bool preallocate) {
	return new EncodedEltEnigmaBNArray(t, len, preallocate);
}


void EncodingEnigmaBN::del_elt(elt_t t, EncodedElt* elt) {
#ifdef PERF_DEBUG
  EncodingEnigmaBN::num_live_enc_elts[t]--;
#endif 
	delete [] elt;
}

EncodedProduct* EncodingEnigmaBN::new_prod() {
	return new EncodedProductEnigmaBN();
}

uint32_t EncodingEnigmaBN::size_of_prod() {
  return sizeof(EncodedProductEnigmaBN);
}

void print_digits(digit_t* d, int count) {
  for (int i=count-1; i>= 0; i--) {
	  printf ("%08X%08X", (uint32_t) (d[i]>>(uint64_t)32), (uint32_t) d[i]);	  
  }
}

void EncodingEnigmaBN::print(elt_t t, EncodedElt* e) {  
  print(t, e, false);
}

void EncodingEnigmaBN::print(elt_t t, EncodedElt* e, bool affine) {
	assert(e);
	EncodedEltEnigmaBN* ee = (EncodedEltEnigmaBN*)e;
	
  if (!affine) {
    affinize(t, e, e);
  }
  int num_digits = t == L ? EncodingEnigmaBN::num_L_digits : EncodingEnigmaBN::num_R_digits;
  int num_affine_digits = (num_digits / 3) * 2;
  int num_elt_digits = (num_digits / 3) * 1;
  digit_t buffer[4];
  memset(((LEncodedEltEnigmaBN*)e)->get_digits() + num_affine_digits, 0, num_elt_digits * sizeof(digit_t));  // Zero out the Z-coordinate.

	if (t == L) {
    printf("(");
    for (int i = 0; i < 3; i++) {
      from_mont((FieldEltDigits*)ee->get_digits() + i * 4, (FieldEltDigits*)buffer);
      print_digits(buffer, 4);
      if (i == 0 || i == 1) {
        printf(",\n ");
      } else {
        printf(")");
      }
    }
	} else {
		assert(t == R);
    printf("(");
    for (int i = 0; i < 6; i++) {
      from_mont((FieldEltDigits*)ee->get_digits() + i * 4, (FieldEltDigits*)buffer);
      print_digits(buffer, 4);
      if (i == 1 || i == 3) {
        printf(",\n ");
      } else if (i == 5) {
        printf(")");
      } else {
        printf(" ");
      }
    }   
	}	
  if (!affine) {
    projectivize(t, e);
  }
}

void EncodingEnigmaBN::print(EncodedProduct* p) {
	assert(p);
  digit_t buffer[4];
	EncodedProductEnigmaBN* ep = (EncodedProductEnigmaBN*)p;
  for (int i = 0; i < 12; i++) {
    from_mont((FieldEltDigits*)ep->digits + 4 * i, (FieldEltDigits*)buffer);
    print_digits(buffer, 4);
    printf("\n");
  }
}

void EncodingEnigmaBN::copy(elt_t t, EncodedElt* src, EncodedElt* dst) {
	if (t == L) {
		memcpy(dst, src, num_L_digits*sizeof(digit_t));
	} else {
		assert(t == R);
		memcpy(dst, src, num_R_digits*sizeof(digit_t));
	}	
}

bool EncodingEnigmaBN::equal(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	EncodedEltEnigmaBN* ae = (EncodedEltEnigmaBN*)a;
	EncodedEltEnigmaBN* be = (EncodedEltEnigmaBN*)b;
	testOnCurve(t, a);
	testOnCurve(t, b);

	return ecprojective_equivalent(ae->get_digits(), be->get_digits(), members->curve[t], members->arithTempSpace, members->PBIGCTX_PASS) != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool	
}

void EncodingEnigmaBN::zero(elt_t t, EncodedElt* a) { 
	copy(t, Zero[t], a);	// Use our precomputed value	
}

void EncodingEnigmaBN::one(elt_t t, EncodedElt* a) {
	copy(t, gen[t], a);	// Use our precomputed value	
}

// Adds two encoded elements.  Result goes in r.
void EncodingEnigmaBN::add(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	addsub(t, a, b, r, Add);
}

// Adds two encoded elements.  Result goes in r.
void EncodingEnigmaBN::sub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r) {
	addsub(t, a, b, r, Sub);
}

// Adds two encoded elements.  Result goes in r.
void EncodingEnigmaBN::addsub(elt_t t, EncodedElt* a, EncodedElt* b, EncodedElt* r, addsub_t op) {
	assert(a && b && r);
	testOnCurve(t, a);
	testOnCurve(t, b);

	EncodedEltEnigmaBN* ae = (EncodedEltEnigmaBN*)a;
	EncodedEltEnigmaBN* be = (EncodedEltEnigmaBN*)b;
	EncodedEltEnigmaBN* re = (EncodedEltEnigmaBN*)r;

	assert(sizeof(members->arithTempSpace) >= members->curve[t]->fdesc->ndigtemps_arith * 2 * sizeof(digit_t));

	int the_op = op == Add ? 1 : -1;
	BOOL OK = ecprojective_addsub_projective(ae->get_digits(), be->get_digits(), re->get_digits(), the_op, members->curve[t], members->arithTempSpace, members->PBIGCTX_PASS);		
	assert(OK);
	testOnCurve(t, re);
  if (collecting_stats) { this->num_addsub[t]++; }
}

void EncodingEnigmaBN::doubleIt(elt_t t, EncodedElt* a, EncodedElt* b) {
	assert(a && b);
	EncodedEltEnigmaBN* ae = (EncodedEltEnigmaBN*)a;
	EncodedEltEnigmaBN* be = (EncodedEltEnigmaBN*)b;
	testOnCurve(t, a);
	BOOL OK = ecprojective_double(ae->get_digits(), be->get_digits(), members->curve[t], members->arithTempSpace, members->PBIGCTX_PASS);
	assert(OK);
	testOnCurve(t, b);
  if (collecting_stats) { this->num_double[t]++; }
}


// Multiplies the value encoded in g by the constant c.  Result goes in r.
// For pairings, we compute r <- g^c
void EncodingEnigmaBN::mul(elt_t t, EncodedElt* g, FieldElt* c, EncodedElt* r) {		
	assert(g && c && r);
	testOnCurve(t, g);
	
	FieldEltPrime* cp = (FieldEltPrime*)c;
	EncodedEltEnigmaBN* ge = (EncodedEltEnigmaBN*)g;
	EncodedEltEnigmaBN* re = (EncodedEltEnigmaBN*)r;

	BOOL OK = ecproj_single_exponentiation(ge->get_digits(), cp->get_digits(), members->bn.pbits, members->curve[t]->lnggorder, re->get_digits(), members->curve[t], members->PBIGCTX_PASS);
	assert(OK);
	testOnCurve(t, r);
  if (collecting_stats) { this->num_constmul[t]++; }
}

void EncodingEnigmaBN::mul(LEncodedElt* a, REncodedElt* b, EncodedProduct* r) {
	pair(a, b, r);
}

void EncodingEnigmaBN::pair(LEncodedElt* L1, REncodedElt* R1, EncodedProduct* result) {
	assert(L1 && R1 && result);
	testOnCurve(L, L1);
	testOnCurve(R, R1);

	EncodedEltEnigmaBN* ae = (EncodedEltEnigmaBN*)a;
	EncodedEltEnigmaBN* be = (EncodedEltEnigmaBN*)b;
	EncodedProductEnigmaBN* resulte = (EncodedProductEnigmaBN*)result;

	BOOL OK;
	// Pairing library doesn't handle g^0 values properly, so check for them
	if (isZero(L, L1) || isZero(R, R1)) {
		OK = Kimmediate(1, resulte->digits, &members->bn.ExtField, members->PBIGCTX_PASS); 
		assert(OK);
	} else {
		// Despite the name, bn_o_ate_proj expects affine coordinates in
		EncodedEltEnigmaBN* Laff = (EncodedEltEnigmaBN*)new_elt(L);
		EncodedEltEnigmaBN* Raff = (EncodedEltEnigmaBN*)new_elt(R);

		affinize(L, L1, Laff);
		affinize(R, R1, Raff);
		OK = bn_o_ate_proj(Raff->get_digits(), Laff->get_digits(), &members->bn, resulte->digits, members->PBIGCTX_PASS);
		assert(OK);
		OK = OK && bn_o_ate_finalexp(resulte->digits, &members->bn, resulte->digits, members->PBIGCTX_PASS);
		assert(OK);
		del_elt(L, Laff);
		del_elt(R, Raff);
	}
  if (collecting_stats) { this->num_pair++; }
}

bool EncodingEnigmaBN::isZero(elt_t t, EncodedElt* elt) {
	assert(elt);
	testOnCurve(t, elt);
	EncodedEltEnigmaBN* elte = (EncodedEltEnigmaBN*)elt;
	
	digit_t temp[260];
	assert(members->bn.BaseCurve.ndigtemps <= sizeof(temp));
	
	BOOL equal = ecprojective_equivalent(elte->get_digits(), Zero[t]->get_digits(), members->curve[t], temp, members->PBIGCTX_PASS);	
	return equal != 0;
}


void EncodingEnigmaBN::add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductEnigmaBN* ae = (EncodedProductEnigmaBN*)a;
	EncodedProductEnigmaBN* be = (EncodedProductEnigmaBN*)b;
	EncodedProductEnigmaBN* re = (EncodedProductEnigmaBN*)r;

	digit_t tmps[500];
	assert(members->bn.ExtField.ndigtemps_arith <= sizeof(tmps));

	BOOL OK = K12mul(ae->digits, be->digits, re->digits, &members->bn, tmps, members->PBIGCTX_PASS);		
	assert(OK);
}

void EncodingEnigmaBN::sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r) {
	assert(a && b && r);
	EncodedProductEnigmaBN* ae = (EncodedProductEnigmaBN*)a;
	EncodedProductEnigmaBN* be = (EncodedProductEnigmaBN*)b;
	EncodedProductEnigmaBN* re = (EncodedProductEnigmaBN*)r;

	EncodedProductEnigmaBN tmp;
	digit_t tmps[500];
	assert(members->bn.ExtField.ndigtemps_arith <= sizeof(tmps));

	BOOL OK = K12inv(be->digits, tmp.digits, &members->bn, tmps, members->PBIGCTX_PASS);
	assert(OK);
	OK = K12mul(ae->digits, tmp.digits, re->digits, &members->bn, tmps, members->PBIGCTX_PASS);		
	assert(OK);
}

bool EncodingEnigmaBN::equals(EncodedProduct* a, EncodedProduct* b) {
	assert(a && b);
	EncodedProductEnigmaBN* ae = (EncodedProductEnigmaBN*)a;
	EncodedProductEnigmaBN* be = (EncodedProductEnigmaBN*)b;

	return 0 != Kequal(ae->digits, be->digits, &members->bn.ExtField, members->PBIGCTX_PASS);
}

// Checks whether the plaintexts inside the encodings obey:
// L1*R1 - L2*R2 == L3*R3
bool EncodingEnigmaBN::mulSubEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2, LEncodedElt* L3, REncodedElt* R3) {
	assert(L1 && R1 && L2 && R2 && L3 && R3);
	EncodedProductEnigmaBN pairing1;
	EncodedProductEnigmaBN pairing2;
	EncodedProductEnigmaBN pairing3;

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

	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}


// Checks whether the plaintexts inside the encodings obey:
// L1*R1 == L2*R2
bool EncodingEnigmaBN::mulEquals(LEncodedElt* L1, REncodedElt* R1, LEncodedElt* L2, REncodedElt* R2) {
	assert(L1 && R1 && L2 && R2);
	EncodedProductEnigmaBN pairing1;
	EncodedProductEnigmaBN pairing2;

	testOnCurve(L, L1); testOnCurve(R, R1);
	testOnCurve(L, L2); testOnCurve(R, R2);

	// Compute the two pairings
	pair(L1, R1, &pairing1);
	pair(L2, R2, &pairing2);

	BOOL result = Kequal(pairing1.digits, pairing2.digits, &members->bn.ExtField, members->PBIGCTX_PASS);
	
	return result != 0;  // The != 0 supresses a compiler warning about casting a BOOL (int) to a bool
}

void setDigitsToZero(digit_t* digits, int lenInDigits) {
  for (int i = 0; i < lenInDigits; i++) {
    digits[i] = 0;
  }
}

// Expects out to point to a buffer with at least condensed_size_in_digits*sizeof(digit_t) space
void EncodingEnigmaBN::compressLaff(bool elt_is_zero, digit_t* aff_in, digit_t* out) {
  digit_t* xCoeff = ((digit_t*)aff_in) + 0;
  digit_t* yCoeff = ((digit_t*)aff_in) + 4;

  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (elt_is_zero) {
    setDigitsToZero(out, condensed_size_in_digits_L);
    return;
  }

  // Compute x^3 + b, where b = 2
  digit_t expon = 3;
  digit_t xCubed[4];
  BOOL OK = Kexpon(xCoeff, &expon, 1, xCubed, &members->bn.BaseField, members->PBIGCTX_PASS);
  assert(OK);

  digit_t two[4];
  OK = Kimmediate(2, two, &members->bn.BaseField, members->PBIGCTX_PASS);
  assert(OK);

  OK = Kadd(xCubed, two, xCubed, &members->bn.BaseField, members->PBIGCTX_PASS);
  assert(OK);

  // Compute positive root
  digit_t sqrtRhs[4];
  BOOL is_square;
  OK = Ksqrt(xCubed, sqrtRhs, &members->bn.BaseField, &is_square, members->PBIGCTX_PASS);
  assert(OK);
  assert(is_square);

#ifdef DEBUG_EC
  // Check the square root is correct
  digit_t tmps[2000];
  digit_t allegedRhs[4];
  OK = Kmul(sqrtRhs, sqrtRhs, allegedRhs, &members->bn.BaseField, tmps, PBIGCTX_PASS);
  assert(OK);
  assert(Kequal(allegedRhs, xCubed, &members->bn.BaseField, PBIGCTX_PASS));
#endif

  BOOL equal = Kequal(sqrtRhs, yCoeff, &members->bn.BaseField, members->PBIGCTX_PASS);

  // We use the top-most bit of X to store the bit indicating whether to use the negative root
  int lastDigit = aff_size_in_digits_L / 2 - 1;	// EncodedElt contains 2 field elements
  if (!equal) {
    // Set the top bit of the last digit of X to indicate that we use the "negative" root		
    xCoeff[lastDigit] |= ((digit_t)1) << 63;
  }
  else {
    // Make sure the top bit is actually clear, indicating the "positive" root		
    assert(0 == (xCoeff[lastDigit] & (((digit_t)1) << 63)));
  }

  memcpy(out, xCoeff, condensed_size_in_digits_L * sizeof(digit_t));
}

// On the twist curve, we use y^2 = x^3 + (1 - i), and each point is represented as: a + i*b
void EncodingEnigmaBN::compressRaff(bool elt_is_zero, digit_t* aff_in, digit_t* out) {
  digit_t* xCoeff = ((digit_t*)aff_in) + 0;
  digit_t* yCoeff = ((digit_t*)aff_in) + 8;

  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (elt_is_zero) {
    setDigitsToZero(out, condensed_size_in_digits_R);
    return;
  }

  // Compute x^3 + (1 - i)
  digit_t expon = 3;
  digit_t xCubed[8];
  BOOL OK = Kexpon(xCoeff, &expon, 1, xCubed, &members->bn.TwistField, members->PBIGCTX_PASS);
  assert(OK);

  digit_t oneNegOne[8];
  OK = Kimmediate(1, oneNegOne, &members->bn.BaseField, members->PBIGCTX_PASS);	// Put 1 into first slot
  assert(OK);
  OK = Knegate(oneNegOne, ((digit_t*)oneNegOne) + 4, &members->bn.BaseField, members->PBIGCTX_PASS);	// Put -1 into the second slot
  assert(OK);

  OK = Kadd(xCubed, oneNegOne, xCubed, &members->bn.TwistField, members->PBIGCTX_PASS);
  assert(OK);

  // Compute a root
  digit_t sqrtRhs[8];
  BOOL is_square;
  OK = Ksqrt(xCubed, sqrtRhs, &members->bn.TwistField, &is_square, members->PBIGCTX_PASS);
  assert(OK);
  assert(is_square);

#ifdef DEBUG_EC
  // Check the square root is correct
  digit_t tmps[2000];
  digit_t allegedRhs[8];
  OK = Kmul(sqrtRhs, sqrtRhs, allegedRhs, &members->bn.TwistField, tmps, members->PBIGCTX_PASS);
  assert(OK);
  assert(Kequal(allegedRhs, xCubed, &members->bn.TwistField, members->PBIGCTX_PASS));
#endif

  // Compute the other root
  digit_t negSqrtRhs[8];
  OK = Knegate(sqrtRhs, negSqrtRhs, &members->bn.TwistField, members->PBIGCTX_PASS);
  assert(OK);

  BOOL equal = Kequal(sqrtRhs, yCoeff, &members->bn.TwistField, members->PBIGCTX_PASS);

  // Since we can't actually test for negative roots, and the library gives us a random one,
  // we distinguish between them based on which last digit is smaller
  bool equalsSmaller = (equal && (sqrtRhs[7] <= negSqrtRhs[7])) || (!equal && (sqrtRhs[7] > negSqrtRhs[7]));

  // We use the top-most bit of X to store the bit indicating whether to use the smaller root
  int lastDigit = aff_size_in_digits_R / 2 - 1;	// EncodedElt contains 2 field elements
  if (equalsSmaller) {
    // Set the top bit of the last digit of X to indicate that we use the smaller root
    xCoeff[lastDigit] |= ((digit_t)1) << 63;
  }
  else {
    // Make sure the top bit is actually clear, indicating the larger root
    assert(!(xCoeff[lastDigit] & (((digit_t)1) << 63)));
  }

  memcpy(out, xCoeff, aff_size_in_digits_R*sizeof(digit_t) / 2);
}


void EncodingEnigmaBN::compress(elt_t t, EncodedElt* e) {
  EncodedEltEnigmaBN* elt = (EncodedEltEnigmaBN*)e;
  bool elt_is_zero = isZero(t, e);
  affinize(t, elt, elt);
  if (t == L) {
    compressLaff(elt_is_zero, elt->get_digits(), elt->get_digits());
  }
  else {
    assert(t == R);
    compressRaff(elt_is_zero, elt->get_digits(), elt->get_digits());
  }
}

bool isAllZero(digit_t* value, int num_digits) {
  for (int i = 0; i < num_digits; i++) {
    if (value[i] != 0) {
      return false;
    }
  }
  return true;
}

void EncodingEnigmaBN::projectivize(elt_t t, EncodedElt* elt) {
  assert(elt);

  if (t == L) {
    BOOL OK = ecaffine_projectivize(((LEncodedEltEnigmaBN*)elt)->get_digits(), ((LEncodedEltEnigmaBN*)elt)->get_digits(), &members->bn.BaseCurve, members->arithTempSpace, members->PBIGCTX_PASS);
    assert(OK);
  } else {
    assert(t == R);
    BOOL OK = ecaffine_projectivize(((REncodedEltEnigmaBN*)elt)->get_digits(), ((REncodedEltEnigmaBN*)elt)->get_digits(), &members->bn.TwistCurve, members->arithTempSpace, members->PBIGCTX_PASS);
    assert(OK);
  }
    
  testOnCurve(t, elt);
}

void EncodingEnigmaBN::decompressL(digit_t* in, LEncodedElt* out) {
  assert(in);
  assert(out);
  
  digit_t* xCoeff = ((LEncodedEltEnigmaBN*)out)->get_digits() + 0;
  digit_t* yCoeff = ((LEncodedEltEnigmaBN*)out)->get_digits() + 4;

  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (isAllZero(in, condensed_size_in_digits_L)) {
    this->zero(L, (LEncodedEltEnigmaBN*)out);
    return;
  }

  // Check whether the top bit of X is set.  Indicates which root to use for Y
  int lastDigit = aff_size_in_digits_L / 2 - 1;	// EncodedElt contains 2 field elements
  bool useNegRoot = (in[lastDigit] & (((digit_t)1) << 63)) > 0;

  // Either way, clear out the top bit and copy the input into the x coefficient
  in[lastDigit] &= ((((digit_t)1) << 63) - 1);
  memcpy(xCoeff, in, condensed_size_in_digits_L*sizeof(digit_t));

  // Compute x^3 + b, where b = 2
  digit_t expon = 3;
  digit_t xCubed[4];
  BOOL OK = Kexpon(xCoeff, &expon, 1, xCubed, &members->bn.BaseField, members->PBIGCTX_PASS);
  assert(OK);

  digit_t two[4];
  OK = Kimmediate(2, two, &members->bn.BaseField, members->PBIGCTX_PASS);
  assert(OK);
  OK = Kadd(xCubed, two, xCubed, &members->bn.BaseField, members->PBIGCTX_PASS);
  assert(OK);

  // Compute "positive" root and, if necessary, use it to compute the negative root
  BOOL is_square;
  OK = Ksqrt(xCubed, yCoeff, &members->bn.BaseField, &is_square, members->PBIGCTX_PASS);
  assert(OK);
  assert(is_square);

  if (useNegRoot) {
    OK = Knegate(yCoeff, yCoeff, &members->bn.BaseField, members->PBIGCTX_PASS);
    assert(OK);
  }

  projectivize(L, out);
  testOnCurve(L, out);
}

// On the twist curve, we use y^2 = x^3 + (1 - i), and each point is represented as a + i*b
void EncodingEnigmaBN::decompressR(digit_t* in, REncodedElt* out) {
  digit_t* xCoeff = ((REncodedEltEnigmaBN*)out)->get_digits() + 0;
  digit_t* yCoeff = ((REncodedEltEnigmaBN*)out)->get_digits() + 8;

  // Special-case the point at infinity (since it isn't a solution to the curve equation)
  if (isAllZero(in, condensed_size_in_digits_R)) {
    this->zero(R, (REncodedEltEnigmaBN*)out);
    return;
  }

  // Check whether the top bit of X is set.  Indicates which root to use for Y
  int lastDigit = aff_size_in_digits_R / 2 - 1;	// EncodedElt contains 2 field elements
  bool useSmallerRoot = (in[lastDigit] & (((digit_t)1) << 63)) > 0;

  // Either way, clear out the top bit and copy the input into the x coefficient
  in[lastDigit] &= ((((digit_t)1) << 63) - 1);
  memcpy(xCoeff, in, condensed_size_in_digits_R*sizeof(digit_t));

  // Compute x^3 + (1 - i)
  digit_t expon = 3;
  digit_t xCubed[8];
  BOOL OK = Kexpon(xCoeff, &expon, 1, xCubed, &members->bn.TwistField, members->PBIGCTX_PASS);
  assert(OK);

  digit_t oneNegOne[8];
  OK = Kimmediate(1, oneNegOne, &members->bn.BaseField, members->PBIGCTX_PASS);	// Put 1 into first slot
  assert(OK);
  OK = Knegate(oneNegOne, ((digit_t*)oneNegOne) + 4, &members->bn.BaseField, members->PBIGCTX_PASS);	// Put -1 into the second slot
  assert(OK);

  OK = Kadd(xCubed, oneNegOne, xCubed, &members->bn.TwistField, members->PBIGCTX_PASS);
  assert(OK);

  // Compute a root and use it to compute the other root
  BOOL is_square;
  OK = Ksqrt(xCubed, yCoeff, &members->bn.TwistField, &is_square, members->PBIGCTX_PASS);
  assert(OK);
  assert(is_square);

  digit_t otherRoot[8];
  OK = Knegate(yCoeff, otherRoot, &members->bn.TwistField, members->PBIGCTX_PASS);
  assert(OK);

#ifdef DEBUG_EC
  // Check the other square root is correct
  digit_t tmps[2000];
  digit_t allegedRhs[8];
  OK = Kmul(otherRoot, otherRoot, allegedRhs, &members->bn.TwistField, tmps, members->PBIGCTX_PASS);
  assert(OK);
  assert(Kequal(allegedRhs, xCubed, &members->bn.TwistField, members->PBIGCTX_PASS));
#endif

  if ((useSmallerRoot && otherRoot[7] <= yCoeff[7]) ||
    (!useSmallerRoot && otherRoot[7] > yCoeff[7])) {
    memcpy(yCoeff, otherRoot, 8 * sizeof(digit_t));
  }

  projectivize(R, out);  
  testOnCurve(R, out);
}

void EncodingEnigmaBN::decompress(elt_t t, digit_t* compressed_elt, EncodedElt* elt) {
  assert(compressed_elt);
  assert(elt);
  if (t == L) {
    decompressL(compressed_elt, elt);
  } else {
    assert(t == R);
    decompressR(compressed_elt, elt);
  }
}

void EncodingEnigmaBN::to_mont(FieldElt* normal, FieldElt* mont) {
  base_field->mul(normal, this->montR, mont);
}

void EncodingEnigmaBN::from_mont(FieldElt* mont, FieldElt* normal) {
  base_field->mul(mont, this->montRinv, normal);
}

void EncodingEnigmaBN::write(Archiver* arc, elt_t t, EncodedElt* e, bool simple) {
	assert(arc && e);

	unsigned char* data = (unsigned char*)e;
#ifdef CONDENSE_ARCHIVES
  int num_digits = t == L ? EncodingEnigmaBN::condensed_size_in_digits_L : EncodingEnigmaBN::condensed_size_in_digits_R;
  compress(t, e);
#else
	int num_digits = t == L ? EncodingEnigmaBN::num_L_digits : EncodingEnigmaBN::num_R_digits;
  int num_affine_digits = (num_digits / 3) * 2;
  int num_elt_digits = (num_digits / 3) * 1;

  if (simple) {
    affinize(t, e, e);  // C code expects affine coordinates
    num_digits = num_affine_digits;  // Only need 2/3 the space for affine.  

    // Convert each remaining Fp elt into Normal form
    for (int i = 0; i < num_affine_digits / 4; i++) {
      digit_t* ptr = ((LEncodedEltEnigmaBN*)e)->get_digits() + i * 4;
      from_mont((FieldEltDigits*)ptr, (FieldEltDigits*)ptr);
    }
  }
#endif // CONDENSE_ARCHIVES	
	arc->write((unsigned char*)data, num_digits * sizeof(digit_t));

#ifndef CONDENSE_ARCHIVES
  if (simple) {
    // Convert each Fp elt back into Montgomery form
    for (int i = 0; i < num_affine_digits / 4; i++) {
      digit_t* ptr = ((LEncodedEltEnigmaBN*)e)->get_digits() + i * 4;
      to_mont((FieldEltDigits*)ptr, (FieldEltDigits*)ptr);
    }
    // Restore the point to projective form, in case we want to keep using it
    projectivize(t, e);
  }
#endif //CONDENSE_ARCHIVES
}

void EncodingEnigmaBN::read (Archiver* arc, elt_t t, EncodedElt* e) {
	assert(arc && e);
  EncodedEltEnigmaBN* ee = (EncodedEltEnigmaBN*)e;

	unsigned char* data = (unsigned char*)e;
#ifdef CONDENSE_ARCHIVES
  int num_digits = t == L ? EncodingEnigmaBN::condensed_size_in_digits_L : EncodingEnigmaBN::condensed_size_in_digits_R;
  arc->read((unsigned char*)data, num_digits * sizeof(digit_t));
  decompress(t, ee->get_digits(), e);
#else
	int num_digits = t == L ? EncodingEnigmaBN::num_L_digits : EncodingEnigmaBN::num_R_digits;
	arc->read((unsigned char*)data, num_digits * sizeof(digit_t));
#endif // CONDENSE_ARCHIVES
}


void EncodingEnigmaBN::affinize(elt_t t, EncodedElt* in, EncodedElt* out) {
	assert(in && out);
	testOnCurve(t, in);
	EncodedEltEnigmaBN* ine  = (EncodedEltEnigmaBN*)in;
	EncodedEltEnigmaBN* oute = (EncodedEltEnigmaBN*)out;

	BOOL OK = ecprojective_affinize(ine->get_digits(), oute->get_digits(), members->curve[t], members->arithTempSpace, members->PBIGCTX_PASS);
	assert(OK);
}

Field* EncodingEnigmaBN::get_tiny_bn_field() {
  digit_t prime = 0x15c9b074c05625;

  // Create a small prime field
  mp_modulus_t* p = new mp_modulus_t;
  BOOL OK = create_modulus(&prime, 1, FROM_LEFT, p, members->PBIGCTX_PASS);
  assert(OK);

  field_desc_t* msExponentField = new field_desc_t;

  OK = Kinitialize_prime(p, msExponentField, members->PBIGCTX_PASS);
  assert(OK);

  FieldDescription* desc = new FieldDescription;
  desc->desc = msExponentField;
  FieldPrime* field = new FieldPrime(desc, 53);
  delete desc;

  return field;
}

void EncodingEnigmaBN::testOnCurve(elt_t t, EncodedElt* a) {
#ifdef DEBUG_EC
	assert(a);
	EncodedEltEnigmaBN* ae = (EncodedEltEnigmaBN*)a;
	digit_t temp[260];
	BOOL OK = ecprojective_on_curve(ae->get_digits(), members->curve[t], "I'm on a members->curve", temp, members->PBIGCTX_PASS);
	assert(OK);
#endif
}

void EncodingEnigmaBN::encodeSlow(elt_t t, FieldElt* in, EncodedElt* out) {
	assert(in && out);
	FieldEltPrime* ine = (FieldEltPrime*)in;
	EncodedEltEnigmaBN* oute = (EncodedEltEnigmaBN*)out;

	BOOL OK = ecproj_single_exponentiation(gen[t]->get_digits(), ine->get_digits(), members->bn.pbits, members->curve[t]->lnggorder, oute->get_digits(), members->curve[t], members->PBIGCTX_PASS);
	assert(OK);
}

void EncodingEnigmaBN::print_stats() {
  printf("addsub: %d (%d L, %d R)\n", this->num_addsub[L] + this->num_addsub[R], this->num_addsub[L], this->num_addsub[R]);
  printf("double: %d (%d L, %d R)\n", this->num_double[L] + this->num_double[R], this->num_double[L], this->num_double[R]);
  printf("constmul: %d (%d L, %d R)\n", this->num_constmul[L] + this->num_constmul[R], this->num_constmul[L], this->num_constmul[R]);
  
  printf("Pair: %d\n", this->num_pair);

  if (cache_attempts > 0) {
    printf("EncodedElt cache hit %d times out of %d attempts (%0.2f%%)\n", cache_hits, cache_attempts, 100 * (cache_hits / (float)cache_attempts));
  }
}

void EncodingEnigmaBN::reset_stats() {
  this->num_addsub[0] = 0;
  this->num_double[0] = 0;
  this->num_constmul[0] = 0;
  this->num_addsub[1] = 0;
  this->num_double[1] = 0;
  this->num_constmul[1] = 0;
  this->num_pair = 0;
}

void EncodingEnigmaBN::collect_stats(bool on) {
  this->collecting_stats = on;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//			Arrays
//////////////////////////////////////////////////////////////////////////////////////////////////////

EncodedEltEnigmaBNArray::EncodedEltEnigmaBNArray(elt_t t, uint32_t len, bool allocate)
  : EncodedDigitArray(t, len, allocate, t == L ? EncodingEnigmaBN::num_L_digits : EncodingEnigmaBN::num_R_digits) {
    ; // Parent does all the work
}

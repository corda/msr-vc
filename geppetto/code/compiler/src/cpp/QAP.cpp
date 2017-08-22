#include "Keys.h"
#include "Proof.h"
#include "Field.h"
#include "SparsePolynomial.h"
#include "Commitment.h"
#include "FileArchiver.h"
#include "Config.h"
#include "PolyGeometricInterpolator.h"

#include <iostream>
#include <cassert>
#include <memory>

//#define PERF_DEBUG
//#define TEST_SERIALIZATION

extern "C" {

// Call to initialize FFLib's settings
//         max_mem_gb: Use at most this many GB for crypto table storage (approximate)
//       precomp_free: Should we consider anything precomputable free?  If so, we'll spend much more time computing large tables
//      verbose_level: How verbose should FFLib be?  0 = silent.  Larger values are more verbose (see verbosity enum above)
//   machine_readable: Should the library's output be in a machine friendly format?  0 = human, >0 = machine

_CrtMemState initial_state;

__declspec(dllexport) void configure(int max_mem_gb, int precomp_free, int verbose_level, int machine_readable) {
  Config::max_mem = max_mem_gb;
  Config::precomp_free = precomp_free > 0;
  Config::verbose_level = verbose_level;
  Config::machine_readable = machine_readable > 0;

#ifdef MEM_LEAK
  _CrtMemCheckpoint(&initial_state);  
  //_CrtSetBreakAlloc(18338);
#endif 
}

extern void print_field_stats();  // Defined in FSharpIfc.c

// Please call this before exiting from F#
__declspec(dllexport) void cleanup() {
  timers.cleanup();
#ifdef MEM_LEAK
  //_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  //_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
  //_CrtDumpMemoryLeaks();  // Check for memory leaks
  _CrtMemDumpAllObjectsSince(&initial_state);
#endif
#ifdef PERF_DEBUG
  print_field_stats();
#endif // PERF_DEBUG
  if (Config::verbose_level >= INFO) {
    Encoding::current_encoding->print_stats();
  }
}

//////////////////////////////////////////////////////////////////////////////
//				Key generation helper routines
//////////////////////////////////////////////////////////////////////////////
FieldEltArray* createLagrangeDenominators(FieldEltArray* targets, int num_roots, Field* field)
{
	unique_ptr<FieldEltArray> denominators(field->newEltArray(num_roots, true));
	unique_ptr<FieldEltArray> u(field->newEltArray(num_roots, true));
  FieldElt* tmp = field->newElt();
  FieldElt* one = field->newElt();
  FieldElt* zero = field->newElt();

	// We assume that targets_i = ratio^i. 
	// This helps us compute the denominators in linear time.
	field->one(one);
	field->zero(zero);

	// u[i] = (ratio^i - 1) * u[i-1]
	field->one(u->elt(0));
	for (int i = 1; i < num_roots; ++i)
	{
		field->sub(targets->elt(i), one, tmp);
		field->mul(u->elt(i - 1), tmp, u->elt(i));
	}

  FieldElt* gpref = field->newElt();
  FieldElt* gsuff = field->newElt();
	FieldElt* suffd = field->newElt();
	field->one(gpref);
	field->one(gsuff);
	if(num_roots>=2) field->copy(targets->elt(num_roots-2), suffd);

	// At iteration i: gpref = ratio^(i*(i-1)/2), gsuff = ratio^(i*(d-i)), where d = num_roots-1
	for (int i = 0; i < num_roots; ++i)
	{
		field->mul(u->elt(i), u->elt(num_roots - 1 - i), tmp);
		if ((num_roots - 1 - i) % 2) field->sub(zero, tmp, tmp);
		field->mul(tmp, gpref, tmp);
		field->mul(tmp, gsuff, tmp);
		field->copy(tmp, denominators->elt(i));

		// Increments
		if (i + 1 < num_roots)
		{
			field->mul(gpref, targets->elt(i), gpref);
			field->mul(gsuff, suffd, gsuff);
			if (i + 2 < num_roots) 
				field->div(suffd, targets->elt(2), suffd);
		}
	}

  field->delElt(tmp);
  field->delElt(one);
  field->delElt(zero);
  field->delElt(gpref);
  field->delElt(gsuff);
  field->delElt(suffd);

	return denominators.release();
}

FieldEltArray* evalLagrange(int degree, Field* field, FieldEltArray* targetRoots, FieldElt* evaluationPoint, FieldEltArray* denominators) {
	FieldEltArray* L = field->newEltArray(degree, true);
	FieldEltArray* diff = field->newEltArray(degree, true);
	
	// First, calculate the "master" numerator, i.e., \prod_k (x - r_k)
	FieldElt* superNumerator = field->newElt();
	field->one(superNumerator);

	for (int k = 0; k < degree; k++) { // For each r_k
		field->sub(evaluationPoint, targetRoots->elt(k), diff->elt(k));
		field->mul(superNumerator, diff->elt(k), superNumerator);	// superNumerator *= (evalPt - root[k])
	}

	// Finally, calculate the individual Lagrange polynomial values
	FieldElt* numerator = field->newElt();
	for (int i = 0; i < degree; i++) {	// For each L_i
		// Remove (x-r_k) if k == i
		
		field->div(superNumerator, diff->elt(i), numerator);	

		// Do the final calculation
		field->div(numerator, denominators->elt(i), L->elt(i));
	}

  field->delElt(numerator);		
	field->delElt(superNumerator);
  delete diff;
	
	return L;
}

void evalTargetPoly(FieldElt* evaluationPoint, FieldElt* result, FieldEltArray* targetRoots, int numTargetRoots, Field* field) {
	field->one(result);

	FieldElt* tmp = field->newElt();
	for (int i = 0; i < numTargetRoots; i++) {
		field->sub(evaluationPoint, targetRoots->elt(i), tmp);
		field->mul(result, tmp, result);	// result += evalPt * root[i]		
	}
  field->delElt(tmp);
}

//////////////////////////////////////////////////////////////////////////////
//				Key-manipulation routines
//////////////////////////////////////////////////////////////////////////////
enum key_selector { EK=1, VK=2, BOTH=3 };
// Set simple to true to save in a format the C version will understand
__declspec(dllexport) void save_key(Encoding* encoding, KeyPair* keys, key_selector which_key, char* file_name, bool simple=false) {
  assert(encoding);
  assert(keys);
  FileArchiver file_archiver(file_name, FileArchiver::Write);
  Archive arc(&file_archiver, encoding);
  
  if (which_key == EK || which_key == BOTH) {
    assert(keys->ek);
    keys->ek->serialize(&arc, simple);
  }

  if (which_key == VK || which_key == BOTH) {
    assert(keys->vk);
    keys->vk->serialize(&arc, simple);
    //keys->vk->print();
  }
}

__declspec(dllexport) KeyPair* load_key(Encoding* encoding, key_selector which_key, char* file_name) {
  assert(encoding);
  KeyPair* keys = new KeyPair();

  FileArchiver file_archiver(file_name, FileArchiver::Read);
  Archive arc(&file_archiver, encoding);

  if (which_key == EK || which_key == BOTH) {
    assert(keys->ek);
    keys->ek->deserialize(encoding, &arc);
  }

  if (which_key == VK || which_key == BOTH) {
    assert(keys->vk);
    keys->vk->deserialize(encoding, &arc);
  }

  return keys;
}

__declspec(dllexport) void print_key(Encoding* encoding, KeyPair* keys, key_selector which_key) {
  assert(encoding);
  assert(keys);

  if (which_key == EK || which_key == BOTH) {
    assert(keys->ek);
    keys->ek->print();
  }

  if (which_key == VK || which_key == BOTH) {
    assert(keys->vk);
    keys->vk->print();
  }
}

//////////////////////////////////////////////////////////////////////////////
//				Key-generation routines
//////////////////////////////////////////////////////////////////////////////

class KeyGenContext {
public:
  KeyGenContext(Encoding* encoding, int num_banks);
  ~KeyGenContext();

  KeyPair* kp;
  Encoding* encoding;
  Field* field;
  int num_banks;

  FieldElt* s;  // Secret evaluation point
  FieldElt *rV, *rW, *rY;   // Create customized versions of the generator

  FieldElt* current_poly;     // Accumulate the current v_k, w_k, or y_k
  FieldElt* beta_accumulator; // Accumulate the v_k*w_k*y_k used for the beta term
  FieldElt* tmp;
  FieldElt* tmp2;

  FieldEltArray* alphaV;
  FieldEltArray* alphaW;
  FieldEltArray* alphaY;
  FieldEltArray* beta;
  FieldEltArray* gamma;

  FieldEltArray* lagrange;  // L_i(s)

  // Track the state machine used when streaming the polynomials
  int current_bank;
  int current_coeff_index;
  poly_selector_t current_poly_type;

};

KeyGenContext::KeyGenContext(Encoding* enc, int num_banks) {
  encoding = enc;
  field = encoding->getSrcField();
  num_banks = num_banks;

  this->s  = field->newElt();
	this->rV = field->newElt();
	this->rW = field->newElt();
	this->rY = field->newElt();

  this->current_poly = field->newElt();
  this->beta_accumulator = field->newElt();
  this->tmp  = field->newElt();
  this->tmp2 = field->newElt();

  this->alphaV = field->newEltArray(num_banks, true);
  this->alphaW = field->newEltArray(num_banks, true);
  this->alphaY = field->newEltArray(num_banks, true);
  this->beta   = field->newEltArray(num_banks, true);
  this->gamma  = field->newEltArray(num_banks, true);

  this->current_bank = 0;
  this->current_coeff_index = 0;
  this->current_poly_type = V;

  field->zero(this->current_poly);
  field->zero(this->beta_accumulator);
}

KeyGenContext::~KeyGenContext() {
  field->delElt(s);
  field->delElt(rV);
  field->delElt(rW);
  field->delElt(rY);

  field->delElt(current_poly);
  field->delElt(beta_accumulator);
  field->delElt(tmp);
  field->delElt(tmp2);

  delete alphaV;
  delete alphaW;
  delete alphaY;
  delete beta;
  delete gamma;

  delete lagrange;
}

__declspec(dllexport) KeyGenContext* key_gen_init(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[]) {
  Timer* timer_init = timers.newTimer("KeyGenInit", NULL);  timer_init->start();
  KeyGenContext* context = new KeyGenContext(encoding, num_banks);
  Field* field = encoding->getSrcField();
  context->kp = new KeyPair(encoding, num_roots, num_banks, bank_sizes, bank_type);

  Timer* timer_rand = timers.newTimer("GenRandom", "KeyGenInit"); timer_rand->start();

  // Select a secret evaluation point (s)	
	field->assignRandomElt(context->s);

  // Select random generators
	field->assignRandomElt(context->rV);
	field->assignRandomElt(context->rW);
  field->mul(context->rV, context->rW, context->rY);

	// Select secret multiples to ensure extractability
  for (int i = 0; i < num_banks; i++) {
    field->assignRandomElt(context->alphaV->elt(i));    
    field->assignRandomElt(context->alphaW->elt(i));
    field->assignRandomElt(context->alphaY->elt(i));
    field->assignRandomElt(context->beta->elt(i));
    field->assignRandomElt(context->gamma->elt(i));
  }
  timer_rand->stop();
  
  // Figure out how many encodings we'll do
  int num_polys = 0;
  for (int i = 0; i < num_banks; i++) {
    num_polys += bank_sizes[i];
  }

  Timer* timer_prep = timers.newTimer("PreKeyPwrs", "KeyGenInit");  timer_prep->start();
  int numEncodings =  num_roots + (3*2+1) * num_polys + 5*num_banks;  // num_roots for the g^{s^i} terms, 3*2 for v,w,y and alpha versions, and 1 for beta; 5 VK terms per bank
	encoding->prepareForManyEncodings(numEncodings, Config::max_mem, Config::precomp_free);
  timer_prep->stop();
	
  // Encode the VK terms (we could do this on a per-bank streaming basis, but I don't think it would be a substantial savings
  Timer* timer_vk = timers.newTimer("PvEncodeKEA", "KeyGenInit");  timer_vk->start();
  FieldElt* betaGamma = field->newElt();
  for (int i = 0; i < num_banks; i++) {
    encoding->encode(R, context->alphaV->elt(i), context->kp->vk->alphaV->elt(i));
	  encoding->encode(L, context->alphaW->elt(i), context->kp->vk->alphaW->elt(i));
	  encoding->encode(R, context->alphaY->elt(i), context->kp->vk->alphaY->elt(i));
	  encoding->encode(R, context->gamma ->elt(i), context->kp->vk->gammaR->elt(i));

    field->mul(context->beta->elt(i), context->gamma->elt(i), betaGamma);
    encoding->encode(L, betaGamma, context->kp->vk->betaGammaL->elt(i));
	  encoding->encode(R, betaGamma, context->kp->vk->betaGammaR->elt(i));
  }
  field->delElt(betaGamma);

  // Update alphas to contain appropriate r_v, r_w, or r_y multipliers
  for (int i = 0; i < num_banks; i++) {
    field->mul(context->rV, context->alphaV->elt(i), context->alphaV->elt(i));
    field->mul(context->rW, context->alphaW->elt(i), context->alphaW->elt(i));
    field->mul(context->rY, context->alphaY->elt(i), context->alphaY->elt(i));
  }
  timer_vk->stop();

  // Encode powers of s^i
  Timer* timer_powers = timers.newTimer("PreKeyPwrsOfS", "KeyGenInit");  timer_powers->start();
	FieldElt* sPower = field->newElt();
	field->one(sPower);
	for (int i = 0; i <= num_roots; i++) {		
		encoding->encode(L, sPower, context->kp->ek->powers->elt(i));
		field->mul(sPower, context->s, sPower);
	}
  field->delElt(sPower);
  timer_powers->stop();

  Timer* timer_eval_lagrange = timers.newTimer("EvalLagrangeAtS", "KeyGenInit");  timer_eval_lagrange->start();
  // Assign roots in a geometric sequence
	FieldEltArray* targets = field->newEltArray(num_roots, true);
	FieldElt* ratio = field->newElt(); 
	field->set(ratio, context->kp->ek->root_ratio);
	if (num_roots > 0) field->one(targets->elt(0));
	for (int i = 1; i < num_roots; i++) {
		field->copy(targets->elt(i - 1), targets->elt(i));
		field->mul(targets->elt(i), ratio, targets->elt(i));
	}
  
  // Calculate L_i(s) for i in num_roots
	FieldEltArray* denominators = createLagrangeDenominators(targets, num_roots, field);
	context->lagrange = evalLagrange(num_roots, field, targets, context->s, denominators);	
  delete denominators;
  denominators = NULL;
  timer_eval_lagrange->stop();

  // Calculate g_y^{t(s)}
  Timer* timer_t_at_s = timers.newTimer("ComputeTatS", "KeyGenInit");  timer_t_at_s->start();
	FieldElt* d_at_s = field->newElt();
	evalTargetPoly(context->s, d_at_s, targets, num_roots, field);
	field->mul(context->rY, d_at_s, d_at_s);
	encoding->encode(R, d_at_s, context->kp->vk->RtAtS);
	field->delElt(d_at_s);
  timer_t_at_s->stop();

  Timer* timer_inverse_t = timers.newTimer("InvertT", "KeyGenInit");  timer_inverse_t->start();
	context->kp->ek->doPrecalcs(field, ratio, num_roots);
  field->delElt(ratio);
  timer_inverse_t->stop();

  delete targets;

  timer_init->stop();
  return context;
}

// Encode the current polynomial value into the appropriate key.  Bool is an optimization for polys that should be set to 0
void complete_current_poly(KeyGenContext* context, bool is_zero) {        
  assert(context);
  assert(context->kp);
  Field* field = context->field;  assert(field);
  EvalKey* ek  = context->kp->ek; assert(ek);
  VerifKey* vk = context->kp->vk; assert(vk);

  if (ek->bank_is_mine(context->current_bank)) {
    if (is_zero) {
      ek->apply_zero(context->current_bank, context->current_coeff_index, context->current_poly_type);
    } else {
      // EK has both poly and alpha*poly, so find the right alpha and the right r term
      FieldElt* alpha = NULL;
      FieldElt* r = NULL;
      switch (context->current_poly_type) {
        case V:
          alpha = context->alphaV->elt(context->current_bank); 
          r = context->rV;
          break;
        case W: 
          alpha = context->alphaW->elt(context->current_bank); 
          r = context->rW;
          break;
        case Y:
          alpha = context->alphaY->elt(context->current_bank); 
          r = context->rY;
          break;
      }

      FieldElt* alphaPoly = context->tmp;
      FieldElt* rPoly     = context->tmp2;
      
      field->mul(context->current_poly, alpha, alphaPoly);
      field->mul(context->current_poly, r, rPoly);  // Do this second, since alpha already has an rV/rW/rY factor in it
      ek->apply(context->current_bank, context->current_coeff_index, context->current_poly_type, rPoly, alphaPoly);

      // Fold the poly into the beta accumulator too
      if (context->current_poly_type == V) {        
        field->copy(rPoly, context->beta_accumulator);   // beta_acc <- v_k(s)   Writes it directly in
      } else if (context->current_poly_type == W) {        
        field->add(rPoly, context->beta_accumulator, context->beta_accumulator);  // beta_acc += w_k(s)
      } else {
        field->add(rPoly, context->beta_accumulator, context->beta_accumulator);  // beta_acc += r_y * y_k(s)     
      } 
    }

    if (context->current_poly_type == Y) { 
      // We reached Y, so go ahead and encode the beta term too
      field->mul(context->beta->elt(context->current_bank), context->beta_accumulator, context->beta_accumulator);
      ek->apply(context->current_bank, context->current_coeff_index, context->beta_accumulator);

      // Reset the accumulator
      field->zero(context->beta_accumulator);
    } 

    if (vk->bank_is_mine(context->current_bank)) {
      // This is a shared bank, so recycle the work we already did on EK to fill in VK
      int ek_bank_id = ek->elt(context->current_bank);
      int vk_bank_id = vk->elt(context->current_bank);
      switch (context->current_poly_type) {
        case V:
          vk->encoding->copy(L, (*ek->Vpolys[ek_bank_id])[context->current_coeff_index], (*vk->Vpolys[vk_bank_id])[context->current_coeff_index]);
          break;
        case W:
          vk->encoding->copy(L, (*ek->Wpolys[ek_bank_id])[context->current_coeff_index], (*vk->Wpolys[vk_bank_id])[context->current_coeff_index]);
          break;
        case Y:
          vk->encoding->copy(L, (*ek->Ypolys[ek_bank_id])[context->current_coeff_index], (*vk->Ypolys[vk_bank_id])[context->current_coeff_index]);
          break;
      }      
    }
  } else if (vk->bank_is_mine(context->current_bank)) {
    if (is_zero) {
      vk->apply_zero(context->current_bank, context->current_coeff_index, context->current_poly_type);
    } else {
      FieldElt* r = NULL;
      switch (context->current_poly_type) {
        case V:          
          r = context->rV;
          break;
        case W: 
          r = context->rW;
          break;
        case Y:
          r = context->rY;
          break;
      }     
      field->mul(context->current_poly, r, context->current_poly);

      vk->apply(context->current_bank, context->current_coeff_index, context->current_poly_type, context->current_poly);
    }
  }

  // Reset the poly
  field->zero(context->current_poly);
}

// Transitions through the sequence of bank_id x coeff_index x which_poly to make sure we cover each point,
// even if the caller skips over some entries
void advance_key_gen_state_machine(KeyGenContext* context, int bank_id, int coeff_index, poly_selector_t which_poly) {
  assert(context);
 
  assert(bank_id >= 0 && coeff_index >= 0);
  assert(bank_id == context->current_bank || bank_id == context->current_bank + 1); // Shouldn't skip banks
  assert(bank_id == context->current_bank || (bank_id != context->current_bank && coeff_index == 0));  // Should cover all coeffs
  if (bank_id == context->current_bank && coeff_index == context->current_coeff_index) {
    assert(which_poly >= context->current_poly_type);  // Should process polys in order
  }

	//assert(coeff_index == context->current_coeff_index || coeff_index == context->current_coeff_index + 1 || bank_id != context->current_bank); // Shouldn't skip coeffs, unless we've swapped banks
	if (! (coeff_index == context->current_coeff_index || coeff_index == context->current_coeff_index + 1 || bank_id != context->current_bank)) {
		WARN("WARNING: Skipped coefficient entirely!  Went from %d to %d\n", context->current_coeff_index, coeff_index);
	}
  
 // for (int b = context->current_bank; b < bank_id; b++) {

  if (context->current_bank != bank_id || context->current_poly_type != which_poly || context->current_coeff_index != coeff_index) {  // We've finished the polynomial we were working on
    // Encode the previous polynomial value into the appropriate key
    complete_current_poly(context, false);
    bool swapped = context->current_bank != bank_id || context->current_coeff_index != coeff_index;

    // If we swapped banks or coeffs, check for missing polynomials at the end of the previous bank
    if (swapped && context->current_poly_type != Y) {      
      swapped = true;
      // We swapped to a new bank without finishing all of the polys at the previous bank      
      if (context->current_poly_type == V) { 
        // We skipped W
        context->current_poly_type = W;
        complete_current_poly(context, true);        
      }
      // We skipped Y
      context->current_poly_type = Y;
      complete_current_poly(context, true);      
    }


		if (context->current_bank == bank_id && coeff_index - context->current_coeff_index > 1) {
			// We skipped one or more entire coefficients in the current bank
			do {
				context->current_coeff_index++;
				for (int type = V; type <= Y; type++) {
					context->current_poly_type = (poly_selector_t) type;
					complete_current_poly(context, true);
				}
			} while (context->current_coeff_index < coeff_index);

		}

    // Done with the old bank/coeff
    context->current_bank = bank_id;
    context->current_coeff_index = coeff_index;        

    // If we just arrived at the current coeffs (either by swapping banks or coeffs), check for missing polys at the beginning of this coeff
    if (swapped && which_poly != V) {
      // We swapped to a new coeff and skipped some of the initial polys (V or (V and W))
      context->current_poly_type = V;
      complete_current_poly(context, true);   

      if (which_poly == Y) {
        context->current_poly_type = W;
        complete_current_poly(context, true);   
      }

      // Done with old poly types
      context->current_poly_type = which_poly;
    } 

    // If we stayed at the current coeff, but skipped the middle coeff, fill it in
    if (!swapped && context->current_poly_type == V && which_poly == Y) {
      context->current_poly_type = W;
      complete_current_poly(context, true);   
    }
  }

  // Update all the trackers
  context->current_bank = bank_id;
  context->current_coeff_index = coeff_index;
  context->current_poly_type = which_poly;
}

// poly_index is relative to the bank_id
// We expect bank_id, coeff_index, which_poly, and root_index to be monotonic across calls
// We shouldn't skip any banks along the way, since that would indicate a redundant bank, which should be eliminated.
// Nor should we skip coeff_index values, since skipping would imply a redundant poly.
// We may skip V, W, or Y, which indicates that they should be the 0 poly.
// The constant (v_0, w_0, y_0) polys are simply included in the IO bank.  F# will ensure they always get a 1 coefficient
__declspec(dllexport) void key_gen_add_poly_elt(KeyGenContext* context, int bank_id, int coeff_index, poly_selector_t which_poly, int root_index, FieldElt* val) {
  assert(context);  
  Field* field = context->field; 
  assert(field);

  //printf("Entering key_gen_add_poly_elt\n");
  //print_field_stats();
  
  // Catch up to the most recent invocation, in case the caller skipped some combinations
  advance_key_gen_state_machine(context, bank_id, coeff_index, which_poly);

  // Fold c_i * L_i(s) into the current polynomial
  field->mul(context->lagrange->elt(root_index), val, context->tmp);
  field->add(context->tmp, context->current_poly, context->current_poly);
  //printf("Leaving key_gen_add_poly_elt\n");
  //print_field_stats();
}

// Version of key_gen_add_poly_elt that can be called with an F# int, rather than a FieldElt pointer
FieldElt* key_gen_add_val = NULL;
__declspec(dllexport) void key_gen_add_poly_elt_by_val(KeyGenContext* context, int bank_id, int coeff_index, poly_selector_t which_poly, int root_index, int val) {
  //printf("Entering key_gen_add_poly_elt_by_val\n");
  //print_field_stats();
  if (key_gen_add_val == NULL) {
    key_gen_add_val = context->field->newElt();
  }
  context->field->set(key_gen_add_val, val);
  key_gen_add_poly_elt(context, bank_id, coeff_index, which_poly, root_index, key_gen_add_val);
  //printf("Leaving key_gen_add_poly_elt_by_val\n");
  //print_field_stats();
}

__declspec(dllexport) KeyPair* key_gen_finalize(KeyGenContext* context) {
  // Need to call one more time to finish off any polynomials we may have skipped over
  advance_key_gen_state_machine(context, context->current_bank, context->current_coeff_index + 1, V); // Fake a call for one more than the max possible

  context->encoding->doneWithManyEncodings();

  KeyPair* kp = context->kp;

#ifdef TEST_SERIALIZATION
	// Test key serialization/deserialization
  save_key(context->encoding, kp, BOTH, "test.keys");
  KeyPair* kp2 = load_key(context->encoding, BOTH, "test.keys");

  assert(kp->ek->equals(context->encoding, kp2->ek));
  assert(kp->vk->equals(context->encoding, kp2->vk));
	delete kp2;
#endif //TEST_SERIALIZATION

#ifdef PERF_DEBUG 
  Encoding* encoding = context->encoding;
  int nonzero_v = 0;
  int nonzero_w = 0;
  int nonzero_y = 0;
  //int total_size = 0;
  LEncodedElt* lzero = encoding->new_elt(L);
  REncodedElt* rzero = encoding->new_elt(R);
  encoding->zero(L, lzero);
  encoding->zero(R, rzero);  
  for (int b = 0; b < kp->vk->num_banks_mine; b++) {
    nonzero_v += kp->vk->Vpolys[b]->size();
    nonzero_w += kp->vk->Wpolys[b]->size();
    nonzero_y += kp->vk->Ypolys[b]->size();
  }

  printf("Non-zero VK polynomials: V=%d, W=%d, Y=%d.\n", nonzero_v, nonzero_w, nonzero_y);
  
  nonzero_v = 0;
  nonzero_w = 0;
  nonzero_y = 0;  
  for (int b = 0; b < kp->ek->num_banks_mine; b++) {
    nonzero_v += kp->ek->Vpolys[b]->size();
    nonzero_w += kp->ek->Wpolys[b]->size();
    nonzero_y += kp->ek->Ypolys[b]->size();
  }

  printf("Non-zero EK polynomials: V=%d, W=%d, Y=%d.\n", nonzero_v, nonzero_w, nonzero_y);
#endif // PERF_TEST

  delete context;

  return kp;
}

__declspec(dllexport) void key_delete(KeyPair* keys) {
  delete keys;
}

//////////////////////////////////////////////////////////////////////////////
//				Commitment routines
//////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) void save_commits(Encoding* encoding, int num_commits, Commitment* commits[], char* file_name, bool simple = false) {
  assert(encoding);
  assert(commits);
  FileArchiver file_archiver(file_name, FileArchiver::Write);
  Archive arc(&file_archiver, encoding);

  if (!simple) {
    arc.write(num_commits);
  }

	for (int i = 0; i < num_commits; i++) {
		commits[i]->serialize(&arc, simple);
	}
}

__declspec(dllexport) void save_commit(Encoding* encoding, Commitment* commit, char* file_name, bool simple = false) {
	Commitment* commits[] = { commit };
	save_commits(encoding, 1, commits, file_name, simple);
}

__declspec(dllexport) Commitment** load_commits(Encoding* encoding, KeyPair* keys, char* file_name) {
  assert(encoding);
	assert(keys);

	FileArchiver file_archiver(file_name, FileArchiver::Read);
  Archive arc(&file_archiver, encoding);

	int num_commits;
	arc.read(num_commits);

	Commitment** commits = new Commitment*[num_commits];

	for (int i = 0; i < num_commits; i++) {
		commits[i] = new Commitment(keys);
		commits[i]->deserialize(encoding, &arc);
	}
  
  return commits;
}

__declspec(dllexport) Commitment* load_commit(Encoding* encoding, KeyPair* keys, char* file_name) {
	Commitment** commits = load_commits(encoding, keys, file_name);
	Commitment* commit = commits[0];
	delete commits;
	return commit;
}


// Assumption: When bank_id changes, we're done with the previous bank, and we're not coming back to it
// coeff_index is relative to the bank_id
__declspec(dllexport) Commitment* commit_init(KeyPair* keys, bool prove, int bank_id, bool recomputed) { 
	Commitment* commit = new Commitment(keys, prove, bank_id, recomputed);
	return commit;
}

__declspec(dllexport) void commit_nonzero_poly(Commitment* commit, poly_selector_t which_poly, int root_index, FieldElt* nonzero_val_at_root, FieldElt* coeff_val) {  
  assert(commit->prover);   // Verifier doesn't bother with polynomials
  commit->timer_nonzero->start();

  // Add it to the polynomials we're accumulating
  FieldEltArray* arr = which_poly == V ? commit->dense_V : which_poly == W ? commit->dense_W : commit->dense_Y;
  Field* field = commit->kp->ek->encoding->getSrcField();

  field->mul(coeff_val, nonzero_val_at_root, commit->tmp);    // Coeff * non-zero val at root
  field->add(commit->tmp, arr->elt(root_index), arr->elt(root_index));
  commit->timer_nonzero->stop();  
}

void mul_add(Commitment* c, elt_t t, EncodedElt* base, FieldElt* coeff, EncodedElt* acc) {
  MultiBankKey* key = c->prover ? (MultiBankKey*)c->kp->ek : (MultiBankKey*)c->kp->vk;
  EncodedElt* tmp = t == L ? c->tmpL : c->tmpR;
  
  c->kp->ek->encoding->mul(t, base, coeff, tmp);
  c->kp->ek->encoding->add(t, tmp, acc, acc);
}

// Performs commit *= g^(a[coeff_index] * coeff_val), where g^a[coeff_index] is already in the evaluation key
__declspec(dllexport) void commit_enc_poly(Commitment* commit, int coeff_index, FieldElt* coeff_val) {
  //if (!commit->prover) { commit->encoding()->collect_stats(true); }
  commit->timer_enc->start();
  // Add it to the various commitment crypto elements we're accumulating
  commit->commit_enc_poly(coeff_index, coeff_val);
  commit->timer_enc->stop();
  //if (!commit->prover) { commit->encoding()->collect_stats(false); }
}

// Same as commit_enc_poly(), but the same coeff_val is shared among many coeff_index (useful in conditional cases)
__declspec(dllexport) void commit_enc_poly_multi(Commitment* commit, const int coeff_index[], int st, int en,FieldElt* coeff_val) {
  //if (!commit->prover) { commit->encoding()->collect_stats(true); }
  commit->timer_enc->start();
	commit->commit_enc_poly_multi(coeff_index+st, en-st, coeff_val);
  commit->timer_enc->stop();
  //if (!commit->prover) { commit->encoding()->collect_stats(false); }
}

__declspec(dllexport) void commit_finalize(Commitment* commit) {  
  //if (!commit->prover) { commit->encoding()->collect_stats(true); }
  commit->timer_enc->start();
  commit->finalize_enc_polys();
  commit->timer_enc->stop();
  //if (!commit->prover) { commit->encoding()->collect_stats(false); }
}

// Call to free-up memory after this commitment is no longer needed for proving, but you still want it around
__declspec(dllexport) void commit_streamline(Commitment* commit) {
  // Clean-up intermediate values
	commit->streamline();

#ifdef TEST_SERIALIZATION
  save_commit(commit->encoding(), commit, "test.com");
  Commitment* c2 = load_commit(commit->encoding(), commit->kp, "test.com");
  assert(commit->equals(commit->encoding(), c2));
#endif // TEST_SERIALIZATION
}

__declspec(dllexport) void commit_delete(Commitment* commit) {
  delete commit;
}

__declspec(dllexport) void commit_print(Commitment* commit) {
  if (!commit) {
    WARN("You asked me to print a NULL commitment!\n");
  } else {
    commit->print();
  }
}

//////////////////////////////////////////////////////////////////////////////
//				Proof routines
//////////////////////////////////////////////////////////////////////////////

void add_dense_polys(Field* field, int len, FieldEltArray* a, FieldEltArray* b, FieldEltArray* r) {
  assert(a && b && r);
  assert(len >= 0);

#ifdef DEBUG
  FieldEltArray* debug_array = field->newEltArray(len, true);
  for (int i = 0; i < len; i++) {
    field->add(a->elt(i), b->elt(i), debug_array->elt(i));
  }
#endif // DEBUG

  if (a->complete() && b->complete()) {
    for (int i = 0; i < len; i++) {
      field->add(a->elt(i), b->elt(i), r->elt(i));
    }
  } else if (a->complete() || b->complete()) {
    FieldEltArray* complete = a->complete() ? a : b;
    FieldEltArray* sparse = a->complete() ? b : a;

    // Assuming spare arrays are truly sparse, I think the following is more efficient than the version below it
    if (complete != r) {  // I think most of the time we're computing add_dense_polys(a, r, r), so we can optimize here
      for (uint32_t i = 0; i < (uint32_t)len; i++) {
        field->copy(complete->elt(i), r->elt(i));
      }
    }

    sparse->begin(); 
    while (!sparse->end()) {
      FieldPair* sparse_pair = sparse->next();
      field->add(sparse_pair->second, r->elt(sparse_pair->first), r->elt(sparse_pair->first));
    }
  } else {
    // Try to do the copies on the array with more values
    FieldEltArray* large = a->length() >= b->length() ? a : b;
    FieldEltArray* small = a == large ? b : a;

    large->begin();
    while (!large->end()) {
      FieldPair* sparse_pair = large->next();
      field->copy(sparse_pair->second, r->elt(sparse_pair->first));
    }

    small->begin();
    while (!small->end()) {
      FieldPair* sparse_pair = small->next();
      field->add(sparse_pair->second, r->elt(sparse_pair->first), r->elt(sparse_pair->first));
    }
  }

#ifdef DEBUG
  for (int i = 0; i < len; i++) {
    assert(field->equal(debug_array->elt(i), r->elt(i)));
  }
  delete debug_array;
#endif // DEBUG
}

__declspec(dllexport) void save_proofs(Encoding* encoding, int num_proofs, Proof* proofs[], char* file_name, bool simple = false) {
	assert(encoding);
	assert(proofs);

  FileArchiver file_archiver(file_name, FileArchiver::Write);
  Archive arc(&file_archiver, encoding);

  if (!simple) {
    arc.write(num_proofs);
  }

	for (int i = 0; i < num_proofs; i++) {
		proofs[i]->serialize(&arc, simple);
	}
}

__declspec(dllexport) Proof** load_proofs(Encoding* encoding, char* file_name) {
  assert(encoding);

	FileArchiver file_archiver(file_name, FileArchiver::Read);
  Archive arc(&file_archiver, encoding);

	int num_proofs;
	arc.read(num_proofs);

	Proof** proofs = new Proof*[num_proofs];

	for (int i = 0; i < num_proofs; i++) {
		proofs[i] = new Proof(encoding);
		proofs[i]->deserialize(&arc);
	}
  
  return proofs;
}

// Cedric: added a variant that does not rely on Proof**
__declspec(dllexport) Proof* load_proof(Encoding* encoding, char* file_name) {
	Proof** proofs = load_proofs(encoding, file_name);
	Proof* proof = proofs[0];
	delete proofs;
	return proof;
}

// Expects commits[i] to be either a real commitment, or NULL to indicate this bank isn't active.  |commits| should be # of banks
__declspec(dllexport) Proof* prove(KeyPair* kp, Commitment* commits[]) {
  assert(kp);
  assert(kp && kp->ek);
  Timer* totalPoly = timers.newTimer("PolyOps", NULL);
  totalPoly->start();

  EvalKey* ek = kp->ek;
  Encoding* encoding = ek->encoding;
  Field* field = encoding->getSrcField();

  // Consolidate the polynomials and compute h(x).
  FieldEltArray* dense_V = field->newEltArray(ek->num_roots, true);
  FieldEltArray* dense_W = field->newEltArray(ek->num_roots, true);
  FieldEltArray* dense_Y = field->newEltArray(ek->num_roots, true);

  field->zero(dense_V, ek->num_roots);
  field->zero(dense_W, ek->num_roots);
  field->zero(dense_Y, ek->num_roots);

  for (int i = 0; i < ek->num_banks_total; i++) {
    if (commits[i] != NULL) {
      assert(!commits[i]->done);  // Otherwise, someone cleaned up too soon!
      add_dense_polys(field, ek->num_roots, commits[i]->dense_V, dense_V, dense_V);
      add_dense_polys(field, ek->num_roots, commits[i]->dense_W, dense_W, dense_W);
      add_dense_polys(field, ek->num_roots, commits[i]->dense_Y, dense_Y, dense_Y);
    }
  }

//#define NOGEOM
  FieldElt* ratio = field->newElt();
  FieldElt* one = field->newElt();
  field->one(one);
  field->set(ratio, ek->root_ratio);
#ifdef NOGEOM
  // TODO: Generate these on demand, rather than use up all of this memory!
  // Assign roots in a geometric progression
  FieldEltArray* roots = field->newEltArray(ek->num_roots, true);  
  if (ek->num_roots > 0) field->one(roots->elt(0));
  for (int i = 1; i < ek->num_roots; i++) {
	  field->copy(roots->elt(i - 1), roots->elt(i));
	  field->mul(roots->elt(i), ratio, roots->elt(i));
  }
#endif 
  Timer* t = timers.newTimer("Interpolation", "PolyOps");
  t->start();
#ifdef NOGEOM
  FieldEltArray* denominators = createLagrangeDenominators(roots, ek->num_roots, field);
  PolyTree* poly_tree = new PolyTree(field, roots, ek->num_roots);	
	Poly* V = Poly::interpolate(field, dense_V, ek->num_roots, poly_tree, denominators);
	Poly* W = Poly::interpolate(field, dense_W, ek->num_roots, poly_tree, denominators);
	Poly* Y = Poly::interpolate(field, dense_Y, ek->num_roots, poly_tree, denominators);
  delete roots;
#else
  PolyGeometricInterpolator pgi(field, one, ratio, ek->num_roots - 1);
  Poly* V = Poly::interpolateGeometric(pgi, dense_V);
  Poly* W = Poly::interpolateGeometric(pgi, dense_W);
  Poly* Y = Poly::interpolateGeometric(pgi, dense_Y);
#endif
  t->stop();

  // reclaim memory aggressively
  delete dense_V;
  delete dense_W;
  delete dense_Y;

#ifdef NOGEOM
  poly_tree->deleteAllButRoot();
  delete denominators;
  denominators = NULL;
#endif

	Poly Ht(field);
	Timer* tmul = timers.newTimer("Mul", "PolyOps");
	tmul->start();
	Poly::mul(*V, *W, Ht);
	tmul->stop();
	Poly::sub(Ht, *Y, Ht);

#ifdef NOGEOM
	Poly* T = poly_tree->polys[poly_tree->height-1][0];
#endif
  field->delElt(ratio);
	Poly H(field), R(field), poly_zero(field);
	Poly::zero(poly_zero, 1);
	Timer* tmod = timers.newTimer("TMod", "PolyOps");
	tmod->start();
//#define IGNORE_PRECOMPUTE
#ifdef IGNORE_PRECOMPUTE
  #ifndef NOGEOM
    #error "IGNORE_PRECOMPUTE requires NOGEOM"
	unique_ptr<Poly> T = PolyTree::polyTreeRootGeom(field, one, ratio, ek->num_roots);
  #endif
	Poly::mod(Ht, *T, H, R);
#else
	Poly::reverse(Ht);
	Poly::mul(Ht, *ek->tInv, H);
	Poly::reduce(H, ek->num_roots-1);
	Poly::reverse(H);
#endif
	tmod->stop();
	totalPoly->stop();

  if (Config::verbose_level >= INFO) {
    if (Config::machine_readable) {
      timers.printRaw();
    } else {
      timers.printStats();
    }
  }

#ifdef NOGEOM
	delete poly_tree;
	poly_tree = NULL;
#endif
	delete V;
	delete W;
	delete Y;
  field->delElt(one);
  Proof* proof = new Proof(encoding);
  int field_elt_max_bitsize = field->eltSize()*sizeof(digit_t) * 8;
  ip_handle_t precomputed_data = encoding->prepareForInnerProduct(ek->powers, ek->num_roots, 1, field_elt_max_bitsize, Config::max_mem, Config::precomp_free);
	encoding->innerProduct(precomputed_data, H.coefficients, H.getDegree() + 1, proof->H);
	encoding->doneWithInnerProduct(precomputed_data);

#ifdef TEST_SERIALIZATION
	Proof* proofs[] = { proof };
	save_proofs(encoding, 1, proofs, "test.proof");
	Proof** restored_proofs = load_proofs(encoding, "test.proof");
	assert(restored_proofs);
	assert(proof->equal(encoding, restored_proofs[0]));
#endif // TEST_SERIALIZATION
  //proof->print();
  return proof;
}

__declspec(dllexport) void proof_delete(Proof* proof) {
  delete proof;
}

__declspec(dllexport) void proof_print(Proof* proof) {
  assert(proof);
  proof->print();
}


//////////////////////////////////////////////////////////////////////////////
//				Simple file interface for F#
//////////////////////////////////////////////////////////////////////////////

typedef struct {
  FileArchiver* file_archiver;
  Archive*      arc;
} FileHandle;

__declspec(dllexport) FileHandle* file_open(Encoding* encoding, char* file_name, bool read) {
  FileHandle* fh = new FileHandle;
  fh->file_archiver = new FileArchiver(file_name, read ? FileArchiver::Read : FileArchiver::Write);
  fh->arc = new Archive(fh->file_archiver, encoding);
  return fh;
}

__declspec(dllexport) int file_read_int(Encoding* encoding, FileHandle* fh) {
  int ret = 0;
  fh->arc->read(ret);
  return ret;
}

__declspec(dllexport) int file_read_unsigned_char(Encoding* encoding, FileHandle* fh) {
  unsigned char ret = 0;
  fh->arc->read(ret);
  return ret;
}

__declspec(dllexport) FieldElt* file_read_field_elt(Encoding* encoding, FileHandle* fh) {
  Field* field = encoding->getSrcField();
  FieldElt* elt = field->newElt();
  field->read(fh->file_archiver, elt);
  return elt;
}

__declspec(dllexport) void file_write_int(Encoding* encoding, FileHandle* fh, int val) {
  fh->arc->write(val);
}

__declspec(dllexport) void file_write_field_elt(Encoding* encoding, FileHandle* fh, FieldElt* elt) {
  Field* field = encoding->getSrcField();
  field->write(fh->file_archiver, elt);
}


__declspec(dllexport) void file_close(Encoding* encoding, FileHandle* fh) {
  delete fh->arc;
  delete fh->file_archiver;
  delete fh;
}


//////////////////////////////////////////////////////////////////////////////
//				Verification routines
//////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) int verify_commit(KeyPair* keys, int bank_id, Commitment* c, bool expect_to_fail) {
  assert(keys);
  assert(keys->vk);
  assert(c);

  VerifKey* vk = keys->vk;
  Encoding* encoding = vk->encoding;
	REncodedElt* encoded_one_R = encoding->new_elt(R);
	encoding->one(R, encoded_one_R);

  c->timer_alpha_check->start();
  bool alpha_check = true;
  alpha_check &= encoding->mulEquals(c->c_vals[V], vk->alphaV->elt(bank_id), c->c_vals[alphaV], encoded_one_R);
  assert(alpha_check || expect_to_fail);
  alpha_check &= encoding->mulEquals(vk->alphaW->elt(bank_id), c->c_vals[W], c->c_vals[alphaW], encoded_one_R);
  assert(alpha_check || expect_to_fail);
  alpha_check &= encoding->mulEquals(c->c_vals[Y], vk->alphaY->elt(bank_id), c->c_vals[alphaY], encoded_one_R);
  assert(alpha_check || expect_to_fail);
		
	if (alpha_check) {
		INFO("\tPassed alpha check!\n");
	} else {
		INFO("\tFailed alpha check!\n");
	}
  c->timer_alpha_check->stop();

  c->timer_beta_check->start();
  bool beta_check = true;
	EncodedProduct* gammaZ = encoding->new_prod();
  encoding->mul(c->c_vals[beta], vk->gammaR->elt(bank_id), gammaZ);
	
	LEncodedElt* encodedVY = encoding->new_elt(L);
  encoding->add(L, c->c_vals[V], c->c_vals[Y], encodedVY);
	EncodedProduct* betaGammaVY = encoding->new_prod();
	encoding->mul(encodedVY, vk->betaGammaR->elt(bank_id), betaGammaVY);		// = Enc(beta*gamma*(V(s) + Y(s))

	// Compute (beta*gamma)(W)
	EncodedProduct* betaGammaW = encoding->new_prod();
  encoding->mul(vk->betaGammaL->elt(bank_id), c->c_vals[W], betaGammaW);   // = Enc(beta*gamma*(W(s))

	// Combine to get (beta*gamma)*(V + Y + W)
	EncodedProduct* sum = encoding->new_prod();
	encoding->add(betaGammaVY, betaGammaW, sum);

	// Compare
	beta_check &= encoding->equals(gammaZ, sum);
  assert(beta_check || expect_to_fail);

	if (beta_check) {
    INFO("\tPassed beta check!\n");
	} else {
    INFO("\tFailed beta check!\n");
	}
  c->timer_beta_check->stop();

  delete gammaZ;
  encoding->del_elt(L, encodedVY);
	encoding->del_elt(R, encoded_one_R);
  delete betaGammaVY;
  delete betaGammaW;
  delete sum;
  return alpha_check && beta_check;
}

__declspec(dllexport) int verify_proof(KeyPair* keys, Proof* proof, int num_banks, Commitment** c /*[num_banks]*/, bool expect_to_fail) {
  assert(keys);
  assert(keys->vk);
  assert(proof);

  Timer* timer_verify_proof = timers.newTimer("VerifyProof", NULL); timer_verify_proof->start();
  Timer* timer_consolidation = timers.newTimer("ConsolidateCommits", "VerifyProof");  timer_consolidation->start();

  VerifKey* vk = keys->vk;
  Encoding* encoding = vk->encoding;
	bool quadratic_check = true;

  // Consolidate all of the active banks
  EncodedElt* encodedV = encoding->new_elt(L);
  EncodedElt* encodedW = encoding->new_elt(R);
  EncodedElt* encodedY = encoding->new_elt(L);

  encoding->zero(L, encodedV);
  encoding->zero(R, encodedW);
  encoding->zero(L, encodedY);

  for (int i = 0; i < num_banks; i++) {
    if (c[i] != NULL) {
      encoding->add(L, c[i]->c_vals[V], encodedV, encodedV);
      encoding->add(R, c[i]->c_vals[W], encodedW, encodedW);
      encoding->add(L, c[i]->c_vals[Y], encodedY, encodedY);
    }
  }

  timer_consolidation->stop();
  
  Timer* timer_quad = timers.newTimer("QuadCheck", "VerifyProof");  
  timer_quad->start(); 
	REncodedElt* encoded_one_R = vk->encoding->new_elt(R);
	vk->encoding->one(R, encoded_one_R);
  quadratic_check = quadratic_check && encoding->mulSubEquals(encodedV, encodedW, encodedY, encoded_one_R, proof->H, vk->RtAtS);
	quadratic_check = quadratic_check && encoding->mulSubEquals(encodedV, encodedW, encodedY, encoded_one_R, proof->H, vk->RtAtS);
	vk->encoding->del_elt(R, encoded_one_R);

	if (quadratic_check) { 
    INFO("\tPassed quadratic check!\n");
	} else {
    INFO("\tFailed quadratic check!\n");
	}
	assert(quadratic_check || expect_to_fail);
  timer_quad->stop();

  encoding->del_elt(L, encodedV);
  encoding->del_elt(R, encodedW);
  encoding->del_elt(L, encodedY);  

  timer_verify_proof->stop();
  return quadratic_check;
}


}	// end extern "C"

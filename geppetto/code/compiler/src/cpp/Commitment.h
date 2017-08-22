#pragma once

#include "Keys.h"

class Commitment {
public:
  static const int MULTI_EXP_BATCH_SIZE = 1000;  

	Commitment(KeyPair* kp_i);
  Commitment(KeyPair* kp_i, bool prover_i, int bank_id, bool recomputed);
  ~Commitment();
  void streamline();    // Free up some of our memory, assuming this commitment will never be used in a call to prove

  elt_t Commitment::encFor(int id) const;
  void commit_enc_poly(int coeff_index, FieldElt* coeff_val);
  void commit_enc_poly_multi(const int coeff_index[], int coeff_count, FieldElt* coeff_val);
  void finalize_enc_polys();

  // Does this commitment care about this type of polynomial?
  bool is_mine(poly_selector_t pid);

	MultiBankKey* key();	// Returns the active key
	Encoding* encoding();	// Returns a pointer to the encoding in use

	virtual void serialize(Archive* arc, bool simple=false);
	virtual void deserialize(Encoding* encoding, Archive* arc);
  virtual void print();

	virtual bool equals(Encoding* encoding, Commitment* other);

	// We only want to compute these once, but we can reuse them across all of our proofs
	static PolyTree* polyTree;					// Precompute polynomial tree, with leaves (x-r), for r in target roots, and internal nodes as products of their children
	static FieldEltArray* denominators;	// Pre-computed Lagrange denominators

  KeyPair* kp;
  Encoding* my_encoding;
  int bank_id;
  bool prover;      // Does the prover generate this commitment?
  bool recomputed;  // Will this commitment be recomputed by the verifier?
  bool done;
  FieldEltArray* dense_V;
  FieldEltArray* dense_W;
  FieldEltArray* dense_Y;
  FieldElt* tmp;

  // Committed values
  EncodedElt* c_vals[NUM_POLY_SELECTORS];
 
  // Accumulators for batching multi-exponentiations
  EncodedEltArray* batch_bases[NUM_POLY_SELECTORS];
    
  // Exponents for batching multi-exponentiations
  FieldEltArray* batch_exp[NUM_POLY_SELECTORS];

  // Counters for batching multi-exponentiations
  int batch_ctr[NUM_POLY_SELECTORS];

  // Track statistics about our batching
  int batch_counts[NUM_POLY_SELECTORS];

  int batch_size;   // Batch size used for this particular commitment

  EncodedElt* tmpL;
  EncodedElt* tmpR;

  Timer* timer_nonzero;
  Timer* timer_enc;
  Timer* timer_alpha_check;
  Timer* timer_beta_check;

private:
  void process_batch(int id);
  void init_timers();
};

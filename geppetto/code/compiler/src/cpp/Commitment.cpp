#include "Commitment.h"
#include "Config.h"

Commitment::Commitment(KeyPair* kp_i) {
  assert(kp_i);

  kp = kp_i;
  bank_id = -1;
  prover = false;
  recomputed = false;
  done = false;

  dense_V = NULL;
  dense_W = NULL;
  dense_Y = NULL;
 
  tmp = NULL;
  
  for (int p = 0; p < NUM_POLY_SELECTORS; p++) {
    batch_ctr[p] = 0;
    batch_exp[p] = NULL;
  }

  tmpL = NULL;
  tmpR = NULL;

  // Initialize batch bases, so the destructor doesn't explode on deserialized commitments
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    batch_bases[id] = NULL;
    batch_counts[id] = 0;
  }
  init_timers();
}

Commitment::Commitment(KeyPair* kp_i, bool prover_i, int bank_id_i, bool recomputed) {
  assert(kp_i);
  assert(kp_i->ek || kp_i->vk);

  this->kp = kp_i;
  this->bank_id = bank_id_i;
  this->prover = prover_i;
  this->recomputed = recomputed; 
  this->done = false;

  MultiBankKey* key = prover_i ? (MultiBankKey*)(kp_i->ek) : (MultiBankKey*)(kp_i->vk);
  this->my_encoding = key->encoding;

  Field* field = my_encoding->getSrcField();

  if (prover) {    
    int bank_size = key->bank_size(bank_id);
    if (bank_size < 0.4 * key->num_roots && bank_size != 1) {
      //printf("Using map instead of array for the dense polynomials because size = %d and roots = %d\n", bank_size, kp->ek->num_roots);
      dense_V = new FieldEltDigitsMap(field, field->num_digits_per_elt(), true);
      dense_W = new FieldEltDigitsMap(field, field->num_digits_per_elt(), true);
      dense_Y = new FieldEltDigitsMap(field, field->num_digits_per_elt(), true);
    } else {
      dense_V = field->newEltArray(key->num_roots, true);
      dense_W = field->newEltArray(key->num_roots, true);
      dense_Y = field->newEltArray(key->num_roots, true);
    }

    field->zero(dense_V, key->num_roots);
    field->zero(dense_W, key->num_roots);
    field->zero(dense_Y, key->num_roots);

    tmp = field->newElt();
  }
  
  tmpL = my_encoding->new_elt(L);
  tmpR = my_encoding->new_elt(R);

  this->batch_size = min(MULTI_EXP_BATCH_SIZE, key->bank_size(bank_id));

  // Initialize all of our elements
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    batch_ctr[id] = 0;
    batch_counts[id] = 0;

    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { c_vals[id] = NULL;  continue; }

    elt_t type = pid == W ? R : L;

    c_vals[id] = my_encoding->new_elt(type);
    my_encoding->zero(type, c_vals[id]);
    batch_bases[id] = my_encoding->new_elt_array(type, batch_size, true);
    batch_exp[id] = field->newEltArray(batch_size, true);
  }

  init_timers();
}

Commitment::~Commitment() {
  if (!done) {
    streamline();
  }
 
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { continue; }

    my_encoding->del_elt(pid == W ? R : L, c_vals[id]);
    c_vals[id] = NULL;

    delete batch_bases[id];   
    delete batch_exp[id];
    batch_bases[id] = NULL;  
    batch_exp[id] = NULL;
  }
}

void Commitment::init_timers() {
  // TODO: Cleanup the leaked names.  The problem is that they need to live as long as the timer does.
#define MAX_TIMER_NAME_SIZE 100 
  char* name = (char*)malloc(MAX_TIMER_NAME_SIZE);
  sprintf_s(name, MAX_TIMER_NAME_SIZE, "commit_enc_%s_%p", prover ? "prover" : "verifier", this);
  timer_enc = timers.newTimer(name, NULL);
  name = (char*)malloc(100);
  sprintf_s(name, MAX_TIMER_NAME_SIZE, "commit_nonzero_%s_%p", prover ? "prover" : "verifier", this);
  timer_nonzero = timers.newTimer(name, NULL);
  name = (char*)malloc(100);
  sprintf_s(name, MAX_TIMER_NAME_SIZE, "commit_alpha_%s_%p", prover ? "prover" : "verifier", this);
  timer_alpha_check = timers.newTimer(name, NULL);
  name = (char*)malloc(100);
  sprintf_s(name, MAX_TIMER_NAME_SIZE, "commit_beta_%s_%p", prover ? "prover" : "verifier", this);
  timer_beta_check = timers.newTimer(name, NULL);
}

void Commitment::streamline() {
  if (prover) {
    delete dense_V;
    delete dense_W;
    delete dense_Y;
    dense_V = NULL;
    dense_W = NULL;
    dense_Y = NULL; 
   
    my_encoding->getSrcField()->delElt(tmp);
    tmp = NULL;
  }

  my_encoding->del_elt(L, tmpL);
  my_encoding->del_elt(R, tmpR);

  tmpL = NULL;
  tmpR = NULL;

  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { continue; }

    delete batch_bases[id];
    delete batch_exp[id];
    batch_bases[id] = NULL;
    batch_exp[id] = NULL;
  }

	done = true;
}

elt_t Commitment::encFor(int id) const
{
	return id == W ? R : L;
}

void Commitment::process_batch(int id) {
  poly_selector_t pid = (poly_selector_t)id;
  EncodedElt* tmp = pid == W ? tmpR : tmpL;

  Encoding* encoding = key()->encoding;

  int batch_len = batch_ctr[id];
  
  int field_elt_max_bit_set = encoding->getSrcField()->highest_bit(batch_exp[id]->elt(batch_len - 1));

  ip_handle_t h = encoding->prepareForInnerProduct(batch_bases[id], batch_len, 1, field_elt_max_bit_set+1, Config::max_mem, Config::precomp_free);
  encoding->innerProduct(h, batch_exp[id], batch_len, tmp);
  encoding->doneWithInnerProduct(h);

  encoding->add(pid == W ? R : L, c_vals[id], tmp, c_vals[id]);

  batch_counts[id] += batch_len;
}

void Commitment::finalize_enc_polys() {
  Encoding* encoding = key()->encoding;
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { continue; }

    if (batch_ctr[id] == 0) {
      // Nothing to do, since we initialize the commitment values to 0
    } else {
      process_batch(id);
    }   
  }  
}

static int enc_count_total = 0;
static int enc_count_nonzero = 0;
extern "C" {
	__declspec(dllexport) void enc_count_reset() { enc_count_total = enc_count_nonzero = 0; }
	__declspec(dllexport) int enc_count_total_get() { return enc_count_total; }
	__declspec(dllexport) int enc_count_nonzero_get() { return enc_count_nonzero; }
}

#define COPY_ADD(arr) \
  EncodedElt* key_elt = arr->find(coeff_index[j]);                           \
  if (key_elt != NULL) { /* Key isn't 0 here */                              \
    /* Save the exponent */                                                  \
    field->copy(coeff_val, batch_exp[id]->elt(batch_ctr[id]));               \
    /* Save/accumulate the base */                                           \
    if (coeff_count == 1) { enc->copy(t, key_elt, base); }                   \
    else { enc->add(t, key_elt, base, base); }                               \
    nonzero = true;                                                          \
  }

void Commitment::commit_enc_poly_multi(const int coeff_index[], int coeff_count, FieldElt* coeff_val) {
	// Take note of the exp and base provided
	Encoding* enc = key()->encoding;
	Field* field = enc->getSrcField();
	enc_count_total++;
	if (field->isZero(coeff_val)) return;
	enc_count_nonzero++;                                                                                         

	for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
		poly_selector_t pid = (poly_selector_t)id;
		if (!is_mine(pid)) { continue; }

		int bank_index = key()->elt(bank_id);
		elt_t t = encFor(id);
		EncodedElt* base = batch_bases[id]->elt(batch_ctr[id]);
		enc->zero(t, base);
    bool nonzero = false;

		for (int j = 0; j < coeff_count; ++j) {
			switch (pid) {
        case V:      { COPY_ADD(key()->Vpolys[bank_index]);                   break; }
				case W:      { COPY_ADD(key()->Wpolys[bank_index]);                   break; }
				case Y:      { COPY_ADD(key()->Ypolys[bank_index]);                   break; }
				case alphaV: { COPY_ADD(((EvalKey*)key())->alphaVpolys[bank_index]);  break; }
				case alphaW: { COPY_ADD(((EvalKey*)key())->alphaWpolys[bank_index]);  break; }
				case alphaY: { COPY_ADD(((EvalKey*)key())->alphaYpolys[bank_index]);  break; }
				case beta:   { COPY_ADD(((EvalKey*)key())->betaVWYpolys[bank_index]); break; }
				default: assert(false); // Should never reach here!
			}
		}
    
    if (nonzero) {
      batch_ctr[id]++;
    }

    if (batch_ctr[id] == batch_size) {
      process_batch(id);
      batch_ctr[id] = 0;
    }
	}
}

void Commitment::commit_enc_poly(int coeff_index, FieldElt* coeff_val) {
  commit_enc_poly_multi(&coeff_index, 1, coeff_val);
}

bool Commitment::is_mine(poly_selector_t pid) {
  return prover || (!prover && (pid == V || pid == W || pid == Y));
}

MultiBankKey* Commitment::key() {
	return prover ? (MultiBankKey*)kp->ek : (MultiBankKey*)kp->vk;
}

Encoding* Commitment::encoding() {
  return my_encoding;
}

void Commitment::serialize(Archive* arc, bool simple) {
  if (!done) { WARN("Warning: You told me to serialize a commitment before telling me it's done.  Are you sure you're ready to serialize?\n"); }
  if (!simple) {
    arc->write(bank_id);
    arc->write(prover);
  }	

  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { continue; }

    my_encoding->write(arc->get_archiver(), pid == W ? R : L, c_vals[id], simple);    
  }	
}

void Commitment::deserialize(Encoding* encoding, Archive* arc) {
	assert(kp);
  my_encoding = encoding;
	arc->read(bank_id);
	arc->read(prover);
	done = true;

	dense_V = NULL;
	dense_W = NULL;
	dense_Y = NULL;
	tmp = NULL;

  // Read all of our elements
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { continue; }

    elt_t type = pid == W ? R : L;
    c_vals[id] = my_encoding->new_elt(type);
    encoding->read(arc->get_archiver(), type, c_vals[id]);
  }  

	tmpL = NULL;
	tmpR = NULL;
}

void Commitment::print() {
  Encoding* encoding = this->encoding();

  printf("Commitment for %s bank: %d", prover ? "prover" : "verifier", bank_id);
  printf("\n  V: ");       if (!is_mine(V))      { printf("--"); } else { printf("\n");  encoding->print(L, c_vals[V]);      }
  printf("\n  W: ");       if (!is_mine(W))      { printf("--"); } else { printf("\n");  encoding->print(R, c_vals[W]);      }
  printf("\n  Y: ");       if (!is_mine(Y))      { printf("--"); } else { printf("\n");  encoding->print(L, c_vals[Y]);      }
  printf("\n  alphaV: ");  if (!is_mine(alphaV)) { printf("--"); } else { printf("\n");  encoding->print(L, c_vals[alphaV]); }
  printf("\n  alphaW: ");  if (!is_mine(alphaW)) { printf("--"); } else { printf("\n");  encoding->print(L, c_vals[alphaW]); }
  printf("\n  alphaY: ");  if (!is_mine(alphaY)) { printf("--"); } else { printf("\n");  encoding->print(L, c_vals[alphaY]); }
  printf("\n  beta: ");    if (!is_mine(beta))   { printf("--"); } else { printf("\n");  encoding->print(L, c_vals[beta]);   }
  printf("\n");
  printf("Batch counts:\n");
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    printf("\tcount[%d] = %d\n", id, batch_counts[id]);
  }
}

bool Commitment::equals(Encoding* encoding, Commitment* other) {
	bool equal = true;

	if (bank_id != other->bank_id) {
		INFO("bank_id differ %d %d\n", bank_id, other->bank_id);
		equal = false;
	}

	if (prover != other->prover) {
    INFO("prover differ %d %d\n", prover, other->prover);
		equal = false;
	}

  // Check all of our elements
  for (int id = 0; id < NUM_POLY_SELECTORS; id++) {
    poly_selector_t pid = (poly_selector_t)id;
    if (!is_mine(pid)) { continue; }

    if (!encoding->equal(pid == W ? R : L, c_vals[id], other->c_vals[id])) {
      INFO("%d committed value differs\n", id);
      equal = false;
    }
  }

	return equal;
}

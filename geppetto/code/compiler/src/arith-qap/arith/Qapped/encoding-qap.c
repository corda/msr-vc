////////////////////////////////////////////////////////
//  Code for loading and manipulating QAP-level
//  structures (e.g., proofs, commitments, and keys)
////////////////////////////////////////////////////////


////////////////////////////////////////////////////////
//  The three QAP-related structs
////////////////////////////////////////////////////////
typedef struct {
  struct s_EncodedEltL V[SUM_BANK_SIZES];
  struct s_EncodedEltR W[SUM_BANK_SIZES];
  struct s_EncodedEltL Y[SUM_BANK_SIZES];
  struct s_EncodedEltR alphaV[NUM_BANKS_TOTAL];
  struct s_EncodedEltL alphaW[NUM_BANKS_TOTAL];
  struct s_EncodedEltR alphaY[NUM_BANKS_TOTAL];
  struct s_EncodedEltR gammaR[NUM_BANKS_TOTAL];
  struct s_EncodedEltL betaGammaL[NUM_BANKS_TOTAL];
  struct s_EncodedEltR betaGammaR[NUM_BANKS_TOTAL];  
  struct s_EncodedEltR D;
} cVerifKey;

typedef struct {
  struct s_EncodedEltL V;
  struct s_EncodedEltR W;
  struct s_EncodedEltL Y;
  struct s_EncodedEltL alphaV;
  struct s_EncodedEltL alphaW;
  struct s_EncodedEltL alphaY;
  struct s_EncodedEltL beta;
} cCommitment;

typedef struct {
  struct s_EncodedEltL H;
} cProof;


void init_commit(cCommitment* commit) {
  _encoding_initL(&commit->V);
  _encoding_initR(&commit->W);
  _encoding_initL(&commit->Y);
  _encoding_initL(&commit->alphaV);
  _encoding_initL(&commit->alphaW);
  _encoding_initL(&commit->alphaY);
  _encoding_initL(&commit->beta);
}

////////////////////////////////////////////////////////
//  Loading structs into memory
////////////////////////////////////////////////////////

// the filename is e.g. "cProof.n.data"
#define RUN_TIME 0
#define COMPILE_TIME 1

void load_int32s     (const char* name, int t, int bank_len, int* dst, bool runtime);
void load_cProof     (const char* name, cProof* pi, int n, bool compile_time);
void load_cCommitment(const char* name, cCommitment* c, int, bool compile_time);
void load_cVerifKey  (const char* name, cVerifKey* c, int n, bool compile_time);
//void zero_commit(cCommitment *commit);

// ---------------------------------------------------------

void print_commitment(cCommitment* commit) {
	enc_print_rawL("\ncommit->V",      &commit->V);
	enc_print_rawR("\ncommit->W",      &commit->W);
	enc_print_rawL("\ncommit->Y",      &commit->Y);
	enc_print_rawL("\ncommit->alphaV", &commit->alphaV);
	enc_print_rawL("\ncommit->alphaW", &commit->alphaW);
	enc_print_rawL("\ncommit->alphaY", &commit->alphaY);
	enc_print_rawL("\ncommit->beta",   &commit->beta);
}

#ifndef MQAP 
#include <stdio.h>
#include "encoding-io.c"

void load_int32s(const char* name, int t, int bank_len, int* dst, bool runtime)
{
	char filename[100];
	sprintf(filename, "../input/build/%s.%d.data", name, t);
	FILE* file = fopen(filename, "rb");
  assert(file);
	buffered_read(dst, bank_len * 4, file);	
	fclose(file);
}

void load_cProof(const char* name, cProof* proof, int id, bool compile_time) {
  char filename[100];
  sprintf(filename, "../input/build/proof-%s.%d.data", name, id);
  FILE* file = fopen(filename, "rb");
  assert(file);
  _encoding_initL(&proof->H);
  load_encodedL(file, &proof->H);
  fclose(file);
}

void load_cCommitment(const char* name, cCommitment* commit, int id, bool compile_time) {
  char filename[100];
  sprintf(filename, "../input/build/commit-%s.%d.data", name, id);
  printf("loading commitment from %s\n", filename);
  FILE* file = fopen(filename, "rb");
  assert(file);
  load_encodedL(file, &commit->V);
  load_encodedR(file, &commit->W);
  load_encodedL(file, &commit->Y);
  load_encodedL(file, &commit->alphaV);
  load_encodedL(file, &commit->alphaW);
  load_encodedL(file, &commit->alphaY);
  load_encodedL(file, &commit->beta);
  fclose(file);
}

void load_cVerifKey(const char* name, cVerifKey* vk, int id, bool compile_time) {
	char filename[100];
	sprintf(filename, "../input/build/%s.%d.data", name, id);
	printf("reading %s\n", filename);
	FILE* file = fopen(filename, "rb");
	if (file == NULL) {
		fputs("File error", stderr);
		assert(0);
	}
	if (ferror(file)) {
		perror("Open error");
	}

	for (int i = 0; i < SUM_BANK_SIZES; i++) {
		_encoding_initL(&vk->V[i]);
		load_encodedL(file, &vk->V[i]);
		//enc_print_rawL("&vk->V[i]", &vk->V[i]);
	}
	for (int i = 0; i < SUM_BANK_SIZES; i++) {
		_encoding_initR(&vk->W[i]);
		load_encodedR(file, &vk->W[i]);
		//enc_print_rawR("&vk->W[i]", &vk->W[i]);
	}
	for (int i = 0; i < SUM_BANK_SIZES; i++) {
		_encoding_initL(&vk->Y[i]);
		load_encodedL(file, &vk->Y[i]);
		//enc_print_rawL("&vk->Y[i]", &vk->Y[i]);
	}

	for (int i = 0; i < NUM_BANKS_TOTAL; i++) {
		_encoding_initR(&vk->alphaV[i]);
		load_encodedR(file, &vk->alphaV[i]);
		//enc_print_rawR("&vk->alphaV[i]", &vk->alphaV[i]);
	}
	for (int i = 0; i < NUM_BANKS_TOTAL; i++) {
		_encoding_initL(&vk->alphaW[i]);
		load_encodedL(file, &vk->alphaW[i]);
		//enc_print_rawL("&vk->alphaW[i]", &vk->alphaW[i]);
	}
	for (int i = 0; i < NUM_BANKS_TOTAL; i++) {
		_encoding_initR(&vk->alphaY[i]);
		load_encodedR(file, &vk->alphaY[i]);
		//enc_print_rawR("&vk->alphaY[i]", &vk->alphaY[i]);
	}
	for (int i = 0; i < NUM_BANKS_TOTAL; i++) {
		_encoding_initR(&vk->gammaR[i]);
		load_encodedR(file, &vk->gammaR[i]);
		//enc_print_rawR("&vk->gammaR[i]", &vk->gammaR[i]);
	}
	for (int i = 0; i < NUM_BANKS_TOTAL; i++) {
		_encoding_initL(&vk->betaGammaL[i]);
		load_encodedL(file, &vk->betaGammaL[i]);
		//enc_print_rawL("&vk->betaGammaL[i]", &vk->betaGammaL[i]);
	}
	for (int i = 0; i < NUM_BANKS_TOTAL; i++) {
		_encoding_initR(&vk->betaGammaR[i]);
		load_encodedR(file, &vk->betaGammaR[i]);
		//enc_print_rawR("&vk->betaGammaR[i]", &vk->betaGammaR[i]);
	}
	_encoding_initR(&vk->D);
	load_encodedR(file, &vk->D);
	//enc_print_rawR("&vk->D", &vk->D);
	fclose(file);
}
#endif // !MQAP


void zero_commit(cCommitment *commit) { 
  set_zeroL(&commit->V);
  set_zeroR(&commit->W);
  set_zeroL(&commit->Y);
  // Following are unnecessary, since we don't anticipate computing on them.  Including to be on the safe side, since they should be free.
  set_zeroL(&commit->alphaV);
  set_zeroL(&commit->alphaW);
  set_zeroL(&commit->alphaY);
  set_zeroL(&commit->beta);
}

// Computes a fresh commitment (which only requires V, W, and Y). 
// WARNING: The bank_id is indexed only relative to the recomputable banks, not the total number of banks
void compute_commit(cVerifKey *vk, cCommitment *commit, int global_id, int bank_id, int *inputs, int offset, int length) {
	zero_commit(commit);
	struct s_EncodedEltL tmpL;
	_encoding_initL(&tmpL);
	struct s_EncodedEltR tmpR;
	_encoding_initR(&tmpR);

	QapFp_t f;
	QapFp_init(f);
	// printf("Computing commitment for bank %d\n", bank_id);
	for (int i = 0; i < length; i++) {


		QapFp_set_ui(f, inputs[i]);
		// printf("encoding input #%d ", i);  QapFp_print(f);  printf("\n");
		// enc_print_rawL("V", &commit->V);

		// enc_printL("vk->V", &vk->V[offset + i]);
		_encoding_mulL(&vk->V[offset + i], f, &tmpL);
		// enc_print_rawL("tmpL", &tmpL);
		_encoding_add_rawL(&tmpL, &commit->V, &commit->V);

		// enc_printR("vk->W", &vk->W[offset + i]);
		_encoding_mulR(&vk->W[offset + i], f, &tmpR);
		_encoding_add_rawR(&tmpR, &commit->W, &commit->W);

		// enc_printL("vk->Y", &vk->Y[offset + i]);
		_encoding_mulL(&vk->Y[offset + i], f, &tmpL);
		_encoding_add_rawL(&tmpL, &commit->Y, &commit->Y);
	}
	printf("Computed commitment %2d [%4d] %9d *:\n", global_id, length, nRoot());
	//enc_print_rawL("V", &commit->V);
	//enc_print_rawR("W", &commit->W);
	//enc_print_rawL("Y", &commit->Y);
}

// Verifies variables on a shared bus.  We only encode these on the Y polynomials, so verification is cheaper
bool verify_commit_bus(cVerifKey *vk, cCommitment *commit, int bank_id) {
  bool result = true;
  bool ok;

  //printf("\nlightly verifying commitment to bank %d:\n", bank_id);
  //print_commitment(commit);

  struct s_EncodedEltR encoded_one_R;
  _encoding_initR(&encoded_one_R);
  copyR(&e_genR, &encoded_one_R);

  EncodedProduct lhs, rhs;
  encoded_product_init(&lhs);
  encoded_product_init(&rhs);

  // Alpha check

  // printf("About to verify the three following Y values for bank %d:\n", bank_id);
  // printf("c_vals[Y]\n");
  // enc_print_rawL("Y", &commit->Y);
  // printf("vk->alphaY\n");
  // enc_print_rawR("Y", &vk->alphaY[bank_id]);
  // printf("c_vals[alphaY]\n");
  // enc_print_rawL("Y", &commit->alphaY);
  encoding_pair(&commit->Y, &vk->alphaY[bank_id], &lhs);
  encoding_pair(&commit->alphaY, &encoded_one_R, &rhs);
  ok = encoding_product_equals(&lhs, &rhs);
  printf("Bus alpha check: %d %12d*\n", ok, nRoot());  // 20 K*
  
  result &= ok;

  // Beta check
  encoding_pair(&commit->beta, &vk->gammaR[bank_id], &lhs);
  encoding_pair(&commit->Y, &vk->betaGammaR[bank_id], &rhs);
  ok = encoding_product_equals(&lhs, &rhs);
  printf("Bus beta check:  %d %12d*\n", ok, nRoot());  // 20 K*
  result &= ok;

  printf("Verified bus commitment %d, result %d\n", bank_id, result);
  return result;
}

// Verifies a full commitment (e.g., for the local variables)
bool verify_commit(cVerifKey* vk, cCommitment* commit, int bank_id) {
  bool result = true;
  bool ok;

  //printf("\nfully verifying commitment to bank %d:\n", bank_id);
  //print_commitment(commit);

  struct s_EncodedEltR encoded_one_R;
  _encoding_initR(&encoded_one_R);
  copyR(&e_genR, &encoded_one_R);

  EncodedProduct lhs, rhs;
  encoded_product_init(&lhs);
  encoded_product_init(&rhs);

  // Alpha checks  
  encoding_pair(&commit->V, &vk->alphaV[bank_id], &lhs);  
  encoding_pair(&commit->alphaV, &encoded_one_R, &rhs);
  ok = encoding_product_equals(&lhs, &rhs);
  printf("Full alpha check: %d %12d*\n", ok, nRoot());  // 20 K*
  result &= ok;

  encoding_pair(&vk->alphaW[bank_id], &commit->W, &lhs);
  encoding_pair(&commit->alphaW, &encoded_one_R, &rhs);
  ok = encoding_product_equals(&lhs, &rhs);
  printf("Full alpha check: %d %12d*\n", ok, nRoot());  // 20 K*
  result &= ok;

  encoding_pair(&commit->Y, &vk->alphaY[bank_id], &lhs);
  encoding_pair(&commit->alphaY, &encoded_one_R, &rhs);
  ok = encoding_product_equals(&lhs, &rhs);
  printf("Full alpha check: %d %12d*\n", ok, nRoot());  // 20 K*
  result &= ok;

  // Beta checks
  struct s_EncodedEltL encodedVY;
  _encoding_initL(&encodedVY);
  EncodedProduct betaGammaVY;
  encoded_product_init(&betaGammaVY);
  EncodedProduct betaGammaW;
  encoded_product_init(&betaGammaW);

  encoding_pair(&commit->beta, &vk->gammaR[bank_id], &lhs);
  _encoding_add_rawL(&commit->V, &commit->Y, &encodedVY);
  encoding_pair(&encodedVY, &vk->betaGammaR[bank_id], &betaGammaVY);  // = Enc(beta*gamma*(V(s) + Y(s))
  encoding_pair(&vk->betaGammaL[bank_id], &commit->W, &betaGammaW);    // = Enc(beta*gamma*(W(s))
  encoding_add(&betaGammaVY, &betaGammaW, &rhs);
  ok = encoding_product_equals(&lhs, &rhs);
  printf("Full beta check:  %d %12d*\n", ok, nRoot()); 
  result &= ok;

  printf("Verified full commitment %d, result %d\n", bank_id, result);
  return result;
}


bool verify_proof(cVerifKey *vk, cProof *proof, int num_commits, cCommitment* commits) {
	// Accumulators
	struct s_EncodedEltL encoded_V;
	struct s_EncodedEltR encoded_W;
	struct s_EncodedEltL encoded_Y;
	_encoding_initL(&encoded_V);
	_encoding_initR(&encoded_W);
	_encoding_initL(&encoded_Y);
	set_zeroL(&encoded_V);
	set_zeroR(&encoded_W);
	set_zeroL(&encoded_Y);
  for (int i = 0; i < num_commits; i++) {
		_encoding_add_rawL(&commits[i].V, &encoded_V, &encoded_V);
		_encoding_add_rawR(&commits[i].W, &encoded_W, &encoded_W);
		_encoding_add_rawL(&commits[i].Y, &encoded_Y, &encoded_Y);
		//printf("Adding commitment %d\n", i);
		//enc_print_rawL("V", &commits[i].V);
		//enc_print_rawR("W", &commits[i].W);
		//enc_print_rawL("Y", &commits[i].Y);
	}

  //printf("Final summation commitment: \n");
  //enc_print_rawL("V", &encoded_V);
  //enc_print_rawR("W", &encoded_W);
  //enc_print_rawL("Y", &encoded_Y);

	struct s_EncodedEltR encoded_one_R;
	_encoding_initR(&encoded_one_R);
	copyR(&e_genR, &encoded_one_R);

	EncodedProduct vw, y1, rhs;
	encoded_product_init(&vw);
	encoded_product_init(&y1);
	encoded_product_init(&rhs);

	encoding_pair(&encoded_V, &encoded_W, &vw);
	encoding_pair(&encoded_Y, &encoded_one_R, &y1);
	encoding_pair(&proof->H, &vk->D, &rhs);

	encoding_sub(&vw, &y1, &vw);

	int ok = encoding_product_equals(&vw, &rhs);
	printf("Verified proof, result %d\n", ok);
	return(ok);
}

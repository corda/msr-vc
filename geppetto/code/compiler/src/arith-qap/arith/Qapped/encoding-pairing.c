typedef struct {
	QapFp12_t e;
} EncodedProduct;

void encoded_product_init(EncodedProduct* r)
{
	QapFp12_init(r->e);
}

int nRoot(); //ugly

void encoding_pair(struct s_EncodedEltL* L1, struct s_EncodedEltR* R1, EncodedProduct* r)
{
	testOnCurveL(L1);
	testOnCurveR(R1);

	//enc_print_rawL("PAIRING L", L1);
	//enc_print_rawR("PAIRING R", R1);

	// Library doesn't handle g^0 values properly, so we special-case it
	// we used to exclude it using
	// zeroAssert(isZeroL(L1));
	// zeroAssert(isZeroR(R1));

	int zL = isZeroL(L1);
	int zR = isZeroR(R1);
	int b = (1 - zL) * (1 - zR);

	QapFp12_t x, z; // we compute x anyway, then select b ? x : z 
	QapFp12_init(x);
	QapFp12_init(z);
	QapFp12_set_one(z);

	{
		int n;
		n = nRoot();

		Qapbn_optimal_ate_miller_aff(x, R1->e, L1->e);
		// ~ 3,961* (Left) or ~ 6,216 (Both)  
		n = nRoot();

		Qapbn_finalexpo_neg(x, x);
		n = nRoot();
		// ~ 4,500* --AFF--> 4,472 (Left or Both) unchanged
		// TOTAL * COST:
		// projective was ~ 11,500 (?)
		// affine is: 8,433* (Left) or 10,688* (Both) 
	}
	QapFp12_select(r->e, b, x, z);
}

typedef struct {
  QapFp2_t premem[PREMEM_SIZE];
} pairing_precomp_cache;

void encoding_pair_prep(struct s_EncodedEltR* R1, pairing_precomp_cache cache)
{
  testOnCurveR(R1);

#ifndef MQAP 
  // Recheck it is Ok to omit this one!
  if (isZeroR(R1)) {
    // Library doesn't handle g^0 values properly, so we special-case it
    assert(0); // for MQAP sync!
    return;
  }
#endif
  int n;
  n = nRoot();

#ifndef AFF
  assert(0);  // We only have this implemented for the affine version at present
#endif

  Qapbn_optimal_ate_miller_precompute_aff(R1->e, cache.premem);
  n = nRoot();
}


void encoding_pair_with_prep(struct s_EncodedEltL* L1, pairing_precomp_cache cache, EncodedProduct* r)
{
  testOnCurveL(L1);

#ifndef MQAP 
  // Recheck it is Ok to omit this one!
  if (isZeroL(L1)) {
    // Library doesn't handle g^0 values properly, so we special-case it
    assert(0); // for MQAP sync!
    QapFp12_set_one(r->e);
    return;
  }
#endif
  int n;
  n = nRoot();

  Qapbn_optimal_ate_miller_useprecomputed_aff(r->e, L1->e, cache.premem);
  n = nRoot();
}

void encoding_print(const char* msg, EncodedProduct* a)
{
	printf("%s:\n", msg);
	QapFp12_print(a->e);
	printf("\n");
}

void encoding_add(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r)
{
	QapFp12_mul(r->e, a->e, b->e);
}


void encoding_sub(EncodedProduct* a, EncodedProduct* b, EncodedProduct* r)
{
  EncodedProduct inv;
  encoded_product_init(&inv);
  QapFp12_inv(inv.e, b->e);
  QapFp12_mul(r->e, a->e, inv.e);
}

bool encoding_product_equals(EncodedProduct* a, EncodedProduct* b)
{
	return QapFp12_cmpeq(a->e, b->e);
}


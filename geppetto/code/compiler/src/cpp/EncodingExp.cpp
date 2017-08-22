#include <assert.h>
#include "EncodingExp.h"
#include <bitset>
#include <intrin.h>
#include "Poly.h"
#include "Config.h"
#include <stdint.h>

#pragma intrinsic(_BitScanReverse)

Field* fequal::field;
Field* fhash::field ;

EncodingExp::EncodingExp() {
	// Clear out the encoding tables used to speed up encoding operations
	encodingTable[0] = encodingTable[1] = NULL;
	precomputedPowers = false;
	h_ctr = 0;
  cache_hits = 0;
  cache_attempts = 0;
  fequal::field = NULL;
  fhash::field = NULL;
}

EncodingExp::~EncodingExp() {
	doneWithManyEncodings();
}

Field* EncodingExp::getSrcField() {
	return srcField;
}

void EncodingExp::doneWithManyEncodings() {
	if (encodingTable[L] != NULL) {
		for (uint32_t i = 0; i < v; i++) {
			delete [] encodingTable[L][i];
			delete [] encodingTable[R][i];
		}

		delete [] encodingTable[L];
		delete [] encodingTable[R];

		encodingTable[L] = NULL;
		encodingTable[R] = NULL;
	}

  elt_t t = L;
  for (EncodedCache::iterator iter = enc_cache[t].begin();
    iter != enc_cache[t].end();
    iter++) {
    delete iter->second;
  }
  t = R;
  for (EncodedCache::iterator iter = enc_cache[t].begin();
    iter != enc_cache[t].end();
    iter++) {
    delete iter->second;
  }

	precomputedPowers = false;
}

unsigned long find_msb(unsigned long x) {
	if (x == 0) {
		return -1;
	} else {
		unsigned long ret;
		_BitScanReverse(&ret, x);
		return ret;
	}

#ifdef GCC
	if (x == 0) {
		return 0;
	} else {
		return sizeof(unsigned long) * 8 - __builtin_clz(x);
	}
#endif
}

/*
	Doesn't always seem to choose the optimal value of h.  

	h  Pre   LEnc  REnc  Total
	6	0.15	0.062	0.258	0.47
	7	0.022	0.051	0.218	0.291
	8	0.032	0.055	0.185	0.272   <- Actual min
	9	0.061	0.041	0.197	0.299
	10	0.118	0.037	0.153	0.308  <- Alg below choses this
	11	0.253	0.036	0.144	0.433
*/

void EncodingExp::allocatePreCalcMem(uint64_t numExps, uint32_t maxMemGb, bool precompFree) {
	double maxMem = (double) maxMemGb * (((long long)1)<<30);

	if (precompFree) {
		// Build the largest table that will fit in memory
		for (h = 1; h < num_exponent_bits(); h++) {
			uint64_t twoToH = ((uint64_t)1)<<h;
			uint64_t storage = (twoToH - 1) * h * (bytes_per(L) + bytes_per(R));
			if (storage > maxMem) {
				h--;
				break;				
			} 
		}
		v = h;
	} else {
		// Numerically minimize the expected number of multiplications, subject to the memory limit defined above
		uint64_t minMults = INT64_MAX;
		int minH = -1;
		for (uint32_t hOption = 1; hOption < num_exponent_bits(); hOption++) {
			uint64_t twoToH = ((uint64_t)1)<<hOption;
			uint64_t storage = (uint64_t)(twoToH - 1)*(uint64_t)hOption * (uint64_t)(bytes_per(L) + bytes_per(R));
			if (storage > maxMem) {
				break;
			}

			// See qsp-practical/notes/optimizing-exp.txt for details
			//uint64_t numMults = (int)((hOption-1)/(float)hOption * num_exponent_bits/2.0 + (1<<(hOption-1))*(hOption+1) + numExps*(num_exponent_bits*hOption + 1.5*num_exponent_bits - 2 * hOption*hOption) / (float)(hOption*hOption));
			uint32_t a = (uint32_t)ceil(num_exponent_bits() / (float) hOption);
			uint32_t b = (uint32_t)ceil(num_exponent_bits() / (float)(hOption*hOption));
		
			uint64_t numMults = (uint64_t) (numExps * (((twoToH-1)/(float)twoToH) * a + 1.5*b - 2)		// For the actual exponentations
													+ (hOption-1)*a/2.0 + twoToH * (hOption+1)/2.0);
			if (numMults < minMults) {
				minH = hOption;
				minMults = numMults;
			}
		}

		h = minH;
		v = minH;
	}

	INFO("Chose h=v=%d for preallocation optimization\n", h);

	a = (int)ceil(num_exponent_bits() / (float) h);
	b = (int) ceil(a / (float) v);

	encodingTable[L] = new EncodedEltArray*[v];
	encodingTable[R] = new EncodedEltArray*[v];

	for (uint32_t i = 0; i < v; i++) {
		encodingTable[L][i] = new_elt_array(L, ((long long)1)<<h, true);
		encodingTable[R][i] = new_elt_array(R, ((long long)1)<<h, true);
	}	
}


void EncodingExp::prepareForManyEncodings(uint64_t numExps, uint32_t maxMemGb, bool precompFree) {
	if (precomputedPowers) { return; } // Already precomputed

	allocatePreCalcMem(numExps, maxMemGb, precompFree);

	//printf("Precomputing powersFast with h = %d, v = %d\n", h, v);
	INFO("Size of precomputed table = %d x %d = %d {L,R}EncodedElts = %.01f MB\n", v, 1 << h, v * (1<<h), v * (1<<h) * (bytes_per(L) + bytes_per(R)) / (1024.0 * 1024));

	prepareForManyEncodingsHelper(L);
	prepareForManyEncodingsHelper(R);

	precomputedPowers = true;
}

// Implements more flexible exponentiation with precomputation, C.H.Lim and P.J.Lee, Advances in Cryptology-CRYPTO'94, http://dasan.sejong.ac.kr/~chlim/pub/crypto94.ps
// which is a special case of Pippenger's algorithm, according to http://cr.yp.to/papers/pippenger.pdf
void EncodingExp::prepareForManyEncodingsHelper(elt_t t) {
	// Memoize: PowersFast[i] = x^{2^{i*a}}
	EncodedEltArray* PowersFast = new_elt_array(t, num_exponent_bits());	

	// Do the first one directly
	one(t, PowersFast->elt(0));	

	// Calculate x_i = x^{2^{i*a}}.  Takes (p-1)*a squarings
	for (uint32_t i = 1; i < h; i++) {
		for (uint32_t j = 0; j < a; j++) {
			if (j == 0) {
				doubleIt(t, PowersFast->elt(i-1), PowersFast->elt(i)); 
			} else {
				doubleIt(t, PowersFast->elt(i), PowersFast->elt(i)); 
			}
		}
	}
	
	// Handle G[0][0] directly
	zero(t, encodingTable[t][0]->elt(0));	

	for (uint32_t i = 1; i < ((uint32_t)1 << h); i++) {
		int upperMostOneBit = find_msb(i);
		this->copy(t, PowersFast->elt(upperMostOneBit), encodingTable[t][0]->elt(i));

		// Toggle the msb, and use the rest as an index into previous G calculations
		int lowerBits = (i ^ (1 << upperMostOneBit));
		if (i != 0 && lowerBits != 0) {
			add(t, encodingTable[t][0]->elt(i), encodingTable[t][0]->elt(lowerBits), encodingTable[t][0]->elt(i));			
		}
	}

	for (uint32_t j = 1; j < v; j++) {
		for (uint32_t i = 0; i < ((uint32_t)1 << h); i++) {
			// G[j][i] = (G[j-1][i])^{2^b}
			for (uint32_t k = 0; k < b; k++) {
				if (k == 0) {
					doubleIt(t, encodingTable[t][j-1]->elt(i), encodingTable[t][j]->elt(i)); 					
				} else {
					doubleIt(t, encodingTable[t][j]->elt(i), encodingTable[t][j]->elt(i)); 								
				}
			}
		}
	}

	delete PowersFast;
}


int EncodingExp::computeIndex(FieldElt* in, int j, int k) {
	assert(in);
	int index = 0;

	for (int i = h - 1; i >= 0; i--) {
		int exp_index = i*a + b*j + k;

		if (b * j + k < a) {	// Otherwise, we've exceeded the length of this a-length block
			index <<= 1;
			index |= srcField->bit(exp_index, in);
		}
	}

	return index;
}

// Raise the generator of the base group to the element provided: out <- g_1^in
//#define CACHE_ENCODE
void EncodingExp::encode(elt_t t, FieldElt* in, EncodedElt* out) {
	assert(in && out);
#ifdef CACHE_ENCODE
  if (fhash::field == NULL) {
    fhash::field = getSrcField();
    fequal::field = getSrcField();
  }
  EncodedCache::iterator cache_result = enc_cache[t].find(in);
  cache_attempts++;

  if (cache_result != enc_cache[t].end()) {
    cache_hits++;
    copy(t, cache_result->second, out);
//#ifdef DEBUG
//    EncodedElt* test_elt = this->new_elt(t);
//    encodeSlow(t, in, test_elt);
//    assert(this->equal(t, cache_result->second, test_elt));
//    this->del_elt(t, test_elt);
//#endif //DEBUG
  } else {
#endif // CACHE_ENCODE
    if (precomputedPowers) {
      copy(t, encodingTable[t][0]->elt(0), out);		// out <- g^0 = 1

      for (int k = b - 1; k >= 0; k--) {
        if (k != b - 1) {
          doubleIt(t, out, out);		// No point the first time around
        }
        for (int j = v - 1; j >= 0; j--) {
          int index = computeIndex(in, j, k);
          if (index != 0) {	// lG[j][0] = g^0 = 1
            add(t, out, encodingTable[t][j]->elt(index), out);
          }
        }
      }
    } else {
      //printf("Warning, encoding via a slow routine\n");
      encodeSlow(t, in, out);
    }
#ifdef CACHE_ENCODE
    Field* field = this->getSrcField();
    FieldElt* cached_key = field->newElt();
    field->copy(in, cached_key);
    EncodedElt* cached_elt = this->new_elt(t);
    this->copy(t, out, cached_elt);
    enc_cache[t].insert(CachedEncoding(cached_key, cached_elt));
  }
#endif // CACHE_ENCODE

}


// TODO: Reduce memory footprint by not storing g^0 in every table!
ip_handle_t EncodingExp::prepareForInnerProduct(EncodedEltArray* bases, int len, int numUses, int max_exponent_size, int maxMemGb, bool precompFree) {
	double maxMem = (double) maxMemGb * (((long long)1)<<30);
	int numTerms = 2;		// Can't handle more than 2
	int numBits = 1;
	elt_t t = bases->t();

	MultiExpPrecompTable* table = new MultiExpPrecompTable;
	table->t = t;

	if (len == 0) { 
		assert(0); // Why did we end up here? 
		table->precompTables = NULL; 
		return (uint32_t)-1; 
	}

	uint64_t numTables = (int)ceil(len / (float)numTerms);
	uint64_t minCost = UINT64_MAX;
  int max_bits = min(max_exponent_size, (int)num_exponent_bits());
	while (true) {
		uint64_t perTableSize = 1 << (numBits * numTerms);		// 2^{numBits*numTerms}
		uint64_t totalSizeInB = perTableSize * numTables * bytes_per(t);

		uint64_t preCompDbls = (uint64_t) (numTables * ((1 << numBits) / 2.0 + (1 << numBits) / 2.0));		// Slight under-counting, since it doesn't include the (j=k and j/2=k/2=0 mod 2) case
		uint64_t preCompMults = numTables * (perTableSize - numTerms - 1) - preCompDbls;	// Worst-case: 1 mult per table entry.  Each base (plus 0), is simply a copy.  We replace some mults with doublings

    uint64_t onlineMults = numUses * (uint64_t)(max_bits * len / (float)(numTerms * numBits) - 1);	// (DW-1)
    uint64_t onlineDbls = numUses * (uint64_t)(max_bits / (float)numBits);	// (W-1)

		uint64_t preCost = preCompMults * 3 + preCompDbls;		// Mults are ~3x dbls
		uint64_t onlineCost = onlineMults * 3 + onlineDbls;
		uint64_t totalCost = preCost + onlineCost;

		if (totalSizeInB > maxMem || (!precompFree && totalCost > minCost)) {  // We're always constrained by memory.  If precomp isn't free, we're also constrained to avoid spending too much time on precomp
			// Too big!  Back up.
			numBits--;
			break;
		} else {
			minCost = totalCost;
			numBits++;
		}
	}

	table->numTerms = numTerms;
  table->numBits = min(max_exponent_size, numBits);  // No need to prep for more bits than will be present in the exponents
  if (max_exponent_size < numBits) { 
    INFO("Saving effort by computing tables for %d bits instead of %d bits\n", max_exponent_size, numBits); 
  }
	table->len = len;
	table->numTables = (int)ceil(len / (float)table->numTerms);
	int tableSize = 1 << (numBits * numTerms);		// 2^{numBits*numTerms}
	table->sizeInMB = (uint64_t) (tableSize * table->numTables * bytes_per(t) / (1024.0 * 1024.0));
	INFO("Precomputing multi-exp tables: %d bits from %d terms over %d bases. Size=%dMB\n", numBits, numTerms, len, table->sizeInMB);
	
	table->precompTables = new EncodedEltArray*[table->numTables];
	for (int i = 0; i < table->numTables; i++) {
		table->precompTables[i] = new_elt_array(bases->t(), tableSize, true);
		int baseIndex = numTerms * i;
		bool oddNumBases = false;

		if (baseIndex >= len - 1) {	// Deal with the case where there are an odd number of bases
			baseIndex--;
			oddNumBases = true;
		}

		zero(t, table->precompTables[i]->elt(0));  // g^0		
		// Compute powers of the first base
		for (int j = 1; j < (1 << numBits); j++) {
			if (j == 1) {
				copy(t, bases->elt(baseIndex + 1), table->precompTables[i]->elt(j));	// g_2^1
			} else if (j % 2 == 0) {	// Double instead of adding, since doubling is cheaper
				doubleIt(t, table->precompTables[i]->elt(j/2), table->precompTables[i]->elt(j));
		  }	else {				
				add(t, bases->elt(baseIndex + 1), table->precompTables[i]->elt(j-1), table->precompTables[i]->elt(j)); // g_2^j
			}
		}

		if (oddNumBases) { // If there are an odd # of bases, then we're done
			break;
		}

		// Multiply by powers of the second base
		for (int j = 1; j < (1 << numBits); j++) {
			int oldTableIndex = (j-1) * (1 << numBits);
			int tableIndex = j * (1 << numBits);
			if (j == 1) {
				copy(t, bases->elt(baseIndex), table->precompTables[i]->elt(tableIndex));	// g_1^1
			} else if (j % 2 == 0) {	// Double instead of adding, since doubling is cheaper
				oldTableIndex = (j/2) * (1 << numBits);
				doubleIt(t, table->precompTables[i]->elt(oldTableIndex), table->precompTables[i]->elt(tableIndex));
			} else {
				add(t, bases->elt(baseIndex), table->precompTables[i]->elt(oldTableIndex), table->precompTables[i]->elt(tableIndex)); // g_1^j
			}

			for (int k = 1; k < (1 << numBits); k++) {
				if (j == k && j % 2 == 0 && k % 2 == 0) {	// We can squeeze out an additional doubling as g_1^j * g_2^k = (g_1^j/2 * g_2^k/2)^2
					int halfIndex = (j/2) * (1 << numBits) + k/2;
					doubleIt(t, table->precompTables[i]->elt(halfIndex), table->precompTables[i]->elt(tableIndex + k));
				} else {
					add(t, table->precompTables[i]->elt(tableIndex), table->precompTables[i]->elt(k), table->precompTables[i]->elt(tableIndex + k));		// g_1^j * g_2^k
				}
			}
		}
	}	

	// Add the table to our collection and return a handle to it
	ip_handle_t ret = h_ctr++;
	PrecompTables.insert(MultiExpTable(ret, table));
	return ret;
}



// Computes r <- g_{a_i}^{b_i} for 0 <= i < len
// Computes "horizontally" in the sense that for each chunk of the exponent space,
// we combine the appropriate table entries for all exponent-base sets,
// and then we use doubling to "move" the values into place.
// This saves a factor of chunkSize doublings compared to the vertical approach

void EncodingExp::innerProduct(ip_handle_t handle, FieldEltArray* exp, int len, EncodedElt* r) {
	assert(exp && r);

	// Uses the precomputed tables to step through each set of bases and each set of exponents faster
	MultiExpTables::iterator search_result = PrecompTables.find(handle);
	assert(search_result != PrecompTables.end());
	MultiExpPrecompTable* precompTables = search_result->second;
	elt_t t = precompTables->t;

	EncodedElt* tmp = new_elt(t);

	int chunkSize = (num_exponent_bits()-1 - 1) / precompTables->numBits + 1; // ceil((num_exponent_bits-1) / numBits)
	zero(t, r);  // g^0
	bool rIsZero = true;		// As long as r remain zero, we can skip many ops
	bool tmpIsZero = true;

	for (int j = 0; j < chunkSize; j++) {	// For each chunk of exponent bits
		zero(t, tmp);  // g^0
		tmpIsZero = true;

		int bitIndex = max(((int)num_exponent_bits()-1) - ((j+1) * precompTables->numBits), 0);
		int prevBitIndex = ((int)num_exponent_bits()-1) - ((j) * precompTables->numBits);
		int numBitsToGet = min(precompTables->numBits, prevBitIndex);

		if (!rIsZero) {
			for (int k = 0; k < numBitsToGet; k++) {
				doubleIt(t, r, r);
			}
		}

		for (int i = 0; i < precompTables->numTables && i*precompTables->numTerms < len ; i++) {	// For each grouping of numTerm terms
			int expIndex = i * precompTables->numTerms;

				// Build an index into the table, based on the values of the exponents
			int tableIndex = 0;	 			
			for (int k = expIndex; k < expIndex + precompTables->numTerms && k < len; k++) {	// For each exponent in this group
				// Grab the appropriate bits of the current exponent				
				int bits = srcField->getBits(exp->elt(k), bitIndex, numBitsToGet);

				tableIndex <<= precompTables->numBits;	// Make room for the new contribution (we still use numBits to shift, since that's how the table is built)
				tableIndex |= bits;
			}

			// How many bases are in this chunk?
			int leftOver = 0;
			if (len <= (precompTables->len / precompTables->numTerms) * precompTables->numTerms) {	// Len is not in the (possible) leftover chunk at the end
				leftOver = precompTables->numTerms - (len - expIndex);
			} else {	// Can only have as many leftovers as there were assigned the extra chunk at the end to begin with
				leftOver = (precompTables->len % precompTables->numTerms) - (len - expIndex);
			}

			// Adjust tableIndex to account for missing terms in this chunk
			for (int k = 0; k < leftOver; k++) {
				tableIndex <<= precompTables->numBits;
			}

			//printf("Table index = %d\n", tableIndex);
			if (tableIndex != 0) {	// No point adding a 0 term
				if (tmpIsZero) {	// Copy instead of adding 0
					copy(t, precompTables->precompTables[i]->elt(tableIndex), tmp);
					tmpIsZero = false;
				} else {
					add(t, tmp, precompTables->precompTables[i]->elt(tableIndex), tmp);			
				}
				rIsZero = false;
			} 
		}

		if (!tmpIsZero) {  // No point adding a 0 term
			if (rIsZero) { // Copy instead of adding 0
				copy(t, tmp, r);
			} else {
				add(t, r, tmp, r);
			}
		}		
	}

	del_elt(t, tmp);
}

void EncodingExp::doneWithInnerProduct(ip_handle_t handle) {
	MultiExpTables::iterator search_result = PrecompTables.find(handle);
	assert(search_result != PrecompTables.end());
	MultiExpPrecompTable* precompTables = search_result->second;

	for (int i = 0; i < precompTables->numTables; i++) {
		delete [] precompTables->precompTables[i];
	}

	delete [] precompTables->precompTables;

	delete precompTables;
}

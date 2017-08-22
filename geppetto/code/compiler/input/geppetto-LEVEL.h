#pragma once 

/* for bootstrapping, the plan is to substitute _LEVEL with "", then "_BOOT", then "_2" ...
   we just use replace/regexp, as macros are not recursive */

/* We natively run C programs with generated definitions for banks and
outsourced calls in these modes:

DEFAULT the "reference code", as simple as possible, used 
    as a basis for interpreted runs. Banks are just indexed buffers; 
	outsourcing simply updates function and proof indexes before
	just calling the function. 

SCHEDULE collects additional information about outsourced calls and banks,
	so that we can run keygen & emit custom data for the verifier.
	This mode is lightweight, hopefully just macro-expansions and printfs.

VERIFY verifies an outsourced run, given keys, and evidence (currently
	many small files). This mode requires linking to the crypto
	libraries, key information issued by Gepetto during the
	compilation, and schedule information produced by (G).
	This was "specified by example" in boot-loop.c

Geppetto also interprets C programs in more complex modes. 
To this end, we "clang -DMQAP" the source program that includes this file. 
The main difference from DEFAULT is that the functions 

	load_  store_  __enter  __exit  verify

become primitives managed by the interpreter. We currently have 3 interpreted modes:

COMPILE to compile the outsourced functions to QAPs (QAPGen).

PROVE to prove an outsourced run and generate evidence (QAPEval)

VERIFY to verify an outsourced run (QAPEval too, a bit redundant, but good sanity check)

We try to keep all these modes as close as possible. 
Note that Geppetto depends on the details of this file!

*/

#ifdef MQAP_LEVEL
#ifndef BOOTSTRAP_LEVEL
#define INTERPRET_LEVEL
#endif
int printf(const char *, ...);
void* malloc(unsigned int);
#else // !MQAP_LEVEL
#include <stdio.h>
#include <stdlib.h>
#endif

#ifdef SCHEDULE
#define trace(...) fprintf(STATE.trace, __VA_ARGS__); printf(__VA_ARGS__)
#else
#define trace(...) printf(__VA_ARGS__)
#endif 

#ifndef MAX_BANK_TYPES // this must be manually set for "wide" mQAPs
#define MAX_BANK_TYPES 20
#endif
#ifndef MAX_OUTSOURCE_CALLS // this must be manually set for larger numbers of proofs
#define MAX_OUTSOURCE_CALLS 100 
#endif 
// NUM_BANKS_LEVEL must be precisely set (otherwise the function indexes are off)
// #define NUM_BANKS_LEVEL 3

#define GEPPETTO_DATA_DIR "./build/"
#define TOPLEVEL (-1)

#ifdef VERIFY_LEVEL
#include "common.h"
#include GEPPETTO_VERIFY_LEVEL_H 
#include "encoding-qap.c"   // Must come after definition of NUM_BANKS_TOTAL and SUM_BANK_SIZES above

// TODO: Sync indexing below with compiler-written values
// TODO: Compiler needs to output a .h file with OFFSET_##BT for each bank to point to bank's location in VK
// TODO: Compiler-generated .h file needs to specify NUM_BANKS
// TODO: Need to obtain reader/writer data from G
// TODO: Need to obtain recompute data from G? Compiler?


// Finally, we need reader, writer, and recompute arrays, each NUM_BANKS x MAX_OUTSOURCE_CALLS, filled in based on G and/or compiler

// When generating commitments for non-local banks, the compiler should index each commitment using:
//    b * MAX_OUTSOURCE_CALLS + tBank[b]
// where b is the bank's ID and tBank[b] is the ID for this instance of the bank
// This is true for the IDs of files containing saved bank values (e.g., for recomputation)

/////////////////////////////////////////////////////
#endif  

typedef struct {
	int function_id;  // function being outsourced, if any, or TOPLEVEL otherwise 
	int outsource_id; // sequence number for complete proofs, starting from 0
	int tBank[MAX_BANK_TYPES];                    // last instance of t allocated (recall we start at 1) 
	// as a sanity check, we could also enforce the bank discipline: a single store followed by loads
#ifdef VERIFY_LEVEL
	cVerifKey vk; // verification key
	int ok; // cumulative checks so far (could be eliminated)
#endif
#ifdef SCHEDULE_LEVEL	
	FILE* trace;
#endif
} GeppettoState_LEVEL;

GeppettoState_LEVEL STATE_LEVEL;

// primitive functions to control prover/verifier interpretation 
// the integer parameters are function_id and outsource_id
int  __Enter_LEVEL(int, int);
void __Exit_LEVEL(int, int);
void* Verify_LEVEL(int, int, ...);

#ifdef INTERPRET_LEVEL
// the only difference is that we leave load_, store_ and Verify undefined (they have a primitive interpretation)
#define BANK_LEVEL(BT,BR)                                                                                         \
	const int index_##BT = __COUNTER__;                                                                           \
	typedef struct { int t; BR value; } Bank_##BT;                                                          \
	typedef Bank_##BT* BT;                                                                                  \
	BT alloc_##BT(void) { BT b = malloc(sizeof(Bank_##BT)); b->t = ++STATE_LEVEL.tBank[index_##BT]; return b; }  \
	void load_##BT(BT b, BR* p);                                                                            \
	void store_##BT(BT b, BR* p);                                                                           \
	BT save_##BT(BR* p) { BT b = alloc_##BT(); store_##BT(b, p); return b; }
#else
// when interpreted by the prover, initialize the QAP state for the outsourced call
// when interpreted by the verifier, verify the proof instead of running f
int __Enter_LEVEL(int ctr, int idx) { return 1; }
void __Exit_LEVEL(int ctr, int idx) {}
void* Verify_LEVEL(int i, int j, ...) {
	printf("verification placeholder for %d. %d\n",i,j);
	exit(-1);
}
#ifdef VERIFY_LEVEL
// this variant differs by maintaining an extra commitment in each bus
// and managing them when loading/saving at toplevel
#define BANK_LEVEL(BT,BR)  \
	const int index_##BT = __COUNTER__;                                   \
	typedef struct { int t; BR value; cCommitment c; } Bank_##BT;          \
	typedef Bank_##BT* BT;                                                \
	BT alloc_##BT(void) { BT b = malloc(sizeof(Bank_##BT)); init_commit(&(b->c)); b->t = ++STATE_LEVEL.tBank[index_##BT]; return b; }  \
	void load_##BT(BT b, BR* p);                                                                            \
	void store_##BT(BT b, BR* p);                                                                           \
	BT save_##BT(BR* p) { BT b = alloc_##BT(); store_##BT(b, p); return b; }\
	\
	void load_##BT(BT b, BR* p) {\
		trace_LEVEL("[%2d.%2d ]_LEVEL load %s %d\n", STATE_LEVEL.function_id, STATE_LEVEL.outsource_id, #BT, b->t);     \
		*p = b->value;                                                                                      \
	};                                                                    \
	\
	void store_##BT(BT b, BR* p) {\
		compute_commit(&STATE_LEVEL.vk, &b->c, C_##BT, R_##BT, (int*)p, O_##BT, L_##BT); \
		b->value = *p;                                                                                      \
	}; \
	BT recommit_verify_##BT() { \
		BT b = alloc_##BT();\
		if (recommit_##BT[b->t]) { \
			load_int32s(#BT, b->t, sizeof(BR)/4, (int*)(&b->value), RUN_TIME); \
			compute_commit(&STATE_LEVEL.vk, &b->c, C_##BT, R_##BT, (int*)(&b->value), O_##BT, L_##BT); \
		} else { \
			load_cCommitment(#BT, &b->c, b->t, RUN_TIME); \
			STATE_LEVEL.ok *= verify_commit_bus(&STATE_LEVEL.vk, &b->c, index_##BT ); \
		} \
		return b; \
	}
#else 
// neither INTERPRET nor VERIFY
#define BANK_LEVEL(BT,BR)                                                                                     \
	const int index_##BT = __COUNTER__;                                                                   \
	typedef struct { int t; BR value; } Bank_##BT;                                                        \
	typedef Bank_##BT* BT;                                                                                \
	BT alloc_##BT(void) { BT b = malloc(sizeof(Bank_##BT)); b->t = ++STATE_LEVEL.tBank[index_##BT]; return b; }  \
	void load_##BT(BT b, BR* p);                                                                            \
	void store_##BT(BT b, BR* p);                                                                           \
	BT save_##BT(BR* p) { BT b = alloc_##BT(); store_##BT(b, p); return b; }\
	\
	void load_##BT(BT b, BR* p) {\
		trace_LEVEL("[%2d.%2d ]_LEVEL load %s %d\n", STATE_LEVEL.function_id, STATE_LEVEL.outsource_id, #BT, b->t);     \
		*p = b->value;                                                                                      \
	};                                                                                                    \
	void store_##BT(BT b, BR* p) {         \
		trace_LEVEL("[%2d.%2d ]_LEVEL save %s %d\n", STATE_LEVEL.function_id, STATE_LEVEL.outsource_id, #BT, b->t);     \
		b->value = *p;                                                                                      \
	};
#endif // VERIFY
#endif // INTERPRET

void init_LEVEL() {
	STATE_LEVEL.function_id = TOPLEVEL;
	STATE_LEVEL.outsource_id = 0;
	for (int t = 0; t < MAX_BANK_TYPES; t++) STATE_LEVEL.tBank[t] = 0;

#ifdef SCHEDULE_LEVEL
	char filename[200];
	sprintf(filename, "%s/%s.trace", GEPPETTO_DATA_DIR, GEPPETTO_NAME_LEVEL);
	STATE_LEVEL.trace = fopen(filename, "w");
	if (!STATE_LEVEL.trace) {
		printf("Failed to open trace file!\n");
		exit(2);
	}
#endif 
#ifdef VERIFY_LEVEL
	STATE_LEVEL.ok = true;    // Start by assuming everything is correct
	common_init();
	load_cVerifKey(&STATE_LEVEL.vk, 0, COMPILE_TIME);
#endif
}

int Enter_LEVEL(int ctr, const char* f, char* args) {
	STATE_LEVEL.function_id = ctr - NUM_BANKS;
	trace_LEVEL("[%2d.%2d ]_LEVEL call %s %d\n", STATE_LEVEL.function_id, STATE_LEVEL.outsource_id, f, STATE_LEVEL.function_id);
	return __Enter_LEVEL(STATE_LEVEL.function_id, STATE_LEVEL.outsource_id);
}
void* Exit_LEVEL(void* p) {
	trace_LEVEL("[%2d.%2d ]_LEVEL return\n", STATE_LEVEL.function_id, STATE_LEVEL.outsource_id);
	__Exit_LEVEL(STATE_LEVEL.function_id, STATE_LEVEL.outsource_id);
	STATE_LEVEL.function_id = TOPLEVEL;
	STATE_LEVEL.outsource_id++;
	return p;
}

// We need a conditional to prevent any inlining of the outsourced call
// TODO simplify VERIFY, as we always call the custom verifier
#ifdef VERIFY_LEVEL
#define OUTSOURCE_LEVEL(F,...) \
	Exit_LEVEL(\
	Enter_LEVEL(__COUNTER__, #F, #__VA_ARGS__) ? \
		Verify_##F(STATE_LEVEL.outsource_id, __VA_ARGS__) : \
		Verify_##F(STATE_LEVEL.outsource_id, __VA_ARGS__))
#else
// this definition covers all other modes
#define OUTSOURCE_LEVEL(F,...) \
	Exit_LEVEL(\
		Enter_LEVEL(__COUNTER__,#F,#__VA_ARGS__) ? \
		F(__VA_ARGS__) : \
		Verify_LEVEL(STATE_LEVEL.function_id, STATE_LEVEL.outsource_id, __VA_ARGS__))
#endif







/* Given a dump of outsource, and of the commitment schedule generated
   by the bank definitions below, we can easily reconstruct all the
   information we need (and perform sanity checks) in F#. 

   Ideally, the information issued by F# should *not* depend on the actual
   schedule though.
*/

/*
                        BANKS AND COMMITMENTS

BANK(bank,t) declares & implements a bank type, the only kind of types
enabled as arguments and results of outsourced calls.

  bank is a fresh type identifier for the resulting type
  t is an existing C type identifier, consisting of integers, i.e.
  t must support typecasts from/to int[sizeof(t)/4].

Geppetto provides 2 functions to operate on bank types [hiding store/alloc for now]

void load_bank(bank b, t* p); // copy b's content into *p
bank save_bank(t* content);   // creates a bank with contents *p
 
For each bank type, we internally keep track of 
  the number of instances allocated so far (t, starting from 1) and,
  for each instance, 
    its index t, 
    its saved content, 
    its writer, and 
    its readers (ranging over TOPLEVEL,0,1,...STATE.num_outsource_calls-1)

As the program completes, we know what the verifier must do with each
bank; we may record that information e.g. as constant arrays for the C
verifier. (For instance, we can't select between OUTPUT and SHARED yet
as the bank is saved.)

INPUT  if written by TOPLEVEL, the verifier recomputes the commitment (on the fly)

OUTPUT if written by i and read by TOPLEVEL, the verifier loads the contents
       from a file in the evidence and recomputes the commitment
       during the i-th call verification.

SHARED if written by i and accessed by other outsourced calls, the
       verifier loads and verifies its commiment from a file in the
       evidence during the i-th call verification.

       (if written by i and not accessed elsewhere, issue a warning.)

More generally, Geppetto/KeyGen can use all reader-writer constraints
to select the kind of banks to use (we sometimes have several
reasonable choices; SHARED tells us which QAPs are sharing the bank;
we need a bus if there are at least two; otherwise a LOCAL bank may be
enough. I/O banks may be buses, or not.) Assume Geppetto --keygen
issues all these details as a header for the C verifier.

*/


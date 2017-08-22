#include "pinocchio.h" 
#ifndef MQAP
#include <stdio.h>
#include <string.h>
#include "pinocchio.c"
#endif

#define GEPPETTO_NAME "whole"
#define GEPPETTO_VERIFY_H "build/whole.verify.h"
#define GEPPETTO_VERIFY_C "build/whole.verify.c"
#define NUM_BANKS 2
#include "geppetto.h"

struct Sample {
	int a;
	int b;
	int c;
};

BANK(SampleIN, struct Sample)
BANK(SampleOUT, struct Sample)
//BANK(SampleEXTRA, struct Sample)

#ifdef VERIFY
#include GEPPETTO_VERIFY_C
#endif

//#define BOUND(X) bound(X, 0, 1 << 30)  // this triggers an annoying corner case: an empty W for the locals.
#define BOUNDED(X) X 

SampleOUT test(SampleIN x) {
	struct Sample s; load_SampleIN(x, &s);
	s.a = BOUNDED(s.a * s.b);
	s.b = BOUNDED(s.a + s.b);
	s.c = BOUNDED(s.b + s.c);  
	return(save_SampleOUT(&s));
}

//#define VKEY_TIME RUN_TIME // meaning the verication keys are unknown at compile-time!
#define BOOTSTRAP 1 // meaning the program is bootstrapped, not necessarily compiling in bootstrap mode

#ifdef BOOTSTRAP
// bootstrapping involves several successive runs, with different flags.
// in particular, we need to compile the source program twice.
// - schedule (both at the same time)
// - for the inner outsourced calls, we just go through outer calls (for now a single one)
//   we do that for keygen, prove, and verify (a sanity check) following the inner schedule
// - for bootstrapping, we only consider the inner code in VERIFY mode. We need to complete the inner steps before. 
//   we do that for keygen, prove, and verify, following the outer schedule.
//
// Conversely, the MQAP flag is global, and we should only be using one level of crypto at a time.

#define GEPPETTO_NAME_BOOT "whole-verify"
#define GEPPETTO_VERIFY_BOOT_H "build/whole-verify.verify.h"
#define GEPPETTO_VERIFY_BOOT_C "build/whole-verify.verify.c"
#define NUM_BANKS_BOOT 2
#include "geppetto-boot.h"

BANK_BOOT(0,SampleIN_BOOT, struct Sample)
BANK_BOOT(1,SampleOUT_BOOT, struct Sample)

SampleOUT_BOOT test_BOOT(SampleIN_BOOT x0) {
	init();
	struct Sample s;
	load_SampleIN_BOOT(x0, &s);
	SampleIN x = save_SampleIN(&s);
	SampleOUT r = OUTSOURCE(test, x);
	load_SampleOUT(r, &s);
	return(save_SampleOUT_BOOT(&s));
}
#endif

int main() {
#ifdef BOOTSTRAP
	init_BOOT();
#else
	init();
#endif
	struct Sample s;
	s.a = 2;
	s.b = 3;
	s.c = 5;
	printf(" a=%d b=%d c=%d\n", s.a, s.b, s.c);
#ifdef BOOTSTRAP
	SampleIN_BOOT x = save_SampleIN_BOOT(&s);
	SampleOUT_BOOT r = OUTSOURCE_BOOT(0,test_BOOT, x);
	load_SampleOUT_BOOT(r, &s);
#else
	SampleIN x = save_SampleIN(&s);
	SampleOUT r = OUTSOURCE(test, x);
	load_SampleOUT(r, &s);
#endif
	printf(" a=%d b=%d c=%d\n", s.a, s.b, s.c);
}

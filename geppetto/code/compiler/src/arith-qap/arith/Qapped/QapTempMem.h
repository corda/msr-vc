#pragma once

// Let's just dynamically allocate them and leak them; that's good enough
// for small runs to validate correctness.

// TODO these probably all need corresponding frees! Yikes.

typedef struct s_Dummy { int dummy; } Dummy;
#define SOME_OFFSET					/*want to reuse these in !PINOCCHIO_MODE? fill in these offsets */
#define TEMP_FORMAL						Dummy* dummy
#define TEMP_ALLOC_FP(lval, offset)		QapFp_init(lval);
#define TEMP_ALLOC_FPL(lval, offset)	QapFpl_init(lval);
#define TEMP_ALLOC(lval, offset)		TEMP_ALLOC_FP(lval, );	/* deprecated term; no distinction between Fp/Fpl */
#define TEMP_ALLOC_FP2(lval, offset)	\
	TEMP_ALLOC_FP(lval->a0, offset); TEMP_ALLOC_FP(lval->a1, offset+SOME_OFFSET);
#define TEMP_ALLOC_FP2L(lval, offset)	\
	TEMP_ALLOC_FPL(lval->A0, offset); TEMP_ALLOC_FPL(lval->A1, offset+SOME_OFFSET);
#define TEMP_ALLOC_FP4(lval)			\
	TEMP_ALLOC_FP2(lval->a0, offset); TEMP_ALLOC_FP2(lval->a1, offset+SOME_OFFSET);
#define TEMP_ALLOC_FP6(lval, offset)	\
	TEMP_ALLOC_FP2(lval->a0, offset); \
	TEMP_ALLOC_FP2(lval->a1, offset+SOME_OFFSET); \
	TEMP_ALLOC_FP2(lval->a2, offset+SOME_OFFSET);
#define TEMP_ALLOC_FP6L(lval, offset)	\
	TEMP_ALLOC_FP2L(lval->A0, offset); \
	TEMP_ALLOC_FP2L(lval->A1, offset+SOME_OFFSET); \
	TEMP_ALLOC_FP2L(lval->A2, offset+SOME_OFFSET);
#define TEMP_FREE_FP(fp)				QapFp_free(fp);
#define TEMP_FREE_FPL(fpl)				QapFp_free(fpl);
#define TEMP_FREE_FP2(fp)				TEMP_FREE_FP(fp->a0); TEMP_FREE_FP(fp->a1);
#define TEMP_FREE_FP2L(fp)				TEMP_FREE_FP(fp->A0); TEMP_FREE_FP(fp->A1);
#define TEMP_FREE_FP4(fp)				TEMP_FREE_FP2(fp->a0); TEMP_FREE_FP2(fp->a1);
#define TEMP_FREE_FP6(fp)				TEMP_FREE_FP2(fp->a0); TEMP_FREE_FP2(fp->a1); TEMP_FREE_FP2(fp->a2);
#define TEMP_FREE_FP6L(fp)				TEMP_FREE_FP2L(fp->A0); TEMP_FREE_FP2L(fp->A1); TEMP_FREE_FP2L(fp->A2);
#define TEMP_PASS(offset)				dummy
extern Dummy* g_dummy;
#define TEMP_PASS_MEM()					g_dummy
#define TEMP_PASS_MEM2()				g_dummy
#define TEMP_PASS_MEM4()				g_dummy
#define TEMP_PASS_MEM6()				g_dummy
#define TEMP_PASS_MEM12()				g_dummy
#define TEMP_ALLOC_MEMFP4(lval, offset)	TEMP_ALLOC_FP4(lval)
#define TEMP_FREE_MEMFP4(lval)			TEMP_FREE_FP4(lval)


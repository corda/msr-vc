#pragma once

// Microsoft's bignum library expects NDEBUG to be defined (or badness occurs!)
#ifndef NDEBUG
  #define UNDEFINENDEBUG
  #define NDEBUG
#endif
#include <msbignum.h>
#include <field.h>
#include "pfc.h"
#include "modfft.h"
#ifdef UNDEFINENDEBUG
  #undef NDEBUG
  #undef UNDEFINENDEBUG
#endif

// These classes are needed to keep Enigma header files from colliding with ARITH's header files

class FFT_info {
public:
	modfft_info_t info;
};

class FieldPrimeMembers {
public:
	bigctx_t BignumCtx;
	bigctx_t* pbigctx;
	digit_t mulTempSpace[16];
	digit_t divTempSpace[100];
	field_desc_tc* msfield;
	FFT_info fft_info;
};

class Field32Members {
public:
	mp_modulus_t* mod;
	FFT_info fft_info;
};

class FieldDescription {
public:
	field_desc_tc* desc;
};

class EncodingEnigmaBNMembers {
public:
	pairing_friendly_curve_t bn;  
	bigctx_t BignumCtx;
	bigctx_t* pbigctx;
	mp_modulus_t* p;
	field_desc_t* msExponentField;

	digit_t arithTempSpace[260];

  FieldDescription* base_desc;
  mp_modulus_t* base_p;

	// Pointers to the base and twist curves
	const ecurve_t* curve[2];
};
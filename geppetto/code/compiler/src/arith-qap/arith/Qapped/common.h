#pragma once

/////////////////////////////////////////////////////////////////
//  Set the common parameters
//  Typically, these should be overridden on the command line
//  (e.g., from the makefile)

// Define TINY to use EC curves with artifically small parameters, e.g., for bootstrap debugging
//#define TINY  

// Define SELECTMUL to use explicit selects instead of branches, as required for QAP compilation
#define SELECTMUL 1


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//
//   DO NOT CHANGE CODE BELOW HERE
//
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////
//  Basic definitions and includes

#include "QapInterface.h"
static QapInterface qap_ifc;

// tweaks for testing the compiler 
#ifdef MQAP
#include "../../../../input/pinocchio.h"
#include "PrimitiveIfc.h"
void assert (int i) {
	zeroAssert(i == 0); //$ inefficient
}
void elem_assert(int i) {
	assert(i); 
}
#else
#define NOASSERT 1
#include "../../../../input/pinocchio.c"
#include "PrimitiveIfc-Arith.c"
#endif

typedef int bool;

// Pull in QapArith machinery.
#include "arith-qap-aggregator.c"

////////////////////////////////////////////////////////////////////////////////
// Instantiate macro-templated code for encodings to get Left & Right versions

#define SIDEL 10
#define SIDER 11

#define SIDE(f) f##L
#define SIDENUM SIDEL
#include "encoding-curves.c"
#undef EncodedElt
#undef SIDENUM
#undef SIDE

#define SIDE(f) f##R
#define SIDENUM SIDER
#include "encoding-curves.c"
#undef EncodedElt
#undef SIDENUM
#undef SIDE

//////////////////////////////////////////////////////////////////////////////
// Pull in pairing code

#include "encoding-pairing.c"

//////////////////////////////////////////////////////////////////////////////
// Required initialization steps

void common_init() {
  elem_static_init();
  QapInterface_init(&qap_ifc);
  encoding_static_initL();
  encoding_static_initR();
}

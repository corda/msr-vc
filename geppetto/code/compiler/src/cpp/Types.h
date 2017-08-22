#pragma once

#define WIN32_LEAN_AND_MEAN  // Prevent conflicts in Windows.h and WinSock.h
#pragma warning( disable : 4005 )	// Ignore warning about intsafe.h macro redefinitions

// Uncomment to compress keys & proofs before writing them to an archive
//#define CONDENSE_ARCHIVES
#ifdef _DEBUG
// Uncomment this line to check for memory leaks using VS.NET
//#define MEM_LEAK

// Uncomment this line to use sequential, rather than random, field elements
//#define DEBUG_RANDOM

// Perform extra checks on elliptic curve operations
//#define DEBUG_EC

//#define VERBOSE

#endif	// _DEBUG

// Hunting for memory leaks
#ifdef MEM_LEAK
#define CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
// Print file names for the locations of memory leaks
#define DEBUG_NEW new(_NORMAL_BLOCK ,__FILE__, __LINE__)
#define new DEBUG_NEW
// See https://msdn.microsoft.com/en-us/library/x98tx3cf.aspx for useful tips
#endif
// Done with code for hunting memory leaks

// Enable additional performance counters and checks
//#define PERF_DEBUG

#include "timer.h"

extern TimerList timers;

#ifdef __cplusplus
extern "C" {
#endif
extern BOOL WINAPI random_bytes(BYTE* buffer, const size_t count, PVOID pContext);
#ifdef __cplusplus
}
#endif

typedef unsigned int uint32_t;
typedef unsigned short uint16_t;

enum poly_selector_t { V = 0, W, Y, alphaV, alphaW, alphaY, beta };
static const int NUM_POLY_SELECTORS = 7;
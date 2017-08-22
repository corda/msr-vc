// Interface for F# code

#pragma warning( disable : 4005 )	// Ignore warning about intsafe.h macro redefinitions
#include "Types.h"
#include <iostream>
#include <fstream>
#include "timer.h"
#include <objbase.h>
#include <comutil.h>

#pragma comment(lib, "comsuppw.lib")
#pragma comment(lib, "kernel32.lib")

using namespace std;

#include "Encoding.h"
#include "Field.h"
#include "SparsePolynomial.h"
#include "EncodingEnigmaBN.h"
#include "FieldPrime.h"
#define MEM_LEAK
extern "C" {
#ifdef MEM_LEAK
// Memory tracking
int num_field_arrays = 0;
int num_field_elts = 0;
int num_field_ops = 0;
#endif // MEM_LEAK

#ifdef MEM_LEAK
#define field_elt_created()   { num_field_elts++; }
#define field_elt_destroyed() { num_field_elts--; }
#define field_array_created()   { num_field_arrays++; }
#define field_array_destroyed() { num_field_arrays--; }
#define field_op() { num_field_ops++; }

char* make_field_stats_str() {
  char* ret = new char[500];
  int len = 0;
  len += sprintf_s(ret + len, 500 - len, "Number of F# allocated field elements still alive: %d.  Field arrays: %d.  Number of field ops: %d\n", num_field_elts, num_field_arrays, num_field_ops);
  len += sprintf_s(ret + len, 500 - len, "Total number of allocated field elements still alive: %d.  Total number of encoded elts still alive: L = %d, R = %d\n", FieldPrime::num_live_field_elts, EncodingEnigmaBN::num_live_enc_elts[0], EncodingEnigmaBN::num_live_enc_elts[1]);
  return ret;
}

void print_field_stats() {
  char* stats = make_field_stats_str();
  printf("%s", stats);
  delete[] stats;
}

#else
#define field_elt_created()   { ; }
#define field_elt_destroyed() { ; }
#define field_array_created()   { ; }
#define field_array_destroyed() { ; }
#define field_op() { ; }
void print_field_stats() {;}
#endif // MEM_LEAK

__declspec(dllexport) BSTR field_stats() {
  char* stats = make_field_stats_str();
  BSTR ret = _com_util::ConvertStringToBSTR(stats);  
  delete[] stats;
  return ret;
}

__declspec(dllexport) void print_current_field_stats()
{
  print_field_stats();
}

__declspec(dllexport) int num_live_field_elts()
{
  return num_field_elts;
}

__declspec(dllexport) int num_live_field_arrays()
{
  return num_field_arrays;
}

__declspec(dllexport) int num_live_field_elts_total()
{
  return FieldPrime::num_live_field_elts;
}

__declspec(dllexport) int num_live_enc_elts()
{
  return EncodingEnigmaBN::num_live_enc_elts[0] + EncodingEnigmaBN::num_live_enc_elts[1];
}


__declspec(dllexport) Encoding* makeEncoding(encoding_t type)
{
  return NewEncoding(type);
}

__declspec(dllexport) Field* getField(Encoding* enc)
{	
	return enc->getSrcField();
}

__declspec(dllexport) FieldElt* getRandomFieldElt(Field* f)
{
  field_elt_created();
	return f->newRandomElt();
}

// Number of int64s in a field element
__declspec(dllexport) int eltSize(Field* f)
{
  assert(f);
  return f->eltSize();
}

__declspec(dllexport) void getFull(const FieldElt* f1, _int64 w[4], Field* f)
{
  assert(f1 && f);
  f->getFull(f1, w);
}

__declspec(dllexport) void setFull(FieldElt* f1, const _int64 w[4], Field* f)
{
  assert(f1 && f);
  f->setFull(f1, w);
}

__declspec(dllexport) void getMod(_int64 w[4], Field* f) 
{
  assert(f);
  f->getMod(w);
}


__declspec(dllexport) FieldElt* field_add(const FieldElt* f1, const FieldElt* f2, Field* f)
{
	assert(f1 && f2 && f);
  FieldElt* ret = f->newElt(); field_elt_created(); field_op();
	f->add(f1,f2,ret);
	return ret;
}

__declspec(dllexport) FieldElt* field_sub(const FieldElt* f1, const FieldElt* f2, Field* f)
{
	assert(f1 && f2 && f);
  FieldElt* ret = f->newElt(); field_elt_created(); field_op();
	f->sub(f1, f2, ret);
	return ret;
}

__declspec(dllexport) FieldElt* field_mul(const FieldElt* f1, const FieldElt* f2, Field* f)
{
	assert(f1 && f2 && f);
  FieldElt* ret = f->newElt(); field_elt_created(); field_op();
	f->mul(f1,f2,ret);
	return ret;
}

__declspec(dllexport) void field_mul_alt(const FieldElt* f1, const FieldElt* f2, FieldElt* r, Field* f)
{
	assert(f1 && f2 && f);
	f->mul(f1,f2,r);
}

__declspec(dllexport) FieldElt* field_div(const FieldElt* f1, const FieldElt* f2, Field* f)
{
	assert(f1 && f2 && f);
  FieldElt* ret = f->newElt(); field_elt_created(); field_op();
	f->div(f1, f2, ret);
	return ret;
}

__declspec(dllexport) FieldElt* field_exp(const FieldElt* f1, int ex, Field* f)
{
	assert(f1 && f);
  FieldElt* ret = f->newElt(); field_elt_created(); field_op();
	f->exp(f1, ex, ret);
	return ret;
}

__declspec(dllexport) int fieldEltEql(const FieldElt* f1, const FieldElt* f2, Field* f)
{
	assert(f1 && f2 && f);
	return f->equal(f1,f2);
}

__declspec(dllexport) int field_compare(const FieldElt* f1, const FieldElt* f2, Field* f)
{
	assert(f1 && f2 && f);
	return f->compare(f1, f2);
}

__declspec(dllexport) FieldElt* makeFieldElt(int i, Field* f)
{
	assert(f);
  FieldElt* ret = f->newElt(); field_elt_created();
	f->set(ret, i);
	return ret;
}

__declspec(dllexport) void printFieldElt(FieldElt* f1, Field * f)
{
	assert(f1 && f);
  f->print(f1);
}


__declspec(dllexport) char* fieldEltToString(FieldElt* f1, Field * f)
{
  assert(f1 && f);
  char* c_str = f->to_string(f1);

  // Allocate a string that F# can safely deallocate
  size_t len = strlen(c_str) + sizeof(char);
  char* ret = (char*)::CoTaskMemAlloc(len); 
  strcpy_s(ret, len, c_str);
  delete[] c_str;

  return ret;
}


_declspec(dllexport) int fieldEltBit(const FieldElt* f1, const int n, Field* f)
{
	assert(f1 && f);
	return f->bit(n, f1);
}

__declspec(dllexport) void deleteFieldEltArray(FieldEltArray* arr)
{
  delete arr;  field_array_destroyed();
}

__declspec(dllexport) void deleteFieldElt(FieldElt* elt, Field* f)
{
  f->delElt(elt);  field_elt_destroyed();
}

__declspec(dllexport) FieldEltArray* newFieldEltArray(Field* field, int size)
{
  assert(field);
  return field->newEltArray(size, true);  field_array_created();
}

__declspec(dllexport) void setElem(FieldEltArray* arr, FieldElt* elem, int idx, Field* f)
{
  assert(arr && elem && f);
  assert(idx >= 0);
  f->copy(elem, arr->elt(idx));
}

__declspec(dllexport) uint64_t cGetRDTSC() {
  return GetRDTSC();
}

}

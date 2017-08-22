#pragma once

#include "Encoding.h"
#include "Poly.h"
#include "Archive.h"
#include <unordered_map>	

enum bank_type_t { IO=0, NON_IO, BOTH_BANK_TYPES };

typedef unordered_map<unsigned long, EncodedElt*> KeyElts;
typedef pair<unsigned long, EncodedElt*> KeyElt;

class KeyElements {
public:
  KeyElements();
  KeyElements(Encoding* encoding, elt_t t, int max_size);
  ~KeyElements();

  void set_max_size(int max) { max_size = max; }
  int get_max_size() { return max_size; }

  elt_t get_type() { return type; }
  KeyElts* get_elts() { return &elts; }

  int size();

  // Returns the elt requested or NULL if absent
  EncodedElt* find(unsigned long index);  

  bool operator==(const KeyElements& other);
  bool operator!=(const KeyElements& other);

  // Returns the elt requested or allocates a fresh elt for this index
  EncodedElt* operator[](unsigned long index);

  void serialize(Archive* arc, bool simple = false);
  void deserialize(Encoding* encoding, Archive* arc);

  void print();
  

private:
  Encoding* encoding;
  KeyElts elts;
  elt_t type;
  int max_size;
};

class MultiBankKey {
public:
  MultiBankKey();
  MultiBankKey(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[], bank_type_t my_bank_type);
  virtual ~MultiBankKey();
  
  virtual bool bank_is_mine(int bank_id);
  EncodedEltArray** new_array(elt_t t, int bank_sizes[]);
  void del_array(EncodedEltArray** arr);

  // Optimized version for setting polys to 0
  virtual void apply_zero(int bank_id, int coeff_index, poly_selector_t which_poly) = 0;

	virtual void serialize(Archive* arc, bool simple=false);
	virtual void deserialize(Encoding* encoding, Archive* arc);
  virtual void print();  

  // Lookup this bank's ID
  int elt(int bank_id);

  // Size of this bank
  int bank_size(int bank_id);

  virtual bool equals(Encoding* encoding, MultiBankKey* other);

  Encoding* encoding;
  int num_banks_total;
  int num_banks_mine;
  int* bank_sizes;
  bank_type_t* bank_type;
  int num_roots;

  int* bank_index_map;

	// Lots of polys evaluated at s. One set for each bank.
  KeyElements** Vpolys;
  KeyElements** Wpolys;
  KeyElements** Ypolys;

protected:
  int count_my_banks();
  bool is_my_type(bank_type_t bank_type);
  virtual void print(EncodedEltArray* a);
  bank_type_t my_bank_type;
};

class EvalKey : public MultiBankKey {
public:
  EvalKey();
  EvalKey(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[]);
  ~EvalKey();
  void apply(int bank_id, int coeff_index, poly_selector_t which_poly, FieldElt* val, FieldElt* alpha_val);
  void apply(int bank_id, int coeff_index, FieldElt* beta_val);
  virtual void apply_zero(int bank_id, int coeff_index, poly_selector_t which_poly);

  virtual void serialize(Archive* arc, bool simple = false);
	virtual void deserialize(Encoding* encoding, Archive* arc);
  virtual void print();

  virtual bool equals(Encoding* encoding, EvalKey* other);

  static const uint64_t root_ratio = 2;
  unique_ptr<Poly> tInv;
  void doPrecalcs(Field* field, FieldElt* root_ratio, int num_roots);
	// Powers of s^i
	LEncodedEltArray* powers;

	KeyElements** alphaVpolys;
	KeyElements** alphaWpolys;
	KeyElements** alphaYpolys;

  KeyElements** betaVWYpolys;
};

class VerifKey : public MultiBankKey {
public:
  VerifKey();
  VerifKey(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[]);
  ~VerifKey();
  void apply(int bank_id, int coeff_index, poly_selector_t which_poly, FieldElt* val);
  virtual void apply_zero(int bank_id, int coeff_index, poly_selector_t which_poly);

  virtual void serialize(Archive* arc, bool simple = false);
	virtual void deserialize(Encoding* encoding, Archive* arc);
  virtual void print();

  virtual bool equals(Encoding* encoding, VerifKey* other);

	// KEA values
	REncodedEltArray* alphaV;
	LEncodedEltArray* alphaW;
	REncodedEltArray* alphaY;

	REncodedEltArray* gammaR;
	LEncodedEltArray* betaGammaL;
	REncodedEltArray* betaGammaR;

	REncodedElt* RtAtS;	
};

class KeyPair {
public:
  KeyPair();
  KeyPair(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[]);
  ~KeyPair();
  
  EvalKey*  ek;
  VerifKey* vk;
};
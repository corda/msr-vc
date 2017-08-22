#include "Keys.h"
#include <stdio.h>
#include <sstream>
#include <assert.h>

KeyElements::KeyElements() {
  this->encoding = NULL;
  this->type = L;
  this->max_size = 0;
}

KeyElements::KeyElements(Encoding* encoding, elt_t t, int max_size) {
  this->encoding = encoding;
  this->type = t;
  this->max_size = max_size;
}

KeyElements::~KeyElements() {
  for (KeyElts::iterator iter = elts.begin();
    iter != elts.end();
    iter++) {
    delete iter->second;
  }
}

int KeyElements::size() {
  return (int)elts.size();
}

// Returns the elt requested or NULL if absent
EncodedElt* KeyElements::find(unsigned long index) {
  KeyElts::iterator result = elts.find(index);
  if (result == elts.end()) {
    return NULL;
  } else {
    return result->second;
  }
}

bool KeyElements::operator==(const KeyElements& other) {
  bool equal = true;

  if (elts.size() != other.elts.size()) {
    printf("Size differs: %d vs %d\n", elts.size(), other.elts.size());
    equal = false;
  } else {
    for (KeyElts::const_iterator iter = elts.begin();
      iter != elts.end();
      iter++) {
      KeyElts::const_iterator result = other.elts.find(iter->first);
      if (result == other.elts.end()) {
        printf("Could not find element %d\n", iter->first);
        equal = false;
      } else if (!encoding->equal(type, iter->second, result->second)) {
        printf("Differ at elt %d\n", iter->first);
        equal = false;
      }
    }
  }
  return equal;  
}

bool KeyElements::operator!=(const KeyElements& other) {
  return !(*this == other);
}

EncodedElt* KeyElements::operator[](unsigned long index) {
  KeyElts::iterator result = elts.find(index);
  if (result == elts.end()) {
    EncodedElt* ret = encoding->new_elt(type);
    elts.insert(KeyElt(index, ret));
    encoding->zero(type, ret);
    return ret;
  } else {
    return result->second;
  }
}

void KeyElements::print() {
  for (KeyElts::const_iterator iter = elts.begin();
    iter != elts.end();
    iter++) {
    printf("%0.2d: \n", iter->first);
    encoding->print(type, iter->second);
    printf("\n");
  }
}

void KeyElements::serialize(Archive* arc, bool simple) {
  if (!simple) {
    unsigned int tmp = type;
    arc->write(tmp);
    tmp = this->get_max_size();
    arc->write(tmp);
    tmp = this->size();
    arc->write(tmp);

    // Write out the dense version, pairs of (index, EncodedElt) 
    for (KeyElts::const_iterator iter = elts.begin();
      iter != elts.end();
      iter++) {
      tmp = iter->first;
      arc->write(tmp);
      encoding->write(arc->get_archiver(), type, iter->second, simple);
    }
  } else {
    // Write out the sparse version
    EncodedElt* zero = encoding->new_elt(type);
    encoding->zero(type, zero);
    for (int i = 0; i < max_size; i++) {
      EncodedElt* elt = this->find(i);
      if (elt == NULL) {
        encoding->write(arc->get_archiver(), type, zero, simple);
      } else {
        encoding->write(arc->get_archiver(), type, elt, simple);
      }
    }
  }
}

void KeyElements::deserialize(Encoding* encoding, Archive* arc) {
  int max = 0;
  elt_t t = L;
  int size = 0;

  arc->read(t);
  arc->read(max);
  arc->read(size); 
  
  this->type = t;
  this->max_size = max;
  this->encoding = encoding;
  
  for (int i = 0; i < size; i++) {
    int index = 0;
    arc->read(index);
    encoding->read(arc->get_archiver(), t, (*this)[index]);
  }
}

MultiBankKey::MultiBankKey() {
  encoding = NULL;
  num_roots = -1;
  num_banks_total = -1;
  num_banks_mine = -1;
  bank_sizes = NULL;
  bank_type = NULL;
  bank_index_map = NULL;

  Vpolys = NULL;
  Wpolys = NULL;
  Ypolys = NULL;
}

MultiBankKey::MultiBankKey(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_types[], bank_type_t my_bank_type) {
  assert(encoding);
  assert(num_banks >= 1);

  this->encoding = encoding;
  this->num_roots = num_roots;
  this->num_banks_total = num_banks;
  this->my_bank_type = my_bank_type;
  this->bank_sizes = new int[num_banks];
  this->bank_type = new bank_type_t[num_banks];

  for (int i = 0; i < num_banks; i++) {
    this->bank_sizes[i] = bank_sizes[i];
    this->bank_type[i] = bank_types[i];
  }

  this->num_banks_mine = count_my_banks();

  // Set up the index map from the global bank index to my local index
  bank_index_map = new int[num_banks];
  int my_bank_index = 0;
  for (int i = 0; i < num_banks; i++) {
    bank_index_map[i] = is_my_type(bank_type[i]) ? my_bank_index++ : -1;
  }

  // Initialize all of the containers for our key elements
  this->Vpolys = new KeyElements*[num_banks_mine];
  this->Wpolys = new KeyElements*[num_banks_mine];
  this->Ypolys = new KeyElements*[num_banks_mine];

  my_bank_index = 0;
  for (int i = 0; i < num_banks; i++) {
    if (is_my_type(bank_type[i])) {
      this->Vpolys[my_bank_index] = new KeyElements(encoding, L, bank_sizes[i]);
      this->Wpolys[my_bank_index] = new KeyElements(encoding, R, bank_sizes[i]);
      this->Ypolys[my_bank_index] = new KeyElements(encoding, L, bank_sizes[i]);
      my_bank_index++;
    }
  }
}

MultiBankKey::~MultiBankKey() {
  for (int i = 0; i < num_banks_mine; i++) {
    delete this->Vpolys[i];
    delete this->Wpolys[i];
    delete this->Ypolys[i];
  }

  delete[] this->Vpolys;
  delete[] this->Wpolys;
  delete[] this->Ypolys;

  delete[] bank_index_map;
  delete[] bank_sizes;
  delete[] bank_type;
}


bool MultiBankKey::bank_is_mine(int bank_id) {
  assert(0 <= bank_id && bank_id < num_banks_total);
  return bank_index_map[bank_id] != -1;
}

EncodedEltArray** MultiBankKey::new_array(elt_t t, int bank_sizes[]) {
  EncodedEltArray** ret = NULL;

  if (t == L) {
    ret = new LEncodedEltArray*[num_banks_mine];
  } else {
    assert(t == R);
    ret = new REncodedEltArray*[num_banks_mine];
  }

  int my_bank_index = 0;
  for (int i = 0; i < num_banks_total; i++) {
    if (is_my_type(bank_type[i])) {
      ret[my_bank_index++] = encoding->new_elt_array(t, bank_sizes[i], true);
    }
  }

  return ret;
}

void MultiBankKey::del_array(EncodedEltArray** arr) {
  for (int i = 0; i < num_banks_mine; i++) {
    delete arr[i];
  }
  delete [] arr;
}

int MultiBankKey::elt(int bank_id) {
  assert(bank_index_map[bank_id] != -1);
  assert(0 <= bank_id && bank_id < num_banks_total);
  return bank_index_map[bank_id];
}

int MultiBankKey::bank_size(int bank_id) {
  return bank_sizes[bank_id];
}

bool MultiBankKey::is_my_type(bank_type_t t) {
  return t == my_bank_type || t == BOTH_BANK_TYPES;
}

int MultiBankKey::count_my_banks() {
  int sum = 0;
  for (int i = 0; i < num_banks_total; i++) {
    if (is_my_type(bank_type[i])) {
      sum++;
    }
  }
  return sum;
}

void MultiBankKey::serialize(Archive* arc, bool simple) {
  if (!simple) {    
    arc->write(num_roots);    
    arc->write(num_banks_total);
    arc->write(num_banks_mine);
    arc->write(bank_sizes, num_banks_total);
    arc->write(bank_type, num_banks_total);
    arc->write(bank_index_map, num_banks_total);
    arc->write(my_bank_type);
  }

  // Write all the Vs, then Ws, then Ys to stay in sync with the C code (encoding-qap.c)
  for (int i = 0; i < num_banks_mine; i++) {
    Vpolys[i]->serialize(arc, simple);
  }
  for (int i = 0; i < num_banks_mine; i++) {
    Wpolys[i]->serialize(arc, simple);
  }
  for (int i = 0; i < num_banks_mine; i++) {
    Ypolys[i]->serialize(arc, simple);
  }
}

void MultiBankKey::deserialize(Encoding* encoding, Archive* arc) {
  this->encoding = encoding;
  arc->read(num_roots);
  arc->read(num_banks_total);
  arc->read(num_banks_mine);

  bank_sizes = new int[num_banks_total];
  arc->read(bank_sizes, num_banks_total);
  bank_type = new bank_type_t[num_banks_total];
  arc->read(bank_type, num_banks_total);

  bank_index_map = new int[num_banks_total];
  arc->read(bank_index_map, num_banks_total);

  arc->read(my_bank_type);

  Vpolys = new KeyElements*[num_banks_mine];
  Wpolys = new KeyElements*[num_banks_mine];
  Ypolys = new KeyElements*[num_banks_mine];

  for (int i = 0; i < num_banks_mine; i++) {
    Vpolys[i] = new KeyElements();
    Vpolys[i]->deserialize(encoding, arc);
  }
  for (int i = 0; i < num_banks_mine; i++) {
    Wpolys[i] = new KeyElements();
    Wpolys[i]->deserialize(encoding, arc);
  }
  for (int i = 0; i < num_banks_mine; i++) {
    Ypolys[i] = new KeyElements();
    Ypolys[i]->deserialize(encoding, arc);
  }
}

void MultiBankKey::print(EncodedEltArray* a) {
  for (int i = 0; i < a->length(); i++) {
    printf("%0.2d: \n", i);
    encoding->print(a->t(), a->elt(i));
    printf("\n");
  }
}

void MultiBankKey::print() {
  int my_bank_index = 0;
  for (int i = 0; i < num_banks_total; i++) {
    if (bank_is_mine(i)) {
      printf("Bank %d\n", i);
      printf("\nE(V(s)): \n");  Vpolys[my_bank_index]->print();
      printf("\nE(W(s)): \n");  Wpolys[my_bank_index]->print();
      printf("\nE(Y(s)): \n");  Ypolys[my_bank_index]->print();
      my_bank_index++;
    }
  }
}

bool compare_array(Encoding* encoding, EncodedEltArray* a, EncodedEltArray* b, string name) {
  bool equal = true;

  if (a->length() != b->length()) {
      printf("%s->length() differs %d vs %d\n", name.c_str(), a->length(), b->length());       
      equal = false;
  } else {
    for (int p = 0; p < a->length(); p++) {
		  if (!encoding->equal(a->t(), a->elt(p), b->elt(p))) { 
			  printf("%s differs at elt %d\n", name.c_str(), p);
			  equal = false;
		  }
    }
  }
  return equal;
}

bool compare_multi_array(Encoding* encoding, EncodedEltArray** a, EncodedEltArray** b, int len, string name) {
  bool equal = true;
  for (int i = 0; i < len; i++) {
    string name_with_index = name;
    name_with_index.append(static_cast<ostringstream*>( &(ostringstream() << i) )->str());
    equal == equal && compare_array(encoding, a[i], b[i], name_with_index);
	}
  return equal;
}

bool compare_KeyElts_array(KeyElements** a, KeyElements** b, int len) {
  bool equal = true;
  for (int i = 0; i < len; i++) {
    equal == equal && (*a[i] == *b[i]);  
  }
  return equal;
}

bool MultiBankKey::equals(Encoding* encoding, MultiBankKey* other) {
	bool equal = true;

  if (num_roots != other->num_roots) {
    printf("num_roots differ %d %d\n", num_roots, other->num_roots);
    equal = false;
  }

	if (num_banks_total != other->num_banks_total) {
		printf("num_banks_total differ %d %d\n", num_banks_total, other->num_banks_total);
		equal = false;
	}
	if (num_banks_mine != other->num_banks_mine) {
		printf("numPolys differ %d %d\n", num_banks_mine, other->num_banks_mine);
		equal = false;
	}

	for (int i = 0; i < num_banks_total; i++) {
		if (bank_type[i] != other->bank_type[i]) { 
			printf("bank_type[%d] differs\n", i);
			equal = false;
		}
	}
  for (int i = 0; i < num_banks_total; i++) {
		if (bank_index_map[i] != other->bank_index_map[i]) { 
			printf("bank_index_map[%d] differs\n", i);
			equal = false;
		}
	}

  equal = equal && compare_KeyElts_array(Vpolys, other->Vpolys, num_banks_mine);
  equal = equal && compare_KeyElts_array(Wpolys, other->Wpolys, num_banks_mine);
  equal = equal && compare_KeyElts_array(Ypolys, other->Ypolys, num_banks_mine);
	
	return equal;
}

/////////////////// New Eval Key //////////////////////////////////////

EvalKey::EvalKey() : MultiBankKey() {
  this->powers = NULL;

  this->alphaVpolys  = NULL;
  this->alphaWpolys  = NULL;
  this->alphaYpolys  = NULL;
  this->betaVWYpolys = NULL;
}

EvalKey::EvalKey(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[])
  : MultiBankKey(encoding, num_roots, num_banks, bank_sizes, bank_type, NON_IO) {
  
  powers = encoding->new_elt_array(L, num_roots + 1);

  this->alphaVpolys  = new KeyElements*[num_banks_mine];
  this->alphaWpolys  = new KeyElements*[num_banks_mine];
  this->alphaYpolys  = new KeyElements*[num_banks_mine];
  this->betaVWYpolys = new KeyElements*[num_banks_mine];

  int my_bank_index = 0;
  for (int i = 0; i < num_banks; i++) {
    if (is_my_type(bank_type[i])) {
      this->alphaVpolys [my_bank_index] = new KeyElements(encoding, L, bank_sizes[i]);
      this->alphaWpolys [my_bank_index] = new KeyElements(encoding, L, bank_sizes[i]);
      this->alphaYpolys [my_bank_index] = new KeyElements(encoding, L, bank_sizes[i]);
      this->betaVWYpolys[my_bank_index] = new KeyElements(encoding, L, bank_sizes[i]);
      my_bank_index++;
    }    
  }
}

void EvalKey::serialize(Archive* arc, bool simple) {
  MultiBankKey::serialize(arc, simple);

  arc->write_elt_array(encoding, powers, simple);

  for (int i = 0; i < num_banks_mine; i++) {
    alphaVpolys[i]->serialize(arc, simple);
    alphaWpolys[i]->serialize(arc, simple);
    alphaYpolys[i]->serialize(arc, simple);
                  
    betaVWYpolys[i]->serialize(arc, simple);
  }
}

void EvalKey::deserialize(Encoding* encoding, Archive* arc) {
  MultiBankKey::deserialize(encoding, arc);

  arc->read_elt_array(encoding, &powers);

  alphaVpolys  = new KeyElements*[num_banks_mine];
  alphaWpolys  = new KeyElements*[num_banks_mine];
  alphaYpolys  = new KeyElements*[num_banks_mine];
  betaVWYpolys = new KeyElements*[num_banks_mine];

  for (int i = 0; i < num_banks_mine; i++) {
    alphaVpolys[i]  = new KeyElements();
    alphaWpolys[i]  = new KeyElements();
    alphaYpolys[i]  = new KeyElements();
    betaVWYpolys[i] = new KeyElements();

    alphaVpolys[i]->deserialize(encoding, arc);
    alphaWpolys[i]->deserialize(encoding, arc);
    alphaYpolys[i]->deserialize(encoding, arc);
    betaVWYpolys[i]->deserialize(encoding, arc);
  }
  // Regenerate tInv coefficients
  Field* field = encoding->getSrcField();
  FieldElt* ratio = field->newElt();
  field->set(ratio, root_ratio);
  doPrecalcs(field,ratio,this->num_roots);
  field->delElt(ratio);
}

bool EvalKey::equals(Encoding* encoding, EvalKey* other) {
  bool equal = MultiBankKey::equals(encoding, other);

  equal = equal && compare_array(encoding, powers, other->powers, "powers");
  equal = equal && compare_KeyElts_array(alphaVpolys, other->alphaVpolys, num_banks_mine);
  equal = equal && compare_KeyElts_array(alphaWpolys, other->alphaWpolys, num_banks_mine);
  equal = equal && compare_KeyElts_array(alphaYpolys, other->alphaYpolys, num_banks_mine);
  equal = equal && compare_KeyElts_array(betaVWYpolys, other->betaVWYpolys, num_banks_mine);
  equal = equal && tInv->equals(*other->tInv);

  return equal;
}

// TODO: This ignores the optimization where the verifier applies W's IO to alphaW to minimize work on the twist curve :(
void EvalKey::apply(int bank_id, int coeff_index, poly_selector_t which_poly, FieldElt* val, FieldElt* alpha_val) {
  assert(bank_is_mine(bank_id));
  int my_bank_index = bank_index_map[bank_id];
  assert(my_bank_index != -1);

  switch (which_poly) {
  case V:
    encoding->encode(L,       val,      (*Vpolys[my_bank_index])[coeff_index]);
    encoding->encode(L, alpha_val, (*alphaVpolys[my_bank_index])[coeff_index]);
    break;
  case W:
    encoding->encode(R,       val,      (*Wpolys[my_bank_index])[coeff_index]);
    encoding->encode(L, alpha_val, (*alphaWpolys[my_bank_index])[coeff_index]);
    break;
  case Y:
    encoding->encode(L,       val,      (*Ypolys[my_bank_index])[coeff_index]);
    encoding->encode(L, alpha_val, (*alphaYpolys[my_bank_index])[coeff_index]);
    break;
  } 
}

void EvalKey::apply(int bank_id, int coeff_index, FieldElt* beta_val) {
  assert(bank_is_mine(bank_id));
  int my_index = bank_index_map[bank_id];
  assert(my_index != -1);

  encoding->encode(L, beta_val, (*betaVWYpolys[my_index])[coeff_index]);
}

void EvalKey::apply_zero(int bank_id, int coeff_index, poly_selector_t which_poly) {
  // We can safely ignore this, now that we've swapped to densely encoded keys
 }

EvalKey::~EvalKey() {
  delete powers;

  for (int b = 0; b < num_banks_mine; b++) {
    delete this->alphaVpolys[b];
    delete this->alphaWpolys[b];
    delete this->alphaYpolys[b];
    delete this->betaVWYpolys[b];
  }
}

void EvalKey::print() {
  MultiBankKey::print();
  int my_bank_index = 0;
  for (int i = 0; i < num_banks_total; i++) {
    if (bank_is_mine(i)) {
      printf("Bank %d", i);
      printf("\nE(alphaV(s)): \n");  alphaVpolys[my_bank_index]->print();
      printf("\nE(alphaW(s)): \n");  alphaWpolys[my_bank_index]->print();
      printf("\nE(alphaY(s)): \n");  alphaYpolys[my_bank_index]->print();
      printf("\nE(beta(s)): \n");    betaVWYpolys[my_bank_index]->print();
      my_bank_index++;
    }
  }
  printf("\n");
}

void EvalKey::doPrecalcs(Field* field, FieldElt* ratio, int num_roots)
{
	// Calculate inverse of t(s)
	FieldElt* one = field->newElt(); field->one(one);
	unique_ptr<Poly> T = PolyTree::polyTreeRootGeom(field, one, ratio, num_roots);
  int len = T->coefficients->length();
	tInv.reset(new Poly(field));
	Poly::reverse(*T); // Why am I doing this?
	Poly::invert(*T, num_roots - 1, *tInv);
  field->delElt(one);
}


/////////////////// New Verif Key //////////////////////////////////////

VerifKey::VerifKey() : MultiBankKey() {
	alphaV = NULL;
  alphaW = NULL;
	alphaY = NULL;

	gammaR     = NULL;
	betaGammaL = NULL;
	betaGammaR = NULL;
		
  RtAtS = NULL;	
}

VerifKey::VerifKey(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[])
  : MultiBankKey(encoding, num_roots, num_banks, bank_sizes, bank_type, IO) {

	alphaV = encoding->new_elt_array(R, num_banks);
  alphaW = encoding->new_elt_array(L, num_banks);
	alphaY = encoding->new_elt_array(R, num_banks);

	gammaR     = encoding->new_elt_array(R, num_banks);
	betaGammaL = encoding->new_elt_array(L, num_banks);
	betaGammaR = encoding->new_elt_array(R, num_banks);
		
  RtAtS = encoding->new_elt(R);	
}

VerifKey::~VerifKey() {
  delete alphaV;
  delete alphaW;
  delete alphaY;
  
  delete gammaR;
  delete betaGammaL;
  delete betaGammaR;

  encoding->del_elt(R, RtAtS);	
}

void VerifKey::apply(int bank_id, int coeff_index, poly_selector_t which_poly, FieldElt* val) {
  assert(bank_is_mine(bank_id));
  int my_index = bank_index_map[bank_id];
  assert(my_index != -1);

  switch (which_poly) {
  case V:
    encoding->encode(L, val, (*Vpolys[my_index])[coeff_index]);    
    break;
  case W:
    encoding->encode(R, val, (*Wpolys[my_index])[coeff_index]);
    break;
  case Y:
    encoding->encode(L, val, (*Ypolys[my_index])[coeff_index]);
    break;
  }  
}

void VerifKey::apply_zero(int bank_id, int coeff_index, poly_selector_t which_poly) {
  // We can safely ignore this, now that we've swapped to densely encoded keys  
}

void VerifKey::serialize(Archive* arc, bool simple) {
  MultiBankKey::serialize(arc, simple);

  arc->write_elt_array(encoding, alphaV, simple);
  arc->write_elt_array(encoding, alphaW, simple);
  arc->write_elt_array(encoding, alphaY, simple);

  arc->write_elt_array(encoding, gammaR, simple);
  arc->write_elt_array(encoding, betaGammaL, simple);
  arc->write_elt_array(encoding, betaGammaR, simple);

	encoding->write(arc->get_archiver(), R, RtAtS, simple);
}

void VerifKey::deserialize(Encoding* encoding, Archive* arc) {
  MultiBankKey::deserialize(encoding, arc);

  arc->read_elt_array(encoding, &alphaV);
  arc->read_elt_array(encoding, &alphaW);
  arc->read_elt_array(encoding, &alphaY);

  arc->read_elt_array(encoding, &gammaR);
  arc->read_elt_array(encoding, &betaGammaL);
  arc->read_elt_array(encoding, &betaGammaR);

  RtAtS = encoding->new_elt(R);
	encoding->read(arc->get_archiver(), R, RtAtS);
}


void VerifKey::print() {
  MultiBankKey::print();
  
  printf("\nE(alphaV):\n");       MultiBankKey::print(alphaV);
  printf("\nE(alphaW):\n");       MultiBankKey::print(alphaW);
  printf("\nE(alphaY):\n");       MultiBankKey::print(alphaY);
  printf("\nE(gamma): \n");       MultiBankKey::print(gammaR);
  printf("\nE(betaGammaL):\n");   MultiBankKey::print(betaGammaL);
  printf("\nE(betaGammaR):\n");   MultiBankKey::print(betaGammaR);
  printf("\nE(RtAtS):\n");        encoding->print(R, RtAtS);
  printf("\n");
}


bool VerifKey::equals(Encoding* encoding, VerifKey* other) {
  bool equal = MultiBankKey::equals(encoding, other);

  equal = equal && compare_array(encoding, alphaV, other->alphaV, "alphaV");
  equal = equal && compare_array(encoding, alphaW, other->alphaW, "alphaW");
  equal = equal && compare_array(encoding, alphaY, other->alphaY, "alphaY");

  equal = equal && compare_array(encoding, gammaR, other->gammaR, "gammaR");
  equal = equal && compare_array(encoding, betaGammaL, other->betaGammaL, "betaGammaL");
  equal = equal && compare_array(encoding, betaGammaR, other->betaGammaR, "betaGammaR");

  if (!encoding->equal(R, RtAtS, other->RtAtS)) {
    printf("RtAtS differs\n");
    equal = false;
  }

  return equal;
}

/////////////////// Container for EK and VK //////////////////////////////////////
KeyPair::KeyPair() {
  this->ek = new EvalKey();
  this->vk = new VerifKey();
}

KeyPair::KeyPair(Encoding* encoding, int num_roots, int num_banks, int bank_sizes[], bank_type_t bank_type[]) {
  this->ek = new EvalKey (encoding, num_roots, num_banks, bank_sizes, bank_type);
  this->vk = new VerifKey(encoding, num_roots, num_banks, bank_sizes, bank_type);
}
 
KeyPair::~KeyPair() {
  delete ek;
  delete vk;
}
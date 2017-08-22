#include "Archive.h"

void Archive::write_elt_array(Encoding* encoding, EncodedEltArray* arr, bool simple) {
  if (!simple) {
    int tmp = arr->length();
    write(tmp);
    tmp = arr->t();
    write(tmp);
    bool pre = arr->pre();
    write(pre);
  }

  for (int i = 0; i < arr->length(); i++) {
		encoding->write(arc, arr->t(), arr->elt(i), simple);		
  }
}


void Archive::read_elt_array(Encoding* encoding, EncodedEltArray** arr) {
  int len = 0;
  elt_t t = L;
  bool pre = false;

  read(len);
  read(t);
  read(pre);
  *arr = encoding->new_elt_array(t, len, pre);

  for (int i = 0; i < len; i++) {
		encoding->read(arc, t, (*arr)->elt(i));
  }
}

#define WIN32_LEAN_AND_MEAN  // Prevent conflicts in Window.h and WinSock.h
#include "Proof.h"
#include <string.h>


Proof::Proof() {
	H = NULL;
}

Proof::Proof(Encoding* encoding)  { 	
  this->encoding = encoding;

	H = encoding->new_elt(L);  
}


Proof::~Proof() { 
	encoding->del_elt(L, H);
}

void Proof::serialize(Archive* arc, bool simple) {
	encoding->write(arc->get_archiver(), L, H, simple);
}

void Proof::deserialize(Archive* arc) {
	assert(H);
	encoding->read(arc->get_archiver(), L, H);
}

#define CHECK(t, entry) \
if (!encoding->equal(t, entry, other->entry)) { printf("%s differs\n", #entry); equal = false; }

bool Proof::equal(Encoding* encoding, Proof* other) {
	bool equal = true;
	CHECK(L, H);

	return equal;
}

void Proof::print() {
  this->encoding->print(L, this->H);
}

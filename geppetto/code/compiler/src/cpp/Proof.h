#pragma once

#include "Encoding.h"
#include "Archive.h"

class Proof {
public:
  Encoding* encoding;

	LEncodedElt* H;

	Proof();
	Proof(Encoding* encoding);
	~Proof();

	bool equal(Encoding* encoding, Proof* other);
  void serialize(Archive* arc, bool simple = false);
	void deserialize(Archive* arc);
  void print();
};


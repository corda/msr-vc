#pragma once

#include "Types.h"
#include "Field.h"
#include "Poly.h"
#include <memory>

class Field;
class Poly;

class PolyTree {
public:
	PolyTree(Field* field, FieldEltArray* x, int len);
	~PolyTree();

  void deleteAllButRoot();

	Poly*** polys;
	int height;
	int numBasePolys;
	static unique_ptr<Poly> polyTreeRootGeom(Field* field, FieldElt* x0, FieldElt* xr, int len);
};

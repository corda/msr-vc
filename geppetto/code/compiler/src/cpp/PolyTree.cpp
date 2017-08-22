#include "PolyTree.h"
#include "Config.h"
#include <stdio.h>
#include <memory>

#include "Util.h"


PolyTree::PolyTree(Field* field, FieldEltArray* x, int len) {
	// Precompute P_{ij}, where P_{0j} = (x - x_j) and P_{i+1,j} = P_{i,2j} * P_{i,2j+1}
	// All of the these polynomials are monic, so we omit the x^degree coefficient (since it's always 1)
	int numPolys = len;
	
	this->height = my_log2(numPolys) + 1;
	this->polys = new Poly**[this->height];
	this->numBasePolys = len;

	FieldElt* nada = field->newElt();
	field->zero(nada);	

	this->polys[0] = new Poly*[numPolys];
	for (int i = 0; i < numPolys; i++) {
		this->polys[0][i] = new Poly(field, 1);
		field->one(this->polys[0][i]->coefficients->elt(1));
		field->sub(nada, x->elt(i), this->polys[0][i]->coefficients->elt(0));	 // Poly_{0,i} = (x - x_i)
		//printf("\t*");
	}

	// Create the tree structure
	for (int i = 1; i < this->height; i++) {
		numPolys /= 2;
		this->polys[i] = new Poly*[numPolys];		
	}

	// Now fill it in.  Each interior node is the product of its children
	numPolys = len;
	for (int i = 0; i < this->height - 1; i++) {	// For each level of the tree
		//printf("\n");
		for (int j = 0; j < numPolys / 2; j++) {	// For each node in the next level up
			this->polys[i+1][j] = new Poly(field);
			Poly::mul(*this->polys[i][2*j], *this->polys[i][2*j + 1], *this->polys[i+1][j]);
			//printf("\t*");
		}

		if (numPolys % 2 == 1) { // We have an extra poly to deal with
			// Multiply it into the last poly
			Poly::mul(*this->polys[i+1][numPolys/2-1], *this->polys[i][numPolys-1], *this->polys[i+1][numPolys/2-1]);
		}

		numPolys /= 2;
	}
  field->delElt(nada);
	//printf("\n");
}

void PolyTree::deleteAllButRoot() {
	int numPolys = numBasePolys;
	for (int i = 0; i < height-1; i++) {			
		for (int j = 0; j < numPolys; j++) {
			delete polys[i][j];
		}
		delete [] polys[i];
    polys[i] = NULL;
		numPolys /= 2;
	}
}

PolyTree::~PolyTree() {
	int numPolys = numBasePolys;
	for (int i = 0; i < height; i++) {			
    if (polys[i] != NULL) {
      for (int j = 0; j < numPolys; j++) {
        delete polys[i][j];
      }
      delete [] polys[i];
    }
    numPolys /= 2;
	}
	delete [] polys;
}
Timer* total = timers.newTimer("Tree total", "PolyOps");
Timer* odd = timers.newTimer("odd", "Tree total");
Timer* copyTimer = timers.newTimer("copy", "Tree total");
Timer* mul = timers.newTimer("mul", "Tree total");
Timer* recur2Timer = timers.newTimer("recur2", "Tree total");

unique_ptr<Poly> PolyTree::polyTreeRootGeom(Field* field, FieldElt* x0, FieldElt* xr, int len)
{
	if (len == 0)
	{	// return f(x) = 1, the constant function
		unique_ptr<Poly> one(new Poly(field, 0));
		INFO("Lengths: %d\n", one->coefficients->length());
		field->one(one->coefficients->elt(0));
		return one;
	} 

	unique_ptr<Poly> recur1 = polyTreeRootGeom(field, x0, xr, len / 2);
	total->start();
	unique_ptr<Poly> recur2(new Poly(field, len / 2));
	//printf("Lengths: %d %d %d\n", len/2, recur1->coefficients->length(), recur2->coefficients->length());
	copyTimer->start();
	field->copyArray(recur1->coefficients, recur2->coefficients);
	copyTimer->stop();

	recur2Timer->start();
	// recur2(x) = xr^(n^2) * recur1(x/(xr^n)) where n = len/2
  FieldElt* rn  = field->newElt();
  FieldElt* rn2 = field->newElt();	
	field->exp(xr, len / 2, rn);        // Doing this in two steps so that 32-bit
	field->exp(rn, len / 2, rn2); // integer parameter does not overflow
	for (int i = 0; i <= len / 2; ++i)
	{
		field->mul(rn2, recur1->coefficients->elt(i), recur2->coefficients->elt(i));
		field->div(rn2, rn, rn2);
	}
  field->delElt(rn2);
	recur2Timer->stop();
	unique_ptr<Poly> res(new Poly(field, len));
	mul->start();
	Poly::mul(*recur1, *recur2, *res);
	mul->stop();

	if (len % 2 == 0) {
		total->stop(); 
    field->delElt(rn);
    return res;
	}

	odd->start();
	FieldElt* zero = field->newElt();
	field->zero(zero);

	// Else multiply res by (x - rn), where rn = x0 * xr^(len-1))
	unique_ptr<Poly> shifted(new Poly(field, len));
	field->exp(xr, len - 1, rn);
	field->mul(rn, x0, rn);
	// shifted[0] = -rn*res[0]
	field->mul(rn, res->coefficients->elt(0), shifted->coefficients->elt(0));
	field->sub(zero, shifted->coefficients->elt(0), shifted->coefficients->elt(0));
  field->delElt(zero);
	// shifted[len] = res[len-1]
	field->copy(res->coefficients->elt(len - 1), shifted->coefficients->elt(len));
	for (int i = 1; i <= len - 1; ++i)
	{
		// shifted[i] = res[i-1] - rn*res[i]
		FieldElt* dest = shifted->coefficients->elt(i);
		field->mul(rn, res->coefficients->elt(i), dest);
		field->sub(res->coefficients->elt(i - 1), dest, dest);
	}
  field->delElt(rn);
	odd->stop();
	total->stop();
	return shifted;
}
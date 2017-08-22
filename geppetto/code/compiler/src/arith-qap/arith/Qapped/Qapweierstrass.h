static void VARIANT(weierstrass_affdbl_template)(POINT c, POINT a, MEM mem) {
	//lambda: = x1 ^ 2;
	SQR(mem->t0, a->X);
	//c: = lambda + lambda;
	MUL2(mem->t1, mem->t0);
	//lambda: = c + lambda;
	ADD(mem->t0, mem->t0, mem->t1);
	//c: = 2 * y1;
	MUL2(mem->t1, a->Y);
	//c: = 1 / c;
	INV(mem->t1, mem->t1);
	//lambda: = lambda*c;
	MUL(mem->t0, mem->t0, mem->t1);
	//c: = lambda ^ 2;
	SQR(mem->t1, mem->t0);
	//c: = c - x1;
	SUB(mem->t1, mem->t1, a->X);
	//x3: = c - x1;
	SUB(mem->t2, mem->t1, a->X);
	//y3: = lambda*x3;
	MUL(mem->t3, mem->t0, mem->t2);
	//lambda: = lambda*x1;
	MUL(mem->t0, mem->t0, a->X);
	//y3: = y3 + y1;
	ADD(c->Y, mem->t3, a->Y);
	//y3: = lambda - y3;
    SUB(c->Y, mem->t0, c->Y);
	SET(c->X, mem->t2);
}

static void VARIANT(weierstrass_affadd_template)(POINT c, POINT a, POINT b, MEM mem) {
	//lambda: = y2 - y1;
	SUB(mem->t0, b->Y, a->Y);
	//c: = x2 - x1;
	SUB(mem->t1, b->X, a->X);
	//c: = 1 / c;
	INV(mem->t1, mem->t1);
	//lambda: = lambda*c;
	MUL(mem->t0, mem->t0, mem->t1);
	//c: = lambda ^ 2;
	SQR(mem->t1, mem->t0);
	//x3: = c - x1;
	SUB(mem->t2, mem->t1, a->X);
	//x3: = x3 - x2;
	SUB(mem->t2, mem->t2, b->X);
	//y3: = lambda*x3;
	MUL(mem->t3, mem->t0, mem->t2);
	//lambda: = x1*lambda;
	MUL(mem->t0, mem->t0, a->X);
	//y3: = y3 + y1;
	ADD(c->Y, mem->t3, a->Y);
	//y3: = lambda - y3;
	SUB(c->Y, mem->t0, c->Y);
	SET(c->X, mem->t2);
}

static int VARIANT(weierstrass_oncurve_aff_template) (POINT x, MEM mem, CURVE curve) {
  SQR (mem->t0, x->Y); // left = t0
  SQR (mem->t1, x->X); // right = t1
  MUL (mem->t1, mem->t1, x->X);
  ADD (mem->t1, mem->t1, curve->b);
  return (CMPEQ (mem->t0, mem->t1));
}

static int VARIANT(weierstrass_equal_aff_template) (POINT x, POINT y, CURVE curve) {
  return ((CMPEQ (x->X, y->X)) && (CMPEQ (x->Y, y->Y) ));
}




struct elm {
  long repr[4];
};

typedef float /*struct elm* */ field;

field zero;
field  add(field x, field y);
field mlt(field x, field  y);

/* etc, with arith as an implementation, and custom support to compile calls into primitive operations. */ 

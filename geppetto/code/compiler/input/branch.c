#include "pinocchio.h"
#ifndef MQAP
#include <stdio.h>
#endif
// small examples for compiling branches. 

struct Matrix { int M[100*100]; };
typedef struct Matrix *t;

struct tB; 
typedef struct tB *t_bank;

struct bank_In { int b; int x; }; 
struct bank_Out { int r; }; 

void g(struct bank_In* in, struct bank_Out* out)
{
  out->r  = in->x && in->b;
  // more efficient, but hard to get from clang:
  // out->r = (in->x * in->b) != 0; 
}
// interestingly compiled to a branch, then x is converted, then x*b is converted:
// 
// 0:  (x2)(x4) = (1 - x5)
// 1:  (x2)(x5) = (0)
// 2:  (1 - x5)(x1) = (x6)
// 3:  (x6)(x7) = (1 - x8)
// 4:  (x6)(x8) = (0)
// 5:  (1 - x8)(1) = (x3)


void f(int* s, struct bank_In* in, struct bank_Out* out)
{
  int x = /* in->x; */ bound(in->x,0,100);
  int b = /* in->b; */ bound(in->b,0,1); // ensuring b is a boolean (TODO: check why we don't use it)
  int t = 2;

  // a static branch
  if (*s) x = x * x;

  // a dynamic branch
  if (b) t = x * x;
  else t = 3 * x + 5;

  printf("f(%4d,%4d,%4d) = %4d\n",*s,b,x,t);
  out->r = t;
}

#ifndef MQAP
int main() 
{
int s;
struct bank_In in;
struct bank_Out out;
in.x = 4;
for(int i = 0; i < 4; i++) {
in.b = i % 2; s = i / 2;
f(&s,&in,&out);
//printf("f(%4d,%4d,%4d) = %4d\n",s,in.b,in.x,out.r);
}
}
#endif

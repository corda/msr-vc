#pragma once

#if PARAM==0 
#define SIZE 2
#elif PARAM==1
#define SIZE 3
#elif PARAM==2
#define SIZE 10
#elif PARAM==3
#define SIZE 20
#elif PARAM==4
#define SIZE 50	
#elif PARAM==5
#define SIZE 110
#else
#error unknown PARAM
#endif

#ifndef TRUNCATE
#define TRUNCATE 0
#endif

//#define MAT(m, r, c)	(m)[r][c]

struct M {
  int x[SIZE][SIZE];
};

struct bank_A {
  struct M A;
};

struct bank_R {
  struct M R;
};

void outsource(struct bank_A *a, struct M *b,struct bank_R *r);

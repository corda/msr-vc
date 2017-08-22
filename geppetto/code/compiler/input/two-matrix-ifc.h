#pragma once

#if PARAM==0 
#define SIZE 2
#elif PARAM==1
#define SIZE 5 
#elif PARAM==2
#define SIZE 10
#elif PARAM==3
#define SIZE 20
#elif PARAM==4
#define SIZE 40	
#elif PARAM==5
#define SIZE 80
#elif PARAM==6
#define SIZE 160
#elif PARAM==99
#define SIZE 30
#else
#error unknown PARAM
#endif

#define TRUNCATE 0 
#define MAT(m, r, c)	(m)[(r)*SIZE+(c)]

struct bank_Input {
  int a[SIZE*SIZE];
  int b[SIZE*SIZE];
};

struct bank_Output {
  int r[SIZE*SIZE];
};

void outsource(struct bank_Input *input, struct bank_Output *output);

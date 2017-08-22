
// sample code

#define N 2
typedef struct { int M[N][N]; } matrix;

void square(matrix* a, matrix* s) {
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++) {
      s->M[i][j] = 0;
      for (int k = 0; k < N; k++) 
        s->M[i][j] += a->M[i][k] * a->M[k][j];
    }
}

void print(matrix* a) {
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) 
      printf("%6d ",a->M[i][j]);
    printf("\n");
  };
  printf("\n");
}

BANK(bx, matrix)
BANK(br, matrix)
PROVE(outsource)

br outsource(bx x) {
    matrix in, out;
    load_bx(x, &in);
    square(&in, &out);
    return save_br(&out);
};

void square_proxy(matrix* in, matrix* out) { 
  br r = outsource(save_bx(in));
  load_br(r, out); 
};

int main() {
  matrix A;
  matrix S;
  A.M[0][0] = 1;
  A.M[0][1] = 2;
  A.M[1][0] = 3;
  A.M[1][1] = 4;
  print(&A);
  square_proxy(&A, &S);
  print(&S);
  square(&A, &S);
  print(&S);
  return 0;
}

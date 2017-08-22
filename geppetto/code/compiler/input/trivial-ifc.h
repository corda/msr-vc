#pragma once

struct Config {
  int a;
};

struct bank_In {
  int x;
  int y;
  int z;
};

struct bank_Out {
  int r;
};

void outsource(struct Config *config, struct bank_In *input, struct bank_Out *output);

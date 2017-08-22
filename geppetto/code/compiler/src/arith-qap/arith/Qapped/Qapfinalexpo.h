#ifndef __FINALEXPO_H
#define __FINALEXPO_H

#include "Qaparith_error.h"
#include "QapFp4.h"
#include "QapFp12.h"


/***************** Internal stuff *******************/

typedef struct {
  int lenu;
  int *nafu;
  int signu;
} bn_finalexpo_config_t;

extern bn_finalexpo_config_t bn_finalexpo_config; /* Defined in finalexpo.c */

void Qapbn_finalexpo_neg (QapFp12_t c, QapFp12_t a);

int bn_finalexpo_initialize_config (void);

void bn_finalexpo_free_config (void);

#endif /* __FINALEXPO_H */

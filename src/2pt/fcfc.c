/*******************************************************************************
* 2pt/fcfc.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "eval_cf.h"

int main(int argc, char *argv[]) {
  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return FCFC_ERR_CONF;
  }

  CF *cf;
  if (!(cf = cf_init(conf))) {
    printf(FMT_FAIL);
    P_EXT("failed to initialise correlation function evaluations\n");
    conf_destroy(conf);
    return FCFC_ERR_CF;
  }

  if (eval_cf(conf, cf)) {
    printf(FMT_FAIL);
    P_EXT("failed to evaluate correlation functions\n");
    conf_destroy(conf); cf_destroy(cf);
    return FCFC_ERR_CF;
  }

  conf_destroy(conf);
  cf_destroy(cf);
  return 0;
}

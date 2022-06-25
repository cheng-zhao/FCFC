/*******************************************************************************
* 2pt/fcfc.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "eval_cf.h"
#include "define_para.h"
#ifdef MPI
#include <stdlib.h>
#endif

int main(int argc, char *argv[]) {
  CONF *conf = NULL;
  CF *cf = NULL;

#ifdef WITH_PARA
  /* Initialize parallelisms. */
  PARA para;
  para_init(&para);
#endif

#ifdef MPI
  /* Initialize configurations with the root rank only. */
  if (para.rank == para.root) {
#endif

    if (!(conf = load_conf(argc, argv
#ifdef WITH_PARA
        , &para
#endif
        ))) {
      printf(FMT_FAIL);
      P_EXT("failed to load configuration parameters\n");
      FCFC_QUIT(FCFC_ERR_CONF);
    }

    if (!(cf = cf_setup(conf
#ifdef OMP
        , &para
#endif
        ))) {
      printf(FMT_FAIL);
      P_EXT("failed to initialise correlation function evaluations\n");
      conf_destroy(conf);
      FCFC_QUIT(FCFC_ERR_CF);
    }

#ifdef MPI
  }

  /* Broadcast configurations. */
  cf_setup_worker(&cf, &para);
#endif

  if (eval_cf(conf, cf
#ifdef MPI
      , &para
#endif
      )) {
    printf(FMT_FAIL);
    P_EXT("failed to evaluate correlation functions\n");
    conf_destroy(conf); cf_destroy(cf);
    FCFC_QUIT(FCFC_ERR_CF);
  }

  conf_destroy(conf);
  cf_destroy(cf);
#ifdef MPI
  if (MPI_Finalize()) {
    P_ERR("failed to finalize MPI\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }
#endif
  return 0;
}

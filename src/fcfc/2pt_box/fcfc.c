/*******************************************************************************
* 2pt_box/fcfc.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "eval_cf.h"

#ifdef MPI
#include <stdlib.h>
#include <mpi.h>
#endif
#ifdef OMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
  CONF *conf = NULL;
  CF *cf = NULL;

#ifdef MPI
  /* MPI initialization. */
  if (MPI_Init(NULL, NULL)) {
    P_EXT("failed to initialize MPI\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }
  int rank, ntask;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)) {
    P_EXT("failed to retrieve MPI ranks\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }
  if (MPI_Comm_size(MPI_COMM_WORLD, &ntask)) {
    P_EXT("failed to retrive the number of MPI tasks\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  /* Initialize configurations with the root rank only. */
  if (rank == FCFC_MPI_ROOT) {
#endif
#ifdef OMP
    const int nthread = omp_get_max_threads();
#endif

    if (!(conf = load_conf(argc, argv
#ifdef MPI
        , ntask
#endif
#ifdef OMP
        , nthread
#endif
        ))) {
      printf(FMT_FAIL);
      P_EXT("failed to load configuration parameters\n");
      FCFC_QUIT(FCFC_ERR_CONF);
    }

    if (!(cf = cf_setup(conf
#ifdef OMP
        , nthread
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
  cf_setup_worker(&cf, rank);
#endif

  if (eval_cf(conf, cf
#ifdef MPI
      , ntask, rank
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

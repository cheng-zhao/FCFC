/*******************************************************************************
* parallel.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#if defined(OMP) || defined(MPI)

#include "define_comm.h"
#include "define_para.h"

#ifdef MPI
#include <stdio.h>
#include <stdlib.h>
#endif
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                      Function for setting up parallelisms
\*============================================================================*/

/******************************************************************************
Function `para_init`:
  Initialise the structure for parallelisms.
******************************************************************************/
void para_init(PARA *para) {
#ifdef MPI
  /* MPI initialization. */
  if (MPI_Init(NULL, NULL)) {
    P_EXT("failed to initialize MPI\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  para->root = FCFC_MPI_ROOT;
  para->comm = MPI_COMM_WORLD;

  if (MPI_Comm_rank(para->comm, &para->rank)) {
    P_EXT("failed to retrieve MPI ranks\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }
  if (MPI_Comm_size(para->comm, &para->ntask)) {
    P_EXT("failed to retrive the number of MPI tasks\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  if (para->root < 0 || para->root >= para->ntask) {
    P_EXT("invalid root rank defined by `FCFC_MPI_ROOT`: %d\n", para->root);
    FCFC_QUIT(FCFC_ERR_MPI);
  }
#endif
#ifdef OMP
  para->nthread = omp_get_max_threads();
#endif
}

#endif

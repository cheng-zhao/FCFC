/*******************************************************************************
* 2pt/eval_cf.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "eval_cf.h"
#include "count_func.h"
#include "read_file.h"
#include "read_res.h"
#include "save_res.h"
#include "legpoly.h"
#include "build_tree.h"
#include <stdlib.h>
#include <math.h>

#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
#define CLEAN_TREE                                                      \
  for (int _i = 0; _i < cf->ncat; _i++)                                 \
    tree_destroy(tree[_i], cf->treetype);                               \
  free(tree);


/*============================================================================*\
            Functions for steps of correlation function evaluations
\*============================================================================*/

/******************************************************************************
Function `eval_pairs`:
  Evaluate pair counts for the input catalogs.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations;
  * `ntask`:    number of MPI tasks;
  * `rank`:     ID of MPI task.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int eval_pairs(const CONF *conf, CF *cf
#ifdef MPI
    , const int ntask, const int rank
#endif
    ) {
  /* Allocate memory for trees. */
  void **tree;
  if (!(tree = malloc(sizeof(void *) * cf->ncat))) {
    P_ERR("failed to allocate memory for trees\n");
    FCFC_QUIT(FCFC_ERR_MEMORY);
  }
  for (int i = 0; i < cf->ncat; i++) tree[i] = NULL;

  for (int i = 0; i < cf->npc; i++) {
    /* Read pair counts from file if applicable. */
    if (!cf->comp_pc[i]) {
#ifdef MPI
      if (rank == FCFC_MPI_ROOT) {
#endif
        printf("Reading %s pairs ...", conf->pc[i]);
        if (conf->verbose) printf("\n  Filename: %s\n", conf->pcout[i]);
        fflush(stdout);

        int e = read_pair_count(conf->pcout[i], cf, i);
        if (e) {
          CLEAN_TREE;
          FCFC_QUIT(e);
        }

        printf(FMT_DONE);
#ifdef MPI
        fflush(stdout);
      }
#endif
      continue;
    }

    int cat[2];
    cat[0] = cf->pc_idx[0][i];
    cat[1] = cf->pc_idx[1][i];

    /* Read catalogue and construct the tree if necessary. */
    for (int j = 0; j < 2; j++) {
      int idx = cat[j];
      if (!tree[idx]) {
        if (!(tree[idx] = tree_create(conf, cf, idx
#ifdef MPI
            , rank
#endif
            ))) {
          CLEAN_TREE;
          FCFC_QUIT(FCFC_ERR_TREE);
        }
      }
    }

    /* Count pairs. */
#ifdef MPI
    if (rank == FCFC_MPI_ROOT) {
#endif
      printf("Counting %c%c pairs ...", cf->label[cat[0]], cf->label[cat[1]]);
      if (conf->verbose) printf("\n");
      fflush(stdout);
#ifdef MPI
    }
#endif

    if (cat[0] == cat[1]) {             /* auto pairs */
#ifndef MPI
      count_pairs(tree[cat[0]], tree[cat[0]], cf, cf->cnt[i], true, cf->wt[i]);
#else
      if (rank == FCFC_MPI_ROOT) {
        count_pairs(tree[cat[0]], tree[cat[0]], cf, cf->cnt[i], true,
            cf->wt[i], ntask, rank);
#endif
        /* Double auto pairs. */
        if (cf->wt[i]) {
          for (size_t k = 0; k < cf->ntot; k++) cf->cnt[i][k].d *= 2;
          cf->norm[i] = cf->data[cat[0]].wt * (cf->data[cat[0]].wt - 1);
        }
        else {
          for (size_t k = 0; k < cf->ntot; k++) cf->cnt[i][k].i *= 2;
          cf->norm[i] = (double) cf->data[cat[0]].n * (cf->data[cat[0]].n - 1);
        }
#ifdef MPI
      }
      else count_pairs(tree[cat[0]], tree[cat[0]], cf, cf->cnt[0], true,
          cf->wt[i], ntask, rank);
#endif
    }
    else {                              /* cross counts */
#ifndef MPI
      count_pairs(tree[cat[0]], tree[cat[1]], cf, cf->cnt[i], false, cf->wt[i]);
#else
      if (rank == FCFC_MPI_ROOT) {
        count_pairs(tree[cat[0]], tree[cat[1]], cf, cf->cnt[i], false,
            cf->wt[i], ntask, rank);
#endif
        if (cf->wt[i]) cf->norm[i] = cf->data[cat[0]].wt * cf->data[cat[1]].wt;
        else cf->norm[i] = (double) cf->data[cat[0]].n * cf->data[cat[1]].n;
#ifdef MPI
      }
      else count_pairs(tree[cat[0]], tree[cat[1]], cf, cf->cnt[0], false,
          cf->wt[i], ntask, rank);
#endif
    }
    /* Normalise pair counts. */
#ifdef MPI
    if (rank == FCFC_MPI_ROOT) {
#endif
      if (cf->wt[i]) {
        for (size_t k = 0; k < cf->ntot; k++)
          cf->ncnt[i][k] = cf->cnt[i][k].d / cf->norm[i];
      }
      else {
        for (size_t k = 0; k < cf->ntot; k++)
          cf->ncnt[i][k] = cf->cnt[i][k].i / cf->norm[i];
      }

      /* Save pair counts. */
      int e = save_res(conf, cf, i, FCFC_OUTPUT_PAIR_COUNT);
      if (e) {
        CLEAN_TREE;
        FCFC_QUIT(e);
      }
#ifdef MPI
    }
#endif

    /* Release memory if necessary. */
    for (int j = 0; j < 2; j++) {
      bool free_data = true;
      for (int k = i + 1; k < cf->npc; k++) {
        if (cat[j] == cf->pc_idx[0][k] || cat[j] == cf->pc_idx[1][k]) {
          free_data = false;
          break;
        }
      }
      if (free_data) {
        DATA *data = cf->data + cat[j];
        for (int k = 0; k < FCFC_XDIM; k++) {
          if (data->x[k]) {
            free(data->x[k]);
            data->x[k] = NULL;
          }
        }
        if (data->w) {
          free(data->w);
          data->w = NULL;
        }
        tree_destroy(tree[cat[j]], cf->treetype);
        tree[cat[j]] = NULL;
      }
    }
#ifdef MPI
    if (rank == FCFC_MPI_ROOT) {
#endif
      printf(FMT_DONE);
#ifdef MPI
      fflush(stdout);
    }
#endif
  }

  CLEAN_TREE;

  return 0;
}

/******************************************************************************
Function `eval_cf_exp`:
  Evaluate correlation functions given the expressions for estimators.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int eval_cf_exp(const CONF *conf, CF *cf) {
  printf("Evaluate correlation function estimators ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* Prepare the array of normalised pair counts, for libast evaluations. */
#ifdef OMP
  double *pc = calloc((size_t) cf->nthread * cf->npc, sizeof(double));
#else
  double *pc = calloc(cf->npc, sizeof(double));
#endif
  if (!pc) {
    P_ERR("failed to allocate memory for 2PCF evaluation\n");
    return FCFC_ERR_MEMORY;
  }

  for (int i = 0; i < cf->ncf; i++) {
#ifdef OMP
  #pragma omp parallel num_threads(cf->nthread) default(none)   \
    firstprivate(pc) shared(cf,i,stderr)
    {
      const int tid = omp_get_thread_num();
      pc += (size_t) tid * cf->npc;
  #pragma omp for
#endif
      for (size_t j = 0; j < cf->ntot; j++) {
        /* Set pair counts to be passed to libast. */
        for (int k = 0; k < cf->npc; k++) pc[k] = cf->ncnt[k][j];
        /* Evaluate the 2PCF. */
        if (ast_eval_num(cf->ast_cf[i], cf->cf[i] + j, pc, cf->npc)) {
          ast_perror(cf->ast_cf[i], stderr,
              FMT_ERR " failed to evaluate 2PCF:");
#ifdef OMP
          exit(FCFC_ERR_AST);
#else
          free(pc); return FCFC_ERR_AST;
#endif
        }
      }
#ifdef OMP
    }
#endif

    /* Save the correlation function. */
    int e = save_res(conf, cf, i, FCFC_OUTPUT_2PCF_RAW);
    if (e) {
      free(pc); return e;
    }
  }

  free(pc);
  printf(FMT_DONE);
#ifdef MPI
  fflush(stdout);
#endif
  return 0;
}

/******************************************************************************
Function `eval_cf_mp`:
  Evaluate correlation function multipoles given xi(s,mu).
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int eval_cf_mp(const CONF *conf, CF *cf) {
  printf("Compute correlation function multipoles ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  for (int i = 0; i < cf->ncf; i++) {
    for (int l = 0; l < cf->nl; l++) {
      int ell = cf->poles[l];
      double fac = (2 * ell + 1) / (double) cf->nmu;
      for (int j = 0; j < cf->nmu; j++) {
        double mu = (j + 0.5) / (double) cf->nmu;
        for (int k = 0; k < cf->ns; k++) {
          cf->mp[i][k + l * cf->ns] +=
              cf->cf[i][k + j * cf->ns] * fac * legpoly(ell, mu);
        }
      }
    }

    /* Save the correlation function multipoles. */
    int e = save_res(conf, cf, i, FCFC_OUTPUT_2PCF_INTEG);
    if (e) return e;
  }

  printf(FMT_DONE);
#ifdef MPI
  fflush(stdout);
#endif
  return 0;
}

/******************************************************************************
Function `eval_cf_wp`:
  Evaluate projected correlation functions given xi(s_perp,pi).
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int eval_cf_wp(const CONF *conf, CF *cf) {
  printf("Compute projected correlation functions ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  for (int i = 0; i < cf->ncf; i++) {
    for (int j = 0; j < cf->np; j++) {
      double dpi = cf->pbin[j + 1] - cf->pbin[j];
      for (int k = 0; k < cf->ns; k++) {
        cf->wp[i][k] += 2 * cf->cf[i][k + j * cf->ns] * dpi;
      }
    }

    /* Save the projected correlation functions. */
    int e = save_res(conf, cf, i, FCFC_OUTPUT_2PCF_INTEG);
    if (e) return e;
  }

  printf(FMT_DONE);
#ifdef MPI
  fflush(stdout);
#endif
  return 0;
}


/*============================================================================*\
                 Interface for correlation function evaluations
\*============================================================================*/

/******************************************************************************
Function `eval_cf`:
  Evaluate correlation functions.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations;
  * `ntask`:    number of MPI tasks;
  * `rank`:     ID of MPI task.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int eval_cf(const CONF *conf, CF *cf
#ifdef MPI
    , const int ntask, const int rank
#endif
    ) {
  int e;
  if ((e = eval_pairs(conf, cf
#ifdef MPI
      , ntask, rank
#endif
      ))) return e;
#ifdef MPI
  if (rank == FCFC_MPI_ROOT) {
#endif
    if (cf->ncf) {
      if ((e = eval_cf_exp(conf, cf))) return e;
      if (cf->nl) if ((e = eval_cf_mp(conf, cf))) return e;
      if (cf->comp_wp && cf->bintype == FCFC_BIN_SPI) {
        if ((e = eval_cf_wp(conf, cf))) return e;
      }
    }
#ifdef MPI
  }
#endif
  return 0;
}

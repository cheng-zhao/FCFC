/*******************************************************************************
* 2pt_box/init_cf.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "eval_cf.h"
#include "read_file.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
              Function for parsing correlation function estimators
\*============================================================================*/

/******************************************************************************
Function `cf_expression`:
  Convert the strings of correlation function estimators to libast expressions.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cf_expression(const CONF *conf, CF *cf) {
  if (cf->npc >= 999) {       /* ensure that 3 digits are enough */
    P_ERR("unexpected number of pairs: %d\n", conf->npc);
    return FCFC_ERR_CFG;
  }
  for (int i = 0; i < cf->ncf; i++) {
    size_t len = strlen(conf->cf[i]);
    if (len >= (SIZE_MAX - 1) / 3) {
      P_ERR("unexpected length of the %d-th " FMT_KEY(CF_ESTIMATOR) "\n",
          i + 1);
      return FCFC_ERR_CFG;
    }
    /* There are at most 26 catalogs, so (26^2 + 1) different pairs, including
     * the analytical RR. Hence at most 3 digits are needed for representing
     * the pairs. Given the format of libast variables: ${NUM}, 6 characters
     * are enough for every pair combination, which is 3 times the length of
     * two letters. An extra character is requested for the terminate null. */
    if (!(cf->cf_exp[i] = calloc(len * 3 + 1, sizeof(char)))) {
      P_ERR("failed to allocate memory for correlation function expressions\n");
      return FCFC_ERR_MEMORY;
    }

    size_t k = 0;
    for (char *c = conf->cf[i]; *c; c++) {
      /* Copy the character directly if it is not a catalog label. */
      if (*c != FCFC_SYM_ANA_RR && (*c < 'A' || *c > 'Z'))
        cf->cf_exp[i][k++] = *c;
      else {                    /* replace pairs by libast variables */
        int num = 0;
        /* Replace analytical RR by ${cf->npc + 1}. */
        if (*c == FCFC_SYM_ANA_RR && *c == c[1]) {
          num = cf->npc + 1;
          /* Allocate memory for analytical RR. */
          if (!cf->rr && !(cf->rr = malloc(sizeof(double) * cf->ntot))) {
            P_ERR("failed to allocate memory for analytical RR\n");
            return FCFC_ERR_MEMORY;
          }
        }
        /* Invalid expression if there is only one captial letter. */
        else if (*c == FCFC_SYM_ANA_RR || c[1] < 'A' || c[1] > 'Z') {
          P_ERR("unbalanced label for pairs:\n");
          fprintf(stderr, "%s\n", conf->cf[i]);
          for (size_t ii = 0; conf->cf[i] + ii < c; ii++) fprintf(stderr, " ");
          fprintf(stderr, "^\n");
          return FCFC_ERR_CFG;
        }
        /* Find the index of the pair. */
        else {
          bool found = false;
          int ii;
          for (ii = 0; ii < conf->npc; ii++) {
            if (c[0] == conf->pc[ii][0] && c[1] == conf->pc[ii][1]) {
              found = true;
              break;
            }
          }
          if (!found) {
            P_ERR("pair not found in " FMT_KEY(PAIR_COUNT) ":\n");
            fprintf(stderr, "%s\n", conf->cf[i]);
            for (size_t jj = 0; conf->cf[i] + jj < c; jj++)
              fprintf(stderr, " ");
            fprintf(stderr, "^\n");
            return FCFC_ERR_CFG;
          }
          num = ii + 1;
        }

        /* Add the variable indicater to the expression. */
        cf->cf_exp[i][k++] = AST_VAR_FLAG;
        cf->cf_exp[i][k++] = AST_VAR_START;
        int val = num / 100;
        if (val) {
          cf->cf_exp[i][k++] = val + '0';
          num -= val * 100;
        }
        val = num / 10;
        cf->cf_exp[i][k++] = val + '0';
        num -= val * 10;
        cf->cf_exp[i][k++] = num + '0';
        cf->cf_exp[i][k++] = AST_VAR_END;
        c += 1;
      }
    }
  }
  return 0;
}


/*============================================================================*\
        Interfaces for correlation function initialisation and clean-up
\*============================================================================*/

/******************************************************************************
Function `cf_init`:
  Initialise the structure for correlation function evaluations.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
CF *cf_init(const CONF *conf) {
  printf("Initialising the correlation function calculation ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  CF *cf = calloc(1, sizeof(CF));
  if (!cf) {
    P_ERR("failed to allocate memory for the intialisation\n");
    return NULL;
  }

  /* Initialise all pointers for better error handling. */
  cf->sbin = cf->s2bin = cf->pbin = NULL;
  cf->stab = cf->mutab = cf->ptab = cf->ndata = NULL;
  cf->data = NULL;
  cf->pc_idx[0] = cf->pc_idx[1] = NULL;
  cf->cnt = NULL;
  cf->norm = cf->rr = NULL;
  cf->pcnt = NULL;
  cf->ncnt = cf->cf = cf->mp = cf->wp = NULL;
  cf->cf_exp = NULL;
  cf->ast_cf = NULL;

  /* Initialisation. */
#ifdef OMP
  cf->nthread = omp_get_max_threads();
#else
  cf->nthread = 1;
#endif
  cf->bintype = conf->bintype;
  cf->bsize = conf->bsize;
  cf->prec = REAL_NAN;
  if (conf->dprec != INT_MAX) {
    cf->prec = 1;
    int prec = conf->dprec;
    while (prec > 0) {
      cf->prec *= 0.1;
      prec -= 1;
    }
    while (prec < 0) {
      cf->prec *= 10;
      prec += 1;
    }
  }

  cf->ns = conf->nsbin;
  cf->nmu = conf->nmu;
  cf->np = conf->npbin;

  cf->ncat = conf->ninput;
  cf->label = conf->label;
  cf->npc = conf->npc;
  cf->comp_pc = conf->comp_pc;
  cf->ncf = conf->ncf;
  cf->nl = conf->npole;
  cf->poles = conf->poles;
  cf->comp_wp = conf->wp;

  /* Define separation bins. */
  real *tx, *ty;        /* arrays for storing tables from file */
  tx = ty = NULL;
  if (conf->fsbin) {    /* read separation bins from file */
    size_t num;
    if (read_ascii_table(conf->fsbin, &tx, &ty, &num)) {
      P_ERR("failed to read separation bins from file: `%s'\n", conf->fsbin);
      cf_destroy(cf); return NULL;
    }
    if (num > FCFC_MAX_BIN_NUM) {
      P_ERR("too many separation bins in file: `%s'\n", conf->fsbin);
      free(tx); free(ty); cf_destroy(cf); return NULL;
    }
    for (size_t i = 1; i < num; i++) {
      if (tx[i] - ty[i - 1] > REAL_TOL || tx[i] - ty[i - 1] < -REAL_TOL) {
        P_ERR("discontinuous distance bin (" REAL_OFMT "," REAL_OFMT ") and ("
            REAL_OFMT "," REAL_OFMT ") in file: `%s'\n",
            tx[i - 1], ty[i - 1], tx[i], ty[i], conf->fsbin);
        free(tx); free(ty); cf_destroy(cf); return NULL;
      }
    }
    cf->ns = (int) num;
  }

  if (!(cf->sbin = malloc(sizeof(real) * (cf->ns + 1))) ||
      !(cf->s2bin = malloc(sizeof(real) * (cf->ns + 1)))) {
    P_ERR("failed to allocate memory for separation bins\n");
    if (tx) free(tx);
    if (ty) free(ty);
    free(cf); return NULL;
  }

  /* Compute edges of separation bins. */
  if (conf->fsbin) {
    cf->sbin[0] = tx[0];
    for (int i = 1; i < cf->ns; i++) {
      if (tx[i] != ty[i - 1]) cf->sbin[i] = (tx[i] + ty[i - 1]) * 0.5;
      else cf->sbin[i] = tx[i];
    }
    cf->sbin[cf->ns] = ty[cf->ns - 1];
    free(tx); free(ty);
    tx = ty = NULL;
    if (conf->verbose) printf("  %d separation bins loaded from file `%s'\n",
        cf->ns, conf->fsbin);
  }
  else {
    for (int i = 0; i <= cf->ns; i++) cf->sbin[i] = conf->smin + conf->ds * i;
  }
  for (int i = 0; i <= cf->ns; i++) cf->s2bin[i] = cf->sbin[i] * cf->sbin[i];
  cf->s2min = cf->s2bin[0];
  cf->s2max = cf->s2bin[cf->ns];

  /* Define pi bins */
  if (cf->bintype == FCFC_BIN_SPI) {
    if (conf->fpbin) {  /* read pi bins from file */
      size_t num;
      if (read_ascii_table(conf->fpbin, &tx, &ty, &num)) {
        P_ERR("faiedl to read pi bins from file: `%s'\n", conf->fpbin);
        cf_destroy(cf); return NULL;
      }
      if (num > FCFC_MAX_BIN_NUM) {
        P_ERR("too many pi bins in file: `%s'\n", conf->fpbin);
        free(tx); free(ty); cf_destroy(cf); return NULL;
      }
      for (size_t i = 1; i < num; i++) {
        if (tx[i] - ty[i - 1] > REAL_TOL || tx[i] - ty[i - 1] < -REAL_TOL) {
          P_ERR("discontinuous distance bin (" REAL_OFMT "," REAL_OFMT ") and ("
              REAL_OFMT "," REAL_OFMT ") in file: `%s'\n",
              tx[i - 1], ty[i - 1], tx[i], ty[i], conf->fsbin);
          free(tx); free(ty); cf_destroy(cf); return NULL;
        }
      }
      cf->np = (int) num;
    }

    if (!(cf->pbin = malloc(sizeof(real) * (cf->np + 1)))) {
      P_ERR("failed to allocate memory for pi bins\n");
      if (tx) free(tx);
      if (ty) free(ty);
      cf_destroy(cf); return NULL;
    }

    /* Compute edges of pi bins */
    if (conf->fpbin) {
      cf->pbin[0] = tx[0];
      for (int i = 1; i < cf->np; i++) {
        if (tx[i] != ty[i - 1]) cf->pbin[i] = (tx[i] + ty[i - 1]) * 0.5;
        else cf->pbin[i] = tx[i];
      }
      cf->pbin[cf->np] = ty[cf->np - 1];
      free(tx); free(ty);
      tx = ty = NULL;
      if (conf->verbose) printf("  %d pi bins loaded from file `%s'\n",
          cf->np, conf->fpbin);
    }
    else {
      for (int i = 0; i <= cf->np; i++)
        cf->pbin[i] = conf->pmin + conf->dpi * i;
    }
    cf->pmin = cf->pbin[0];
    cf->pmax = cf->pbin[cf->np];
  }

  if (cf->bintype == FCFC_BIN_SMU) cf->ntot = (size_t) cf->ns * cf->nmu;
  else if (cf->bintype == FCFC_BIN_SPI) cf->ntot = (size_t) cf->ns * cf->np;
  else  cf->ntot = cf->ns;      /* cf->bintype == FCFC_BIN_ISO */

  /* Setup lookup tables. */
  if (cf->prec != REAL_NAN) {
    real min = cf->s2bin[0];
    real max = cf->s2bin[cf->ns];
    size_t offset = (size_t) (min * cf->prec);
    /* Number of elements in the lookup table. */
    size_t ntab = (size_t) (max * cf->prec) - offset;

    if (!(cf->stab = malloc(sizeof(size_t) * ntab))) {
      P_ERR("failed to allocate memory for the lookup table of separations\n");
      cf_destroy(cf); return NULL;
    }

    /* Set values for the lookup table of squared separations. */
    int j = 1;
    for (size_t i = 0; i < ntab; i++) {
      if (i + offset < (size_t) (cf->s2bin[j] * cf->prec)) cf->stab[i] = j - 1;
      else {
        if (++j > cf->ns) {
          P_ERR("failed to create the lookup table of separations\n");
          cf_destroy(cf); return NULL;
        }
        i -= 1;
      }
    }

    /* Setup the table for pi bins. */
    if (cf->bintype == FCFC_BIN_SPI) {
      min = cf->pbin[0];
      max = cf->pbin[cf->np];
      offset = (size_t) (min * cf->prec);
      ntab = (size_t) (max * cf->prec) - offset;

      if (!(cf->ptab = malloc(sizeof(size_t) * ntab))) {
        P_ERR("failed to allocate memory for the lookup table of pi bins\n");
        cf_destroy(cf); return NULL;
      }

      int j = 1;
      for (size_t i = 0; i < ntab; i++) {
        if (i + offset < (size_t) (cf->pbin[j] * cf->prec)) cf->ptab[i] = j - 1;
        else {
          if (++j > cf->np) {
            P_ERR("failed to create the lookup table of pi bins\n");
            cf_destroy(cf); return NULL;
          }
          i -= 1;
        }
      }
    }
  }

  /* Setup the lookup table for mu bins. */
  if (cf->bintype == FCFC_BIN_SMU) {
    cf->nmu2 = (size_t) cf->nmu * cf->nmu;
    if (!(cf->mutab = malloc(sizeof(size_t) * cf->nmu2))) {
      P_ERR("failed to allocate memory for the lookup table of mu bins\n");
      cf_destroy(cf); return NULL;
    }

    int j = 1;
    for (size_t i = 0; i < cf->nmu2; i++) {
      if (i < (size_t) j * j) cf->mutab[i] = j - 1;
      else {
        if (++j > cf->nmu) {
          P_ERR("failed to create the lookup table for mu bins\n");
          cf_destroy(cf); return NULL;
        }
        i -= 1;
      }
    }
  }
  if (conf->verbose) printf("  Separation bins initialised successfully\n");

  /* Initialise pair counts. */
  if (!(cf->pc_idx[0] = malloc(sizeof(int) * cf->npc)) ||
      !(cf->pc_idx[1] = malloc(sizeof(int) * cf->npc))) {
    P_ERR("failed to allocate memory for initialising pair counts\n");
    cf_destroy(cf); return NULL;
  }

  for (int i = 0; i < cf->npc; i++) {
    cf->pc_idx[0][i] = cf->pc_idx[1][i] = -1;
    if (!cf->comp_pc[i]) continue;
    for (int j = 0; j < cf->ncat; j++) {
      if (cf->label[j] == conf->pc[i][0]) cf->pc_idx[0][i] = j;
      if (cf->label[j] == conf->pc[i][1]) cf->pc_idx[1][i] = j;
    }
    if (cf->pc_idx[0][i] == -1 || cf->pc_idx[1][i] == -1) {
      P_ERR("catalog not found for pair count: %s\n", conf->pc[i]);
      cf_destroy(cf); return NULL;
    }
  }

  /* Check if any catalogue is not used. */
  for (int i = 0; i < cf->ncat; i++) {
    bool found = false;
    for (int j = 0; j < cf->npc; j++) {
      if (cf->comp_pc[j] && (i == cf->pc_idx[0][j] || i == cf->pc_idx[1][j])) {
        found = true;
        break;
      }
    }
    if (!found)
      P_WRN("catalog <%c> is not required for pair counting\n", cf->label[i]);
  }

  /* Allocate memory for the catalogues, pair counts, and 2PCFs. */
  if (!(cf->data = malloc(sizeof(DATA *) * cf->ncat))) {
    P_ERR("failed to allocate memory for the input catalogs\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 0; i < cf->ncat; i++) cf->data[i] = NULL;
  if (!(cf->ndata = calloc(cf->ncat, sizeof(size_t)))) {
    P_ERR("failed to allocate memory for the number of input objects\n");
    cf_destroy(cf); return NULL;
  }

  if (!(cf->cnt = malloc(sizeof(size_t *) * cf->npc))) {
    P_ERR("failed to allocate memory for pair counts\n");
    cf_destroy(cf); return NULL;
  }
  cf->cnt[0] = NULL;    /* memory will be allocated only at the first element */
  if (!(cf->norm = calloc(cf->npc, sizeof(double)))) {
    P_ERR("failed to allocate memory for the normalisations of pair counts\n");
    cf_destroy(cf); return NULL;
  }
  if (!(cf->ncnt = malloc(sizeof(double *) * cf->npc))) {
    P_ERR("failed to allocate memory for normalised pair counts\n");
    cf_destroy(cf); return NULL;
  }
  cf->ncnt[0] = NULL;
  /* Allocate memory only for the first elements of arrays. */
  if (!(cf->cnt[0] = calloc(cf->ntot * cf->npc, sizeof(size_t))) ||
      !(cf->ncnt[0] = malloc(sizeof(double) * cf->ntot * cf->npc))) {
    P_ERR("failed to allocate memory for pair counts\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 1; i < cf->npc; i++) {
    cf->cnt[i] = cf->cnt[0] + cf->ntot * i;
    cf->ncnt[i] = cf->ncnt[0] + cf->ntot * i;
  }
#ifdef OMP
  /* Thread-private pair counting pool. */
  if (!(cf->pcnt = malloc(sizeof(size_t) * cf->ntot * cf->nthread))) {
    P_ERR("failed to allocate memory for thread-private counting array\n");
    cf_destroy(cf); return NULL;
  }
#endif

  if (!(cf->cf_exp = malloc(sizeof(char *) * cf->ncf))) {
    P_ERR("failed to allocate memory for correlation function estimators\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 0; i < cf->ncf; i++) cf->cf_exp[i] = NULL;
  if (!(cf->ast_cf = malloc(sizeof(ast_t *) * cf->ncf))) {
    P_ERR("failed to allocate memory for correlation function estimators\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 0; i < cf->ncf; i++) cf->ast_cf[i] = NULL;

  if (!(cf->cf = malloc(sizeof(double *) * cf->ncf))) {
    P_ERR("failed to allocate memory for correlation functions\n");
    cf_destroy(cf); return NULL;
  }
  cf->cf[0] = NULL;
  if (!(cf->cf[0] = malloc(sizeof(double) * cf->ntot * cf->ncf))) {
    P_ERR("failed to allocate memory for correlation functions\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 1; i < cf->ncf; i++) cf->cf[i] = cf->cf[0] + cf->ntot * i;

  if (cf->nl) {                 /* multipoles are required */
    if (!(cf->mp = malloc(sizeof(double *) * cf->ncf))) {
      P_ERR("failed to allocate memory for correlation function multipoles\n");
      cf_destroy(cf); return NULL;
    }
    cf->mp[0] = NULL;
    size_t ntot = (size_t) cf->nl * cf->ns;
    if (!(cf->mp[0] = calloc(ntot * cf->ncf, sizeof(double)))) {
      P_ERR("failed to allocate memory for correlation function multipoles\n");
      cf_destroy(cf); return NULL;
    }
    for (int i = 1; i < cf->ncf; i++) cf->mp[i] = cf->mp[0] + ntot * i;
  }
  else if (cf->comp_wp) {       /* projected CFs are required */
    if (!(cf->wp = malloc(sizeof(double *) * cf->ncf))) {
      P_ERR("failed to allocate memory for projected correlation functions\n");
      cf_destroy(cf); return NULL;
    }
    cf->wp[0] = NULL;
    if (!(cf->wp[0] = calloc(cf->ns * cf->ncf, sizeof(double)))) {
      P_ERR("failed to allocate memory for projected correlation functions\n");
      cf_destroy(cf); return NULL;
    }
    for (int i = 1; i < cf->ncf; i++)
      cf->wp[i] = cf->wp[0] + (size_t) cf->ns * i;
  }

  if (conf->verbose)
    printf("  Memory allocated for pair counts and correlation functions\n");

  /* Construct the abstract syntax tree for CF estimator. */
  if (cf_expression(conf, cf)) {
    cf_destroy(cf); return NULL;
  }
  for (int i = 0; i < cf->ncf; i++) {
    if (!(cf->ast_cf[i] = ast_init())) {
      ast_perror(cf->ast_cf[i], stderr,
          FMT_ERR " cannot initialise the 2PCF estimator:");
      cf_destroy(cf); return NULL;
    }
    if (ast_build(cf->ast_cf[i], cf->cf_exp[i], AST_DTYPE_DOUBLE, true)) {
      ast_perror(cf->ast_cf[i], stderr,
          FMT_ERR " cannot build the 2PCF estimator:");
      cf_destroy(cf); return NULL;
    }
  }

  printf(FMT_DONE);
  return cf;
}

/******************************************************************************
Function `cf_destroy`:
  Release memory allocated for the correlation function calculator.
Arguments:
  * `cf`:       structure for correlation function evaluations.
******************************************************************************/
void cf_destroy(CF *cf) {
  if (!cf) return;
  if (cf->sbin) free(cf->sbin);
  if (cf->s2bin) free(cf->s2bin);
  if (cf->stab) free(cf->stab);
  if (cf->mutab) free(cf->mutab);
  if (cf->pbin) free(cf->pbin);
  if (cf->ptab) free(cf->ptab);
  if (cf->data) {
    for (int i = 0; i < cf->ncat; i++) if (cf->data[i]) free(cf->data[i]);
    free(cf->data);
  }
  if (cf->ndata) free(cf->ndata);
  if (cf->pc_idx[0]) free(cf->pc_idx[0]);
  if (cf->pc_idx[1]) free(cf->pc_idx[1]);
  if (cf->cnt) {
    if (*cf->cnt) free(*cf->cnt);       /* only first element allocated */
    free(cf->cnt);
  }
  if (cf->norm) free(cf->norm);
  if (cf->ncnt) {
    if (*cf->ncnt) free(*cf->ncnt);
    free(cf->ncnt);
  }
  if (cf->rr) free(cf->rr);
  if (cf->pcnt) free(cf->pcnt);
  if (cf->cf_exp) {
    for (int i = 0; i < cf->ncf; i++) if (cf->cf_exp[i]) free(cf->cf_exp[i]);
    free(cf->cf_exp);
  }
  if (cf->ast_cf) {
    for (int i = 0; i < cf->ncf; i++)
      if (cf->ast_cf[i]) ast_destroy(cf->ast_cf[i]);
    free(cf->ast_cf);
  }
  if (cf->cf) {
    if (*cf->cf) free(*cf->cf);
    free(cf->cf);
  }
  if (cf->mp) {
    if (*cf->mp) free(*cf->mp);
    free(cf->mp);
  }
  if (cf->wp) {
    if (*cf->wp) free(*cf->wp);
    free(cf->wp);
  }
  free(cf);
}

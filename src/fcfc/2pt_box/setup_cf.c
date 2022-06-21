/*******************************************************************************
* 2pt_box/init_cf.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "eval_cf.h"
#include "read_file.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif

/*============================================================================*\
                 Functions for setting up correlation functions
\*============================================================================*/

/******************************************************************************
Function `cf_init`:
  Initialise the structure for correlation function evaluations.
Arguments:
  * `nthread`:  number of OpenMP threads.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
static CF *cf_init(void) {
  CF *cf = calloc(1, sizeof(CF));
  if (!cf) return NULL;

  /* Initialize all pointers. */
  cf->sbin = cf->s2bin = cf->pbin = cf->sbin_raw = cf->pbin_raw = NULL;
  cf->stab = cf->ptab = NULL;
  cf->mutab = NULL;
  cf->data = NULL;
  cf->pc_idx[0] = cf->pc_idx[1] = NULL;
  cf->wt = cf->cat_wt = NULL;
#ifdef MPI
  cf->comp_pc = NULL;
#endif
  cf->cnt = NULL;
  cf->norm = cf->rr = NULL;
  cf->pcnt = NULL;
  cf->ncnt = cf->cf = cf->mp = cf->wp = NULL;
  cf->cf_exp = NULL;
  cf->ast_cf = NULL;

  return cf;
}

/******************************************************************************
Function `data_init`:
  Initialise the structure for the input catalog.
Arguments:
  * `data`:     address of the structure for the input catalog.
******************************************************************************/
static inline void data_init(DATA *data) {
  data->x[0] = data->x[1] = data->x[2] = data->w = NULL;
}

/******************************************************************************
Function `read_bins`:
  Read separation bins from a file.
Arguments:
  * `fname`:    name of the file to be read;
  * `num`:      number of separation bins read from file;
  * `bins`:     address of the array for separation bin edges.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_bins(const char *fname, int *num, real **bins) {
  double *tx, *ty;
  size_t nread;
  /* Read bin edges from file. */
  if (read_ascii_table(fname, &tx, &ty, &nread)) {
    P_ERR("failed to read separation bins from files: `%s'\n", fname);
    return FCFC_ERR_FILE;
  }
  /* Validate separation bins. */
  if (nread > FCFC_MAX_NSBIN) {
    P_ERR("too many bins (%zu) read from file: `%s'\n", nread, fname);
    free(tx); free(ty);
    return FCFC_ERR_FILE;
  }
  for (size_t i = 1; i < nread; i++) {
    if (fabs(tx[i] - ty[i - 1]) > DOUBLE_TOL) {
      P_ERR("discontinuous separation bin (" OFMT_DBL "," OFMT_DBL ") and ("
          OFMT_DBL "," OFMT_DBL ") in file: `%s'\n",
          tx[i - 1], ty[i - 1], tx[i], ty[i], fname);
      free(tx); free(ty);
      return FCFC_ERR_FILE;
    }
  }

  /* Allocate memory and store bins. */
  real *sbin = malloc(sizeof(real) * (nread + 1));
  if (!sbin) {
    P_ERR("failed to allocate memory for separation bins\n");
    free(tx); free(ty);
    return FCFC_ERR_MEMORY;
  }

  sbin[0] = tx[0];
  for (size_t i = 1; i < nread; i++) {
    if (tx[i] != ty[i - 1]) sbin[i] = (tx[i] + ty[i - 1]) * 0.5;
    else sbin[i] = tx[i];
  }
  sbin[nread] = ty[nread - 1];

  *num = (int) nread;
  *bins = sbin;

  free(tx);
  free(ty);
  return 0;
}

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
                      Functions for creating lookup tables
\*============================================================================*/

/* Template functions for creating lookup tables. */
#ifdef FCFC_LOOKUP_TABLE_WIDTH
#undef FCFC_LOOKUP_TABLE_WIDTH
#endif

#define FCFC_LOOKUP_TABLE_WIDTH FCFC_LOOKUP_TABLE_W8
#include "create_lut.c"

#define FCFC_LOOKUP_TABLE_WIDTH FCFC_LOOKUP_TABLE_W16
#include "create_lut.c"

/******************************************************************************
Function `create_lut_int`:
  Create the lookup table for integer bin edges with the optimal data type.
Arguments:
  * `bins`:     array for the sorted bin edges;
  * `num`:      the number of bins;
  * `w`:        width of the table entries (in bit).
Return:
  Address of the lookup table on success; NULL on error.
******************************************************************************/
static inline void *create_lut_int(const real *bins, const int num, int *w) {
  if (num <= UINT8_MAX) {
    *w = FCFC_LOOKUP_TABLE_W8;
    return create_lut_int_uint8_t(bins, num);
  }
  else if (num <= UINT16_MAX) {
    *w = FCFC_LOOKUP_TABLE_W16;
    return create_lut_int_uint16_t(bins, num);
  }
  else P_ERR("unexpected number of separation bins: %d\n", num);
  return NULL;
}

/******************************************************************************
Function `create_lut_hybrid`:
  Create the hybrid lookup table with the optimal data type.
Arguments:
  * `bins`:     array for the sorted bin edges;
  * `num`:      the number of bins;
  * `w`:        width of the table entries (in bit).
Return:
  Address of the lookup table on success; NULL on error.
******************************************************************************/
static inline void *create_lut_hybrid(const real *bins, const int num, int *w) {
  if (num <= UINT8_MAX / 2) {
    *w = FCFC_LOOKUP_TABLE_W8;
    return create_lut_hybrid_uint8_t(bins, num);
  }
  else if (num <= UINT16_MAX / 2) {
    *w = FCFC_LOOKUP_TABLE_W16;
    return create_lut_hybrid_uint16_t(bins, num);
  }
  else P_ERR("unexpected number of separation bins: %d\n", num);
  return NULL;
}

/******************************************************************************
Function `gcd_size_t`:
  Compute the greatest common divisor of two size_t integers,
  using the Euclidean algorithm.
Arguments:
  * `a`:        the first integer;
  * `b`:        the second integer.
Return:
  The greatest common divisor.
******************************************************************************/
static inline size_t gcd_size_t(size_t a, size_t b) {
  while (b) {
    size_t tmp = b;
    b = a % b;
    a = tmp;
  }
  return a;
}

/******************************************************************************
Function `least_fac2`:
  Compute the least factor for two floating point numbers, that turns them
  into non-negative integers.
Arguments:
  * `a`:        the first floating point number;
  * `b`:        the second floating point number.
Return:
  The factor on success; zero if no factor exists with the given precision.
******************************************************************************/
static real least_fac2(real a, real b) {
  if (a < 0 || a > REAL_MAXINT || b < 0 || b > REAL_MAXINT) return 0;
  if (a == 0) return 1 / b;
  if (b == 0) return 1 / a;

  /* Convert both numbers to integers with the same rescaling factor. */
  size_t i = 0, ifac = 1;
  do {
    if (REAL_ABS(REAL_ROUND(a) - a) < REAL_TOL &&
        REAL_ABS(REAL_ROUND(b) - b) < REAL_TOL) break;
    a *= 10;
    b *= 10;
    ifac *= 10;
    i++;
  }
  while (i <= FCFC_BIN_MAX_DIGIT_TO_INT);
  if (i == FCFC_BIN_MAX_DIGIT_TO_INT + 1 || a > REAL_MAXINT || b > REAL_MAXINT)
    return 0;                   /* conversion failed */

  size_t ig = gcd_size_t((size_t) REAL_ROUND(a), (size_t) REAL_ROUND(b));
  return ifac / (real) ig;
}

/******************************************************************************
Function `least_fac4`:
  Compute the least factor for four floating point numbers, that turns them
  into non-negative integers.
Arguments:
  * `a`:        the first floating point number;
  * `b`:        the second floating point number;
  * `c`:        the third floating point number;
  * `d`:        the fourth floating point number.
Return:
  The factor on success; zero if no factor exists with the given precision.
******************************************************************************/
static real least_fac4(real a, real b, real c, real d) {
  if (a < 0 || a > REAL_MAXINT || b < 0 || b > REAL_MAXINT ||
      c < 0 || c > REAL_MAXINT || d < 0 || d > REAL_MAXINT) return 0;

  /* Convert both numbers to integers with the same rescaling factor. */
  size_t i = 0, ifac = 1;
  do {
    if (REAL_ABS(REAL_ROUND(a) - a) < REAL_TOL &&
        REAL_ABS(REAL_ROUND(b) - b) < REAL_TOL &&
        REAL_ABS(REAL_ROUND(c) - c) < REAL_TOL &&
        REAL_ABS(REAL_ROUND(d) - d) < REAL_TOL) break;
    a *= 10;
    b *= 10;
    c *= 10;
    d *= 10;
    ifac *= 10;
    i++;
  }
  while (i <= FCFC_BIN_MAX_DIGIT_TO_INT);
  if (i == FCFC_BIN_MAX_DIGIT_TO_INT + 1 || a > REAL_MAXINT ||
      b > REAL_MAXINT || c > REAL_MAXINT || d > REAL_MAXINT) return 0;

  size_t ig = 0;
  if (a != 0) ig = (size_t) REAL_ROUND(a);
  if (b != 0) ig = (ig == 0) ?
      (size_t) REAL_ROUND(b) : gcd_size_t(ig, (size_t) b);
  if (c != 0) ig = (ig == 0) ?
      (size_t) REAL_ROUND(c) : gcd_size_t(ig, (size_t) c);
  if (d != 0) ig = (ig == 0) ?
      (size_t) REAL_ROUND(d) : gcd_size_t(ig, (size_t) d);
  if (ig == 0) ig = 1;

  return ifac / (real) ig;
}

/******************************************************************************
Function `create_tab_sbin`:
  Rescale separation bins and create the table for index lookup.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function settings.
******************************************************************************/
static void create_tab_sbin(const CONF *conf, CF *cf) {
  if (!conf->fsbin) {   /* linear separation bins */
    const real fac = least_fac2(conf->smin, conf->ds);
    if (fac != 0) {
      /* Compute the length of the lookup table. */
      real smin = cf->sbin[0] * fac;
      smin = REAL_ROUND(smin * smin);
      real smax = cf->sbin[cf->ns] * fac;
      smax = REAL_ROUND(smax * smax);
      if (smin <= REAL_MAXINT && smax <= REAL_MAXINT) {
        size_t ntab = (size_t) smax - (size_t) smin;
        if (ntab <= FCFC_LOOKUP_INT_MAX_SIZE) {
          cf->rescale = fac;
          /* Rescale the bin edges. */
          for (int i = 0; i <= cf->ns; i++) {
            cf->sbin[i] = REAL_ROUND(cf->sbin[i] * fac);
            cf->s2bin[i] = cf->sbin[i] * cf->sbin[i];
          }
          /* Create the lookup table for integer bin edges. */
          cf->tabtype = FCFC_LOOKUP_TYPE_INT;
          cf->stab = create_lut_int(cf->s2bin, cf->ns, &cf->swidth);
          return;
        }
      }
    }
  }

  /* Hybrid lookup: nonlinear bin, or the table is too long for integer bins. */
  real smax = cf->sbin[cf->ns];
  smax *= smax;         /* the maximum squared distance */

  /* Compute the rescaling factor. */
  real fac = FCFC_LOOKUP_HYBRID_MAX_SIZE / smax;
  fac = REAL_TRUNC_FRAC(fac);   /* truncate fraction to reduce round error */
  while (fac * smax > FCFC_LOOKUP_HYBRID_MIN_SIZE * 2) fac *= 0.5;
  cf->rescale = fac;

  /* Rescale the bin edges. */
  for (int i = 0; i <= cf->ns; i++) {
    cf->sbin[i] *= fac;
    cf->s2bin[i] = cf->sbin[i] * cf->sbin[i];
  }
  /* Create the lookup table for arbitrary bins. */
  cf->tabtype = FCFC_LOOKUP_TYPE_HYBRID;
  cf->stab = create_lut_hybrid(cf->s2bin, cf->ns, &cf->swidth);
}

/******************************************************************************
Function `create_tab_sp_pi`:
  Rescale both s_perp and pi bins, and create tables for index lookup.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function settings.
******************************************************************************/
static void create_tab_sp_pi(const CONF *conf, CF *cf) {
  if (!conf->fsbin && !conf->fpbin) {   /* linear s_perp and pi bins */
    const real fac = least_fac4(conf->smin, conf->ds, conf->pmin, conf->dpi);
    if (fac != 0) {
      /* Compute the lengths of the lookup table. */
      real s1 = cf->sbin[0] * fac;
      s1 = REAL_ROUND(s1 * s1);
      real s2 = cf->sbin[cf->ns] * fac;
      s2 = REAL_ROUND(s2 * s2);
      if (s1 <= REAL_MAXINT && s2 <= REAL_MAXINT) {
        real p1 = cf->pbin[0] * fac;
        real p2 = cf->pbin[cf->np] * fac;
        size_t ntabs = (size_t) s2 - (size_t) s1;
        size_t ntabp = (size_t) p2 - (size_t) p1;
        if (ntabs <= FCFC_LOOKUP_INT_MAX_SIZE &&
            ntabp <= FCFC_LOOKUP_INT_MAX_SIZE) {
          /* Rescale the bin edhes. */
          cf->rescale = fac;
          for (int i = 0; i <= cf->ns; i++) {
            cf->sbin[i] = REAL_ROUND(cf->sbin[i] * fac);
            cf->s2bin[i] = cf->sbin[i] * cf->sbin[i];
          }
          for (int i = 0; i <= cf->np; i++)
            cf->pbin[i] = REAL_ROUND(cf->pbin[i] * fac);
          /* Create lookup tables for integer bin edges. */
          cf->tabtype = FCFC_LOOKUP_TYPE_INT;
          cf->stab = create_lut_int(cf->s2bin, cf->ns, &cf->swidth);
          cf->ptab = create_lut_int(cf->pbin, cf->np, &cf->pwidth);
          return;
        }
      }
    }
  }

  /* Hybrid lookup: nonlinear bin, or the table is too long for integer bins. */
  real smax = cf->sbin[cf->ns];
  smax *= smax;                         /* the maximum squared s_perp */
  real pmax = cf->pbin[cf->np];         /* the maximum pi */

  /* Compute the rescaling factor. */
  real larger, smaller;
  if (smax >= pmax) {
    larger = smax;
    smaller = pmax;
  }
  else {
    larger = pmax;
    smaller = smax;
  }
  real fac = FCFC_LOOKUP_HYBRID_MAX_SIZE / larger;
  fac = REAL_TRUNC_FRAC(fac);   /* truncate fraction to reduce round error */
  while (fac * smaller > FCFC_LOOKUP_HYBRID_MIN_SIZE * 2) fac *= 0.5;

  /* Rescale the bin edges. */
  cf->rescale = fac;
  for (int i = 0; i <= cf->ns; i++) {
    cf->sbin[i] *= fac;
    cf->s2bin[i] = cf->sbin[i] * cf->sbin[i];
  }
  for (int i = 0; i <= cf->np; i++) cf->pbin[i] *= fac;
  /* Create hybrid lookup tables for arbitrary bin edges. */
  cf->tabtype = FCFC_LOOKUP_TYPE_HYBRID;
  cf->stab = create_lut_hybrid(cf->s2bin, cf->ns, &cf->swidth);
  cf->ptab = create_lut_hybrid(cf->pbin, cf->np, &cf->pwidth);
}

/******************************************************************************
Function `create_tab_mu`:
  Create the lookup table for squared mu bins.
Arguments:
  * `cf`:       structure for correlation function settings.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int create_tab_mu(CF *cf) {
  int nmu2 = cf->nmu * cf->nmu;
  if (!(cf->mutab = malloc(nmu2 * sizeof(uint8_t)))) {
    P_ERR("failed to allocate memory for the lookup table of mu bins\n");
    return FCFC_ERR_MEMORY;
  }
  int n = 1;
  for (int i = 0; i < nmu2; i++) {
    if (i < n * n) cf->mutab[i] = n - 1;
    else {
      if (++n > cf->nmu) {
        P_ERR("failed to create the lookup table for mu bins\n");
        return FCFC_ERR_UNKNOWN;
      }
      i--;
    }
  }
  return 0;
}


/*============================================================================*\
             Interfaces for correlation function setup and clean-up
\*============================================================================*/

/******************************************************************************
Function `cf_setup`:
  Setup the configurations for correlation function evaluations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `nthread`:  number of OpenMP threads.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
CF *cf_setup(const CONF *conf
#ifdef OMP
    , const int nthread
#endif
    ) {
  printf("Initializing the correlation function calculator ...");
  if (conf->verbose) printf("\n");
  fflush(stdout);

  CF *cf = cf_init();
  if (!cf) {
    P_ERR("failed to allocate memory for the intialization\n");
    return NULL;
  }

  /* Initialize configuration parameters. */
#ifdef OMP
  cf->nthread = nthread;
#else
  cf->nthread = 1;
#endif

  for (int i = 0; i < 3; i++) cf->bsize[i] = conf->bsize[i];
  cf->bintype = conf->bintype;
  cf->ns = conf->nsbin;
  cf->np = conf->npbin;
  cf->nmu = conf->nmu;
  cf->treetype = conf->dstruct;

  cf->ncat = conf->ninput;
  cf->label = conf->label;

  cf->npc = conf->npc;
  cf->ncf = conf->ncf;
  cf->nl = conf->npole;
  cf->poles = conf->poles;
  cf->comp_wp = conf->wp;
  cf->verbose = conf->verbose;
#ifdef MPI
  if (!(cf->comp_pc = malloc(sizeof(bool) * cf->npc))) {
    P_ERR("failed to allocate memory for the initialization\n");
    cf_destroy(cf); return NULL;
  }
  memcpy(cf->comp_pc, conf->comp_pc, sizeof(bool) * cf->npc);
#else
  cf->comp_pc = conf->comp_pc;
#endif

  /* Define separation bins. */
  if (conf->fsbin) {    /* read separation bins from file */
    if (read_bins(conf->fsbin, &cf->ns, &cf->sbin)) {
      cf_destroy(cf); return NULL;
    }
    if (conf->verbose) printf("  %d separation bins loaded from file `%s'\n",
        cf->ns, conf->fsbin);
    if (cf->sbin[cf->ns] >= conf->bsize[0] * 0.5 ||
        cf->sbin[cf->ns] >= conf->bsize[1] * 0.5 ||
        cf->sbin[cf->ns] >= conf->bsize[2] * 0.5) {
      P_ERR("the maximum separation must be smaller than half the box size\n");
      cf_destroy(cf); return NULL;
    }
  }
  else {                /* linear bins from configurations */
    if (!(cf->sbin = malloc(sizeof(real) * (cf->ns + 1)))) {
      P_ERR("failed to allocate memory for separation bins\n");
      cf_destroy(cf); return NULL;
    }
    for (int i = 0; i <= cf->ns; i++) cf->sbin[i] = conf->smin + conf->ds * i;
  }
  /* Allocate memory for unrescaled and squared bin edges. */
  if (!(cf->sbin_raw = malloc(sizeof(real) * (cf->ns + 1))) ||
      !(cf->s2bin = malloc(sizeof(real) * (cf->ns + 1)))) {
    P_ERR("failed to allocate memory for separation bins\n");
    cf_destroy(cf); return NULL;
  }
  memcpy(cf->sbin_raw, cf->sbin, sizeof(real) * (cf->ns + 1));

  /* Define pi bins */
  if (cf->bintype == FCFC_BIN_SPI) {
    if (conf->fpbin) {  /* read pi bins from file */
      if (read_bins(conf->fpbin, &cf->np, &cf->pbin)) {
        cf_destroy(cf); return NULL;
      }
      if (conf->verbose) printf("  %d pi bins loaded from file `%s'\n",
          cf->np, conf->fpbin);
      if (cf->pbin[cf->np] >= conf->bsize[0] * 0.5 ||
          cf->pbin[cf->np] >= conf->bsize[1] * 0.5 ||
          cf->pbin[cf->np] >= conf->bsize[2] * 0.5) {
        P_ERR("the maximum pi must be smaller than half the box size\n");
        cf_destroy(cf); return NULL;
      }
    }
    else {              /* linear bins from configurations */
      if (!(cf->pbin = malloc(sizeof(real) * (cf->np + 1)))) {
        P_ERR("failed to allocate memory for pi bins\n");
        cf_destroy(cf); return NULL;
      }
      for (int i = 0; i <= cf->np; i++)
        cf->pbin[i] = conf->pmin + conf->dpi * i;
    }
    /* Allocate memory for unrescaled pi bin edges. */
    if (!(cf->pbin_raw = malloc(sizeof(real) * (cf->np + 1)))) {
      P_ERR("failed to allocate memory for pi bins\n");
      cf_destroy(cf); return NULL;
    }
    memcpy(cf->pbin_raw, cf->pbin, sizeof(real) * (cf->np + 1));
  }

  /* Compute the total number of bins, and create lookup tables. */
  cf->ntot = cf->ns;
  if (cf->bintype == FCFC_BIN_SPI) {
    create_tab_sp_pi(conf, cf);
    if (!cf->stab || !cf->ptab) {
      cf_destroy(cf); return NULL;
    }
    cf->ntot *= cf->np;
  }
  else {
    create_tab_sbin(conf, cf);
    if (!cf->stab) {
      cf_destroy(cf); return NULL;
    }
    if (cf->bintype == FCFC_BIN_SMU) {
      cf->ntot *= cf->nmu;
      /* Create the lookup table for mu bins. */
      if (create_tab_mu(cf)) {
        P_ERR("failed to create the lookup table for mu bins\n");
        cf_destroy(cf); return NULL;
      }
    }
  }
  /* Rescale also the box sizes. */
  for (int i = 0; i < 3; i++) cf->bsize[i] *= cf->rescale;
  if (conf->verbose) printf("  Separation bins initialized successfully\n");

  /* Initialise pair counts. */
  if (!(cf->pc_idx[0] = malloc(sizeof(int) * cf->npc)) ||
      !(cf->pc_idx[1] = malloc(sizeof(int) * cf->npc))) {
    P_ERR("failed to allocate memory for initializing pair counts\n");
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
  if (!(cf->data = calloc((unsigned int) cf->ncat, sizeof(DATA)))) {
    P_ERR("failed to allocate memory for the input catalogs\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 0; i < cf->ncat; i++) data_init(cf->data + i);

  if (!(cf->wt = malloc(sizeof(bool) * cf->npc)) ||
      !(cf->cat_wt = calloc(cf->ncat, sizeof(bool)))) {
    P_ERR("failed to allocate memory for weight indicators\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 0; i < cf->npc; i++) {
    cf->cat_wt[cf->pc_idx[0][i]] = cf->cat_wt[cf->pc_idx[1][i]] = cf->wt[i] =
        conf->has_wt[cf->pc_idx[0][i]] || conf->has_wt[cf->pc_idx[1][i]];
  }

  if (!(cf->cnt = malloc(sizeof(COUNT *) * cf->npc))) {
    P_ERR("failed to allocate memory for pair counts\n");
    cf_destroy(cf); return NULL;
  }
  cf->cnt[0] = NULL;    /* memory will be allocated only at the first element */
  if (!(cf->norm = calloc(cf->npc, sizeof(double)))) {
    P_ERR("failed to allocate memory for the normalizations of pair counts\n");
    cf_destroy(cf); return NULL;
  }
  if (!(cf->ncnt = malloc(sizeof(double *) * cf->npc))) {
    P_ERR("failed to allocate memory for normalized pair counts\n");
    cf_destroy(cf); return NULL;
  }
  cf->ncnt[0] = NULL;
  /* Allocate memory only for the first elements of arrays. */
  if (!(cf->cnt[0] = calloc(cf->ntot * cf->npc, sizeof(COUNT))) ||
      !(cf->ncnt[0] = malloc(sizeof(double) * cf->ntot * cf->npc))) {
    P_ERR("failed to allocate memory for pair counts\n");
    cf_destroy(cf); return NULL;
  }
  for (int i = 1; i < cf->npc; i++) {
    cf->cnt[i] = cf->cnt[0] + cf->ntot * i;
    cf->ncnt[i] = cf->ncnt[0] + cf->ntot * i;
  }
#if     defined(OMP)  ||  FCFC_SIMD  >=  FCFC_SIMD_AVX512
  /* Thread-private pair counting pool. */
  #if   FCFC_SIMD  >=  FCFC_SIMD_AVX512
  if (posix_memalign((void **) &cf->pcnt, FCFC_MEMALIGN_BYTE,
                     cf->ntot * cf->nthread * FCFC_SIMD_BYTES)) {
    P_ERR("failed to allocate memory for private counting array\n");
    cf_destroy(cf); return NULL;
  }
  #else
  if (!(cf->pcnt = malloc(sizeof(COUNT) * cf->ntot * cf->nthread))) {
    P_ERR("failed to allocate memory for thread-private counting array\n");
    cf_destroy(cf); return NULL;
  }
  #endif
#endif

  if (cf->ncf) {
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
        P_ERR("failed to allocate memory for 2PCF multipoles\n");
        cf_destroy(cf); return NULL;
      }
      cf->mp[0] = NULL;
      size_t ntot = (size_t) cf->nl * cf->ns;
      if (!(cf->mp[0] = calloc(ntot * cf->ncf, sizeof(double)))) {
        P_ERR("failed to allocate memory for 2PCF multipoles\n");
        cf_destroy(cf); return NULL;
      }
      for (int i = 1; i < cf->ncf; i++) cf->mp[i] = cf->mp[0] + ntot * i;
    }
    else if (cf->bintype == FCFC_BIN_SPI && cf->comp_wp) {
      /* Projected 2PCFs are required. */
      if (!(cf->wp = malloc(sizeof(double *) * cf->ncf))) {
        P_ERR("failed to allocate memory for projected 2PCFs\n");
        cf_destroy(cf); return NULL;
      }
      cf->wp[0] = NULL;
      if (!(cf->wp[0] = calloc(cf->ns * cf->ncf, sizeof(double)))) {
        P_ERR("failed to allocate memory for projected 2PCFs\n");
        cf_destroy(cf); return NULL;
      }
      for (int i = 1; i < cf->ncf; i++)
        cf->wp[i] = cf->wp[0] + (size_t) cf->ns * i;
    }
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
          FMT_ERR " cannot initialize the 2PCF estimator:");
      cf_destroy(cf); return NULL;
    }
    if (ast_build(cf->ast_cf[i], cf->cf_exp[i], AST_DTYPE_DOUBLE, true)) {
      ast_perror(cf->ast_cf[i], stderr,
          FMT_ERR " cannot build the 2PCF estimator:");
      cf_destroy(cf); return NULL;
    }
  }

#ifndef MPI
  printf(FMT_DONE);
#endif
  return cf;
}

#ifdef MPI
/******************************************************************************
Function `cf_setup_worker`:
  Setup the correlation function configurations for MPI workers.
Arguments:
  * `cf`:       structure for correlation function configurations;
  * `rank`:     ID of MPI task.
Return:
  Address of the structure for correlation function evaluations.
******************************************************************************/
void cf_setup_worker(CF **cf, const int rank) {
  if (!cf) {
    P_ERR("the correlation function settings are not initialized\n");
    FCFC_QUIT(FCFC_ERR_ARG);
  }

  /* Initialize the structure for MPI workers (non-root tasks). */
  if (rank != FCFC_MPI_ROOT && !(*cf = cf_init())) {
    P_ERR("failed to allocate memory for initializing MPI tasks\n");
    FCFC_QUIT(FCFC_ERR_MEMORY);
  }
  CF *c = *cf;

  /* Broadcast variables. */
  MPI_Request req[11];
  if (MPI_Ibcast(c->bsize, 3, FCFC_MPI_REAL, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req) ||
      MPI_Ibcast(&c->bintype, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 1) ||
      MPI_Ibcast(&c->tabtype, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 2) ||
      MPI_Ibcast(&c->ns, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD, req + 3) ||
      MPI_Ibcast(&c->np, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD, req + 4) ||
      MPI_Ibcast(&c->nmu, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD, req + 5) ||
      MPI_Ibcast(&c->nthread, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 6) ||
      MPI_Ibcast(&c->ntot, 1, FCFC_MPI_SIZE_T, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 7) ||
      MPI_Ibcast(&c->treetype, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 8) ||
      MPI_Ibcast(&c->ncat, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 9) ||
      MPI_Ibcast(&c->npc, 1, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + 10) ||
      MPI_Waitall(11, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast correlation function settings\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  /* Allocate memory. */
  if (rank != FCFC_MPI_ROOT) {
    if (!(c->s2bin = malloc(sizeof(real) * (c->ns + 1))) ||
        (c->treetype == FCFC_STRUCT_BALLTREE &&
         !(c->sbin = malloc(sizeof(real) * (c->ns + 1)))) ||
        (c->bintype == FCFC_BIN_SPI &&
        !(c->pbin = malloc(sizeof(real) * (c->np + 1)))) ||
#if     defined(OMP)  ||  FCFC_SIMD  >=  FCFC_SIMD_AVX512
  #if     FCFC_SIMD  >=  FCFC_SIMD_AVX512
        !(c->pcnt = malloc(c->ntot * c->nthread * FCFC_SIMD_BYTES)) ||
  #else
        !(c->pcnt = malloc(sizeof(COUNT) * c->ntot * c->nthread)) ||
  #endif
#endif
        !(c->cnt = malloc(sizeof(COUNT *))) ||
        !(c->cnt[0] = malloc(sizeof(COUNT) * c->ntot)) ||
        !(c->wt = malloc(sizeof(bool) * c->npc)) ||
        !(c->pc_idx[0] = malloc(sizeof(int) * c->npc)) ||
        !(c->pc_idx[1] = malloc(sizeof(int) * c->npc)) ||
        !(c->comp_pc = malloc(sizeof(bool) * c->npc)) ||
        !(c->data = calloc(c->ncat, sizeof(DATA)))) {
      P_ERR("failed to allocate memory for correlation function settings\n");
      FCFC_QUIT(FCFC_ERR_MEMORY);
    }
    for (int i = 0; i < c->ncat; i++) data_init(c->data + i);
  }

  /* Broadcast pair count settings and bin edges. */
  int nreq = 0;
  if (MPI_Ibcast(c->wt, c->npc, MPI_C_BOOL, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + nreq++) ||
      MPI_Ibcast(c->pc_idx[0], c->npc, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + nreq++) ||
      MPI_Ibcast(c->pc_idx[1], c->npc, MPI_INT, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + nreq++) ||
      MPI_Ibcast(c->comp_pc, c->npc, MPI_C_BOOL, FCFC_MPI_ROOT, MPI_COMM_WORLD,
          req + nreq++) ||
      MPI_Ibcast(c->s2bin, c->ns + 1, FCFC_MPI_REAL, FCFC_MPI_ROOT,
          MPI_COMM_WORLD, req + nreq++) ||
      (c->treetype == FCFC_STRUCT_BALLTREE && MPI_Ibcast(c->sbin, c->ns + 1,
          FCFC_MPI_REAL, FCFC_MPI_ROOT, MPI_COMM_WORLD, req + nreq++)) ||
      (c->bintype == FCFC_BIN_SPI && MPI_Ibcast(c->pbin, c->np + 1,
          FCFC_MPI_REAL, FCFC_MPI_ROOT, MPI_COMM_WORLD, req + nreq++)) ||
      MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast pair count settings and separation bins\n");
    FCFC_QUIT(FCFC_ERR_MPI);
  }

  /* Create lookup tables. */
  if (rank != FCFC_MPI_ROOT) {
    if (c->tabtype == FCFC_LOOKUP_TYPE_INT) {
      if (!(c->stab = create_lut_int(c->s2bin, c->ns, &c->swidth)) ||
          (c->bintype == FCFC_BIN_SPI &&
          !(c->ptab = create_lut_int(c->pbin, c->np, &c->pwidth))) ||
          (c->bintype == FCFC_BIN_SMU && create_tab_mu(c))) {
        P_ERR("failed to create lookup table for MPI workers\n");
        FCFC_QUIT(FCFC_ERR_MPI);
      }
    }
    else if (c->tabtype == FCFC_LOOKUP_TYPE_HYBRID) {
      if (!(c->stab = create_lut_hybrid(c->s2bin, c->ns, &c->swidth)) ||
          (c->bintype == FCFC_BIN_SPI &&
          !(c->ptab = create_lut_hybrid(c->pbin, c->np, &c->pwidth))) ||
          (c->bintype == FCFC_BIN_SMU && create_tab_mu(c))) {
        P_ERR("failed to create lookup table for MPI workers\n");
        FCFC_QUIT(FCFC_ERR_MPI);
      }
    }
    else {
      P_ERR("invalid lookup table type: %d\n", c->tabtype);
      FCFC_QUIT(FCFC_ERR_UNKNOWN);
    }
  }

  if (rank == FCFC_MPI_ROOT) {
    if (c->verbose)
      printf("  Correlation function settings initialized for MPI workers\n");
    printf(FMT_DONE);
    fflush(stdout);
  }
}
#endif

/******************************************************************************
Function `cf_destroy`:
  Release memory allocated for the correlation function calculator.
Arguments:
  * `cf`:       structure for correlation function evaluations.
******************************************************************************/
void cf_destroy(CF *cf) {
  if (!cf) return;
  if (cf->s2bin) free(cf->s2bin);
  if (cf->pbin) free(cf->pbin);
  if (cf->stab) free(cf->stab);
  if (cf->ptab) free(cf->ptab);
  if (cf->mutab) free(cf->mutab);
  if (cf->data) {
    for (int i = 0; i < cf->ncat; i++) {
      for (int k = 0; k < FCFC_XDIM; k++) {
        if (cf->data[i].x[k]) free(cf->data[i].x[k]);
      }
      if (cf->data[i].w) free(cf->data[i].w);
    }
    free(cf->data);
  }
  if (cf->sbin) free(cf->sbin);
  if (cf->sbin_raw) free(cf->sbin_raw);
  if (cf->pbin_raw) free(cf->pbin_raw);
  if (cf->pc_idx[0]) free(cf->pc_idx[0]);
  if (cf->pc_idx[1]) free(cf->pc_idx[1]);
  if (cf->wt) free(cf->wt);
  if (cf->cat_wt) free(cf->cat_wt);
#ifdef MPI
  if (cf->comp_pc) free(cf->comp_pc);
#endif
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
    for (int i = 0; i < cf->ncf; i++) ast_destroy(cf->ast_cf[i]);
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

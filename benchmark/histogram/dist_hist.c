/*******************************************************************************
* benchmark/histogram/dist_hist.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "dist_hist.h"
#include "sample_dist.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

/*============================================================================*\
                     Data structure for distance histogram
\*============================================================================*/

typedef struct {
  int type;             /* type of distance bins */
  int n;                /* number of distance bins */
  real dmin;            /* the minimum distance of interest */
  real dmax;            /* the maximum distance of interest */
  real *d2;             /* squared distances of bin edges */
  real *dc;             /* center of distance bins */
  void *tab;            /* distance lookup table */
  size_t *cnt;          /* the distance histogram */
} HIST;

/*============================================================================*\
                 Functions for the distance histogram structure
\*============================================================================*/

/******************************************************************************
Function `destroy_hist`:
  Deconstruct the distance histogram.
Arguments:
  * `hist`:     pointer for the distance histogram.
******************************************************************************/
static void destroy_hist(HIST *hist) {
  if (!hist) return;
  if (hist->d2) free(hist->d2);
  if (hist->cnt) free(hist->cnt);
  if (hist->tab) free(hist->tab);
  free(hist);
}

/******************************************************************************
Function `init_hist`:
  Initialise the distance histogram.
Arguments:
  * `dmin`:     the minimum distance of interest;
  * `dmax`:     the maximum distance of interest;
  * `nbin`:     number of distance bins;
  * `type`:     type of distance bins.
Return:
  Address of the structure for distance histogram on success; NULL on error.
******************************************************************************/
static HIST *init_hist(const real dmin, const real dmax, const int nbin,
    const int type) {
  if (dmin >= dmax || nbin <= 0 || (type != BENCHMARK_HIST_BIN_LIN &&
      type != BENCHMARK_HIST_BIN_LOG)) {
    P_ERR("invalid distance bin settings\n");
    return NULL;
  }

  HIST *hist = calloc(1, sizeof(HIST));
  if (!hist) {
    P_ERR("failed to allocate memory for the distance histogram\n");
    return NULL;
  }

  hist->n = nbin;
  hist->type = type;
  hist->dmin = dmin;
  hist->dmax = dmax;
  hist->d2 = NULL;
  hist->cnt = NULL;
  hist->tab = NULL;

  if (!(hist->d2 = malloc((nbin + 1) * sizeof(real))) ||
      !(hist->dc = malloc(nbin * sizeof(real)))) {
    P_ERR("failed to allocate memory for the distance histogram\n");
    destroy_hist(hist);
    return NULL;
  }
#if BENCHMARK_SIMD >= BENCHMARK_SIMD_AVX512
  if (posix_memalign((void **) &hist->cnt, BENCHMARK_MEMALIGN_BYTE,
      nbin * BENCHMARK_SIMD_BYTES)) {
    P_ERR("failed to allocate aligned memory for the distance histogram\n");
    destroy_hist(hist);
    return NULL;
  }
  memset(hist->cnt, 0, nbin * BENCHMARK_SIMD_BYTES);
#else
  if (!(hist->cnt = calloc(nbin, sizeof(size_t)))) {
    P_ERR("failed to allocate memory for the distance histogram\n");
    destroy_hist(hist);
    return NULL;
  }
#endif

  /* Compute bin edges. */
  real rmin = dmin;
  real rmax = dmax;
  hist->d2[0] = rmin * rmin;
  hist->d2[nbin] = rmax * rmax;
  if (hist->type == BENCHMARK_HIST_BIN_LOG) {
    rmin = log(rmin);
    rmax = log(rmax);
  }
  const real db = (rmax - rmin) / nbin;
  for (int i = 0; i < nbin; i++)
    hist->dc[i] = rmin + (i + 0.5) * db;

  if (hist->type == BENCHMARK_HIST_BIN_LIN) {
    for (int i = 1; i < nbin; i++) {
      hist->d2[i] = rmin + i * db;
      hist->d2[i] *= hist->d2[i];
    }
  }
  else {
    hist->dc[0] = exp(hist->dc[0]);
    for (int i = 1; i < nbin; i++) {
      hist->dc[i] = exp(hist->dc[i]);
      hist->d2[i] = 2 * (rmin + i * db);
      hist->d2[i] = exp(hist->d2[i]);
    }
  }

  return hist;
}

#ifdef PRINT_HISTOGRAM
/******************************************************************************
Function `print_hist`:
  Print the distance histogram.
Arguments:
  * `hist`:     pointer for the distance histogram.
******************************************************************************/
static void print_hist(HIST *hist) {
  printf("The resulting distance histogram:\n");
  for (int i = 0; i < hist->n; i++) {
    printf(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " %zu\n",
        hist->d2[i], hist->d2[i + 1], hist->dc[i], hist->cnt[i]);
  }
}
#endif


/*============================================================================*\
                       Functions for rescaling distances
\*============================================================================*/

/******************************************************************************
Function `least_fac`:
  Compute the least factor for two floating point numbers, that turns them
  into positive integers.
Arguments:
  * `a`:        the first floating point number;
  * `b`:        the second floating point number;
  * `hist`:     the structure for distance histogram.
Return:
  The factor on success; zero if no factor exists with the given precision.
******************************************************************************/
static real least_fac(real a, real b) {
  if (a < 0 || b < 0) return 0;         /* both numbers must be non-negative */
  if (a == 0) return 1 / b;
  if (b == 0) return 1 / a;

  /* Convert both numbers to integers with the same rescaling factor. */
  size_t i = 0, ifac = 1;
  do {
    if (REAL_ABS(round(a) - a) < REAL_TOL &&
        REAL_ABS(round(b) - b) < REAL_TOL) break;
    a *= 10;
    b *= 10;
    ifac *= 10;
    i++;
  }
  while (i <= BENCHMARK_HIST_MAX_DIGIT_TO_INT);
  if (i == BENCHMARK_HIST_MAX_DIGIT_TO_INT + 1) return 0;

  size_t ia = round(a);
  size_t ib = round(b);

  /* Compute the greatest common divisor of the two integers. */
  while (ib) {
    i = ia % ib;
    ia = ib;
    ib = i;
  }
  return ifac / (real) ia;
}

/******************************************************************************
Function `rescale_dist_int`:
  Rescale distances to have integer edges of the bins.
Arguments:
  * `hist`:     the structure for distance histogram.
Return:
  The rescaling factor for squared distances on success; zero on error.
******************************************************************************/
static real rescale_dist_int(HIST *hist) {
  const real step = (hist->dmax - hist->dmin) / hist->n;
  real fac = least_fac(hist->dmin, step);
  if (fac == 0) {
    P_ERR("failed to rescale distances to have integer bin edges\n");
    return 0;
  }
  hist->dmin = round(hist->dmin * fac);
  hist->dmax = round(hist->dmax * fac);

  fac *= fac;
  for (int i = 0; i <= hist->n; i++) hist->d2[i] = round(hist->d2[i] * fac);
  if (hist->d2[hist->n] > BENCHMARK_HIST_MAX_INT_DIST) {
    P_ERR("the rescaled distance is too large\n");
    return 0;
  }
  return fac;
}

/******************************************************************************
Function `rescale_dist_interval`:
  Rescale distances to have the bin width of 1.
Arguments:
  * `hist`:     the structure for distance histogram.
Return:
  The rescaling factor for squared distances on success; zero on error.
******************************************************************************/
static real rescale_dist_interval(HIST *hist) {
  real fac = hist->n / (hist->dmax - hist->dmin);
  fac *= fac;
  for (int i = 0; i <= hist->n; i++) hist->d2[i] *= fac;
  return fac;
}

/******************************************************************************
Function `rescale_dist_hybrid`:
  Rescale distances for the hybrid lookup table.
Arguments:
  * `ntab`:     size of the lookup table;
  * `hist`:     the structure for distance histogram.
Return:
  The rescaling factor for squared distances on success; zero on error.
******************************************************************************/
static real rescale_dist_hybrid(const int ntab, HIST *hist) {
  const real d2max = hist->d2[hist->n];
  real fac = ntab / d2max - REAL_TOL;
  for (int i = 0; i <= hist->n; i++) hist->d2[i] *= fac;
  if (hist->d2[hist->n] > BENCHMARK_HIST_MAX_INT_DIST) {
    P_ERR("the rescaled distance is too large\n");
    return 0;
  }
  return fac;
}

/******************************************************************************
Function `rescale_dist_log`:
  Rescale distances to have the first distance bin starting from 1.
Arguments:
  * `hist`:     the structure for distance histogram.
Return:
  The rescaling factor for squared distances on success; zero on error.
******************************************************************************/
static inline real rescale_dist_log(HIST *hist) {
  real fac = 1 / hist->dmin;
  fac *= fac;
  for (int i = 0; i <= hist->n; i++) hist->d2[i] *= fac;
  return fac;
}


/*============================================================================*\
                Functions for updating the distance histogram
\*============================================================================*/

/******************************************************************************
Function `hist_bsearch`:
  Update the distance histogram using binary search.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void hist_bsearch(const real *d2, const size_t nsp, HIST *hist) {
  const real *dbin = hist->d2;
  const real d2min = dbin[0];
  const real d2max = dbin[hist->n];
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] >= d2max || d2[n] < d2min) continue;
    unsigned int l, u;
    l = 0;
    u = hist->n;
    while (l <= u) {
      unsigned int i = (l + u) >> 1;
      if (dbin[i + 1] <= d2[n]) l = i + 1;
      else if (dbin[i] > d2[n]) u = i - 1;
      else {
        hist->cnt[i]++;
        break;
      }
    }
  }
}

/******************************************************************************
Function `hist_rev_traversal`:
  Update the distance histogram with reversed traversal.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void hist_rev_traversal(const real *d2, const size_t nsp, HIST *hist) {
  const real *dbin = hist->d2;
  const real d2min = dbin[0];
  const real d2max = dbin[hist->n];
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      for (int i = hist->n - 1; i >= 0; i--) {
        if (d2[n] >= dbin[i]) {
          hist->cnt[i]++;
          break;
        }
      }
    }
  }
}

/******************************************************************************
Function `hist_sqrt_trunc`:
  Update the distance histogram by truncating the distance.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void hist_sqrt_trunc(const real *d2, const size_t nsp, HIST *hist) {
  const real d2min = hist->d2[0];
  const real d2max = hist->d2[hist->n];
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      const int i = REAL_SQRT(d2[n]) - hist->dmin;
      if (i < hist->n) hist->cnt[i]++;
    }
  }
}

/******************************************************************************
Function `hist_log_trunc`:
  Update the distance histogram by truncating the logarithm squared distance.
Arguments:
  * `d2`:       array for the sampled squared distances;
  * `nsp`:      number of sampled squared distances;
  * `hist`:     the structure for distance histogram.
******************************************************************************/
static void hist_log_trunc(const real *d2, const size_t nsp, HIST *hist) {
  const real fac = 0.5 * hist->n / REAL_LOG(hist->dmax / hist->dmin);
  const real d2min = hist->d2[0];
  const real d2max = hist->d2[hist->n];
  for (size_t n = 0; n < nsp; n++) {
    if (d2[n] < d2max && d2[n] >= d2min) {
      const int i = REAL_LOG(d2[n]) * fac;
      if (i < hist->n) hist->cnt[i]++;
    }
  }
}


/*============================================================================*\
      Template functions for lookup tables and approximate math functions
\*============================================================================*/

#ifdef BENCHMARK_HIST_TABLE_WIDTH
  #undef BENCHMARK_HIST_TABLE_WIDTH
#endif

#define BENCHMARK_HIST_TABLE_WIDTH      BENCHMARK_HIST_TABLE_W8
#include "hist_lookup.c"

#define BENCHMARK_HIST_TABLE_WIDTH      BENCHMARK_HIST_TABLE_W16
#include "hist_lookup.c"

#define BENCHMARK_HIST_TABLE_WIDTH      BENCHMARK_HIST_TABLE_W32
#include "hist_lookup.c"


#ifdef BENCHMARK_HIST_APPROX_POLY
  #undef BENCHMARK_HIST_APPROX_POLY
#endif

#define BENCHMARK_HIST_APPROX_POLY      1
#include "fast_math.c"

#define BENCHMARK_HIST_APPROX_POLY      2
#include "fast_math.c"

#define BENCHMARK_HIST_APPROX_POLY      3
#include "fast_math.c"

#define BENCHMARK_HIST_APPROX_POLY      4
#include "fast_math.c"

#define BENCHMARK_HIST_APPROX_POLY      5
#include "fast_math.c"


/*============================================================================*\
                  Interface for updating distance histogram
\*============================================================================*/

/******************************************************************************
Function `eval_hist`:
  Evaluate the distance histogram.
Arguments:
  * `smax`:     the maximum squared distance to be sampled;
  * `nsp`:      number of sampled squared distances;
  * `method`:   method for evaluating the distance histogram;
  * `dmin`:     the minimum distance of interest;
  * `dmax`:     the maximum distance of interest;
  * `nbin`:     number of distance bins;
  * `type`:     type of distance bins.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int eval_hist(const real smax, const size_t nsp, const int method,
    const real dmin, const real dmax, const int nbin, const int type) {
  /* Create and rescale distance bins. */
  HIST *hist = init_hist(dmin, dmax, nbin, type);
  if (!hist) return EXIT_FAILURE;
  real rescale = 1;

  if (method < 0) {
    rescale = rescale_dist_hybrid(-method, hist);
  }
  else if (method == BENCHMARK_HIST_SQRT_TRUNC || (method >=
      BENCHMARK_HIST_FASTSQRT_START + BENCHMARK_HIST_APPROX_POLY_MIN && method
      <= BENCHMARK_HIST_FASTSQRT_START + BENCHMARK_HIST_APPROX_POLY_MAX)) {
    rescale = rescale_dist_interval(hist);
  }
  else if (method == BENCHMARK_HIST_LOG_TRUNC || (method >=
      BENCHMARK_HIST_FASTLOG2_START + BENCHMARK_HIST_APPROX_POLY_MIN && method
      <= BENCHMARK_HIST_FASTLOG2_START + BENCHMARK_HIST_APPROX_POLY_MAX)) {
    rescale = rescale_dist_log(hist);
  }
  else if (method == BENCHMARK_HIST_INT_TABLE) {
    rescale = rescale_dist_int(hist);
  }
  if (rescale == 0) {
    destroy_hist(hist);
    return EXIT_FAILURE;
  }

  /* Choose the function for distance histogram update. */
  void (*hist_func) (const real *, const size_t, HIST *) = NULL;
  if (method < 0) {
    if (hist->n <= UINT8_MAX / 2) {
      if (create_table_hybrid_uint8_t(hist)) {
        destroy_hist(hist);
        return EXIT_FAILURE;
      }
      hist_func = hist_table_hybrid_uint8_t;
    }
    else if (hist->n <= UINT16_MAX / 2) {
      if (create_table_hybrid_uint16_t(hist)) {
        destroy_hist(hist);
        return EXIT_FAILURE;
      }
      hist_func = hist_table_hybrid_uint16_t;
    }
    else {
      if (create_table_hybrid_uint32_t(hist)) {
        destroy_hist(hist);
        return EXIT_FAILURE;
      }
      hist_func = hist_table_hybrid_uint32_t;
    }
  }
  else {
    switch (method) {
      case BENCHMARK_HIST_BSEARCH: hist_func = hist_bsearch; break;
      case BENCHMARK_HIST_REV_TRAV: hist_func = hist_rev_traversal; break;
      case BENCHMARK_HIST_SQRT_TRUNC: hist_func = hist_sqrt_trunc; break;
      case BENCHMARK_HIST_FASTSQRT_START + 1:
        hist_func = hist_fast_sqrt_1; break;
      case BENCHMARK_HIST_FASTSQRT_START + 2:
        hist_func = hist_fast_sqrt_2; break;
      case BENCHMARK_HIST_FASTSQRT_START + 3:
        hist_func = hist_fast_sqrt_3; break;
      case BENCHMARK_HIST_FASTSQRT_START + 4:
        hist_func = hist_fast_sqrt_4; break;
      case BENCHMARK_HIST_FASTSQRT_START + 5:
        hist_func = hist_fast_sqrt_5; break;
      case BENCHMARK_HIST_INT_TABLE:
        if (hist->n <= UINT8_MAX) {
          if (create_table_int_uint8_t(hist)) {
            destroy_hist(hist);
            return EXIT_FAILURE;
          }
          hist_func = hist_table_int_uint8_t;
        }
        else if (hist->n <= UINT16_MAX) {
          if (create_table_int_uint16_t(hist)) {
            destroy_hist(hist);
            return EXIT_FAILURE;
          }
          hist_func = hist_table_int_uint16_t;
        }
        else {
          if (create_table_int_uint32_t(hist)) {
            destroy_hist(hist);
            return EXIT_FAILURE;
          }
          hist_func = hist_table_int_uint32_t;
        }
        break;
      case BENCHMARK_HIST_LOG_TRUNC: hist_func = hist_log_trunc; break;
      case BENCHMARK_HIST_FASTLOG2_START + 1:
        hist_func = hist_fast_log2_1; break;
      case BENCHMARK_HIST_FASTLOG2_START + 2:
        hist_func = hist_fast_log2_2; break;
      case BENCHMARK_HIST_FASTLOG2_START + 3:
        hist_func = hist_fast_log2_3; break;
      case BENCHMARK_HIST_FASTLOG2_START + 4:
        hist_func = hist_fast_log2_4; break;
      case BENCHMARK_HIST_FASTLOG2_START + 5:
        hist_func = hist_fast_log2_5; break;
      default:
        P_ERR("unexpected distance histogram update method: %d\n", method);
        destroy_hist(hist);
        return EXIT_FAILURE;
    }
  }

  /* Allocate memory for the squared distance samples. */
  size_t csize = (nsp < BENCHMARK_HIST_DIST_CHUNK) ? nsp :
      BENCHMARK_HIST_DIST_CHUNK;
#if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE
  /* Round up the size of the sample array for memory alignment. */
  size_t tsize = csize * sizeof(real);
  tsize = (tsize + BENCHMARK_MEMALIGN_BYTE - 1) &
      (~(BENCHMARK_MEMALIGN_BYTE - 1));
  real *d2;
  if (posix_memalign((void **) &d2, BENCHMARK_MEMALIGN_BYTE, tsize)) {
    P_ERR("failed to allocate aligned memory for squared distance samples\n");
    destroy_hist(hist);
    return EXIT_FAILURE;
  }
  /* Initialise the extra sample points. */
  size_t num = tsize / sizeof(real);
  for (size_t i = csize; i < num; i++) d2[i] = 0;
#else
  real *d2 = malloc(csize * sizeof(real));
  if (!d2) {
    P_ERR("failed to allocate memory for squared distance samples\n");
    destroy_hist(hist);
    return EXIT_FAILURE;
  }
#endif

  /* Generate samples and measure the computing time. */
  clock_t elapsed = 0;
  size_t nchunk = nsp / BENCHMARK_HIST_DIST_CHUNK;
  if (nsp % BENCHMARK_HIST_DIST_CHUNK != 0) nchunk++;
  for (size_t i = 0; i < nchunk; i++) {
    const size_t num = (i != nchunk - 1) ? BENCHMARK_HIST_DIST_CHUNK :
        nsp - i * BENCHMARK_HIST_DIST_CHUNK;
    sample_dist(d2, smax, num, BENCHMARK_RAND_SEED + i);

    /* Preprocess the input squared distance samples if necessary. */
    if (rescale != 1) {
      for (size_t j = 0; j < num; j++) d2[j] = d2[j] * rescale;
    }

    clock_t start = clock();
    hist_func(d2, num, hist);
    elapsed += clock() - start;
  }

#if     BENCHMARK_SIMD  ==  BENCHMARK_SIMD_AVX512
  /* Reduce private histograms. */
  const size_t npar = hist->n & 0xFFFFFFFFFFFFFFF8ULL;
  for (size_t i = 0; i != npar; i += 8) {
    size_t *cnt = hist->cnt + i;
    __m512i v = _mm512_load_si512(cnt);
    for (int j = hist->n; j < hist->n * 8; j += hist->n) {
      __m512i t = _mm512_loadu_si512(cnt + j);
      v = _mm512_add_epi64(v, t);
    }
    _mm512_store_si512(cnt, v);
  }
  for (size_t i = npar; i < hist->n; i++) {
    for (int j = hist->n; j < hist->n * 8; j += hist->n) {
      hist->cnt[i] += hist->cnt[i + j];
    }
  }
#endif

  double sec = (double) elapsed / CLOCKS_PER_SEC;

  printf(FMT_DONE);

  printf("Time for distance histogram: " OFMT_DBL " seconds\n", sec);
#ifdef PRINT_HISTOGRAM
  print_hist(hist);
#endif

  free(d2);
  destroy_hist(hist);
  return 0;
}

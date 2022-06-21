/*******************************************************************************
* benchmark/histogram/benchmark.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "dist_hist.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>     /* for getopt */

#if _POSIX_C_SOURCE < 200809L
  #error Obsolete POSIX version detected. \
    Please include `-D_POSIX_C_SOURCE=200809` in your compilation flags.
#endif

void help(void) {
  printf("Usage: " CODE_NAME " [OPTION]\n"
      CODE_NAME ": Benchmark program for distance histogram update\n\
  -h    Display this message and exit\n\
  -n    " FMT_KEY(NUM_DIST) "\n\
        Set the total number of distances to be sampled\n\
  -d    " FMT_KEY(MAX_DIST) "\n\
        Set the largest squared distance to be sampled\n\
  -l    " FMT_KEY(R_MIN) "\n\
        Set the minimum separation of interest (default: 0)\n\
  -r    " FMT_KEY(R_MAX) "\n\
        Set the maximum separation of interest\n\
  -b    " FMT_KEY(NUM_BIN) "\n\
        Set the number of separation bins\n\
  -t    " FMT_KEY(BIN_TYPE) "\n\
        Specify the type of distance bins, allowed values are\n\
        %d: linear bin\n\
        %d: logarithm bin\n\
  -m    " FMT_KEY(METHOD) "\n\
        Specify the histogram update method, allowed values are\n\
        %d: binary search (both bin types)\n\
        %d: reverse traversal (both bin types)\n\
        %d: preprocess data and lookup bins with integer edges (linear bin)\n\
        %d: preprocess data and map indices with sqrt (linear bin)\n\
        %d: preprocess data and map indices with log (logarithm bin)\n\
        %d+n: preprocess data and map indices with approximate sqrt\n\
            using a polynomial of order n (linear bin)\n\
        %d+n: preprocess data and map indices with approximate log2\n\
            using a polynomial of order n (logarithm bin)\n\
        negative: hybrid lookup with the given table size (both bin types)\n",
    BENCHMARK_HIST_BIN_LIN, BENCHMARK_HIST_BIN_LOG, BENCHMARK_HIST_BSEARCH,
    BENCHMARK_HIST_REV_TRAV, BENCHMARK_HIST_INT_TABLE,
    BENCHMARK_HIST_SQRT_TRUNC, BENCHMARK_HIST_LOG_TRUNC,
    BENCHMARK_HIST_FASTSQRT_START, BENCHMARK_HIST_FASTLOG2_START);
}

int main(int argc, char *argv[]) {
  int opt, nbin, btype, method;
  size_t nsp;
  double smax, dmin, dmax;

  /* Initialise configurations and read settings from command line options. */
  nsp = 0;
  nbin = btype = method = -1;
  smax = dmax = -1;
  dmin = 0;
  while ((opt = getopt(argc, argv, "hn:d:l:r:b:t:m:")) != -1) {
    switch (opt) {
      case 'h': help(); return 0;
      case 'n':
        if (sscanf(optarg, "%zu", &nsp) != 1) {
          P_ERR("failed to parse " FMT_KEY(NUM_DIST) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'd':
        if (sscanf(optarg, "%lf", &smax) != 1) {
          P_ERR("failed to parse " FMT_KEY(MAX_DIST) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'l':
        if (sscanf(optarg, "%lf", &dmin) != 1) {
          P_ERR("failed to parse " FMT_KEY(R_MIN) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'r':
        if (sscanf(optarg, "%lf", &dmax) != 1) {
          P_ERR("failed to parse " FMT_KEY(R_MAX) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'b':
        if (sscanf(optarg, "%d", &nbin) != 1) {
          P_ERR("failed to parse " FMT_KEY(NUM_BIN) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 't':
        if (sscanf(optarg, "%d", &btype) != 1) {
          P_ERR("failed to parse " FMT_KEY(BIN_TYPE) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'm':
        if (sscanf(optarg, "%d", &method) != 1) {
          P_ERR("failed to parse " FMT_KEY(METHOD) "\n");
          return EXIT_FAILURE;
        }
        break;
      default: break;
    }
  }
  if (optind < argc) {
    P_WRN("unknown command line options:\n ");
    while (optind < argc) printf(" %s", argv[optind++]);
    printf("\n");
  }

  /* Verify configurations. */
  if (nsp == 0 || nsp >= BENCHMARK_HIST_MAX_NSP) {
    P_ERR(FMT_KEY(NUM_DIST) " must be positive and smaller than %zu\n",
        (size_t) BENCHMARK_HIST_MAX_NSP);
    return EXIT_FAILURE;
  }
  if (smax <= 0) {
    P_ERR(FMT_KEY(MAX_DIST) " must be positive\n");
    return EXIT_FAILURE;
  }
  if (dmin < 0 || dmax < 0) {
    P_ERR(FMT_KEY(R_MIN) " and " FMT_KEY(R_MAX) " must be non-negative\n");
    return EXIT_FAILURE;
  }
  if (dmax <= dmin) {
    P_ERR(FMT_KEY(R_MAX) " must be larger than " FMT_KEY(R_MIN) "\n");
    return EXIT_FAILURE;
  }
  if (nbin <= 0 || nbin >= BENCHMARK_HIST_MAX_NBIN) {
    P_ERR(FMT_KEY(NUM_BIN) " must be positive and smaller than %d\n",
        BENCHMARK_HIST_MAX_NBIN);
    return EXIT_FAILURE;
  }
  switch (btype) {
    case BENCHMARK_HIST_BIN_LIN:
      break;
    case BENCHMARK_HIST_BIN_LOG:
      if (dmin == 0) {
        P_ERR(FMT_KEY(R_MIN) " must be non-zero for logarithm bins\n");
        return EXIT_FAILURE;
      }
      break;
    default:
      P_ERR("invalid " FMT_KEY(BIN_TYPE) ": %d\n", btype);
      return EXIT_FAILURE;
  }
  if (method < 0) {
    if (method > - BENCHMARK_HIST_TABLE_MIN_SIZE ||
        method < - BENCHMARK_HIST_TABLE_MAX_SIZE) {
      P_ERR("invalid number of lookup table entries: %d\n", -method);
      return EXIT_FAILURE;
    }
  }
  else if (method >= BENCHMARK_HIST_FASTSQRT_START +
      BENCHMARK_HIST_APPROX_POLY_MIN && method <= BENCHMARK_HIST_FASTSQRT_START
      + BENCHMARK_HIST_APPROX_POLY_MAX) {
    if (btype != BENCHMARK_HIST_BIN_LIN) {
      P_ERR("the histogram update method works only for linear bins\n");
      return EXIT_FAILURE;
    }
  }
  else if (method >= BENCHMARK_HIST_FASTLOG2_START +
      BENCHMARK_HIST_APPROX_POLY_MIN && method <= BENCHMARK_HIST_FASTLOG2_START
      + BENCHMARK_HIST_APPROX_POLY_MAX) {
    if (btype != BENCHMARK_HIST_BIN_LOG) {
      P_ERR("the histogram update method works only for logarithmic bins\n");
      return EXIT_FAILURE;
    }
  }
  else {
    switch (method) {
      case BENCHMARK_HIST_BSEARCH:
      case BENCHMARK_HIST_REV_TRAV: break;
      case BENCHMARK_HIST_SQRT_TRUNC:
      case BENCHMARK_HIST_INT_TABLE:
        if (btype != BENCHMARK_HIST_BIN_LIN) {
          P_ERR("the histogram update method works only for linear bins\n");
          return EXIT_FAILURE;
        }
        break;
      case BENCHMARK_HIST_LOG_TRUNC:
        if (btype != BENCHMARK_HIST_BIN_LOG) {
          P_ERR("the histogram update method works only for logarithm bins\n");
          return EXIT_FAILURE;
        }
        break;
      default:
        P_ERR("invalid " FMT_KEY(METHOD) ": %d\n", method);
        return EXIT_FAILURE;
    }
  }

  /* Print configurations. */
  const char *benchmark_btype_name[] = {"linear", "log"};
  const char *benchmark_method_name[] = {"bsearch", "reverse traversal",
      "lookup table", "index mapping with sqrt", "index mapping with log"};
  printf(FMT_KEY(NUM_DIST) ": %zu    " FMT_KEY(MAX_DIST) ": " OFMT_DBL "    "
      FMT_KEY(R_MIN) ": " OFMT_DBL "    " FMT_KEY(R_MAX) ": " OFMT_DBL "\n"
      FMT_KEY(NUM_BIN) ": %d    " FMT_KEY(BIN_TYPE) ": %s    \n",
      nsp, smax, dmin, dmax, nbin, benchmark_btype_name[btype]);
  if (method >= BENCHMARK_HIST_FASTSQRT_START + BENCHMARK_HIST_APPROX_POLY_MIN
      && method <= BENCHMARK_HIST_FASTSQRT_START +
      BENCHMARK_HIST_APPROX_POLY_MAX) {
    printf(FMT_KEY(METHOD) ": FastSqrt (order %d)\n\n",
        method - BENCHMARK_HIST_FASTSQRT_START);
  }
  else if (method >= BENCHMARK_HIST_FASTLOG2_START +
      BENCHMARK_HIST_APPROX_POLY_MIN && method <= BENCHMARK_HIST_FASTLOG2_START
      + BENCHMARK_HIST_APPROX_POLY_MAX) {
    printf(FMT_KEY(METHOD) ": FastLog2 (order %d)\n\n",
        method - BENCHMARK_HIST_FASTLOG2_START);
  }
  else printf(FMT_KEY(METHOD) ": %s\n\n",
      method >= 0 ? benchmark_method_name[method] : "hybrid lookup table");

  /* Evaluate the distance histogram. */
  if (eval_hist(smax, nsp, method, dmin, dmax, nbin, btype))
    return EXIT_FAILURE;

  return 0;
}

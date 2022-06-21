/*******************************************************************************
* benchmark/struct/benchmark.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "create_data.h"
#include "data_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>     /* for getopt */

#if _POSIX_C_SOURCE < 200809L
  #error Obsolete POSIX version detected. \
    Please include `-D_POSIX_C_SOURCE=200809` in your compilation flags.
#endif

void help(void) {
  printf("Usage: " CODE_NAME " [OPTION]\n"
      CODE_NAME ": Benchmark program for pair counting data structures\n\
  -h    Display this message and exit\n\
  -b    " FMT_KEY(BOXSIZE) "\n\
        Set the side length of the periodic box\n\
  -n    " FMT_KEY(NUM_PART) "\n\
        Set the total number of particles to be sampled\n\
  -s    Use data only in a shell, with inner and outer radii being "
        FMT_KEY(BOXSIZE) "/2\n\
        and " FMT_KEY(BOXSIZE) " respectively, and perform non-periodic "
        "distance evaluations\n\
  -l    " FMT_KEY(R_MIN) "\n\
        Set the minimum separation of interest (default: 0)\n\
  -r    " FMT_KEY(R_MAX) "\n\
        Set the maximum separation of interest, which must be smaller than\n\
        " FMT_KEY(BOXSIZE) "/2\n\
  -d    " FMT_KEY(DATA_STRUCT) "\n\
        Specify the pair counting data structure, allowed values are\n\
        %d: regular grids\n\
        %d: k-d tree\n\
        %d: ball tree\n\
  -N    " FMT_KEY(NUM_CELL) "\n\
        Number of cells per box size for regular grids,\n\
        or the maximum number of particles of leaf nodes for the trees\n",
    BENCHMARK_STRUCT_GRID, BENCHMARK_STRUCT_KDTREE, BENCHMARK_STRUCT_BALLTREE);
}

int main(int argc, char *argv[]) {
  int opt, dtype, csize;
  size_t npar;
  bool box = true;
  double bsize, rmin, rmax;

  /* Initialise configurations and read settings from command line options. */
  npar = 0;
  dtype = csize = -1;
  bsize = rmax = -1;
  rmin = 0;
  while ((opt = getopt(argc, argv, "hsb:n:l:r:d:N:")) != -1) {
    switch (opt) {
      case 'h': help(); return 0;
      case 'b':
        if (sscanf(optarg, "%lf", &bsize) != 1) {
          P_ERR("failed to parse " FMT_KEY(BOXSIZE) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'n':
        if (sscanf(optarg, "%zu", &npar) != 1) {
          P_ERR("failed to parse " FMT_KEY(NUM_PART) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 's': box = false; break;
      case 'l':
        if (sscanf(optarg, "%lf", &rmin) != 1) {
          P_ERR("failed to parse " FMT_KEY(R_MIN) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'r':
        if (sscanf(optarg, "%lf", &rmax) != 1) {
          P_ERR("failed to parse " FMT_KEY(R_MAX) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'd':
        if (sscanf(optarg, "%d", &dtype) != 1) {
          P_ERR("failed to parse " FMT_KEY(DATA_STRUCT) "\n");
          return EXIT_FAILURE;
        }
        break;
      case 'N':
        if (sscanf(optarg, "%d", &csize) != 1) {
          P_ERR("failed to parse " FMT_KEY(NUM_CELL) "\n");
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
  if (npar < BENCHMARK_MIN_NPAR || npar > BENCHMARK_MAX_NPAR) {
    P_ERR(FMT_KEY(NUM_PART) " must be between %d and %zu\n",
        BENCHMARK_MIN_NPAR, (size_t) BENCHMARK_MAX_NPAR);
    return EXIT_FAILURE;
  }
  if (bsize <= 0) {
    P_ERR(FMT_KEY(BOXSIZE) " must be positive\n");
    return EXIT_FAILURE;
  }
  if (rmin < 0 || rmax < 0) {
    P_ERR(FMT_KEY(R_MIN) " and " FMT_KEY(R_MAX) " must be non-negative\n");
    return EXIT_FAILURE;
  }
  if (rmax <= rmin || rmax * 2 >= bsize) {
    P_ERR(FMT_KEY(R_MAX) " must be larger than " FMT_KEY(R_MIN) " and "
        " smaller than " FMT_KEY(BOXSIZE) "/2\n");
    return EXIT_FAILURE;
  }
  if (csize < BENCHMARK_MIN_CSIZE || csize > BENCHMARK_MAX_CSIZE) {
    P_ERR(FMT_KEY(NUM_CELL) " must be between %d and %d\n",
        BENCHMARK_MIN_CSIZE, BENCHMARK_MAX_CSIZE);
    return EXIT_FAILURE;
  }
  if (dtype < 0) {
    P_ERR(FMT_KEY(DATA_STRUCT) " must be non-negative\n");
    return EXIT_FAILURE;
  }
  switch (dtype) {
    case BENCHMARK_STRUCT_GRID:
    case BENCHMARK_STRUCT_KDTREE:
    case BENCHMARK_STRUCT_BALLTREE:
      break;
    default:
      P_ERR("invalid " FMT_KEY(DATA_STRUCT) ": %d\n", dtype);
      return EXIT_FAILURE;
  }

  /* Print configurations. */
  const char *benchmark_struct_name[] = {"grid", "kd-tree", "ball tree"};
  printf(FMT_KEY(BOXSIZE) ": " OFMT_DBL "    " FMT_KEY(NUM_PART) ": %zu    "
      FMT_KEY(BOX_TYPE) ": %s\n" FMT_KEY(R_MIN) ": " OFMT_DBL "    "
      FMT_KEY(R_MAX) ": " OFMT_DBL "\n" FMT_KEY(DATA_STRUCT) ": %s    "
      FMT_KEY(NUM_CELL) ": %d\n\n",
      bsize, npar, box ? "box" : "shell", rmin, rmax,
      benchmark_struct_name[dtype], csize);

  /* Generate data catalogue */
  DATA *data;
  if (!(data = create_data(bsize, npar, box)))  return EXIT_FAILURE;

#ifdef DEBUG
  double r2min = rmin * rmin;
  double r2max = rmax * rmax;
  size_t cnt = 0;
  const double halfb = bsize * 0.5;
  if (box) {
    for (size_t i = 0; i < data->n; i++) {
      for (size_t j = 0; j < data->n; j++) {
        double dx = data->x[i] - data->x[j];
        if (dx > halfb) dx -= bsize;
        else if (dx < -halfb) dx += bsize;
        double dy = data->y[i] - data->y[j];
        if (dy > halfb) dy -= bsize;
        else if (dy < -halfb) dy += bsize;
        double dz = data->z[i] - data->z[j];
        if (dz > halfb) dz -= bsize;
        else if (dz < -halfb) dz += bsize;
        double dist = dx * dx + dy * dy + dz * dz;
        if (dist >= r2min && dist < r2max) cnt++;
      }
    }
  }
  else {
    for (size_t i = 0; i < data->n; i++) {
      for (size_t j = 0; j < data->n; j++) {
        double dx = data->x[i] - data->x[j];
        double dy = data->y[i] - data->y[j];
        double dz = data->z[i] - data->z[j];
        double dist = dx * dx + dy * dy + dz * dz;
        if (dist >= r2min && dist < r2max) cnt++;
      }
    }
  }
  printf("Number of pairs from brute-force counting: %zu\n\n", cnt);
#endif

  /* Evaluate pair counts. */
#ifdef BENCHMARK_TIMING
  size_t npair;
#else
  size_t npair, nnode, ndist;
#endif
  switch (dtype) {
    case BENCHMARK_STRUCT_GRID:
#ifdef BENCHMARK_TIMING
      if (paircount_grid(data, rmin, rmax, csize, &npair))
#else
      if (paircount_grid(data, rmin, rmax, csize, &nnode, &ndist, &npair))
#endif
      {
        destroy_data(data);
        return EXIT_FAILURE;
      }
      break;
    case BENCHMARK_STRUCT_KDTREE:
#ifdef BENCHMARK_TIMING
      if (paircount_kdtree(data, rmin, rmax, csize, &npair))
#else
      if (paircount_kdtree(data, rmin, rmax, csize, &nnode, &ndist, &npair))
#endif
      {
        destroy_data(data);
        return EXIT_FAILURE;
      }
      break;
    case BENCHMARK_STRUCT_BALLTREE:
#ifdef BENCHMARK_TIMING
      if (paircount_balltree(data, rmin, rmax, csize, &npair))
#else
      if (paircount_balltree(data, rmin, rmax, csize, &nnode, &ndist, &npair))
#endif
      {
        destroy_data(data);
        return EXIT_FAILURE;
      }
      break;
    default:
      P_ERR("invalid " FMT_KEY(DATA_STRUCT) ": %d\n", dtype);
      destroy_data(data);
      return EXIT_FAILURE;
  }

  printf("Number of pairs found: %zu\n", npair);
#ifndef BENCHMARK_TIMING
  printf("Number of nodes visited: %zu\n"
      "Number of pair distance evaluations: %zu\n", nnode, ndist);
#endif

  destroy_data(data);
  return 0;
}

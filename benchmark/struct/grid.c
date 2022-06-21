/*******************************************************************************
* benchmark/struct/grid.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "define.h"
#include "data_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/*============================================================================*\
                       Definition of the data structures
\*============================================================================*/

/* Data structure for grid cells. */
typedef struct {
  real *x;                      /* pointer to x coordinates */
  real *y;                      /* pointer to y coordinates */
  real *z;                      /* pointer to z coordinates */
  size_t n;                     /* number of data points    */
} GRID;

/* Data structure for offsets of cells to be visited. */
typedef struct {
  int *i;                       /* pointer to offsets of x index */
  int *j;                       /* pointer to offsets of y index */
  int *k;                       /* pointer to offsets of z index */
  size_t n;                     /* number of index offsets       */
} OFFSET;


/*============================================================================*\
                   Functions for data structure construction
\*============================================================================*/

/******************************************************************************
Function `create_grid`:
  Construct grids for a dataset.
Arguments:
  * `data`:     the input dataset;
  * `ngrid`:    number of grids per box side.
Return:
  Address of the structure for grids on  success; NULL on error.
******************************************************************************/
static GRID *create_grid(DATA *data, const int ngrid) {
  /* Allocate memory. */
  const size_t size = (size_t) ngrid * ngrid * ngrid;
  GRID *grid = calloc(size, sizeof(GRID));
  if (!grid) {
    P_ERR("failed to allocate memory for the grids\n");
    return NULL;
  }
  size_t *gidx = calloc(size, sizeof(size_t));
  if (!gidx) {
    P_ERR("failed to allocate memory for assigning data to grids\n");
    free(grid);
    return NULL;
  }

  /* Count the number of points for each cell. */
  const double fac = ngrid / data->bsize;
  for (size_t i = 0; i < data->n; i++) {
    const int ix = data->x[i] * fac;
    const int iy = data->y[i] * fac;
    const int iz = data->z[i] * fac;
    if (ix < 0 || ix >= ngrid || iy < 0 || iy >= ngrid ||
        iz < 0 || iz >= ngrid) {
      P_ERR("input coordinates outside the periodic box: (" OFMT_DBL ", "
          OFMT_DBL ", " OFMT_DBL ")\n", data->x[i], data->y[i], data->z[i]);
      free(grid); free(gidx);
      return NULL;
    }

    const size_t idx = ((size_t) ix * ngrid + iy) * ngrid + iz;
    grid[idx].n++;
  }

  /* Set pointers to the input dataset for grid cells. */
  size_t offset = 0;
  for (size_t i = 0; i < size; i++) {
    grid[i].x = data->x + offset;
    grid[i].y = data->y + offset;
    grid[i].z = data->z + offset;
    offset += grid[i].n;
  }

  /* Sort the input dataset. */
  for (size_t i = 0; i < size; i++) {
    for (size_t n = gidx[i]; n < grid[i].n; n++) {
      /* Find the right grid cell for the current data point. */
      const int ix = grid[i].x[n] * fac;
      const int iy = grid[i].y[n] * fac;
      const int iz = grid[i].z[n] * fac;
      const size_t j = ((size_t) ix * ngrid + iy) * ngrid + iz;
      /* Swap data points if necessary. */
      if (i != j) {
        real tmp = grid[i].x[n];
        grid[i].x[n] = grid[j].x[gidx[j]];
        grid[j].x[gidx[j]] = tmp;
        tmp = grid[i].y[n];
        grid[i].y[n] = grid[j].y[gidx[j]];
        grid[j].y[gidx[j]] = tmp;
        tmp = grid[i].z[n];
        grid[i].z[n] = grid[j].z[gidx[j]];
        grid[j].z[gidx[j]] = tmp;
        gidx[j]++;
        n--;
      }
    }
  }

  free(gidx);
  return grid;
}

/******************************************************************************
Function `destroy_offset`:
  Deconstruct pre-computed offsets of cells to be visited.
Arguments:
  * `off`:      structure for storing offsets.
******************************************************************************/
static void destroy_offset(OFFSET *off) {
  if (!off) return;
  if (off->i) free(off->i);
  if (off->j) free(off->j);
  if (off->k) free(off->k);
  free(off);
}

/******************************************************************************
Function `offset_resize`:
  Enlarge the memory allocated for index offsets.
Arguments:
  * `off`:      structure for storing offsets;
  * `capacity`: capacity of the structure.
Return:
  Zero on success; EXIT_FAILURE on error.
******************************************************************************/
static int offset_resize(OFFSET *off, size_t *capacity) {
  if (*capacity > off->n) return 0;
  else if (!(*capacity)) {
    *capacity = BENCHMARK_STACK_INIT_SIZE;
    if (!(off->i = malloc(*capacity * sizeof(int))) ||
        !(off->j = malloc(*capacity * sizeof(int))) ||
        !(off->k = malloc(*capacity * sizeof(int)))) {
      P_ERR("failed to allocate memory for grid index offsets\n");
      destroy_offset(off);
      return EXIT_FAILURE;
    }
  }
  else {
    if (*capacity >= SIZE_MAX / 2) {
      P_ERR("not enough memory for grid index offsets\n");
      destroy_offset(off);
      return EXIT_FAILURE;
    }
    *capacity <<= 1;
    int *tmp;
    if ((tmp = realloc(off->i, *capacity * sizeof(int)))) {
      off->i = tmp;
      if ((tmp = realloc(off->j, *capacity * sizeof(int)))) {
        off->j = tmp;
        if ((tmp = realloc(off->k, *capacity * sizeof(int)))) {
          off->k = tmp;
          return 0;
        }
      }
    }
    P_ERR("failed to enlarge memory allocated for grid index offsets\n");
    destroy_offset(off);
    return EXIT_FAILURE;
  }
  return 0;
}

/******************************************************************************
Function `compare_offset`:
  Compare two index offsets for sorting, based on their location in memory.
Arguments:
  * `a`:        pointer to the first offset;
  * `b`:        pointer to the second offset.
Return:
  positive if a < b, negative if a > b, 0 if a = b.
******************************************************************************/
static int compare_offset(const void *a, const void *b) {
  return *((long *) a) - *((long *) b);
}

/******************************************************************************
Function `compute_offset`:
  Compute offsets of the indices of cells to be visited.
Arguments:
  * `bsize`:    side length of the box;
  * `ngrid`:    number of grids per box side;
  * `rmin`:     minimum separation of interest;
  * `rmax`:     maximum separation of interest;
  * `ndist`:    number of distance evaluations while computing offsets.
Return:
  Address of the structure for index offsets on success; NULL on error.
******************************************************************************/
static OFFSET *compute_offset(const double bsize, const int ngrid,
    const double rmin, const double rmax
#ifndef BENCHMARK_TIMING
    , size_t *ndist
#endif
    ) {
  /* Allocate memory. */
  OFFSET *off = calloc(1, sizeof(OFFSET));
  off->i = off->j = off->k = NULL;
  size_t capacity = 0;
  if (offset_resize(off, &capacity)) return NULL;

  const double fac = ngrid / bsize;
  const double csize = bsize / ngrid;
  /* Minimum offset: floor(rmin / sqrt(3) * ngrid / bsize) */
  const double rmini = rmin * 0x1.279a74590331cp-1;
  const int offmin = (rmini <= bsize) ? 0 : (rmini - bsize) * fac;
  /* Maximum offset: ceil(rmax * ngrid / bsize) */
  int offmax = ceil(rmax * fac);        /* no overflow due to limit of ngrid*/
  if (offmax >= ngrid) offmax = ngrid - 1;

  const double r2min = rmin * rmin;
  const double r2max = rmax * rmax;
#ifndef BENCHMARK_TIMING
  size_t cnt = 0;
#endif

  /* Checking only the first octant is enough. */
  for (int i = offmin; i <= offmax; i++) {
    double d2min_i = (i == 0) ? 0 : (i - 1) * csize;
    d2min_i *= d2min_i;
    double d2max_i = (i + 1) * csize;
    d2max_i *= d2max_i;
#ifndef BENCHMARK_TIMING
    cnt++;
#endif
    for (int j = offmin; j <= offmax; j++) {
      double d2min_ij = (j == 0) ? 0 : (j - 1) * csize;
      d2min_ij *= d2min_ij;
      d2min_ij += d2min_i;
      double d2max_ij = (j + 1) * csize;
      d2max_ij *= d2max_ij;
      d2max_ij += d2max_i;
#ifndef BENCHMARK_TIMING
      cnt++;
#endif
      for (int k = offmin; k <= offmax; k++) {
        double d2min = (k == 0) ? 0 : (k - 1) * csize;
        d2min *= d2min;
        d2min += d2min_ij;
        double d2max = (k + 1) * csize;
        d2max *= d2max;
        d2max += d2max_ij;
#ifndef BENCHMARK_TIMING
        cnt++;
#endif
        if (d2max < r2min) continue;
        else if (d2min >= r2max) break; /* d2min only increases in k loop */

        /* Save the current offsets and counterparts in the other octants. */
        if (i) {
          if (off->n >= capacity && offset_resize(off, &capacity)) return NULL;
          off->i[off->n] = -i;
          off->j[off->n] = j;
          off->k[off->n] = k;
          off->n++;
          if (j) {
            if (off->n >= capacity && offset_resize(off, &capacity))
              return NULL;
            off->i[off->n] = -i;
            off->j[off->n] = -j;
            off->k[off->n] = k;
            off->n++;
            if (k) {
              if (off->n >= capacity && offset_resize(off, &capacity))
                return NULL;
              off->i[off->n] = -i;
              off->j[off->n] = -j;
              off->k[off->n] = -k;
              off->n++;
            }
          }
          if (k) {
            if (off->n >= capacity && offset_resize(off, &capacity))
              return NULL;
            off->i[off->n] = -i;
            off->j[off->n] = j;
            off->k[off->n] = -k;
            off->n++;
          }
        }
        if (j) {
          if (off->n >= capacity && offset_resize(off, &capacity)) return NULL;
          off->i[off->n] = i;
          off->j[off->n] = -j;
          off->k[off->n] = k;
          off->n++;
          if (k) {
            if (off->n >= capacity && offset_resize(off, &capacity))
              return NULL;
            off->i[off->n] = i;
            off->j[off->n] = -j;
            off->k[off->n] = -k;
            off->n++;
          }
        }
        if (k) {
          if (off->n >= capacity && offset_resize(off, &capacity)) return NULL;
          off->i[off->n] = i;
          off->j[off->n] = j;
          off->k[off->n] = -k;
          off->n++;
        }
        if (off->n >= capacity && offset_resize(off, &capacity)) return NULL;
        off->i[off->n] = i;
        off->j[off->n] = j;
        off->k[off->n] = k;
        off->n++;
      }
    }
  }
  if (!off->n) {
    P_ERR("no cell is found to be inside the separation range\n");
    destroy_offset(off);
    return NULL;
  }
#ifndef BENCHMARK_TIMING
  *ndist = cnt / 3;
#endif

  /* Sort offsets based on the location in memory, to reduce cache misses. */
  long *idx = malloc(off->n * sizeof(long));
  if (!idx) {
    P_ERR("failed to allocate memory for sorting index offsets\n");
    destroy_offset(off);
    return NULL;
  }

  /* Convert indices to memory locations, sort, and convert back. */
  const int shift = ngrid >> 1;
  for (size_t n = 0; n < off->n; n++) {
    const int i = off->i[n] + shift;
    const int j = off->j[n] + shift;
    const int k = off->k[n] + shift;
    idx[n] = ((long) i * ngrid + j) * ngrid + k;
  }

  qsort(idx, off->n, sizeof(long), compare_offset);

  for (size_t n = 0; n < off->n; n++) {
    const int k = idx[n] % ngrid;
    idx[n] = (idx[n] - k) / ngrid;
    const int j = idx[n] % ngrid;
    const int i = (idx[n] - j) / ngrid;
    off->i[n] = i - shift;
    off->j[n] = j - shift;
    off->k[n] = k - shift;
  }
  free(idx);

  /* Reduce the memory cost if applicable. */
  if (off->n < capacity) {
    int *tmp = realloc(off->i, off->n * sizeof(int));
    if (tmp) off->i = tmp;
    tmp = realloc(off->j, off->n * sizeof(int));
    if (tmp) off->j = tmp;
    tmp = realloc(off->k, off->n * sizeof(int));
    if (tmp) off->k = tmp;
  }

  return off;
}

/*============================================================================*\
                     Pair counting functions from templates
\*============================================================================*/

/* Clean all the relevant macros first */
#ifdef BENCHMARK_PAIRCOUNT_TYPE
  #undef BENCHMARK_PAIRCOUNT_TYPE
#endif
#ifdef BENCHMARK_BIN_SMIN
  #undef BENCHMARK_BIN_SMIN
#endif

/* paircount_box_smin0 */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_BOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_ZERO
#include "count_cell.c"

/* paircount_box */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_BOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_NONZERO
#include "count_cell.c"

/* paircount_smin0 */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_NOBOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_ZERO
#include "count_cell.c"

/* paircount */
#define BENCHMARK_PAIRCOUNT_TYPE        BENCHMARK_PAIRCOUNT_NOBOX
#define BENCHMARK_BIN_SMIN              BENCHMARK_BIN_MIN_NONZERO
#include "count_cell.c"

/*============================================================================*\
                          Interface for pair counting
\*============================================================================*/

/******************************************************************************
Function `paircount_grid`:
  Count pairs using grids.
Arguments:
  * `data`:     structure for the data catalogue;
  * `rmin`:     minimum separation of interest;
  * `rmax`:     maximum separation of interest;
  * `csize`:    number of grids per side, or number of partilces in leaf nodes;
  * `nnode`:    number of distance evaluations between nodes;
  * `ndist`:    number of distance evaluations performed in total;
  * `npair`:    number of pairs in the separation range of interest.
Return:
  Zero on success; EXIT_FAILURE on error.
******************************************************************************/
int paircount_grid(DATA *data, const double rmin, const double rmax,
    const int csize,
#ifndef BENCHMARK_TIMING
    size_t *nnode, size_t *ndist,
#endif
    size_t *npair) {
  if(data->isbox && csize <= 2) {
    P_ERR(FMT_KEY(NUM_CELL) " must be larger than " FMT_KEY(BOXSIZE)
        "/2 for grid based pair counting with periodic boundaries\n");
    return EXIT_FAILURE;
  }

  printf("Constructing the data structures ...");
  fflush(stdout);

  GRID *grid = create_grid(data, csize);
  if (!grid) return EXIT_FAILURE;

  const double r2min = rmin * rmin;
  const double r2max = rmax * rmax;
  *npair = 0;

#ifdef BENCHMARK_TIMING
  OFFSET *off = compute_offset(data->bsize, csize, rmin, rmax);
#else
  *nnode = *ndist = 0;
  OFFSET *off = compute_offset(data->bsize, csize, rmin, rmax, ndist);
  if (rmin != 0) *ndist *= 2;
#endif
  if (!off) {
    free(grid);
    return EXIT_FAILURE;
  }

  printf(FMT_DONE "Start pair counting ...");
  fflush(stdout);

#ifdef BENCHMARK_TIMING
  clock_t start = clock(), diff;
  if (data->isbox) {
    if (rmin == 0)
      paircount_box_smin0(grid, off, csize, data->bsize, r2max, npair);
    else paircount_box(grid, off, csize, data->bsize, r2min, r2max, npair);
  }
  else {
    if (rmin == 0) paircount_smin0(grid, off, csize, r2max, npair);
    else paircount(grid, off, csize, r2min, r2max, npair);
  }
#else
  if (data->isbox) {
    if (rmin == 0) {
      paircount_box_smin0(grid, off, csize, data->bsize, r2max,
          nnode, ndist, npair);
    }
    else {
      paircount_box(grid, off, csize, data->bsize, r2min, r2max,
          nnode, ndist, npair);
    }
  }
  else {
    if (rmin == 0)
      paircount_smin0(grid, off, csize, r2max, nnode, ndist, npair);
    else paircount(grid, off, csize, r2min, r2max, nnode, ndist, npair);
  }
#endif

#ifdef BENCHMARK_TIMING
  diff = clock() - start;
  double sec = (double) diff / CLOCKS_PER_SEC;
#endif

  destroy_offset(off);
  free(grid);

  printf(FMT_DONE);

#ifdef BENCHMARK_TIMING
  printf("Time for pair counting: " OFMT_DBL " seconds\n", sec);
#endif
  return 0;
}

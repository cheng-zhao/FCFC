/*******************************************************************************
* benchmark/struct/count_cell.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/*******************************************************************************
  Implementation of grid-based pair counting.
*******************************************************************************/

/* Macros for the template functions. */
#if defined(BENCHMARK_PAIRCOUNT_TYPE) && defined(BENCHMARK_BIN_SMIN)

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef BENCHMARK_PAIRCOUNT_NAME
  #undef BENCHMARK_PAIRCOUNT_NAME
#endif
#ifdef BENCHMARK_BIN_SMIN_NAME
  #undef BENCHMARK_BIN_SMIN_NAME
#endif

#if     BENCHMARK_PAIRCOUNT_TYPE == BENCHMARK_PAIRCOUNT_BOX
  #define BENCHMARK_PAIRCOUNT_NAME      _box
#elif   BENCHMARK_PAIRCOUNT_TYPE == BENCHMARK_PAIRCOUNT_NOBOX
  #define BENCHMARK_PAIRCOUNT_NAME
#else
  #error unexpected definition of `BENCHMARK_PAIRCOUNT_TYPE`
#endif

#if     BENCHMARK_BIN_SMIN == BENCHMARK_BIN_MIN_ZERO
  #define BENCHMARK_BIN_SMIN_NAME       _smin0
#elif   BENCHMARK_BIN_SMIN == BENCHMARK_BIN_MIN_NONZERO
  #define BENCHMARK_BIN_SMIN_NAME
#else
  #error unexpected definition of `BENCHMARK_BIN_SMIN`
#endif

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c)           a##b##c
#endif

#ifndef PAIR_COUNT_FUNC
  #define PAIR_COUNT_FUNC(a,b,c)        CONCAT_FNAME(a,b,c)
#endif


/*============================================================================*\
                       Function for dual cell pair counts
\*============================================================================*/

/******************************************************************************
Function `paircount<BENCHMARK_PAIRCOUNT_NAME><BENCHMARK_BIN_SMIN_NAME>`:
  Count pairs of objects given the structures for grids.
Arguments:
  * `grid`:     structure for grid cells;
  * `off`:      pre-computed index offsets of cells to be visited;
  * `ngrid`:    number of grids per side;
  * `bsize`:    side length of the periodic box;
  * `dmin`:     minimum squared distance of interest;
  * `dmax`:     maximum squared distance of interest;
  * `nnode`:    number of nodes visited in total;
  * `ndist`:    number of distance evaluations performed in total;
  * `npair`:    number of pairs in the separation range of interest.
******************************************************************************/
static void PAIR_COUNT_FUNC(paircount, BENCHMARK_PAIRCOUNT_NAME,
    BENCHMARK_BIN_SMIN_NAME)
    (const GRID *grid, const OFFSET *off, const int ngrid,
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
    const real bsize,
#endif
#if     BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
    const real dmin,
#endif
    const real dmax,
#ifndef BENCHMARK_TIMING
    size_t *nnode, size_t *ndist,
#endif
    size_t *npair) {

  /* Visit all grid cells. */
  for (int i1 = 0; i1 < ngrid; i1++) {
    const size_t idx_i1 = i1 * ngrid;
    for (int j1 = 0; j1 < ngrid; j1++) {
      const size_t idx_ij1 = (idx_i1 + j1) * ngrid;
      for (int k1 = 0; k1 < ngrid; k1++) {
        const size_t idx1 = idx_ij1 + k1;
        if (!grid[idx1].n) continue;

        /* Visit second grids based on the offsets. */
        for (size_t n = 0; n < off->n; n++) {
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
          real sx, sy, sz;      /* coordinate offsets for periodic conditions */
          sx = sy = sz = 0;
#endif
          int i2 = i1 + off->i[n];
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
          if (i2 < 0) {
            i2 += ngrid;
            sx = bsize;
          }
          else if (i2 >= ngrid) {
            i2 -= ngrid;
            sx = -bsize;
          }
#else
          if (i2 < 0) continue;
          else if (i2 >= ngrid) continue;
#endif
          int j2 = j1 + off->j[n];
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
          if (j2 < 0) {
            j2 += ngrid;
            sy = bsize;
          }
          else if (j2 >= ngrid) {
            j2 -= ngrid;
            sy = -bsize;
          }
#else
          if (j2 < 0) continue;
          else if (j2 >= ngrid) continue;
#endif
          int k2 = k1 + off->k[n];
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
          if (k2 < 0) {
            k2 += ngrid;
            sz = bsize;
          }
          else if (k2 >= ngrid) {
            k2 -= ngrid;
            sz = -bsize;
          }
#else
          if (k2 < 0) continue;
          else if (k2 >= ngrid) continue;
#endif

          const size_t idx2 = ((size_t) i2 * ngrid + j2) * ngrid + k2;
          if (!grid[idx2].n) continue;
          /* Visit all pairs associated with the two grids. */
#ifndef BENCHMARK_TIMING
          (*nnode)++;
#endif
          for (size_t n1 = 0; n1 < grid[idx1].n; n1++) {
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
            const real nx = grid[idx1].x[n1] + sx;
            const real ny = grid[idx1].y[n1] + sy;
            const real nz = grid[idx1].z[n1] + sz;
#endif
            for (size_t n2 = 0; n2 < grid[idx2].n; n2++) {
#if     BENCHMARK_PAIRCOUNT_TYPE  ==  BENCHMARK_PAIRCOUNT_BOX
              const real dx = nx - grid[idx2].x[n2];
              const real dy = ny - grid[idx2].y[n2];
              const real dz = nz - grid[idx2].z[n2];
#else
              const real dx = grid[idx1].x[n1] - grid[idx2].x[n2];
              const real dy = grid[idx1].y[n1] - grid[idx2].y[n2];
              const real dz = grid[idx1].z[n1] - grid[idx2].z[n2];
#endif
              const real dist = dx * dx + dy * dy + dz * dz;
#ifndef BENCHMARK_TIMING
              (*ndist)++;
#endif
#if     BENCHMARK_BIN_SMIN  ==  BENCHMARK_BIN_MIN_NONZERO
              if (dist >= dmin && dist < dmax) (*npair)++;
#else
              if (dist < dmax) (*npair)++;
#endif
            }
          }
        }
      }
    }
  }
}

#undef BENCHMARK_PAIRCOUNT_NAME
#undef BENCHMARK_BIN_SMIN_NAME

#undef BENCHMARK_PAIRCOUNT_TYPE
#undef BENCHMARK_BIN_SMIN

#endif


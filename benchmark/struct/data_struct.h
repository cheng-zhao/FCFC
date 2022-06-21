/*******************************************************************************
* benchmark/struct/data_struct.h: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __DATA_STRUCT_H__
#define __DATA_STRUCT_H__

#include "create_data.h"
#ifdef BENCHMARK_TIMING
#include <time.h>
#endif

/*============================================================================*\
                     Data structure for storing dual nodes
\*============================================================================*/

#ifdef BENCHMARK_TREE_AS_ARRAY
typedef struct {
  size_t a;
  size_t b;
} DUAL_NODE;
#else
typedef struct {
  const void *a;
  const void *b;
} DUAL_NODE;
#endif

typedef struct {
  DUAL_NODE *nodes;
  size_t size;
  size_t capacity;
#ifdef BENCHMARK_TREE_AS_ARRAY
  void *root;
  size_t nleaf;
#endif
} STACK_DUAL_NODE;


/*============================================================================*\
                        Functions for stack manipulation
\*============================================================================*/

/******************************************************************************
Function `stack_push`:
  Push an element to the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack;
  * `a`:        the first node to be pushed to the stack;
  * `b`:        the second node to be pushed to the stack.
******************************************************************************/
#ifdef BENCHMARK_TREE_AS_ARRAY
void stack_push(STACK_DUAL_NODE *s, const size_t a, const size_t b);
#else
void stack_push(STACK_DUAL_NODE *s, const void *a, const void *b);
#endif

/******************************************************************************
Function `stack_pop`:
  Pop an element from the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack.
Return:
  Address of the dual nodes.
******************************************************************************/
static inline DUAL_NODE *stack_pop(STACK_DUAL_NODE *s) {
  if (!s->size) return NULL;
  s->size -= 1;
  return s->nodes + s->size;
}

/******************************************************************************
Function `stack_destroy`:
  Deconstruct the stack for dual nodes.
Arguments:
  * `s`:        pointer to the stack.
******************************************************************************/
void stack_destroy(STACK_DUAL_NODE *s);


/*============================================================================*\
          Interfaces for pair counting with different data structures
\*============================================================================*/

/******************************************************************************
Function `paircount_<DATA_STRUCT>`:
  Count pairs using different data structures,
  and report the number of pairs and distance evaluations.
Arguments:
  * `data`:     structure for the data catalogue;
  * `rmin`:     minimum separation of interest;
  * `rmax`:     maximum separation of interest;
  * `csize`:    number of grids per side, or number of partilces in leaf nodes;
  * `nnode`:    number of nodes visited in total;
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
    size_t *npair);

int paircount_kdtree(DATA *data, const double rmin, const double rmax,
    const int csize,
#ifndef BENCHMARK_TIMING
    size_t *nnode, size_t *ndist,
#endif
    size_t *npair);

int paircount_balltree(DATA *data, const double rmin, const double rmax,
    const int csize,
#ifndef BENCHMARK_TIMING
    size_t *nnode, size_t *ndist,
#endif
    size_t *npair);

#endif

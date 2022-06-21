/*******************************************************************************
* benchmark/struct/define.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_STRUCT_H__
#define __DEFINE_STRUCT_H__

#include "define_comm.h"

#if defined(BENCHMARK_TREE_AS_ARRAY) && defined(BENCHMARK_TREE_PREALLOC)
#error Conflict between `BENCHMARK_TREE_AS_ARRAY` & `BENCHMARK_TREE_PREALLOC`
#endif

#if BENCHMARK_SIMD != BENCHMARK_SIMD_NONE && !defined(BENCHMARK_TREE_PREALLOC)
#error `BENCHMARK_TREE_PREALLOC` must be set for vectorization
#endif

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
#define BENCHMARK_MIN_NPAR      8
#define BENCHMARK_MAX_NPAR      INT_MAX
#define BENCHMARK_MIN_CSIZE     1
#define BENCHMARK_MAX_CSIZE     1024

#define BENCHMARK_STRUCT_GRID           0
#define BENCHMARK_STRUCT_KDTREE         1
#define BENCHMARK_STRUCT_BALLTREE       2

/*============================================================================*\
                        Definitions of runtime constants
\*============================================================================*/
#define CODE_NAME       "BENCHMARK_struct"

#define BENCHMARK_STACK_INIT_SIZE       128      /* initial size of stacks */
#define BENCHMARK_STACK_MAX_SIZE        SIZE_MAX /* maximum size of stacks */
#define BENCHMARK_TREE_MAX_SIZE         INT_MAX  /* maximum number of nodes */

#define BENCHMARK_PAIRCOUNT_BOX         0
#define BENCHMARK_PAIRCOUNT_NOBOX       1
#define BENCHMARK_BIN_MIN_ZERO          0
#define BENCHMARK_BIN_MIN_NONZERO       1

#define BENCHMARK_PCA_EPSILON           1e-8   /* precision for PCA */

#endif

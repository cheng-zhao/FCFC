/*******************************************************************************
* benchmark/histogram/define.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_HIST_H__
#define __DEFINE_HIST_H__

#include "define_comm.h"

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
#define BENCHMARK_HIST_MAX_NSP          SIZE_MAX
#define BENCHMARK_HIST_MAX_NBIN         65536

#define BENCHMARK_HIST_BIN_LIN          0
#define BENCHMARK_HIST_BIN_LOG          1

/* Methods for updating the distance histogram. */
#define BENCHMARK_HIST_BSEARCH          0
#define BENCHMARK_HIST_REV_TRAV         1
#define BENCHMARK_HIST_INT_TABLE        2
#define BENCHMARK_HIST_SQRT_TRUNC       3
#define BENCHMARK_HIST_LOG_TRUNC        4

#define BENCHMARK_HIST_FASTSQRT_START   100
#define BENCHMARK_HIST_FASTLOG2_START   200
#define BENCHMARK_HIST_APPROX_POLY_MIN  1
#define BENCHMARK_HIST_APPROX_POLY_MAX  5

/* Widths of lookup tables */
#define BENCHMARK_HIST_TABLE_W8         1
#define BENCHMARK_HIST_TABLE_W16        2
#define BENCHMARK_HIST_TABLE_W32        4

/* Limits of the lookup table size. */
#define BENCHMARK_HIST_TABLE_MIN_SIZE   1
#define BENCHMARK_HIST_TABLE_MAX_SIZE   1073741824

/*============================================================================*\
                        Definitions of runtime constants
\*============================================================================*/
#define CODE_NAME       "BENCHMARK_hist"

/* Number of random squared distances generated at once. */
#define BENCHMARK_HIST_DIST_CHUNK       67108864

/* Maximum number of fraction digits for integerizing distance bins. */
#define BENCHMARK_HIST_MAX_DIGIT_TO_INT 5
/* Maximum allowed integer squared distance edge:
   2^53, largest consecutive integer */
#define BENCHMARK_HIST_MAX_INT_DIST     9007199254740992ULL

#endif

/*******************************************************************************
* benchmark/common/define_comm.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_COMM_H__
#define __DEFINE_COMM_H__

#include "define_simd.h"

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
#define BENCHMARK_RAND_SEED     1

/*============================================================================*\
                        Definitions of runtime constants
\*============================================================================*/
#ifndef EXIT_FAILURE
#define EXIT_FAILURE    1
#endif

#ifdef SINGLE_PREC
  typedef float real;
  #define REAL_TOL      1e-6    /* tolerance for float number comparison */
  #ifdef FLT_EPSILON
    #define REAL_EPS    REAL_EPS
  #else
    #define REAL_EPS    1e-7
  #endif
  #define REAL_ABS(x)   fabsf(x)
  #define REAL_SQRT(x)  sqrtf(x)
  #define REAL_LOG(x)   logf(x)
#else
  typedef double real;
  #define REAL_TOL      1e-10   /* tolerance for double number comparison */
  #ifdef DBL_EPSILON
    #define REAL_EPS    DBL_EPSILON
  #else
    #define REAL_EPS    1e-16
  #endif
  #define REAL_ABS(x)   fabs(x)
  #define REAL_SQRT(x)  sqrt(x)
  #define REAL_LOG(x)   log(x)
#endif

#define DBL_TOL         1e-10

/*============================================================================*\
                         Definitions for output formats
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL        "%g"             /* Output format for double numbers */

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)

#endif

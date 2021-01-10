/*******************************************************************************
* define_comm.h: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>
 
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

#include <float.h>

/*============================================================================*\
                           Definitions of data types
\*============================================================================*/
#ifdef SINGLE_PREC
  typedef float real;
  #define REAL_TOL      1e-6    /* tolerance for float number comparison */
  #ifdef FLT_EPSILON
    #define REAL_EPS    FLT_EPSILON
  #else
    #define REAL_EPS    1e-7
  #endif
  #define REAL_FMT      "%f"
  #define REAL_OFMT     "%.6g"
  #define REAL_NAN      HUGE_VALF
  #define FCFC_DTYPE_DATA       AST_DTYPE_FLOAT
#else
  typedef double real;
  #define REAL_TOL      1e-10   /* tolerance for double number comparison */
  #ifdef DBL_EPSILON
    #define REAL_EPS    DBL_EPSILON
  #else
    #define REAL_EPS    1e-16
  #endif
  #define REAL_FMT      "%lf"
  #define REAL_OFMT     "%.10g"
  #define REAL_NAN      HUGE_VAL
  #define FCFC_DTYPE_DATA       AST_DTYPE_DOUBLE
#endif

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/
#define SPEED_OF_LIGHT  299792.458
#ifndef M_PI
  #define M_PI            0x1.921fb54442d18p+1    /* PI */
#endif
#ifndef DBL_EPSILON
  #define DBL_EPSILON   1e-16   /* ~ machine epsilon for double numbers */
#endif

#define DEGREE_2_RAD    0x1.1df46a2529d39p-6    /* M_PI / 180 */
#define DBL_TOL         1e-10   /* tolerance for double number comparison */

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
/* Priority of parameters from different sources. */
#define FCFC_PRIOR_CMD               5
#define FCFC_PRIOR_FILE              1

/*============================================================================*\
                            Definitions for file IO
\*============================================================================*/
#define FCFC_PATH_SEP        '/'     /* separator for file paths     */
#define FCFC_FILE_CHUNK      1048576 /* chunk size for ASCII file IO */
#define FCFC_MAX_CHUNK       INT_MAX /* maximum allowed chunk size   */
/* Initial number of objects allocated for the catalogs.        */
#define FCFC_DATA_INIT_NUM   256
/* Comment symbol for the input files (apart from the catalog). */
#define FCFC_READ_COMMENT    '#'
/* Comment symbol for the output files.                         */
#define FCFC_SAVE_COMMENT    '#'

/*============================================================================*\
                            Other runtime constants
\*============================================================================*/
#define FCFC_MAX_ELL            6       /* maximum ell for multipoles        */
#define FCFC_MAX_BIN_NUM        1024    /* maximum number of bins            */
#define FCFC_SYM_ANA_RR         '@'     /* character for analytical RR count */

#define FCFC_STACK_INIT_SIZE    128       /* initial stack size              */
#define FCFC_STACK_MAX_SIZE     SIZE_MAX  /* maximum allowed stack size      */

/*============================================================================*\
                     Definitions for the format of outputs
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m"               /* Red "Exit"        */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"  /* Red "FAIL"        */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL "%.10lg"             /* Output format for double parameters */

/* Green logo of FCFC. */
#define FCFC_LOGO "\x1B[32;1m\
\x1B[30C       FCFC     FC \n\
\x1B[30C     FC     FC FC  \n\
\x1B[30C    FC       FC    \n\
\x1B[30C  FC FC     FC     \n\
\x1B[30C FC     FCFC       \x1B[0m\n\n\
\x1B[20C \x1B[92;1mF\x1B[34;1mast \x1B[92;1mC\x1B[34;1morrelation \
\x1B[92;1mF\x1B[34;1munction \x1B[92;1mC\x1B[34;1malculator\x1B[0m\n"

/*============================================================================*\
                          Definitions for error codes
\*============================================================================*/
#define FCFC_ERR_MEMORY         (-1)
#define FCFC_ERR_ARG            (-2)
#define FCFC_ERR_FILE           (-3)
#define FCFC_ERR_CFG            (-4)
#define FCFC_ERR_AST            (-5)
#define FCFC_ERR_ASCII          (-6)
#define FCFC_ERR_CONF           (-10)
#define FCFC_ERR_DATA           (-11)
#define FCFC_ERR_CNVT           (-12)
#define FCFC_ERR_TREE           (-13)
#define FCFC_ERR_CF             (-14)
#define FCFC_ERR_SAVE           (-15)
#define FCFC_ERR_UNKNOWN        (-99)

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_EXIT " " __VA_ARGS__)

#endif


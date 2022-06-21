/*******************************************************************************
* 2pt/define.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_H__
#define __DEFINE_H__

#include "define_comm.h"

/*============================================================================*\
               Definitions for different pair counting functions
\*============================================================================*/
#define FCFC_STRUCT_KDTREE      0
#define FCFC_STRUCT_BALLTREE    1

#define FCFC_COUNT_AUTO         0
#define FCFC_COUNT_CROSS        1

#define FCFC_BIN_ISO            0
#define FCFC_BIN_SMU            1
#define FCFC_BIN_SPI            2

#define FCFC_LOOKUP_TYPE_INT    0
#define FCFC_LOOKUP_TYPE_HYBRID 1

#define FCFC_BIN_MIN_ZERO       0
#define FCFC_BIN_MIN_NONZERO    1

#define FCFC_COUNT_NO_WT        0
#define FCFC_COUNT_WITH_WT      1

/*============================================================================*\
                    Settings for pair counting optimizations
\*============================================================================*/

/* Maximum number of objects on leaf nodes. */
#if     FCFC_SIMD  ==  FCFC_SIMD_NONE
  #define FCFC_KDTREE_LEAF_SIZE         16
  #define FCFC_BALLTREE_LEAF_SIZE       8
#elif   FCFC_SIMD  ==  FCFC_SIMD_AVX  ||  FCFC_SIMD  ==  FCFC_SIMD_AVX2
  #define FCFC_KDTREE_LEAF_SIZE         32
  #define FCFC_BALLTREE_LEAF_SIZE       16
#else
  #define FCFC_KDTREE_LEAF_SIZE         32
  #define FCFC_BALLTREE_LEAF_SIZE       32
#endif

/* Entry widths of lookup tables. */
#define FCFC_LOOKUP_TABLE_W8            0       /* table with 8-bit entires */
#define FCFC_LOOKUP_TABLE_W16           1       /* table with 16-bit entries */

/* Number of entries for lookup tables. */
#define FCFC_LOOKUP_INT_MAX_SIZE        40960
#define FCFC_LOOKUP_HYBRID_MAX_SIZE     32768
#define FCFC_LOOKUP_HYBRID_MIN_SIZE     8192

/* Number of digits examined for common rescaling factors. */
#define FCFC_BIN_MAX_DIGIT_TO_INT       4

/* Least number of remainders to be processed with SIMD.
 * Process remainders always with scalar codes if it is set to 0. */
#define FCFC_SIMD_MIN_REM_SIZE          2

/*============================================================================*\
                        Definition of runtime constants
\*============================================================================*/
#define FCFC_CODE_NAME                  "FCFC_2PT"

/* Dimension of coordinates (with the sum of coordinate squares). */
#define FCFC_XDIM                       4

/* Number of redshifts for convergency tests of integrations. */
#define FCFC_INT_NUM_ZSP                128

/* Minimum stack size to be assigned to each MPI task / OpenMP thread. */
#define FCFC_STACK_SIZE_PER_TASK        4
#define FCFC_STACK_SIZE_PER_THREAD      4

#ifdef FCFC_METRIC_PERIODIC
  #undef FCFC_METRIC_PERIODIC
#endif

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
#define DEFAULT_CONF_FILE               "fcfc_2pt.conf"
#define DEFAULT_FILE_TYPE               FCFC_FFMT_ASCII
#define DEFAULT_ASCII_SKIP              0
#define DEFAULT_ASCII_COMMENT           '\0'
#define DEFAULT_COORD_CNVT              false
#define DEFAULT_DE_EOS_W                -1
#define DEFAULT_CNVT_ERR                1e-8
#define DEFAULT_STRUCT                  FCFC_STRUCT_KDTREE
#define DEFAULT_BINNING                 FCFC_BIN_ISO
#define DEFAULT_PROJECTED_CF            false
#define DEFAULT_OUTPUT_FORMAT           FCFC_OFMT_BIN
#define DEFAULT_OVERWRITE               FCFC_OVERWRITE_NONE
#define DEFAULT_VERBOSE                 true

#endif

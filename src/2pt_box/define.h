/*******************************************************************************
* 2pt_box/define.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_H__
#define __DEFINE_H__

#include "define_comm.h"

/*============================================================================*\
                         Definition of data structures
\*============================================================================*/
typedef struct {
  real x[3];
} DATA;

typedef enum {
  FCFC_OUT_LIST         = 0,
  FCFC_OUT_MATRIX       = 1
} fcfc_ostyle_t;

/*============================================================================*\
               Definitions for different pair counting functions
\*============================================================================*/
#define FCFC_TREE_TYPE_KDTREE   0

#define FCFC_PAIR_COUNT_AUTO    0
#define FCFC_PAIR_COUNT_CROSS   1

#define FCFC_BIN_ISO            0
#define FCFC_BIN_SMU            1
#define FCFC_BIN_SPI            2

#define FCFC_BIN_EXACT          0
#define FCFC_BIN_INTEG          1
#define FCFC_BIN_TRUNC          2

#define FCFC_BIN_MIN_ZERO       0
#define FCFC_BIN_MIN_NONZERO    1

/*============================================================================*\
                            Configurations of trees
\*============================================================================*/
#define KDTREE_LEAF_SIZE        4

/*============================================================================*\
                        Definition of runtime constants
\*============================================================================*/
#define FCFC_CODE_NAME          "FCFC_2PT_BOX"
#define FCFC_MAX_DIST_PREC      4

/* Minimum stack size to be assigned to each thread. */
#define FCFC_STACK_SIZE_PER_THREAD      8

#ifdef FCFC_DATA_WEIGHT
  #undef FCFC_DATA_WEIGHT
#endif

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
#define DEFAULT_CONF_FILE               "fcfc_2pt_box.conf"
#define DEFAULT_FILE_TYPE               FCFC_FFMT_ASCII
#define DEFAULT_ASCII_SKIP              0
#define DEFAULT_ASCII_COMMENT           '\0'
#define DEFAULT_BINNING                 FCFC_BIN_ISO
#define DEFAULT_PROJECTED_CF            false
#define DEFAULT_OUTPUT_STYLE            0
#define DEFAULT_OVERWRITE               0
#define DEFAULT_VERBOSE                 true

#endif

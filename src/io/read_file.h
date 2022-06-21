/*******************************************************************************
* read_file.h: this file is part of the FCFC program.

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

#ifndef __READ_FILE_H__
#define __READ_FILE_H__

#include "define.h"
#include <stdio.h>
#include <stdbool.h>

#define P_AST_ERR(ast)  ast_perror(ast, stderr, FMT_ERR);

/*============================================================================*\
                           Interfaces for ASCII files
\*============================================================================*/

/******************************************************************************
Function `read_ascii_table`:
  Read the first two columns of an ASCII file as double arrays.
Arguments:
  * `fname`:    filename of the input catalog;
  * `x`:        array for the first column;
  * `y`:        array for the second column;
  * `num`:      number of lines read successfully.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_table(const char *fname, double **x, double **y, size_t *num);

/******************************************************************************
Function `read_ascii_data`, `read_fits_data`, `read_hdf5_data`:
  Read an ASCII file for the positions and weights.
Arguments:
  * `fname`:    filename of the input catalog;
  * `skip`:     number of lines to be skipped before reading positions;
  * `comment`:  character indicating the beginning of a comment line;
  * `fmtr`:     formatter string for `sscanf`;
  * `rcol_ids`: identifiers of columns for floating-point numbers;
  * `nrcol`:    number of columns for floating-point numbers;
  * `sel`:      data selection criteria;
  * `rout`:     address for storing the output floating-point columns;
  * `num`:      number of objects read in total;
  * `verb`:     indicate whether to show detailed standard outputs
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_data(const char *fname, const size_t skip, const char comment,
    const char *fmtr, char *const *rcol_ids, const int nrcol, const char *sel,
    real ***rout, size_t *num, const int verb);

#ifdef WITH_CFITSIO
int read_fits_data(const char *fname, char *const *rcol_ids, const int nrcol,
    const char *sel, real ***rout, size_t *num, const int verb);
#endif

#ifdef WITH_HDF5
int read_hdf5_data(const char *fname, char *const *rcol_ids, const int nrcol,
    const char *sel, real ***rout, size_t *num, const int verb);
#endif

#endif

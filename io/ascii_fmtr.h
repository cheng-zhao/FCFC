/*******************************************************************************
* ascii_fmtr.h: this file is part of the FCFC program.

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

#ifndef __ASCII_FMTR_H__
#define __ASCII_FMTR_H__

#include "define_comm.h"
#include <stdio.h>

/*============================================================================*\
                           Definitions of data types
\*============================================================================*/
typedef enum {
  ASCII_DTYPE_SKIP,             /* not parsed           */
  ASCII_DTYPE_CHAR,             /* char                 */
  ASCII_DTYPE_SCHAR,            /* signed char          */
  ASCII_DTYPE_UCHAR,            /* unsigned char        */
  ASCII_DTYPE_WCHAR,            /* wchar_t              */
  ASCII_DTYPE_INT,              /* int                  */
  ASCII_DTYPE_UINT,             /* unsigned int         */
  ASCII_DTYPE_SHRT,             /* short                */
  ASCII_DTYPE_USHRT,            /* unsigned short       */
  ASCII_DTYPE_LONG,             /* long                 */
  ASCII_DTYPE_ULONG,            /* unsigned long        */
  ASCII_DTYPE_LLNG,             /* long long            */
  ASCII_DTYPE_ULLNG,            /* unsigned long long   */
  ASCII_DTYPE_INTMX,            /* intmax_t             */
  ASCII_DTYPE_UINTMX,           /* uintmax_t            */
  ASCII_DTYPE_SIZE,             /* size_t               */
  ASCII_DTYPE_PTR,              /* ptrdiff_t            */
  ASCII_DTYPE_FLT,              /* float                */
  ASCII_DTYPE_DBL,              /* double               */
  ASCII_DTYPE_LDBL,             /* long double          */
  ASCII_DTYPE_VOID,             /* void *               */
  ASCII_DTYPE_STR,              /* char *               */
  ASCII_DTYPE_WSTR              /* wchar_t *            */
} asc_dtype_t;

/*============================================================================*\
                   Data structure for the parsed information
\*============================================================================*/
typedef struct {
  asc_dtype_t dtype;    /* data type for the specifier */
  char *fmtr;           /* position in the formatter string */
} asc_arg_t;

/*============================================================================*\
                       Interface for the formatter parser
\*============================================================================*/

/******************************************************************************
Function `parse_ascii_fmtr`:
  Parse the formatter string for reading ASCII files.
  The formatter is compliant with fscanf, see C99 standard: section 7.19.6.2.
Arguments:
  * `fmtr`:     the formatter string;
  * `num`:      number of parsed arguments;
  * `rnum`:     number of arguments that are not suppressed.
Return:
  Arguments parsed from the formatter on success; NULL on error.
******************************************************************************/
asc_arg_t *parse_ascii_fmtr(const char *fmtr, int *num, int *rnum);

/******************************************************************************
Function `ascii_arg_destroy`:
  Deconstruct the structure array for parsed arguments.
Arguments:
  * `arg`:      parsed arguments from a formatter string;
  * `num`:      number of parsed elements in the array.
******************************************************************************/
void ascii_arg_destroy(asc_arg_t *arg, const int num);

#endif

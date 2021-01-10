/*******************************************************************************
* 2pt_box/save_res.h: this file is part of the FCFC program.

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

#ifndef __SAVE_RES_H__
#define __SAVE_RES_H__

#include "load_conf.h"
#include "eval_cf.h"

/* Identifiers of the results. */
typedef enum {
  FCFC_OUTPUT_PAIR_COUNT        = 0,    /* pair counts                   */
  FCFC_OUTPUT_2PCF_RAW          = 1,    /* raw 2PCF from the pair counts */
  FCFC_OUTPUT_2PCF_INTEG        = 2     /* integrated 2PCF               */
} fcfc_out_t;

/*============================================================================*\
                          Interface for saving results
\*============================================================================*/

/******************************************************************************
Function `save_res`:
  Write the pair counts or correlation functions to a text file.
Arguments:
  * `conf`:     structure for storing all configurations;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the pair counts or correlation function to be saved;
  * `flag`:     identifier of the results to be saved.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_res(const CONF *conf, const CF *cf, const int idx,
    const fcfc_out_t flag);

#endif

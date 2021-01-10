/*******************************************************************************
* 2pt/build_tree.h: this file is part of the FCFC program.

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

#ifndef __BUILD_TREE_H__
#define __BUILD_TREE_H__

#include "load_conf.h"
#include "eval_cf.h"

/******************************************************************************
Function `tree_create`:
  Construct the tree from an input catalogue for pair counting.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the catalogue to be processed;
  * `type`:     type of the tree.
Return:
  Address of the tree on success; NULL on error.
******************************************************************************/
void *tree_create(const CONF *conf, CF *cf, const int idx, const int type);

/******************************************************************************
Function `tree_destroy`:
  Deconstruct a tree used for pair counting.
Arguments:
  * `tree`:     address of the tree;
  * `type`:     type of the tree.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
void tree_destroy(void *tree, const int type);

#endif

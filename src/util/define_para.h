/*******************************************************************************
* define_para.h: this file is part of the FCFC program.

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

#ifndef __DEFINE_PARA_H__
#define __DEFINE_PARA_H__

#if defined(OMP) || defined(MPI)

#ifndef WITH_PARA
  #define WITH_PARA
#endif

#ifdef MPI
#include <mpi.h>
#endif

/*============================================================================*\
               Definition of the data structure for parallelisms
\*============================================================================*/

typedef struct {
#ifdef MPI
  int ntask;
  int root;
  int rank;
  MPI_Comm comm;
#endif
#ifdef OMP
  int nthread;
#endif
} PARA;

/*============================================================================*\
                      Function for setting up parallelisms
\*============================================================================*/

/******************************************************************************
Function `para_init`:
  Initialise the structure for parallelisms.
******************************************************************************/
void para_init(PARA *para);

#else

#ifdef WITH_PARA
  #undef WITH_PARA
#endif

#endif

#endif

/*******************************************************************************
* libast.h: this file is part of the libast library.

* libast: C library for evaluating expressions with the abstract syntax tree.

* Github repository:
        https://github.com/cheng-zhao/libast

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

#ifndef _LIBAST_H_
#define _LIBAST_H_

#include <stdio.h>
#include <stdbool.h>

/*============================================================================*\
                            Indicators of variables
\*============================================================================*/

#define AST_VAR_FLAG            '$'
#define AST_VAR_START           '{'
#define AST_VAR_END             '}'


/*============================================================================*\
                         Definitions of data structures
\*============================================================================*/

/* Enumeration of supported data types. */
typedef enum {
  AST_DTYPE_BOOL   = 1,
  AST_DTYPE_INT    = 2,
  AST_DTYPE_LONG   = 4,
  AST_DTYPE_FLOAT  = 8,
  AST_DTYPE_DOUBLE = 16,
  AST_DTYPE_STRING = 32
} ast_dtype_t;

/* The interface of the abstract syntax tree. */
typedef struct {
  ast_dtype_t dtype;    /* Data type for the expression.        */
  long nvar;            /* Number of unique variables.          */
  void *var;            /* The list of unique variables.        */
  long *vidx;           /* Unique indices of variables.         */
  char *exp;            /* A copy of the expression string.     */
  void *ast;            /* The root node of the AST.            */
  void *error;          /* Data structure for error handling.   */
} ast_t;


/*============================================================================*\
                            Definitions of functions
\*============================================================================*/

/******************************************************************************
Function `ast_init`:
  Initialise the interface of the abstract syntax tree.
Return:
  The pointer to the interface on success; NULL on error.
******************************************************************************/
ast_t *ast_init(void);

/******************************************************************************
Function `ast_build`:
  Build the abstract syntax tree given the expression and data type.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `str`:      null terminated string for the expression;
  * `dtype`:    data type for the abstract syntax tree;
  * `eval`:     true for pre-evaluating values.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ast_build(ast_t *ast, const char *str, const ast_dtype_t dtype,
    const bool eval);

/******************************************************************************
Function `ast_set_var`:
  Set the value of a variable in the variable array
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `idx`:      index of the variable (starting from 1);
  * `value`:    pointer to a variable holding the value to be set;
  * `size`:     length of the string type variable;
  * `dtype`:    data type of the value.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ast_set_var(ast_t *ast, const long idx, const void *value,
    const size_t size, const ast_dtype_t dtype);

/******************************************************************************
Function `ast_eval`:
  Evaluate the expression given the abstract syntax tree and the variable array.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `value`:    address of the variable holding the evaluated value.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ast_eval(ast_t *ast, void *value);

/******************************************************************************
Function `ast_eval_num`:
  Evaluate the numerical expression given the variable array with the same
  data type.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `value`:    address of the variable holding the evaluated value;
  * `var`:      pointer to the variable array;
  * `size`:     number of elements in the variable array.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ast_eval_num(ast_t *ast, void *value, const void *var, const long size);

/******************************************************************************
Function `ast_perror`:
  Print the error message if there is an error.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `fp`:       output file stream;
  * `msg`:      string to be printed before the error message.
******************************************************************************/
void ast_perror(const ast_t *ast, FILE *fp, const char *msg);

/******************************************************************************
Function `ast_destroy`:
  Release memory allocated for the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree.
******************************************************************************/
void ast_destroy(ast_t *ast);

#endif

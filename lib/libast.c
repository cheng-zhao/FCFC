/*******************************************************************************
* libast.c: this file is part of the libast library.

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

#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "libast.h"

/*============================================================================*\
                                     Macros
\*============================================================================*/

/* Error codes */
#define AST_ERR_MEMORY          (-1)
#define AST_ERR_INIT            (-2)
#define AST_ERR_STRING          (-3)
#define AST_ERR_DTYPE           (-4)
#define AST_ERR_TOKEN           (-5)
#define AST_ERR_VAR             (-6)
#define AST_ERR_EXIST           (-7)
#define AST_ERR_NOEXP           (-8)
#define AST_ERR_VALUE           (-9)
#define AST_ERR_SIZE            (-10)
#define AST_ERR_EVAL            (-11)
#define AST_ERR_NVAR            (-12)
#define AST_ERR_MISMATCH        (-13)
#define AST_ERR_UNKNOWN         (-99)

#define AST_ERRNO(ast)          (((ast_error_t *)ast->error)->errno)
#define AST_IS_ERROR(ast)       (AST_ERRNO(ast) != 0)

/* Mixture data types. */
#define AST_DTYPE_NULL          0
#define AST_DTYPE_INTEGER       (AST_DTYPE_INT | AST_DTYPE_LONG)
#define AST_DTYPE_REAL          (AST_DTYPE_FLOAT | AST_DTYPE_DOUBLE)
#define AST_DTYPE_NUMBER        (AST_DTYPE_INTEGER | AST_DTYPE_REAL)
#define AST_DTYPE_NATIVE        (AST_DTYPE_BOOL | AST_DTYPE_NUMBER)
#define AST_DTYPE_ALL           (AST_DTYPE_NATIVE | AST_DTYPE_STRING)
/* Data type for numerical literals in boolean expressions. */
#define AST_DTYPE_NUM4BOOL      (AST_DTYPE_LONG | AST_DTYPE_DOUBLE)

/*============================================================================*\
                            Internal data structures
\*============================================================================*/

/* Definitions of tokens. */
typedef enum {
  AST_TOK_UNDEF       = 0,      /*  undefined token         */
  AST_TOK_NUM         = 1,      /*  a number                */
  AST_TOK_STRING      = 2,      /*  a string literal        */
  AST_TOK_VAR         = 3,      /*  a variable              */
  AST_TOK_PAREN_LEFT  = 4,      /*  `(` : left parenthesis  */
  AST_TOK_PAREN_RIGHT = 5,      /*  `)` : right parenthesis */
  AST_TOK_ABS         = 6,      /*  abs : absolute value    */
  AST_TOK_SQRT        = 7,      /* sqrt : square root       */
  AST_TOK_LN          = 8,      /*   ln : log_e             */
  AST_TOK_LOG         = 9,      /*  log : log_10            */
  AST_TOK_ISFINITE    = 10,     /*        isfinite          */
  AST_TOK_NEG         = 11,     /*  `-` : negative          */
  AST_TOK_LNOT        = 12,     /*  `!` : logical NOT       */
  AST_TOK_BNOT        = 13,     /*  `~` : bitwise NOT       */
  AST_TOK_EXP         = 14,     /* `**` : Exponent          */
  AST_TOK_MUL         = 15,     /*  `*` : multiplication    */
  AST_TOK_DIV         = 16,     /*  `/` : division          */
  AST_TOK_REM         = 17,     /*  `%` : remainder         */
  AST_TOK_ADD         = 18,     /*  `+` : addition          */
  AST_TOK_MINUS       = 19,     /*  `-` : subtraction       */
  AST_TOK_LEFT        = 20,     /* `<<` : left shift        */
  AST_TOK_RIGHT       = 21,     /* `>>` : right shift       */
  AST_TOK_LT          = 22,     /*  `<` : less than         */
  AST_TOK_LE          = 23,     /* `<=` : less or equal     */
  AST_TOK_GT          = 24,     /*  `>` : greater than      */
  AST_TOK_GE          = 25,     /* `>=` : greater or equal  */
  AST_TOK_EQ          = 26,     /* `==` : equal to          */
  AST_TOK_NEQ         = 27,     /* `!=` : not equal to      */
  AST_TOK_BAND        = 28,     /*  `&` : bitwise AND       */
  AST_TOK_BXOR        = 29,     /*  `^` : bitwise XOR       */
  AST_TOK_BOR         = 30,     /*  `|` : bitwise OR        */
  AST_TOK_LAND        = 31,     /* `&&` : logical AND       */
  AST_TOK_LOR         = 32      /* `||` : logical OR        */
} ast_tok_t;

/* Types of the tokens. */
typedef enum {
  AST_TOKT_NULL,                /* undefined token          */
  AST_TOKT_UOPT,                /* a unary operator         */
  AST_TOKT_BOPT,                /* a binary operator        */
  AST_TOKT_PAREN,               /* parenthesis              */
  AST_TOKT_FUNC,                /* pre-defined function     */
  AST_TOKT_VALUE,               /* number or string literal */
  AST_TOKT_VAR                  /* variable                 */
} ast_tok_type_t;

/* Operator attributes. */
typedef struct {
  ast_tok_type_t type;          /* type of the token        */
  int precedence;               /* precedence of the token  */
  int argc;                     /* number of arguments      */
  int idtype;                   /* accepted argument dtype  */
  int odtype;                   /* returned data type       */
} ast_tok_attr_t;

static const ast_tok_attr_t ast_tok_attr[] = {
  /** token **/
  /*    type    precedence    argc    idtype    odtype    */
  /** AST_TOK_UNDEF                     **/
  {AST_TOKT_NULL,   -1,  1,     AST_DTYPE_NULL,     AST_DTYPE_NULL},
  /** AST_TOK_NUM                       **/
  {AST_TOKT_VALUE,  99,  0,     AST_DTYPE_NULL,   AST_DTYPE_NUMBER},
  /** AST_TOK_STRING                    **/
  {AST_TOKT_VALUE,  99,  0,     AST_DTYPE_NULL,   AST_DTYPE_STRING},
  /** AST_TOK_VAR                       **/
  {AST_TOKT_VAR,    99,  0,     AST_DTYPE_NULL,      AST_DTYPE_ALL},
  /** AST_TOK_PAREN_LEFT  :     `(`     **/
  {AST_TOKT_PAREN,  88,  2,     AST_DTYPE_NULL,      AST_DTYPE_ALL},
  /** AST_TOK_PAREN_RIGHT :     `)`     **/
  {AST_TOKT_PAREN,  88,  2,     AST_DTYPE_NULL,      AST_DTYPE_ALL},
  /** AST_TOK_ABS         :     `abs`   **/
  {AST_TOKT_FUNC,   88,  1,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_SQRT        :     `sqrt`  **/
  {AST_TOKT_FUNC,   88,  1,   AST_DTYPE_NUMBER,     AST_DTYPE_REAL},
  /** AST_TOK_LN          :     `ln`    **/
  {AST_TOKT_FUNC,   88,  1,   AST_DTYPE_NUMBER,     AST_DTYPE_REAL},
  /** AST_TOK_LOG         :     `log`   **/
  {AST_TOKT_FUNC,   88,  1,   AST_DTYPE_NUMBER,     AST_DTYPE_REAL},
  /** AST_TOK_ISFINITE    :    isfinite **/
  {AST_TOKT_FUNC,   88,  1,     AST_DTYPE_REAL,     AST_DTYPE_BOOL},
  /** AST_TOK_NEG         :     `-`     **/
  {AST_TOKT_UOPT,   12,  1,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_LNOT        :     `!`     **/
  {AST_TOKT_UOPT,   12,  1,   AST_DTYPE_NATIVE,     AST_DTYPE_BOOL},
  /** AST_TOK_BNOT        :     `~`     **/
  {AST_TOKT_UOPT,   12,  1,  AST_DTYPE_INTEGER,  AST_DTYPE_INTEGER},
  /** AST_TOK_EXP         :     `**`    **/
  {AST_TOKT_BOPT,   11,  2,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_MUL         :     *`      **/
  {AST_TOKT_BOPT,   10,  2,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_DIV         :     `/`     **/
  {AST_TOKT_BOPT,   10,  2,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_REM         :     `%`     **/
  {AST_TOKT_BOPT,   10,  2,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_ADD         :     `+`     **/
  {AST_TOKT_BOPT,    9,  2,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_MINUS       :     `-`     **/
  {AST_TOKT_BOPT,    9,  2,   AST_DTYPE_NUMBER,   AST_DTYPE_NUMBER},
  /** AST_TOK_LEFT        :     `<<`    **/
  {AST_TOKT_BOPT,    8,  2,  AST_DTYPE_INTEGER,  AST_DTYPE_INTEGER},
  /** AST_TOK_RIGHT       :     `>>`    **/
  {AST_TOKT_BOPT,    8,  2,  AST_DTYPE_INTEGER,  AST_DTYPE_INTEGER},
  /** AST_TOK_LT          :     `<`     **/
  {AST_TOKT_BOPT,    7,  2,   AST_DTYPE_NUMBER,     AST_DTYPE_BOOL},
  /** AST_TOK_LE          :     `<=`    **/
  {AST_TOKT_BOPT,    7,  2,   AST_DTYPE_NUMBER,     AST_DTYPE_BOOL},
  /** AST_TOK_GT          :     `>`     **/
  {AST_TOKT_BOPT,    7,  2,   AST_DTYPE_NUMBER,     AST_DTYPE_BOOL},
  /** AST_TOK_GE          :     `>=`    **/
  {AST_TOKT_BOPT,    7,  2,   AST_DTYPE_NUMBER,     AST_DTYPE_BOOL},
  /** AST_TOK_EQ          :     `==`    **/
  {AST_TOKT_BOPT,    6,  2,      AST_DTYPE_ALL,     AST_DTYPE_BOOL},
  /** AST_TOK_NEQ         :     `!=`    **/
  {AST_TOKT_BOPT,    6,  2,      AST_DTYPE_ALL,     AST_DTYPE_BOOL},
  /** AST_TOK_BAND        :     `&`     **/
  {AST_TOKT_BOPT,    5,  2,  AST_DTYPE_INTEGER,  AST_DTYPE_INTEGER},
  /** AST_TOK_BXOR        :     `^`     **/
  {AST_TOKT_BOPT,    4,  2,  AST_DTYPE_INTEGER,  AST_DTYPE_INTEGER},
  /** AST_TOK_BOR         :     `|`     **/
  {AST_TOKT_BOPT,    3,  2,  AST_DTYPE_INTEGER,  AST_DTYPE_INTEGER},
  /** AST_TOK_LAND        :     `&&`    **/
  {AST_TOKT_BOPT,    2,  2,     AST_DTYPE_BOOL,     AST_DTYPE_BOOL},
  /** AST_TOK_LOR         :     `||`    **/
  {AST_TOKT_BOPT,    1,  2,     AST_DTYPE_BOOL,     AST_DTYPE_BOOL}
};

/* Tagged union for variables with different data types. */
typedef struct {
  int dtype;
  union {
    bool bval; int ival; long lval; float fval; double dval;
    struct ast_string_struct_t { size_t len; char *str; } sval;
  } v;
} ast_var_t;

/* Data structure for error handling. */
typedef struct {
  int errno;                    /* identifier of the error        */
  long vidx;                    /* index of the variable on error */
  const char *tpos;             /* position of the token on error */
  const char *msg;              /* error message                  */
  bool *vset;                   /* check if the variable is set   */
} ast_error_t;

/* The abstract syntax tree (AST). */
typedef struct ast_tree_struct {
  ast_tok_t type;                       /* type of the token           */
  ast_var_t value;                      /* value of the token          */
  const char *ptr;                      /* poisition in the expression */
  struct ast_tree_struct *parent;       /* parent node                 */
  struct ast_tree_struct *left;         /* left child node             */
  struct ast_tree_struct *right;        /* right child node            */
} ast_node_t;


/*============================================================================*\
                        Functions for interface handling
\*============================================================================*/

/******************************************************************************
Function `ast_init`:
  Initialise the interface of the abstract syntax tree.
Return:
  The pointer to the interface on success; NULL on error.
******************************************************************************/
ast_t *ast_init(void) {
  ast_t *ast = malloc(sizeof *ast);
  if (!ast) return NULL;

  ast_error_t *err = malloc(sizeof(ast_error_t));
  if (!err) {
    free(ast);
    return NULL;
  }
  err->errno = 0;
  err->vidx = 0;
  err->tpos = err->msg = NULL;
  err->vset = NULL;
  ast->error = err;

  ast->nvar = 0;
  ast->var = NULL;
  ast->vidx = NULL;
  ast->exp = NULL;
  ast->ast = NULL;
  return ast;
}

/******************************************************************************
Function `ast_msg`:
  Record the error message.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `msg`:      the null terminated error message;
  * `vidx`:     index of the variable on error;
  * `tpos`:     pointer to the token on error.
******************************************************************************/
static void ast_msg(ast_t *ast, const char *msg, const long vidx,
    const char *tpos) {
  if (!ast || !ast->error) return;
  ast_error_t *err = (ast_error_t *) ast->error;
  err->msg = msg;
  err->vidx = vidx;
  err->tpos = tpos;
}


/*============================================================================*\
           Functions for handling variables with different data types
\*============================================================================*/

/******************************************************************************
Function `ast_set_var_value`:
  Set the value of a variable.
Arguments:
  * `var`:      pointer to the variable (assume memory has been allocated);
  * `value`:    pointer to the value;
  * `size`:     length of the string type variable;
  * `dtype`:    data type of the value (no validation performed).
******************************************************************************/
static inline void ast_set_var_value(ast_var_t *var, const void *value,
    const size_t size, const ast_dtype_t dtype) {
  var->dtype = dtype;
  switch (dtype) {
    case AST_DTYPE_BOOL:   var->v.bval = *((bool *) value);   break;
    case AST_DTYPE_INT:    var->v.ival = *((int *) value);    break;
    case AST_DTYPE_LONG:   var->v.lval = *((long *) value);   break;
    case AST_DTYPE_FLOAT:  var->v.fval = *((float *) value);  break;
    case AST_DTYPE_DOUBLE: var->v.dval = *((double *) value); break;
    case AST_DTYPE_STRING:
      var->v.sval.str = (char *) value;
      var->v.sval.len = size;
      break;
    default: return;
  }
}

/******************************************************************************
Function `ast_vidx_pos`:
  Find the insertion position of the variable index in a sorted index array.
Arguments:
  * `arr`:      pointer to the array for variable indices;
  * `num`:      number of elements in the array;
  * `var`:      index of the variable to be inserted.
Return:
  -(index + 1) if the array contains the value;
  the index for the insertion otherwise.
******************************************************************************/
static long ast_vidx_pos(const long *arr, const long num, const long vidx) {
  long l, u;
  l = 0;
  u = num - 1;
  while (l <= u) {
    long i = ((unsigned long) l + (unsigned long) u) >> 1;
    if (arr[i] < vidx) l = i + 1;
    else if (arr[i] > vidx) u = i - 1;
    else return -(i + 1);       /* the element is already there */
  }
  return l;
}

/******************************************************************************
Function `ast_init_var`:
  Initialise the variable array.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `dtype`:    data type of the expression.
Return:
  Pointer to the array.
******************************************************************************/
static void *ast_init_var(ast_t *ast, const ast_dtype_t dtype) {
  ast_error_t *err = (ast_error_t *) ast->error;
  err->vset = calloc(ast->nvar, sizeof(bool));
  if (!err->vset) return NULL;

  if (dtype == AST_DTYPE_BOOL) {
    ast_var_t *var = calloc(ast->nvar, sizeof *var);
    if (!var) return NULL;
    for (long i = 0; i < ast->nvar; i++) var[i].dtype = AST_DTYPE_ALL;
    return var;
  }
  if (dtype == AST_DTYPE_INT) {
    int *var = calloc(ast->nvar, sizeof *var);
    return var;
  }
  if (dtype == AST_DTYPE_LONG) {
    long *var = calloc(ast->nvar, sizeof *var);
    return var;
  }
  if (dtype == AST_DTYPE_FLOAT) {
    float *var = calloc(ast->nvar, sizeof *var);
    return var;
  }
  if (dtype == AST_DTYPE_DOUBLE) {
    double *var = calloc(ast->nvar, sizeof *var);
    return var;
  }
  AST_ERRNO(ast) = AST_ERR_DTYPE;
  return NULL;
}

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
    const size_t size, ast_dtype_t dtype) {
  if (!ast) return AST_ERR_INIT;
  if(AST_IS_ERROR(ast)) return AST_ERRNO(ast);
  if (!ast->nvar) return 0;
  if (!value) return AST_ERRNO(ast) = AST_ERR_VALUE;

  if (idx <= 0) {
    ast_msg(ast, "unexpected variable index", idx, NULL);
    return AST_ERRNO(ast) = AST_ERR_VAR;
  }

  /* Nothing to be done if this variable is not required. */
  const long pos = -ast_vidx_pos(ast->vidx, ast->nvar, idx) - 1;
  if (pos < 0) return 0;
  long lval;
  double dval;
  ast_var_t *var;

  switch (ast->dtype) {
    case AST_DTYPE_BOOL:
      var = (ast_var_t *) ast->var + pos;
      /* Convert int to long. */
      if (dtype == AST_DTYPE_INT) {
        lval = (long) (*((int *) value));
        value = &lval;
        dtype = AST_DTYPE_LONG;
      }
      /* Convert float to double. */
      else if (dtype == AST_DTYPE_FLOAT) {
        dval = (double) (*((float *) value));
        value = &dval;
        dtype = AST_DTYPE_DOUBLE;
      }

      /* Raise an error if the data type is not valid for this variable. */
      if (!(dtype & var->dtype)) {
        ast_msg(ast, "unexpected data type for variable", idx, NULL);
        return AST_ERRNO(ast) = AST_ERR_VAR;
      }
      ast_set_var_value(var, value, size, dtype);
      break;
    case AST_DTYPE_INT:
      if (dtype != ast->dtype) {
        ast_msg(ast, "int type variable expected", idx, NULL);
        return AST_ERRNO(ast) = AST_ERR_VAR;
      }
      *((int *) ast->var + pos) = *((int *) value);
      break;
    case AST_DTYPE_LONG:
      if (dtype == AST_DTYPE_LONG)
        *((long *) ast->var + pos) = *((long *) value);
      /* Convert int to long. */
      else if (dtype == AST_DTYPE_INT)
        *((long *) ast->var + pos) = (long) (*((int *) value));
      else {
        ast_msg(ast, "long type variable expected", idx, NULL);
        return AST_ERRNO(ast) = AST_ERR_VAR;
      }
      break;
    case AST_DTYPE_FLOAT:
      if (dtype == AST_DTYPE_FLOAT)
        *((float *) ast->var + pos) = *((float *) value);
      /* Convert int to float. */
      else if (dtype == AST_DTYPE_INT)
        *((float *) ast->var + pos) = (float) (*((int *) value));
      /* Convert long to float. */
      else if (dtype == AST_DTYPE_LONG)
        *((float *) ast->var + pos) = (float) (*((long *) value));
      else {
        ast_msg(ast, "float type variable expected", idx, NULL);
        return AST_ERRNO(ast) = AST_ERR_VAR;
      }
      break;
    case AST_DTYPE_DOUBLE:
      if (dtype == AST_DTYPE_DOUBLE)
        *((double *) ast->var + pos) = *((double *) value);
      else if (dtype == AST_DTYPE_INT)
        *((double *) ast->var + pos) = (double) (*((int *) value));
      else if (dtype == AST_DTYPE_LONG)
        *((double *) ast->var + pos) = (double) (*((long *) value));
      else if (dtype == AST_DTYPE_FLOAT)
        *((double *) ast->var + pos) = (double) (*((float *) value));
      else {
        ast_msg(ast, "double type variable expected", idx, NULL);
        return AST_ERRNO(ast) = AST_ERR_VAR;
      }
      break;
    default:
      ast_msg(ast, "unknown error for variable", idx, NULL);
      return AST_ERRNO(ast) = AST_ERR_VAR;
  }

  ((ast_error_t *) ast->error)->vset[pos] = true;
  return 0;
}

/******************************************************************************
Function `ast_save_vidx`:
  Record the index of a variable if necessary;
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `idx`:      the index to be recorded.
******************************************************************************/
static void ast_save_vidx(ast_t *ast, const long idx) {
  const long pos = (ast->vidx) ? ast_vidx_pos(ast->vidx, ast->nvar, idx) : 0;
  if (pos < 0) return;          /* the index has already been recorded */

  if (ast->nvar == LONG_MAX) {  /* there is no more space for the insertion */
    AST_ERRNO(ast) = AST_ERR_NVAR;
    return;
  }

  /* Check if the allocated space is enough. */
  if ((ast->nvar & (ast->nvar - 1)) == 0) {     /* nvar is 0 or power of 2 */
    long size = 1;
    if (LONG_MAX / 2 < ast->nvar) size = LONG_MAX;
    else if (ast->nvar) size = ast->nvar << 1;  /* double the size */

    long *tmp = realloc(ast->vidx, size * sizeof(long));
    if (!tmp) {
      AST_ERRNO(ast) = AST_ERR_MEMORY;
      return;
    }
    ast->vidx = tmp;
  }

  /* Right shift the existing elements. */
  if (ast->nvar && pos < ast->nvar)
    memmove(ast->vidx + pos + 1, ast->vidx + pos,
        (ast->nvar - pos) * sizeof(long));
  ast->vidx[pos] = idx;
  ast->nvar += 1;
}


/*============================================================================*\
                       Functions for string manipulation
\*============================================================================*/

/******************************************************************************
Function `ast_skip_space`:
  Skip whitespaces.
Arguments:
  * `src`:      the null-terminated input string.
Return:
  Pointer to the first non-whitespace character of the string.
******************************************************************************/
static inline const char *ast_skip_space(const char *src) {
  const char *dst = src;
  while (isspace(*dst)) dst++;
  return dst;
}

/******************************************************************************
Function `ast_copy_str`:
  Copy a string to a new block of memory.
Arguments:
  * `src`:      the null-terminated input string.
Return:
  Pointer to the copied string.
******************************************************************************/
static inline char *ast_copy_str(const char *src) {
  size_t size = strlen(src);
  char *dst = calloc(size + 1, sizeof(char));
  if (!dst) return NULL;
  memcpy(dst, src, size * sizeof(char));
  return dst;
}

/******************************************************************************
Function `ast_parse_num`:
  Convert a numerical token into a number.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `str`:      the input string;
  * `res`:      the resulting number;
  * `end`:      pointer to the first character that is not interpreted.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int ast_parse_num(ast_t *ast, const char *str, ast_var_t *res,
    char **end) {
  /* Validate arguments. */
  if (!ast) return AST_ERR_INIT;
  if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);

  int ival;
  long lval;
  float fval;
  double dval;

  switch (ast->dtype) {
    case AST_DTYPE_INT:
      lval = strtol(str, end, 10);
      if (lval >= INT_MAX || lval <= INT_MIN) {
        ast_msg(ast, "overflow detected for int type", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      ival = (int) lval;
      ast_set_var_value(res, &ival, 0, AST_DTYPE_INT);
      break;
    case AST_DTYPE_LONG:
      lval = strtol(str, end, 10);
      if (lval == LONG_MAX || lval == LONG_MIN) {
        ast_msg(ast, "overflow detected for long type", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      ast_set_var_value(res, &lval, 0, AST_DTYPE_LONG);
      break;
    case AST_DTYPE_FLOAT:
      fval = strtof(str, end);
      if (fval == HUGE_VALF) {
        ast_msg(ast, "overflow detected for float type", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      ast_set_var_value(res, &fval, 0, AST_DTYPE_FLOAT);
      break;
    case AST_DTYPE_DOUBLE:
      dval = strtod(str, end);
      if (dval == HUGE_VAL) {
        ast_msg(ast, "overflow detected for double type", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      ast_set_var_value(res, &dval, 0, AST_DTYPE_DOUBLE);
      break;
    case AST_DTYPE_BOOL:
      /* Try to convert to double. */
      dval = strtod(str, end);
      if (dval == HUGE_VAL) {
        ast_msg(ast, "overflow detected for double type", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      /* Check if it can also be long. */
      if (dval > LONG_MIN && dval < LONG_MAX && str != *end) {
        bool is_long = true;
        for (const char *tmp = str; tmp < *end; tmp++) {
          if (!isdigit(*tmp)) {
            is_long = false;
            break;
          }
        }
        if (is_long) {
          lval = strtol(str, end, 10);
          ast_set_var_value(res, &lval, 0, AST_DTYPE_LONG);
          break;
        }
      }
      ast_set_var_value(res, &dval, 0, AST_DTYPE_DOUBLE);
      break;
    default:
      return AST_ERRNO(ast) = AST_ERR_DTYPE;
  }

  if (*end - str == 0) {                        /* no character is parsed */
    ast_msg(ast, "unrecognised number", 0, str);
    return AST_ERRNO(ast) = AST_ERR_TOKEN;
  }
  return 0;
}


/******************************************************************************
Function `ast_parse_str`:
  Retrieve a string literal from the expression.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `str`:      the input string;
  * `res`:      the resulting tagged union for different data types;
  * `end`:      pointer to the first character that is not interpreted.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int ast_parse_str(ast_t *ast, const char *str, ast_var_t *res,
    const char **end) {
  const char *c = str + 1;
  while (*c != '\0' && *c != *str) c++;
  if (*c == '\0') {
    ast_msg(ast, "unbalanced quotation mark", 0, str);
    return AST_ERRNO(ast) = AST_ERR_TOKEN;
  }
  ast_set_var_value(res, str + 1, c - str - 1, AST_DTYPE_STRING);
  *end = c + 1;
  return 0;
}

/******************************************************************************
Function `ast_parse_var`:
  Convert the index of a variable into a long number.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `str`:      the input string;
  * `vidx`:     retrieved index of the variable;
  * `end`:      pointer to the first character that is not interpreted.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int ast_parse_var(ast_t *ast, const char *str, long *vidx,
    const char **end) {
  const char *c = str + 1;
  if (*c >= '1' && *c <= '9') {         /* number 1 - 9 */
    *vidx = *c - '0';
    *end = c + 1;
    return 0;
  }
  else if (*c == AST_VAR_START) {
    /* Get the variable index (long integer). */
    long idx = 0;
    for (++c; isdigit(*c); c++) {
      int digit = *c - '0';
      /* Overflow detection. */
      if (LONG_MAX / 10 < idx) {
        ast_msg(ast, "the variable index is too large", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      idx *= 10;
      if (LONG_MAX - digit < idx) {
        ast_msg(ast, "the variable index is too large", 0, str);
        return AST_ERRNO(ast) = AST_ERR_TOKEN;
      }
      idx += digit;
    }
    if (idx > 0 && *c == AST_VAR_END) {
      *vidx = idx;
      *end = c + 1;
      return 0;
    }
  }
  ast_msg(ast, "unrecognised variable", 0, str);
  return AST_ERRNO(ast) = AST_ERR_TOKEN;
}


/*============================================================================*\
              Functions for the abstract syntax tree manipulation
\*============================================================================*/

/******************************************************************************
Function `ast_create`:
  Create a new node of the abstract syntax tree.
Arguments:
  * `type`:     type of the node;
  * `value`:    value of the node.
Return:
  The address of the node.
******************************************************************************/
static ast_node_t *ast_create(const ast_tok_t type, const ast_var_t value) {
  ast_node_t *node = malloc(sizeof *node);
  if (!node) return NULL;
  node->type = type;
  node->value = value;
  node->ptr = NULL;
  node->parent = node->left = node->right = NULL;
  return node;
}

/******************************************************************************
Function `ast_root`:
  Find the root node of the abstract syntax tree.
Arguments:
  * `node`:     an arbitrary node of the tree.
Return:
  Address of the root node.
******************************************************************************/
static inline ast_node_t *ast_root(ast_node_t *node) {
  while (node->parent) node = node->parent;
  return node;
}

/******************************************************************************
Function `ast_delete`:
  Deleta a node from the abstract syntax tree, once a pair of parenthesis
  is parsed.
Arguments:
  * `node`:     the node to be removed.
******************************************************************************/
static void ast_delete(ast_node_t *node) {
  /* This node should have only the left child. */
  ast_node_t *tmp = node->left;

  /* Copy contents from the left child. */
  node->type = tmp->type;
  node->value = tmp->value;
  node->ptr = tmp->ptr;
  node->left = tmp->left;
  node->right = tmp->right;

  if (node->left) node->left->parent = node;
  if (node->right) node->right->parent = node;

  /* Delete the left child. */
  free(tmp);
}

/******************************************************************************
Function `ast_insert`:
  Insert a token to the abstract syntax tree.
Arguments:
  * `node`:     current node of the tree;
  * `tok`:      token for the new node;
  * `value`:    value for the new node;
  * `str`:      position of the token in the expression.
Return:
  The address of the inserted node.
******************************************************************************/
static ast_node_t *ast_insert(ast_node_t *node, const ast_tok_t tok,
    const ast_var_t value, const char *str) {
  /* No need to ceate a new node if this is the first token. */
  if (node->type == AST_TOK_UNDEF) {
    node->type = tok;
    node->value = value;
    node->ptr = str;
    return node;
  }

  /* Create a new node. */
  ast_node_t *new = ast_create(tok, value);
  if (!new) return NULL;
  new->ptr = str;

  /* Insert a new node to the abstract syntax tree. */
  if (ast_tok_attr[tok].type == AST_TOKT_BOPT) {
    /* Find the right ancestor given the precedence and associativity. */
    while (node->parent && node->parent->type != AST_TOK_PAREN_LEFT &&
        ast_tok_attr[node->parent->type].type != AST_TOKT_FUNC &&
        ((ast_tok_attr[node->parent->type].type != AST_TOKT_UOPT &&
        ast_tok_attr[node->parent->type].precedence >=
        ast_tok_attr[tok].precedence) ||
        /* Right-to-left associativity for unary operators. */
        (ast_tok_attr[node->parent->type].type == AST_TOKT_UOPT &&
        ast_tok_attr[node->parent->type].precedence >
        ast_tok_attr[tok].precedence)))
      node = node->parent;
    new->left = node;
    if (node->parent) {                 /* insert between two nodes */
      new->parent = node->parent;
      if (node == node->parent->left) node->parent->left = new;
      else node->parent->right = new;
    }
    node->parent = new;                 /* new root */
    return new;
  }
  /* Insert an effective leaf to the abstract syntax tree. */
  else {        /* AST_TOK_PAREN_LEFT || AST_TOKT_UOPT/FUNC/VALUE */
    if (!node->left) {
      node->left = new;
      new->parent = node;
    }
    else {
      node->right = new;
      new->parent = node;
    }
  }
  return new;
}


/*============================================================================*\
                             Functions for clean-up
\*============================================================================*/

/******************************************************************************
Function `ast_free`:
  Free memory allocated for the abstract syntax tree.
Arguments:
  * `node`:     the root node of the abstract syntax tree.
******************************************************************************/
static void ast_free(ast_node_t *node) {
  if (!node) return;
  ast_free(node->left);
  ast_free(node->right);
  free(node);
}

/******************************************************************************
Function `ast_destroy`:
  Release memory allocated for the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree.
******************************************************************************/
void ast_destroy(ast_t *ast) {
  if (!ast) return;
  ast_error_t *err = (ast_error_t *) ast->error;
  if (err->vset) free(err->vset);
  free(ast->error);
  if (ast->exp) free(ast->exp);
  if (ast->var) free(ast->var);
  if (ast->vidx) free(ast->vidx);
  ast_free((ast_node_t *) ast->ast);
  free(ast);
}


/*============================================================================*\
                            Functions for the parser
\*============================================================================*/

/******************************************************************************
Function `ast_parse_token`:
  Parse the token and push it to the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `src`:      the string on processing.
******************************************************************************/
static void ast_parse_token(ast_t *ast, ast_node_t *node, const char *src) {
  if (!ast || AST_IS_ERROR(ast)) return;
  src = ast_skip_space(src);                    /* skip whitespaces */
  const char *c = src;

  /* Number of leaves. */
  int argc = 0;
  if (node->left) argc++;
  if (node->right) argc++;

  if (!(*c)) {
    /* Number of leaves is smaller than the arguments of this operator. */
    if (argc < ast_tok_attr[node->type].argc) {
      ast_msg(ast, "incomplete expression", 0, src);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
    }
    /* Open parenthesis. */
    do node = node->parent;
    while (node && node->type != AST_TOK_PAREN_LEFT &&
        ast_tok_attr[node->type].type != AST_TOKT_FUNC);
    if (node) {
      ast_msg(ast, "unclosed parenthesis", 0, src);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
    }
    return;
  }

  /* Check the token, so no check is performed in `ast_insert`. */
  ast_tok_t tok = AST_TOK_UNDEF;
  switch (*c) {
    case 'i':                           /* for inf or isfinite */
      if (c[1] == 's' && c[2] == 'f' && c[3] == 'i' && c[4] == 'n'
          && c[5] == 'i' && c[6] == 't' && c[7] == 'e' && c[8] == '(') {
        c += 8;
        tok = AST_TOK_ISFINITE;
      }
      break;
    case '.':
      if (ast->dtype != AST_DTYPE_BOOL && (ast->dtype & AST_DTYPE_REAL) == 0)
        break;
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9': tok = AST_TOK_NUM; break;
    case '\'':
    case '"': tok = AST_TOK_STRING; break;
    case AST_VAR_FLAG: tok = AST_TOK_VAR; break;
    case '(': tok = AST_TOK_PAREN_LEFT; break;
    case ')': tok = AST_TOK_PAREN_RIGHT; break;
    case 'a':
      if (c[1] == 'b' && c[2] == 's' && c[3] == '(') {
        c += 3;
        tok = AST_TOK_ABS;
      }
      break;
    case 's':
      if (c[1] == 'q' && c[2] == 'r' && c[3] == 't' && c[4] == '(') {
        c += 4;
        tok = AST_TOK_SQRT;
      }
      break;
    case 'l':
      if (c[1] == 'n' && c[2] == '(') {
        c += 2;
        tok = AST_TOK_LN;
      }
      else if (c[1] == 'o' && c[2] == 'g' && c[3] == '(') {
        c += 3;
        tok = AST_TOK_LOG;
      }
      break;
    case '+': tok = AST_TOK_ADD; break;
    case '-':
      if (argc >= ast_tok_attr[node->type].argc) tok = AST_TOK_MINUS;
      else tok = AST_TOK_NEG;
      break;
    case '!':
      if (c[1] == '=') {
        c++;
        tok = AST_TOK_NEQ;
      }
      else tok = AST_TOK_LNOT;
      break;
    case '~': tok = AST_TOK_BNOT; break;
    case '*':
      if (c[1] == '*') {
        c++;
        tok = AST_TOK_EXP;
      }
      else tok = AST_TOK_MUL;
      break;
    case '/': tok = AST_TOK_DIV; break;
    case '%': tok = AST_TOK_REM; break;
    case '<':
      if (c[1] == '<') {
        c++;
        tok = AST_TOK_LEFT;
      }
      else if (c[1] == '=') {
        c++;
        tok = AST_TOK_LE;
      }
      else tok = AST_TOK_LT;
      break;
    case '>':
      if (c[1] == '>') {
        c++;
        tok = AST_TOK_RIGHT;
      }
      else if (c[1] == '=') {
        c++;
        tok = AST_TOK_GE;
      }
      else tok = AST_TOK_GT;
      break;
    case '=':
      if (*(++c) == '=') tok = AST_TOK_EQ;
      break;
    case '&':
      if (c[1] == '&') {
        c++;
        tok = AST_TOK_LAND;
      }
      else tok = AST_TOK_BAND;
      break;
    case '^': tok = AST_TOK_BXOR; break;
    case '|':
      if (c[1] == '|') {
        c++;
        tok = AST_TOK_LOR;
      }
      else tok = AST_TOK_BOR;
      break;
    default: break;
  }

  if (tok == AST_TOK_UNDEF) {
    ast_msg(ast, "unrecognised token", 0, src);
    AST_ERRNO(ast) = AST_ERR_TOKEN;
    return;
  }

  /* Validate the data type for this token. */
  if (ast->dtype != AST_DTYPE_BOOL &&
      (ast->dtype & ast_tok_attr[tok].odtype) == 0) {
    ast_msg(ast, "invalid token for the desired data type", 0, src);
    AST_ERRNO(ast) = AST_ERR_TOKEN;
    return;
  }

  /* Validate the number of arguments. */
  if(argc >= ast_tok_attr[node->type].argc) {
    /* This token can only be added as a new node. */
    if (tok == AST_TOK_PAREN_LEFT ||
        ast_tok_attr[tok].type == AST_TOKT_UOPT ||
        ast_tok_attr[tok].type == AST_TOKT_FUNC ||
        ast_tok_attr[tok].type == AST_TOKT_VALUE ||
        ast_tok_attr[tok].type == AST_TOKT_VAR) {
      ast_msg(ast, "missing operator", 0, src);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
      return;
    }
  }
  else {
    /* This token can only be added as an effective leaf. */
    if (tok == AST_TOK_PAREN_RIGHT || ast_tok_attr[tok].type == AST_TOKT_BOPT) {
      ast_msg(ast, "missing value", 0, src);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
      return;
    }
  }

  /* Validate right parenthesis and remove the node if necessary. */
  if (tok == AST_TOK_PAREN_RIGHT) {
    if (node->type == AST_TOK_PAREN_LEFT)
      ast_msg(ast, "empty parenthesis", 0, src);
    else {
      /* Check all the ancestors for left parenthesis or functions. */
      do node = node->parent;
      while (node && node->type != AST_TOK_PAREN_LEFT &&
          ast_tok_attr[node->type].type != AST_TOKT_FUNC);
      if (!node) ast_msg(ast, "unbalanced parenthesis", 0, src);
      else if (node->type == AST_TOK_PAREN_LEFT) ast_delete(node);
    }
  }
  /* The error message is set. */
  if (((ast_error_t *) ast->error)->msg) {
    AST_ERRNO(ast) = AST_ERR_TOKEN;
    return;
  }

  ast_var_t v = {0, .v.ival = 0};
  /* Retrieve the value if this is a number. */
  if (tok == AST_TOK_NUM) {
    char *end;
    if (ast_parse_num(ast, c, &v, &end)) return;
    c = end;
  }
  /* Retrieve the string literal. */
  else if (tok == AST_TOK_STRING) {
    const char *end;
    if (ast_parse_str(ast, c, &v, &end)) return;
    c = end;
  }
  /* Retrieve the index if this is a variable. */
  else if (tok == AST_TOK_VAR) {
    long vidx;
    const char *end;
    if (ast_parse_var(ast, c, &vidx, &end)) return;
    ast_save_vidx(ast, vidx);
    ast_set_var_value(&v, &vidx, 0, AST_DTYPE_LONG);
    c = end;
  }
  else c++;

  /* Insert the token to the abstract syntax tree. */
  if (tok != AST_TOK_PAREN_RIGHT) {
    node = ast_insert(node, tok, v, src);
    if (!node) {
      AST_ERRNO(ast) = AST_ERR_MEMORY;
      return;
    }
  }

  /* Parse the next token. */
  ast_parse_token(ast, node, c);
}

/******************************************************************************
Function `ast_check_dtype`:
  Validata the data type a node in the boolean type abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     pointer to a node of the abstract syntax tree.
Return:
  The output data type of the node.
******************************************************************************/
static int ast_check_dtype(ast_t *ast, ast_node_t *node) {
  if (AST_IS_ERROR(ast) || !node) return AST_DTYPE_NULL;

  /* Return data type of leaf nodes. */
  if (node->type == AST_TOK_NUM) return node->value.dtype;
  else if (node->type == AST_TOK_STRING) return AST_DTYPE_STRING;
  else if (node->type == AST_TOK_VAR) {
    if (node->parent) {         /* limits of the operator input type */
      long i = node->value.v.lval;
      ((ast_var_t *) ast->var + i)->dtype &=
        ast_tok_attr[node->parent->type].idtype;
      if (((ast_var_t *) ast->var + i)->dtype == 0) {
        ast_msg(ast, "conflict operators for variable", ast->vidx[i], NULL);
        AST_ERRNO(ast) = AST_ERR_VAR;
        return AST_DTYPE_NULL;
      }
    }
    return AST_DTYPE_ALL;
  }

  const int ltype = ast_check_dtype(ast, node->left);

  /* Check unary operators. */
  if (ast_tok_attr[node->type].argc == 1) {
    if ((ast_tok_attr[node->type].idtype & ltype) == 0) {
      ast_msg(ast, "unexpected data type for operator", 0, node->ptr);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
      return AST_DTYPE_NULL;
    }
    return ast_tok_attr[node->type].odtype;
  }

  /* Check binary operators. */
  const int rtype = ast_check_dtype(ast, node->right);
  if ((ast_tok_attr[node->type].idtype & ltype & rtype) == 0) {
    /* Allow two different numerical data types. */
    if (!(ltype & AST_DTYPE_NUMBER) || !(rtype & AST_DTYPE_NUMBER)) {
      ast_msg(ast, "unexpected data type for operator", 0, node->ptr);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
      return AST_DTYPE_NULL;
    }
  }
  return ast_tok_attr[node->type].odtype;
}

/******************************************************************************
Function `ast_reset_idx`:
  Reset variable indices of the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree.
******************************************************************************/
static void ast_reset_idx(ast_t *ast, ast_node_t *node) {
  if (AST_IS_ERROR(ast) || !node) return;
  if (node->type == AST_TOK_VAR) {
    const long p = -ast_vidx_pos(ast->vidx, ast->nvar, node->value.v.lval) - 1;
    if (p < 0 || p >= ast->nvar) {
      ast_msg(ast, "unknown error for variable", node->value.v.lval, NULL);
      AST_ERRNO(ast) = AST_ERR_TOKEN;
      return;
    }
    node->value.v.lval = p;
    return;
  }
  if (node->left) ast_reset_idx(ast, node->left);
  if (node->right) ast_reset_idx(ast, node->right);
  return;
}


/*============================================================================*\
                            Functions for evaluation
\*============================================================================*/

/******************************************************************************
Function `ast_ipow`:
  Exponentiation of integers.
  Taken from: https://gist.github.com/orlp/3551590
Arguments:
  * `base`:     the base;
  * `exp`:      the exponent.
Return:
  The resulting integer on success; 0 on overflow/underflow.
******************************************************************************/
static int64_t ipow(int64_t base, uint8_t exp) {
  static const uint8_t highest_bit_set[] = {
    0, 1, 2, 2, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 255, /* 255 for overflow */
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255,
  };

  int64_t result = 1;

  switch (highest_bit_set[exp]) {
  case 255:             /* 255 is the overflow marker */
    if (base == 1) return 1;
    if (base == -1) return 1 - 2 * (exp & 1);
    return 0;
  case 6:
    if (exp & 1) result *= base;
    exp >>= 1;
    base *= base;
  case 5:
    if (exp & 1) result *= base;
    exp >>= 1;
    base *= base;
  case 4:
    if (exp & 1) result *= base;
    exp >>= 1;
    base *= base;
  case 3:
    if (exp & 1) result *= base;
    exp >>= 1;
    base *= base;
  case 2:
    if (exp & 1) result *= base;
    exp >>= 1;
    base *= base;
  case 1:
    if (exp & 1) result *= base;
  default:
    return result;
  }
}

/******************************************************************************
Function `ast_eval_int`:
  Evaluate the value in int type, given the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree;
  * `var`:      the user-supplied variable array.
Return:
  The resulting integer on success; 0 on error.
******************************************************************************/
static int ast_eval_int(ast_t *ast, const ast_node_t *node, const int *var) {
  if (AST_IS_ERROR(ast)) return 0;
  if (node->type == AST_TOK_NUM) return node->value.v.ival;
  else if (node->type == AST_TOK_VAR) {
    if (var) return var[ast->vidx[node->value.v.lval] - 1];
    else return *((int *) ast->var + node->value.v.lval);
  }
  else if (ast_tok_attr[node->type].argc == 1) {
    const int v = ast_eval_int(ast, node->left, var);
    switch (node->type) {
      case AST_TOK_NEG: return -v;
      case AST_TOK_ABS: return (v < 0) ? -v : v;
      case AST_TOK_BNOT: return ~v;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return 0;
    }
  }
  else {
    const int v1 = ast_eval_int(ast, node->left, var);
    const int v2 = ast_eval_int(ast, node->right, var);
    switch (node->type) {
      case AST_TOK_ADD: return v1 + v2;
      case AST_TOK_MINUS: return v1 - v2;
      case AST_TOK_MUL: return v1 * v2;
      case AST_TOK_DIV: return v1 / v2;
      case AST_TOK_REM: return v1 % v2;
      case AST_TOK_EXP:
        if (v2 > UINT8_MAX) {
          if (v1 == 1) return 1;
          if (v1 == -1) return 1 - 2 * (v2 & 1);
          return 0;
        }
        int64_t res = ipow(v1, v2);
        if (res > INT_MAX) return 0;
        return (int) res;
      case AST_TOK_LEFT: return v1 << v2;
      case AST_TOK_RIGHT: return v1 >> v2;
      case AST_TOK_BAND: return v1 & v2;
      case AST_TOK_BXOR: return v1 ^ v2;
      case AST_TOK_BOR: return v1 | v2;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return 0;
    }
  }
}

/******************************************************************************
Function `ast_eval_long`:
  Evaluate the value in long int type, given the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree;
  * `var`:      the user-supplied variable array.
Return:
  The resulting long integer on success; 0 on error.
******************************************************************************/
static long ast_eval_long(ast_t *ast, const ast_node_t *node, const long *var) {
  if (AST_IS_ERROR(ast)) return 0;
  if (node->type == AST_TOK_NUM) return node->value.v.lval;
  else if (node->type == AST_TOK_VAR) {
    if (var) return var[ast->vidx[node->value.v.lval] - 1];
    else return *((long *) ast->var + node->value.v.lval);
  }
  else if (ast_tok_attr[node->type].argc == 1) {
    const long v = ast_eval_long(ast, node->left, var);
    switch (node->type) {
      case AST_TOK_NEG: return -v;
      case AST_TOK_ABS: return (v < 0) ? -v : v;
      case AST_TOK_BNOT: return ~v;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return 0;
    }
  }
  else {
    const long v1 = ast_eval_long(ast, node->left, var);
    const long v2 = ast_eval_long(ast, node->right, var);
    switch (node->type) {
      case AST_TOK_ADD: return v1 + v2;
      case AST_TOK_MINUS: return v1 - v2;
      case AST_TOK_MUL: return v1 * v2;
      case AST_TOK_DIV: return v1 / v2;
      case AST_TOK_REM: return v1 % v2;
      case AST_TOK_EXP:
        if (v2 > UINT8_MAX) {
          if (v1 == 1) return 1;
          if (v1 == -1) return 1 - 2 * (v2 & 1);
          return 0;
        }
        int64_t res = ipow(v1, v2);
        if (res > LONG_MAX) return 0;
        return (long) res;
      case AST_TOK_LEFT: return v1 << v2;
      case AST_TOK_RIGHT: return v1 >> v2;
      case AST_TOK_BAND: return v1 & v2;
      case AST_TOK_BXOR: return v1 ^ v2;
      case AST_TOK_BOR: return v1 | v2;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return 0;
    }
  }
}

/******************************************************************************
Function `ast_eval_float`:
  Evaluate the value in float type, given the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree;
  * `var`:      the user-supplied variable array.
Return:
  The resulting float number on success; HUGE_VALF on error.
******************************************************************************/
static float ast_eval_float(ast_t *ast, const ast_node_t *node,
    const float *var) {
  if (AST_IS_ERROR(ast)) return HUGE_VALF;
  if (node->type == AST_TOK_NUM) return node->value.v.fval;
  else if (node->type == AST_TOK_VAR) {
    if (var) return var[ast->vidx[node->value.v.lval] - 1];
    else return *((float *) ast->var + node->value.v.lval);
  }
  else if (ast_tok_attr[node->type].argc == 1) {
    const float v = ast_eval_float(ast, node->left, var);
    switch (node->type) {
      case AST_TOK_NEG: return -v;
      case AST_TOK_ABS: return fabsf(v);
      case AST_TOK_SQRT: return sqrtf(v);
      case AST_TOK_LN: return logf(v);
      case AST_TOK_LOG: return log10f(v);
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return HUGE_VALF;
    }
  }
  else {
    const float v1 = ast_eval_float(ast, node->left, var);
    const float v2 = ast_eval_float(ast, node->right, var);
    switch (node->type) {
      case AST_TOK_ADD: return v1 + v2;
      case AST_TOK_MINUS: return v1 - v2;
      case AST_TOK_MUL: return v1 * v2;
      case AST_TOK_DIV: return v1 / v2;
      case AST_TOK_EXP: return powf(v1, v2);
      case AST_TOK_REM: return fmodf(v1, v2);
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return HUGE_VALF;
    }
  }
}

/******************************************************************************
Function `ast_eval_double`:
  Evaluate the value in double type, given the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree;
  * `var`:      the user-supplied variable array.
Return:
  The resulting double number on success; HUGE_VAL on error.
******************************************************************************/
static double ast_eval_double(ast_t *ast, const ast_node_t *node,
    const double *var) {
  if (AST_IS_ERROR(ast)) return HUGE_VAL;
  if (node->type == AST_TOK_NUM) return node->value.v.dval;
  else if (node->type == AST_TOK_VAR) {
    if (var) return var[ast->vidx[node->value.v.lval] - 1];
    else return *((double *) ast->var + node->value.v.lval);
  }
  else if (ast_tok_attr[node->type].argc == 1) {
    const double v = ast_eval_double(ast, node->left, var);
    switch (node->type) {
      case AST_TOK_NEG: return -v;
      case AST_TOK_ABS: return fabs(v);
      case AST_TOK_SQRT: return sqrt(v);
      case AST_TOK_LN: return log(v);
      case AST_TOK_LOG: return log10(v);
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return HUGE_VAL;
    }
  }
  else {
    const double v1 = ast_eval_double(ast, node->left, var);
    const double v2 = ast_eval_double(ast, node->right, var);
    switch (node->type) {
      case AST_TOK_ADD: return v1 + v2;
      case AST_TOK_MINUS: return v1 - v2;
      case AST_TOK_MUL: return v1 * v2;
      case AST_TOK_DIV: return v1 / v2;
      case AST_TOK_EXP: return pow(v1, v2);
      case AST_TOK_REM: return fmod(v1, v2);
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return HUGE_VAL;
    }
  }
}

/******************************************************************************
Function `ast_eval_bool`:
  Evaluate the value in bool type, given the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree.
******************************************************************************/
static void ast_eval_bool(ast_t *ast, ast_node_t *node) {
  if (AST_IS_ERROR(ast)) return;
  if (ast_tok_attr[node->type].argc == 0) return;
  /* Unary operators. */
  if (ast_tok_attr[node->type].argc == 1) {
    /* The child node is not evaluated. */
    if (ast_tok_attr[node->left->type].argc != 0)
      ast_eval_bool(ast, node->left);
    if (AST_IS_ERROR(ast)) return;

    ast_var_t *v;
    if (node->left->type == AST_TOK_VAR)
      v = (ast_var_t *) ast->var + node->left->value.v.lval;
    else v = &node->left->value;
    const int dtype = v->dtype;
    bool bres;
    long lres;
    double dres;

    switch (node->type) {
      case AST_TOK_LNOT:
        if (dtype == AST_DTYPE_BOOL) bres = !v->v.bval;
        else if (dtype == AST_DTYPE_LONG) bres = !v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = !v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_NEG:
        if (dtype == AST_DTYPE_LONG) {
          lres = -v->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = -v->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_ABS:
        if (dtype == AST_DTYPE_LONG) {
          lres = (v->v.lval < 0) ? -v->v.lval : v->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = fabs(v->v.dval);
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_SQRT:
        if (dtype == AST_DTYPE_LONG) dres = (double) v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) dres = v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        dres = sqrt(dres);
        ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        return;
      case AST_TOK_LN:
        if (dtype == AST_DTYPE_LONG) dres = (double) v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) dres = v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        dres = log(dres);
        ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        return;
      case AST_TOK_LOG:
        if (dtype == AST_DTYPE_LONG) dres = (double) v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) dres = v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        dres = log10(dres);
        ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        return;
      case AST_TOK_ISFINITE:
        if (dtype == AST_DTYPE_FLOAT) bres = isfinite(v->v.fval) ? true : false;
        else if (dtype == AST_DTYPE_DOUBLE)
          bres = isfinite(v->v.dval) ? true : false;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_BNOT:
        if (dtype == AST_DTYPE_LONG) {
          lres = ~v->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
    }
  }
  /* Binary operators. */
  else {
    /* Evaluate child nodes. */
    if (ast_tok_attr[node->left->type].argc != 0)
      ast_eval_bool(ast, node->left);
    if (AST_IS_ERROR(ast)) return;
    if (ast_tok_attr[node->right->type].argc != 0)
      ast_eval_bool(ast, node->right);
    if (AST_IS_ERROR(ast)) return;

    ast_var_t *v1, *v2;
    if (node->left->type == AST_TOK_VAR)
      v1 = (ast_var_t *) ast->var + node->left->value.v.lval;
    else v1 = &node->left->value;
    if (node->right->type == AST_TOK_VAR)
      v2 = (ast_var_t *) ast->var + node->right->value.v.lval;
    else v2 = &node->right->value;

    int dtype = v1->dtype;
    /* Type cast for numerical types. */
    if (v1->dtype != v2->dtype) {
      if (v1->dtype == AST_DTYPE_LONG) {
        double tmp = (double) v1->v.lval;
        ast_set_var_value(v1, &tmp, 0, AST_DTYPE_DOUBLE);
        dtype = AST_DTYPE_DOUBLE;
      }
      else if (v2->dtype == AST_DTYPE_LONG) {
        double tmp = (double) v2->v.lval;
        ast_set_var_value(v2, &tmp, 0, AST_DTYPE_DOUBLE);
      }
      else {
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      }
    }
    bool bres;
    long lres;
    double dres;

    switch (node->type) {
      case AST_TOK_LAND:
        if (dtype == AST_DTYPE_BOOL) {
          bres = v1->v.bval && v2->v.bval;
          ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_LOR:
        if (dtype == AST_DTYPE_BOOL) {
          bres = v1->v.bval || v2->v.bval;
          ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_LT:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval < v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval < v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_LE:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval <= v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval <= v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_GT:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval > v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval > v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_GE:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval >= v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval >= v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_EQ:
        if (dtype == AST_DTYPE_BOOL) bres = v1->v.bval == v2->v.bval;
        else if (dtype == AST_DTYPE_LONG) bres = v1->v.lval == v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval == v2->v.dval;
        else if (dtype == AST_DTYPE_STRING) {
          if (v1->v.sval.len != v2->v.sval.len) bres = false;
          else bres = !strncmp(v1->v.sval.str, v2->v.sval.str, v1->v.sval.len);
        }
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_NEQ:
        if (dtype == AST_DTYPE_BOOL) bres = v1->v.bval != v2->v.bval;
        else if (dtype == AST_DTYPE_LONG) bres = v1->v.lval != v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval != v2->v.dval;
        else if (dtype == AST_DTYPE_STRING) {
          if (v1->v.sval.len != v2->v.sval.len) bres = true;
          else bres = strncmp(v1->v.sval.str, v2->v.sval.str, v1->v.sval.len);
        }
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        return;
      case AST_TOK_ADD:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval + v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval + v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_MINUS:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval - v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval - v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_MUL:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval * v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval * v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_DIV:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval / v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval / v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_EXP:
        if (dtype == AST_DTYPE_LONG) {
          lres = 0;
          if (v2->v.lval > UINT8_MAX) {
            if (v1->v.lval == 1) lres = 1;
            else if (v1->v.lval == -1) lres = 1 - 2 * (v2->v.lval & 1);
          }
          else {
            int64_t tmp = ipow(v1->v.lval, v2->v.lval);
            if (tmp <= LONG_MAX) lres = (long) tmp;
          }
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = pow(v1->v.dval, v2->v.dval);
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_REM:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval % v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = fmod(v1->v.dval, v2->v.dval);
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_LEFT:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval << v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_RIGHT:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval >> v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_BAND:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval & v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_BXOR:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval ^ v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      case AST_TOK_BOR:
        if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval | v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
    }
  }
}

/******************************************************************************
Function `ast_eval_pre`:
  Pre-evaluate values in the abstract syntax tree.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `node`:     a node of the abstract syntax tree.
******************************************************************************/
static void ast_eval_pre(ast_t *ast, ast_node_t *node) {
  if (AST_IS_ERROR(ast)) return;
  if (ast_tok_attr[node->type].argc == 0) return;
  /* Unary operators. */
  if (ast_tok_attr[node->type].argc == 1) {
    /* The child node is not evaluated. */
    if (ast_tok_attr[node->left->type].argc != 0)
      ast_eval_pre(ast, node->left);
    if (AST_IS_ERROR(ast)) return;

    if (node->left->type != AST_TOK_NUM && node->left->type != AST_TOK_STRING)
      return;

    ast_var_t *v = &node->left->value;
    const int dtype = v->dtype;
    bool bres;
    int ires;
    long lres;
    float fres;
    double dres;

    switch (node->type) {
      case AST_TOK_LNOT:
        if (dtype == AST_DTYPE_BOOL) bres = !v->v.bval;
        else if (dtype == AST_DTYPE_LONG) bres = !v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = !v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_NEG:
        if (dtype == AST_DTYPE_INT) {
          ires = -v->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = -v->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = -v->v.fval;
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = -v->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_ABS:
        if (dtype == AST_DTYPE_INT) {
          ires = (v->v.ival < 0) ? -v->v.ival: v->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = (v->v.lval < 0) ? -v->v.lval: v->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = fabsf(v->v.fval);
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = fabs(v->v.dval);
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_SQRT:
        if (dtype == AST_DTYPE_LONG) dres = (double) v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) dres = v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        dres = sqrt(dres);
        ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_LN:
        if (dtype == AST_DTYPE_LONG) dres = (double) v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) dres = v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        dres = log(dres);
        ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_LOG:
        if (dtype == AST_DTYPE_LONG) dres = (double) v->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) dres = v->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        dres = log10(dres);
        ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_ISFINITE:
        if (dtype == AST_DTYPE_FLOAT) bres = isfinite(v->v.fval) ? true : false;
        else if (dtype == AST_DTYPE_DOUBLE)
          bres = isfinite(v->v.dval) ? true : false;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_BNOT:
        if (dtype == AST_DTYPE_INT) {
          ires = ~v->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = ~v->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
    }
  }
  /* Binary operators. */
  else {
    /* Evaluate child nodes. */
    if (ast_tok_attr[node->left->type].argc != 0)
      ast_eval_pre(ast, node->left);
    if (AST_IS_ERROR(ast)) return;
    if (node->left->type != AST_TOK_NUM && node->left->type != AST_TOK_STRING)
      return;

    if (ast_tok_attr[node->right->type].argc != 0)
      ast_eval_pre(ast, node->right);
    if (AST_IS_ERROR(ast)) return;
    if (node->right->type != AST_TOK_NUM && node->right->type != AST_TOK_STRING)
      return;

    ast_var_t *v1 = &node->left->value;
    ast_var_t *v2 = &node->right->value;
    int dtype = v1->dtype;
    /* Type cast for numerical types. */
    if (v1->dtype != v2->dtype) {
      if (v1->dtype == AST_DTYPE_LONG) {
        double tmp = (double) v1->v.lval;
        ast_set_var_value(v1, &tmp, 0, AST_DTYPE_DOUBLE);
        dtype = AST_DTYPE_DOUBLE;
      }
      else if (v2->dtype == AST_DTYPE_LONG) {
        double tmp = (double) v2->v.lval;
        ast_set_var_value(v2, &tmp, 0, AST_DTYPE_DOUBLE);
      }
      else {
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
      }
    }
    bool bres;
    int ires;
    long lres;
    float fres;
    double dres;

    switch (node->type) {
      case AST_TOK_LAND:
        if (dtype == AST_DTYPE_BOOL) {
          bres = v1->v.bval && v2->v.bval;
          ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_LOR:
        if (dtype == AST_DTYPE_BOOL) {
          bres = v1->v.bval || v2->v.bval;
          ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_LT:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval < v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval < v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_LE:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval <= v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval <= v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_GT:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval > v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval > v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_GE:
        if (dtype == AST_DTYPE_LONG) bres = v1->v.lval >= v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval >= v2->v.dval;
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_EQ:
        if (dtype == AST_DTYPE_BOOL) bres = v1->v.bval == v2->v.bval;
        else if (dtype == AST_DTYPE_LONG) bres = v1->v.lval == v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval == v2->v.dval;
        else if (dtype == AST_DTYPE_STRING) {
          if (v1->v.sval.len != v2->v.sval.len) bres = false;
          else bres = !strncmp(v1->v.sval.str, v2->v.sval.str, v1->v.sval.len);
        }
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_NEQ:
        if (dtype == AST_DTYPE_BOOL) bres = v1->v.bval != v2->v.bval;
        else if (dtype == AST_DTYPE_LONG) bres = v1->v.lval != v2->v.lval;
        else if (dtype == AST_DTYPE_DOUBLE) bres = v1->v.dval != v2->v.dval;
        else if (dtype == AST_DTYPE_STRING) {
          if (v1->v.sval.len != v2->v.sval.len) bres = true;
          else bres = strncmp(v1->v.sval.str, v2->v.sval.str, v1->v.sval.len);
        }
        else {
          AST_ERRNO(ast) = AST_ERR_EVAL;
          return;
        }
        ast_set_var_value(&node->value, &bres, 0, AST_DTYPE_BOOL);
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_ADD:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival + v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval + v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = v1->v.fval + v2->v.fval;
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval + v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_MINUS:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival - v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval - v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = v1->v.fval - v2->v.fval;
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval - v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_MUL:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival * v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval * v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = v1->v.fval * v2->v.fval;
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval * v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_DIV:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival / v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval / v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = v1->v.fval / v2->v.fval;
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = v1->v.dval / v2->v.dval;
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_EXP:
        if (dtype == AST_DTYPE_INT) {
          ires = 0;
          if (v2->v.ival > UINT8_MAX) {
            if (v1->v.ival == 1) ires = 1;
            else if (v1->v.ival == -1) ires = 1 - 2 * (v2->v.ival & 1);
          }
          else {
            int64_t tmp = ipow(v1->v.ival, v2->v.ival);
            if (tmp <= INT_MAX) ires = (int) tmp;
          }
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = 0;
          if (v2->v.lval > UINT8_MAX) {
            if (v1->v.lval == 1) lres = 1;
            else if (v1->v.lval == -1) lres = 1 - 2 * (v2->v.lval & 1);
          }
          else {
            int64_t tmp = ipow(v1->v.lval, v2->v.lval);
            if (tmp <= LONG_MAX) lres = (long) tmp;
          }
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = powf(v1->v.fval, v2->v.fval);
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = pow(v1->v.dval, v2->v.dval);
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_REM:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival % v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval % v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else if (dtype == AST_DTYPE_FLOAT) {
          fres = fmod(v1->v.fval, v2->v.fval);
          ast_set_var_value(&node->value, &fres, 0, AST_DTYPE_FLOAT);
        }
        else if (dtype == AST_DTYPE_DOUBLE) {
          dres = fmod(v1->v.dval, v2->v.dval);
          ast_set_var_value(&node->value, &dres, 0, AST_DTYPE_DOUBLE);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_LEFT:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival << v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval << v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_RIGHT:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival >> v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval >> v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_BAND:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival & v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval & v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_BXOR:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival ^ v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval ^ v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      case AST_TOK_BOR:
        if (dtype == AST_DTYPE_INT) {
          ires = v1->v.ival | v2->v.ival;
          ast_set_var_value(&node->value, &ires, 0, AST_DTYPE_INT);
        }
        else if (dtype == AST_DTYPE_LONG) {
          lres = v1->v.lval | v2->v.lval;
          ast_set_var_value(&node->value, &lres, 0, AST_DTYPE_LONG);
        }
        else AST_ERRNO(ast) = AST_ERR_EVAL;
        node->type = AST_TOK_NUM;
        return;
      default:
        AST_ERRNO(ast) = AST_ERR_EVAL;
        return;
    }
  }
}


/*============================================================================*\
                    Interfaces for the parser and evaluator
\*============================================================================*/

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
    const bool eval) {
  if (!ast) return AST_ERR_INIT;
  if (ast->ast) return AST_ERRNO(ast) = AST_ERR_EXIST;

  if (!str || *(str = ast_skip_space(str)) == '\0')
    return AST_ERRNO(ast) = AST_ERR_STRING;
  if (!(ast->exp = ast_copy_str(str))) return AST_ERRNO(ast) = AST_ERR_MEMORY;

  if (dtype != AST_DTYPE_BOOL && dtype != AST_DTYPE_INT &&
      dtype != AST_DTYPE_LONG && dtype != AST_DTYPE_FLOAT &&
      dtype != AST_DTYPE_DOUBLE)
    return AST_ERRNO(ast) = AST_ERR_DTYPE;
  ast->dtype = dtype;

  /* Initialise the first node. */
  ast_var_t v = {0, .v.ival = 0};
  ast->ast = ast_create(AST_TOK_UNDEF, v);
  if (!ast->ast) return AST_ERRNO(ast) = AST_ERR_MEMORY;

  /* Parse the expression. */
  ast_node_t *node = (ast_node_t *) ast->ast;
  ast_parse_token(ast, node, ast->exp);
  if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);

  /* Redirect the root of the abstract syntax tree. */
  ast->ast = node = ast_root(node);

  /* Reset variable indices. */
  ast_reset_idx(ast, node);
  /* Free pre-allocated memory that is not necessary any more. */
  if (ast->nvar) {
    if (ast->nvar & (ast->nvar - 1)) {          /* nvar is not power of 2 */
      long *tmp = realloc(ast->vidx, ast->nvar * sizeof(long));
      if (tmp) ast->vidx = tmp;         /* do nothing if realloc fails */
    }
  }

  /* Allocate memory for variables. */
  if (ast->nvar) {
    ast->var = ast_init_var(ast, dtype);
    if (!ast->var) return AST_ERRNO(ast) = AST_ERR_MEMORY;
  }

  /* Validate data types for boolean expression. */
  if (dtype == AST_DTYPE_BOOL) {
    if (ast_check_dtype(ast, node) != AST_DTYPE_BOOL) {
      if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
      else return AST_ERRNO(ast) = AST_ERR_MISMATCH;
    }
  }

  /* Pre-evaluate values. */
  if (eval) {
    ast_eval_pre(ast, node);
    if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
  }

  return 0;
}

/******************************************************************************
Function `ast_eval`:
  Evaluate the expression given the abstract syntax tree and the variable array.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `value`:    address of the variable holding the evaluated value.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ast_eval(ast_t *ast, void *value) {
  if (!ast) return AST_ERR_INIT;
  if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
  if (!ast->ast) return AST_ERRNO(ast) = AST_ERR_NOEXP;
  if (!value) return AST_ERRNO(ast) = AST_ERR_VALUE;
  for (long i = 0; i < ast->nvar; i++) {
    if (!((ast_error_t *) ast->error)->vset[i]) {
      ast_msg(ast, "variable not set", ast->vidx[i], NULL);
      return AST_ERRNO(ast) = AST_ERR_VAR;
    }
  }

  ast_node_t *root = (ast_node_t *) ast->ast;
  switch (ast->dtype) {
    case AST_DTYPE_BOOL:
      ast_eval_bool(ast, root);
      if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
      *((bool *) value) = root->value.v.bval;
      break;
    case AST_DTYPE_INT:
      *((int *) value) = ast_eval_int(ast, root, NULL);
      if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
      break;
    case AST_DTYPE_LONG:
      *((long *) value) = ast_eval_long(ast, root, NULL);
      if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
      break;
    case AST_DTYPE_FLOAT:
      *((float *) value) = ast_eval_float(ast, root, NULL);
      if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
      break;
    case AST_DTYPE_DOUBLE:
      *((double *) value) = ast_eval_double(ast, root, NULL);
      if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
      break;
    default:
      return AST_ERRNO(ast) = AST_ERR_DTYPE;
  }
  return 0;
}

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
int ast_eval_num(ast_t *ast, void *value, const void *var, const long size) {
  if (!ast) return AST_ERR_INIT;
  if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
  if (!ast->ast) return AST_ERRNO(ast) = AST_ERR_NOEXP;
  if (!value) return AST_ERRNO(ast) = AST_ERR_VALUE;
  if (!var && size) return AST_ERRNO(ast) = AST_ERR_VAR;
  if (ast->nvar && size < ast->vidx[ast->nvar - 1])
    return AST_ERRNO(ast) = AST_ERR_SIZE;

  switch (ast->dtype) {
    case AST_DTYPE_INT:
      *((int *) value) =
        ast_eval_int(ast, (ast_node_t *) ast->ast, (int *) var);
      break;
    case AST_DTYPE_LONG:
      *((long *) value) =
        ast_eval_long(ast, (ast_node_t *) ast->ast, (long *) var);
      break;
    case AST_DTYPE_FLOAT:
      *((float *) value) =
        ast_eval_float(ast, (ast_node_t *) ast->ast, (float *) var);
      break;
    case AST_DTYPE_DOUBLE:
      *((double *) value) =
        ast_eval_double(ast, (ast_node_t *) ast->ast, (double *) var);
      break;
    default:
      return AST_ERRNO(ast) = AST_ERR_DTYPE;
  }
  if (AST_IS_ERROR(ast)) return AST_ERRNO(ast);
  return 0;
}


/*============================================================================*\
                          Function for error handling
\*============================================================================*/

/******************************************************************************
Function `ast_perror`:
  Print the error message if there is an error.
Arguments:
  * `ast`:      interface of the abstract syntax tree;
  * `fp`:       output file stream;
  * `msg`:      string to be printed before the error message.
******************************************************************************/
void ast_perror(const ast_t *ast, FILE *fp, const char *msg) {
  const char *sep, *errmsg;
  if (!ast) {
    if (!msg || *msg == '\0') msg = sep = "";
    else sep = " ";
    fprintf(fp, "%s%sthe abstract syntax tree is not initialised.\n", msg, sep);
    return;
  }

  if(!(AST_IS_ERROR(ast))) return;
  const ast_error_t *err = (ast_error_t *) ast->error;
  switch (AST_ERRNO(ast)) {
    case AST_ERR_MEMORY:
      errmsg = "failed to allocate memory";
      break;
    case AST_ERR_INIT:
      errmsg = "the abstract syntax tree is not initialised";
      break;
    case AST_ERR_STRING:
      errmsg = "invalid expression string";
      break;
    case AST_ERR_DTYPE:
      errmsg = "invalid data type for the expression";
      break;
    case AST_ERR_TOKEN:
      if (err->msg) errmsg = err->msg;
      else errmsg = "uncaught error of the expression";
      break;
    case AST_ERR_VAR:
      if (err->msg) errmsg = err->msg;
      else errmsg = "uncaught error of the variable";
      break;
    case AST_ERR_EXIST:
      errmsg = "the abstract syntax tree has already been built";
      break;
    case AST_ERR_NOEXP:
      errmsg = "the abstract syntax tree has not been built";
      break;
    case AST_ERR_VALUE:
      errmsg = "value for the evaluation is not set";
      break;
    case AST_ERR_SIZE:
      errmsg = "not enough elements in the variable array";
      break;
    case AST_ERR_EVAL:
      errmsg = "data type error for evaluation";
      break;
    case AST_ERR_NVAR:
      errmsg = "too many number of variables";
      break;
    case AST_ERR_MISMATCH:
      errmsg = "conflict data types in the expression";
      break;
    default:
      errmsg = "unknown error";
      break;
  }

  if (!msg || *msg == '\0') msg = sep = "";
  else sep = " ";
  fprintf(fp, "%s%s%s", msg, sep, errmsg);

  /* Print the specifier of the error location. */
  if (AST_ERRNO(ast) == AST_ERR_TOKEN && ast->exp && err->tpos) {
    fprintf(fp, "\n%s\n", ast->exp);
    for (ptrdiff_t i = 0; i < err->tpos - ast->exp; i++) fprintf(fp, " ");
    fprintf(fp, "^");
  }
  if (AST_ERRNO(ast) == AST_ERR_VAR && err->vidx) {
    if (err->vidx >= 0 && err->vidx < 10)
      fprintf(fp, ": %c%ld", AST_VAR_FLAG, err->vidx);
    else
      fprintf(fp, ": %c%c%ld%c",
          AST_VAR_FLAG, AST_VAR_START, err->vidx, AST_VAR_END);
  }
  fprintf(fp, "\n");
}


/*******************************************************************************
* read_ascii.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "read_file.h"
#include "ascii_fmtr.h"
#include "libast.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                  Definitions for libast compatible data types
\*============================================================================*/
#define AST_DTYPE_NULL          0               /* not supported by libast */
#define AST_DTYPE_STR_SPACE     (-AST_DTYPE_STRING)  /* string with spaces */
#define AST_DTYPE(x)            ((x < 0) ? -(x) : x) /* abs(x) */

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
#define CLEAN_ASCII_PTR                                                 \
  if (ast_real) {                                                       \
    for (int _i = 0; _i < nrcol; _i++) ast_destroy(ast_real[_i]);       \
    free(ast_real);                                                     \
  }                                                                     \
  ast_destroy(ast_sel);                                                 \
  ascii_arg_destroy(arg, nc); free(col);                                \
  if (chunk) free(chunk);                                               \
  if (fp) fclose(fp);                                                   \
  if (res) {                                                            \
    for (int _i = 0; _i < nrcol; _i++) if (res[_i]) free(res[_i]);      \
    free(res);                                                          \
  }

#ifdef OMP
  #define CLEAN_ASCII_PTR_OMP                                           \
    if (nomp > 1) {                                                     \
      for (int _i = 0; _i < nomp - 1; _i++) {                           \
        free(pcol[_i]); ast_destroy(ast_psel[_i]);                      \
      }                                                                 \
      for (int _i = 0; _i < (nomp - 1) * nrcol; _i++)                   \
        ast_destroy(ast_preal[_i]);                                     \
      free(pcol); free(ast_preal); free(ast_psel);                      \
    }                                                                   \
    free(pdata); free(pndata); free(lines);

  #define FCFC_ASCII_QUIT(x) {                                          \
    printf(FMT_FAIL);                                                   \
    P_EXT("failed to read the ASCII file\n");                           \
    exit(x);                                                            \
  }
#else
  #define CLEAN_ASCII_PTR_OMP
#endif

/*============================================================================*\
                        Data structure for ASCII columns
\*============================================================================*/

/* Structure for recording libast compatible data types. */
typedef struct {
  int dtype;            /* AST_DTYPE_NULL or libast compatible data type */
  union {               /* variable for storing the value */
    int ival; long lval; float fval; double dval;
    struct ast_string_struct_t { int len; const char *str; } sval;
  } v;
} ASC_COL;


/*============================================================================*\
                      Functions for parsing ASCII columns
\*============================================================================*/

/******************************************************************************
Function `ascii_col_init`:
  Initialise columns according to the arguments parsed from the formatter.
Arguments:
  * `arg`:      arguments parsed from the formatter;
  * `num`:      number of parsed arguments;
  * `rnum`:     number of arguments that are not suppressed.
Return:
  Pointer to the structure array for ASCII columns on success; NULL on error.
******************************************************************************/
static ASC_COL *ascii_col_init(asc_arg_s *arg, const int num,
    const int rnum) {
  if (rnum <= 0) return NULL;
  ASC_COL *col = malloc(rnum * sizeof(ASC_COL));
  if (!col) return NULL;

  int j = -1;
  for (int i = 0; i < num; i++) {
    if (arg[i].dtype == ASCII_DTYPE_SKIP) continue;
    if (++j >= rnum) {
      P_ERR("unknown error for identifying ASCII columns\n");
      free(col);
      return NULL;
    }
    /* Convert the `fscanf` type to libast type. */
    switch (arg[i].dtype) {
      case ASCII_DTYPE_INT: col[j].dtype = AST_DTYPE_INT; break;
      case ASCII_DTYPE_LONG: col[j].dtype = AST_DTYPE_LONG; break;
      case ASCII_DTYPE_FLT: col[j].dtype = AST_DTYPE_FLOAT; break;
      case ASCII_DTYPE_DBL: col[j].dtype = AST_DTYPE_DOUBLE; break;
      case ASCII_DTYPE_STR: col[j].dtype = AST_DTYPE_STRING; break;
      case ASCII_DTYPE_CHAR: col[j].dtype = AST_DTYPE_STRING; break;
      default:          /* data types that are not supported by libast */
        col[j].dtype = AST_DTYPE_NULL;
        P_WRN("unsupported formatter `%s', "
            "the corresponding column is not used\n", arg[i].fmtr);
        break;
    }

    /* Add a suppressing symbol '*' to the formatter string. */
    if (col[j].dtype == AST_DTYPE_NULL || col[j].dtype == AST_DTYPE_STRING) {
      int len = strlen(arg[i].fmtr) + 1;        /* fmtr is null terminated */
      char *tmp = realloc(arg[i].fmtr, (len + 1) * sizeof(char));
      if (!tmp) {
        free(col);
        return NULL;
      }
      arg[i].fmtr = tmp;

      /* Parse the formatter string to find the right place for '*'. */
      char *fmt = tmp;
      while (*fmt) {
        char c = *fmt++;
        if (c != '%') continue;
        c = *fmt++;
        if (c == '%') continue;
        /* Check if the formatter parses whitespaces. */
        if (col[j].dtype == AST_DTYPE_STRING && fmt - tmp > 2 && !isspace(*tmp))
          col[j].dtype = AST_DTYPE_STR_SPACE;

        /* Now add '*'. */
        memmove(fmt, fmt - 1, tmp + len - fmt + 1);
        *(fmt - 1) = '*';

        /* The formatter can still parse whitespaces with '[]'. */
        if (col[j].dtype == AST_DTYPE_STRING) {
          while (c != '[' && c != '\0') c = *fmt++;
          if (c == '[') {
            if (*fmt == '^') ++fmt;
            if (*fmt == ']') ++fmt;
            while ((c = *fmt++) != '\0' && c != ']' && c != ' ');
            if (c == ' ') col[j].dtype = AST_DTYPE_STR_SPACE;
          }
        }
        break;
      }
    }
  }
  return col;
}

/******************************************************************************
Function `ascii_read_line`:
  Read a string line to column variables.
Arguments:
  * `line`:     the line to be read;
  * `arg`:      arguments parsed from the formatter;
  * `num`:      number of columns to be read;
  * `col`:      variables and their data types for each column;
  * `end`:      pointer to the first character that is not interpreted.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int ascii_read_line(const char *line, const asc_arg_s *arg,
    const int num, ASC_COL *col, const char **end) {
  int j = 0;
  for (int i = 0; i < num; i++) {
    int n = 0;
    /* User-suppressed columns. */
    if (arg[i].dtype == ASCII_DTYPE_SKIP) {
      if (sscanf(line, arg[i].fmtr, &n) != 0 || n == 0) {
        *end = line;
        return FCFC_ERR_ASCII;
      }
    }
    else {
      switch (col[j].dtype) {
        case AST_DTYPE_NULL:
          /* Suppressed columns due to compatibility with libast. */
          if (sscanf(line, arg[i].fmtr, &n) != 0 || n == 0) {
            *end = line;
            return FCFC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_INT:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.ival), &n) != 1 || n == 0) {
            *end = line;
            return FCFC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_LONG:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.lval), &n) != 1 || n == 0) {
            *end = line;
            return FCFC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_FLOAT:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.fval), &n) != 1 || n == 0) {
            *end = line;
            return FCFC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_DOUBLE:
          if (sscanf(line, arg[i].fmtr, &(col[j].v.dval), &n) != 1 || n == 0) {
            *end = line;
            return FCFC_ERR_ASCII;
          }
          break;
        case AST_DTYPE_STRING:
          while (isspace(*line)) line++;        /* omit whitespaces */
        case AST_DTYPE_STR_SPACE:
          /* Do not actually save the string, but compute the position. */
          if (sscanf(line, arg[i].fmtr, &n) != 0 || n == 0) {
            *end = line;
            return FCFC_ERR_ASCII;
          }
          col[j].v.sval.str = line;
          col[j].v.sval.len = n;
          break;
        default:
          P_ERR("unknown data type for column: ${%d}\n", j + 1);
          *end = line;
          return FCFC_ERR_UNKNOWN;
      }
      j++;
    }
    line += n;
  }
  return 0;
}

/******************************************************************************
Function `ascii_apply_sel`:
  Apply the selection criteria based on the AST and ASCII columns.
Arguments:
  * `ast`:      abstract syntax tree for the selection criteria;
  * `col`:      columns to be parsed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int ascii_apply_sel(ast_t *ast, const ASC_COL *col, bool *res) {
  for (int i = 0; i < ast->nvar; i++) {
    long idx = ast->vidx[i];
    const ASC_COL *c = col + (idx - 1);
    switch (AST_DTYPE(c->dtype)) {
      case AST_DTYPE_INT:
        if (ast_set_var(ast, idx, &c->v.ival, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_LONG:
        if (ast_set_var(ast, idx, &c->v.lval, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_FLOAT:
        if (ast_set_var(ast, idx, &c->v.fval, 0, AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_DOUBLE:
        if (ast_set_var(ast, idx, &c->v.dval, 0, AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_STRING:
        if (ast_set_var(ast, idx, c->v.sval.str, c->v.sval.len,
              AST_DTYPE_STRING)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      default:
        P_ERR("column ${%ld} not appropriate for selection\n", idx);
        return FCFC_ERR_AST;
    }
  }

  if (ast_eval(ast, res)) {
    P_AST_ERR(ast);
    return FCFC_ERR_AST;
  }
  return 0;
}

/******************************************************************************
Function `ascii_read_real`:
  Read a floating-point number based on the AST and ASCII columns.
Arguments:
  * `ast`:      abstract syntax tree for the number evaluation;
  * `col`:      columns to be parsed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int ascii_read_real(ast_t *ast, const ASC_COL *col, real *res) {
  for (int i = 0; i < ast->nvar; i++) {
    long idx = ast->vidx[i];
    const ASC_COL *c = col + (idx - 1);
    switch (AST_DTYPE(c->dtype)) {
      case AST_DTYPE_INT:
        if (ast_set_var(ast, idx, &c->v.ival, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast);
          return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_LONG:
        if (ast_set_var(ast, idx, &c->v.lval, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast);
          return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_FLOAT:
        if (ast_set_var(ast, idx, &c->v.fval, 0, AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast);
          return FCFC_ERR_AST;
        }
        break;
      case AST_DTYPE_DOUBLE:
        if (ast_set_var(ast, idx, &c->v.dval, 0, AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast);
          return FCFC_ERR_AST;
        }
        break;
      default:
        P_ERR("column ${%ld} not appropriate for the numerical evaluation\n",
            idx);
        return FCFC_ERR_AST;
    }
  }

  if (ast_eval(ast, res)) {
    P_AST_ERR(ast);
    return FCFC_ERR_AST;
  }
  return 0;
}


/*============================================================================*\
                      Functions for reading file by chunks
\*============================================================================*/

/******************************************************************************
Function `chunk_resize`:
  Enlarge the size of a chunk.
Arguments:
  * `chunk`:    address of the chunk;
  * `size`:     size of the chunk.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int chunk_resize(char **chunk, size_t *size) {
  /* Assume the arguments are not NULL. */
  size_t num;
  if (!(*chunk)) num = FCFC_FILE_CHUNK;
  else {
    if (FCFC_MAX_CHUNK / 2 < *size) return FCFC_ERR_FILE;
    num = *size << 1;
  }

  char *tmp = realloc(*chunk, num * sizeof(char));
  if (!tmp) return FCFC_ERR_MEMORY;

  *chunk = tmp;
  *size = num;
  return 0;
}

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
int read_ascii_table(const char *fname, double **x, double **y, size_t *num) {
  /* Validate arguments. */
  if (!fname) {
    P_ERR("the filename for reading ASCII table is not specified\n");
    return FCFC_ERR_ARG;
  }
  if (!x || !y || !num) {
    P_ERR("ASCII table reader not initialized\n");
    return FCFC_ERR_ARG;
  }

  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    return FCFC_ERR_FILE;
  }

  /* Prepare for the chunk. */
  char *chunk = NULL;
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the file by chunks\n");
    fclose(fp);
    return FCFC_ERR_MEMORY;
  }

  /* Allocate memory for the data. */
  size_t max = FCFC_DATA_INIT_NUM;
  double *nx, *ny;
  nx = ny = NULL;
  if (!(nx = malloc(max * sizeof(double))) ||
      !(ny = malloc(max * sizeof(double)))) {
    P_ERR("failed to allocate memory for the samples\n");
    fclose(fp); free(chunk);
    if (nx) free(nx);
    if (ny) free(ny);
    return FCFC_ERR_MEMORY;
  }

  size_t n, nread, nrest;
  n = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == FCFC_READ_COMMENT || *p == '\0') {   /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      if (sscanf(p, "%lf %lf", nx + n, ny + n) != 2) {
        P_ERR("failed to read line: %s\n", p);
        fclose(fp); free(chunk); free(nx); free(ny);
        return FCFC_ERR_FILE;
      }

      /* Enlarge the memory for the data if necessary. */
      if (++n >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many samples in the file: `%s'\n", fname);
          fclose(fp); free(chunk); free(nx); free(ny);
          return FCFC_ERR_FILE;
        }
        max <<= 1;
        double *tmp = realloc(nx, sizeof(double) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the samples\n");
          fclose(fp); free(chunk); free(nx); free(ny);
          return FCFC_ERR_MEMORY;
        }
        nx = tmp;
        tmp = realloc(ny, sizeof(double) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the samples\n");
          fclose(fp); free(chunk); free(nx); free(ny);
          return FCFC_ERR_MEMORY;
        }
        ny = tmp;
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        fclose(fp); free(chunk); free(nx); free(ny);
        return FCFC_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    fclose(fp); free(chunk); free(nx); free(ny);
    return FCFC_ERR_FILE;
  }

  free(chunk);
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  *x = nx;
  *y = ny;
  *num = n;
  return 0;
}

/******************************************************************************
Function `read_ascii_data`:
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
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_data(const char *fname, const size_t skip, const char comment,
    const char *fmtr, char *const *rcol_ids, const int nrcol, const char *sel,
    real ***rout, size_t *num, const int verb) {
  /* Validate arguments. */
  if (!fname) {
    P_ERR("the filename for ASCII data is not specified\n");
    return FCFC_ERR_ARG;
  }
  if (!fmtr || !rcol_ids || !rout || !num) {
    P_ERR("ASCII data reader not initialized\n");
    return FCFC_ERR_ARG;
  }
  if (nrcol < 3 || nrcol > 4) {
    P_ERR("unexpected number of floating-point columns to be read: %d\n",
        nrcol);
    return FCFC_ERR_ARG;
  }

  /* Parse the formatter. */
  asc_arg_s *arg;
  int nc, nrc;  /* total number of columns and number of columns to be read */
  nc = nrc = 0;
  if (!(arg = parse_ascii_fmtr(fmtr, &nc, &nrc))) return FCFC_ERR_ASCII;
  if (nrc <= 0) {
    P_ERR("no column to be read given the formatter: `%s'\n", fmtr);
    return FCFC_ERR_ASCII;
  }
  if (nrc < 3) P_WRN("reading coordinates from less than 3 columns\n");

  /* Record libast compatible columns. */
  ASC_COL *col;
  if (!(col = ascii_col_init(arg, nc, nrc))) {
    ascii_arg_destroy(arg, nc);
    return FCFC_ERR_ASCII;
  }

  /* Initialise all variables for easy error handling. */
  FILE *fp = NULL;
  char *chunk = NULL;           /* chunk for file reading */
  real **res = NULL;            /* arrays for the outputs */
  ast_t **ast_real = NULL;      /* ASTs for the outputs   */
  ast_t *ast_sel = NULL;        /* AST for the selection  */

  /* Construct ASTs for the columns. */
  if (!(ast_real = malloc(sizeof(ast_t *) * nrcol))) {
    P_ERR("failed to allocate memory for the Abstract Syntax Trees\n");
    CLEAN_ASCII_PTR;
    return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < nrcol; i++) ast_real[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
    if (!(ast_real[i] = ast_init())) {
      P_AST_ERR(ast_real[i]); CLEAN_ASCII_PTR; return FCFC_ERR_AST;
    }
    if (ast_build(ast_real[i], rcol_ids[i], FCFC_AST_REAL, true)) {
      P_AST_ERR(ast_real[i]); CLEAN_ASCII_PTR; return FCFC_ERR_AST;
    }
  }
  if (sel && *sel) {
    if (!(ast_sel = ast_init())) {
      P_AST_ERR(ast_sel); CLEAN_ASCII_PTR; return FCFC_ERR_AST;
    }
    if (ast_build(ast_sel, sel, AST_DTYPE_BOOL, true)) {
      P_AST_ERR(ast_sel); CLEAN_ASCII_PTR; return FCFC_ERR_AST;
    }
  }

  /* Check number of variables for expressions. */
  int max_col = 0;
  for (int i = 0; i < nrcol; i++) {
    if (ast_real[i]->nvar == 0)
      P_WRN("the column expression indicates a constant: `%s'\n", rcol_ids[i]);
    else {
      if (nrc < ast_real[i]->vidx[ast_real[i]->nvar - 1]) {
        P_ERR("not enough columns for expression: `%s'\n", rcol_ids[i]);
        CLEAN_ASCII_PTR; return FCFC_ERR_CFG;
      }
      if (max_col < ast_real[i]->vidx[ast_real[i]->nvar - 1])
        max_col = ast_real[i]->vidx[ast_real[i]->nvar - 1];
    }
  }
  if (ast_sel) {
    if (ast_sel->nvar == 0) {
      P_WRN("the expression for data selection is a constant: `%s'\n", sel);
    }
    if (nrc < ast_sel->vidx[ast_sel->nvar - 1]) {
      P_ERR("not enough columns for data selection: `%s'\n", sel);
      CLEAN_ASCII_PTR; return FCFC_ERR_CFG;
    }
    if (max_col < ast_sel->vidx[ast_sel->nvar - 1])
      max_col = ast_sel->vidx[ast_sel->nvar - 1];
  }

  /* Remove columns that are not needed. */
  if (max_col < nrc) {
    int i, j;
    for (i = j = 0; i < nc; i++) {
      if (arg[i].dtype != ASCII_DTYPE_SKIP) {
        if (++j == max_col) break;
      }
    }
    for (j = i + 1; j < nc; j++) if(arg[j].fmtr) free(arg[j].fmtr);
    nc = i + 1;
    nrc = max_col;
  }
  if (nrc <= 0) {
    P_ERR("no column to be read given the formatter: `%s'\n", fmtr);
    CLEAN_ASCII_PTR; return FCFC_ERR_UNKNOWN;
  }

  /* Prepare for the chunk. */
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the file by chunk\n");
    CLEAN_ASCII_PTR; return FCFC_ERR_MEMORY;
  }

  /* Allocate memory for the data. */
  size_t max = FCFC_DATA_INIT_NUM;
  if (!(res = malloc(sizeof(real *) * nrcol))) {
    P_ERR("failed to allocate memory for the data\n");
    CLEAN_ASCII_PTR; return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < nrcol; i++) res[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
    if (!(res[i] = malloc(max * sizeof(real)))) {
      P_ERR("failed to allocate memory for the data\n");
      CLEAN_ASCII_PTR; return FCFC_ERR_MEMORY;
    }
  }

  /* Dynamic allocations for OpenMP threads. */
#ifdef OMP
  const int nomp = omp_get_max_threads();
  /* Construct the ASCII columns for non-master threads. */
  ASC_COL **pcol = NULL;
  if (nomp > 1) {
    if (!(pcol = malloc(sizeof(ASC_COL *) * (nomp - 1)))) {
      P_ERR("failed to allocate memory for thread-private columns\n");
      FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
    }
    for (int j = 0; j < nomp - 1; j++) {
      if (!(pcol[j] = malloc(sizeof(ASC_COL) * nrc))) {
        P_ERR("failed to allocate memory for thread-private columns\n");
        FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
      }
      memcpy(pcol[j], col, sizeof(ASC_COL) * nrc);
    }
  }
  /* Construct the ASTs again for non-master threads. */
  ast_t **ast_preal, **ast_psel;
  ast_preal = ast_psel = NULL;
  if (nomp > 1) {
    if (!(ast_preal = malloc(sizeof(ast_t *) * (nomp - 1) * nrcol)) ||
        !(ast_psel = malloc(sizeof(ast_t *) * (nomp - 1)))) {
      P_ERR("failed to allocate memory for thread-private ASTs\n");
      FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
    }
  #pragma omp parallel for num_threads(nomp)
    for (int j = 0; j < nomp - 1; j++) {
      /* Positions and weights. */
      for (int i = 0; i < nrcol; i++) {
        int k = j * nrcol + i;
        if (!(ast_preal[k] = ast_init())) {
          P_AST_ERR(ast_preal[k]); FCFC_ASCII_QUIT(FCFC_ERR_AST);
        }
        if (ast_build(ast_preal[k], rcol_ids[i], FCFC_AST_REAL, true)) {
          P_AST_ERR(ast_preal[k]); FCFC_ASCII_QUIT(FCFC_ERR_AST);
        }
      }
      /* Selections. */
      if (sel && *sel) {
        if (!(ast_psel[j] = ast_init())) {
          P_AST_ERR(ast_psel[j]); FCFC_ASCII_QUIT(FCFC_ERR_AST);
        }
        if (ast_build(ast_psel[j], sel, AST_DTYPE_BOOL, true)) {
          P_AST_ERR(ast_psel[j]); FCFC_ASCII_QUIT(FCFC_ERR_AST);
        }
      }
      else ast_psel[j] = NULL;
    }
  }
  /* Construct the private data pool. */
  real *pdata = malloc(sizeof(real) * nomp * nrcol * FCFC_DATA_THREAD_NUM);
  size_t *pndata = calloc(nomp, sizeof(size_t));
  if (!pdata || !pndata) {
    P_ERR("failed to allocate memory for the thread-private data\n");
    FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
  }
  /* Construct the pool for file lines. */
  size_t nlmax = FCFC_DATA_INIT_NUM;
  size_t nl = 0;
  char **lines = malloc(sizeof(char *) * nlmax);
  if (!lines) {
    P_ERR("failed to allocate memory for the thread-private lines\n");
    FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
  }
#endif

  /* Open the file for reading. */
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("failed to open ASCII file for reading: `%s'\n", fname);
    CLEAN_ASCII_PTR; CLEAN_ASCII_PTR_OMP; return FCFC_ERR_FILE;
  }
  if (verb) printf("  Filename: %s\n", fname);

  size_t n, nline, nread, nrest;
  n = nline = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      /* Skip header lines. */
      if (nline++ < skip) {
        p = endl + 1; continue;
      }
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == comment || *p == '\0') {        /* comment or empty */
        p = endl + 1; continue;
      }

#ifdef OMP
      /* Append the line to the pool. */
      lines[nl++] = p;
      /* Enlarge the pool if necessary. */
      if (nl >= nlmax) {
        if (SIZE_MAX / 2 < nlmax) {
          P_ERR("too many lines in the file\n");
          FCFC_ASCII_QUIT(FCFC_ERR_FILE);
        }
        nlmax <<= 1;
        char **tmp = realloc(lines, sizeof(char *) * nlmax);
        if (!tmp) {
          P_ERR("failed to allocate memory for the thread-private lines\n");
          FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
        }
        lines = tmp;
      }
#else
      /* Parse the line. */
      const char *stop = NULL;
      if (ascii_read_line(p, arg, nc, col, &stop)) {
        P_ERR("failed to read the line with format `%s':\n", fmtr);
        fprintf(stderr, "%s\n", p);
        if (stop > p) for (int k = 0; k < stop - p; k++) fprintf(stderr, " ");
        fprintf(stderr, "^\n");
        CLEAN_ASCII_PTR; return FCFC_ERR_ASCII;
      }

      /* Apply selection. */
      if (ast_sel) {
        bool keep = false;
        if (ascii_apply_sel(ast_sel, col, &keep)) {
          CLEAN_ASCII_PTR; return FCFC_ERR_AST;
        }
        if (!keep) {
          p = endl + 1; continue;
        }
      }

      /* Record coordinates and weights. */
      for (int i = 0; i < nrcol; i++) {
        if (ascii_read_real(ast_real[i], col, res[i] + n)) {
          CLEAN_ASCII_PTR; return FCFC_ERR_AST;
        }
      }

      /* Enlarge the memory for the data if necessary. */
      if (++n >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many objects in the file: `%s'\n", fname);
          CLEAN_ASCII_PTR;
          return FCFC_ERR_FILE;
        }
        max <<= 1;
        for (int i = 0; i < nrcol; i++) {
          real *tmp = realloc(res[i], sizeof(real) * max);
          if (!tmp) {
            P_ERR("failed to allocate memory for the data\n");
            CLEAN_ASCII_PTR;
            return FCFC_ERR_MEMORY;
          }
          res[i] = tmp;
        }
      }
#endif
      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        CLEAN_ASCII_PTR; CLEAN_ASCII_PTR_OMP;
        return FCFC_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

#ifdef OMP
    /* Enlarge the memory for the data if necessary. */
    if (nl + n > max) {
      while (nl + n > max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many objects in the file: `%s'.\n", fname);
          CLEAN_ASCII_PTR;
          return FCFC_ERR_FILE;
        }
        max <<= 1;
      }
      for (int i = 0; i < nrcol; i++) {
        real *tmp = realloc(res[i], sizeof(real) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the data.\n");
          CLEAN_ASCII_PTR;
          return FCFC_ERR_MEMORY;
        }
        res[i] = tmp;
      }
    }

  #pragma omp parallel num_threads(nomp) firstprivate(ast_real,ast_sel,col)
    {
      /* Redirect pointers to the private pools. */
      const int tid = omp_get_thread_num();
      if (nomp > 1 && tid > 0) {        /* tid == 0 if nomp == 1 */
        ast_real = ast_preal + (tid - 1) * nrcol;
        ast_sel = ast_psel[tid - 1];
        col = pcol[tid - 1];
      }
      real *pres = pdata + tid * nrcol * FCFC_DATA_THREAD_NUM;
      size_t *pnum = pndata + tid;
      /* Process lines in parallel. */
  #pragma omp for
      for (size_t ii = 0; ii < nl; ii++) {
        /* Parse a line in the pool. */
        const char *pp = lines[ii];
        const char *stop = NULL;
        if (ascii_read_line(pp, arg, nc, col, &stop)) {
  #pragma omp critical
          {
            P_ERR("failed to read the line with format `%s':\n", fmtr);
            fprintf(stderr, "%s\n", pp);
            if (stop > pp) for (int k = 0; k < stop - pp; k++)
              fprintf(stderr, " ");
            fprintf(stderr, "^\n");
          }
          FCFC_ASCII_QUIT(FCFC_ERR_ASCII);
        }

        /* Apply selection. */
        if (ast_sel) {
          bool keep = false;
          if (ascii_apply_sel(ast_sel, col, &keep)) {
            FCFC_ASCII_QUIT(FCFC_ERR_AST);
          }
          if (!keep) continue;
        }

        /* Record coordinates and weights in the private data pool. */
        for (int i = 0; i < nrcol; i++) {
          if (ascii_read_real(ast_real[i], col,
              pres + i * FCFC_DATA_THREAD_NUM + *pnum)) {
            FCFC_ASCII_QUIT(FCFC_ERR_AST);
          }
        }

        /* Record the private data and clear the pool if necessary. */
        if (++(*pnum) >= FCFC_DATA_THREAD_NUM) {
  #pragma omp critical
          {
            /* Enlarge the memory for the data if necessary. */
            if (n + FCFC_DATA_THREAD_NUM >= max) {
              if (SIZE_MAX / 2 < max) {
                P_ERR("too many objects in the file: `%s'.\n", fname);
                FCFC_ASCII_QUIT(FCFC_ERR_FILE);
              }
              max <<= 1;
              if (max < n + FCFC_DATA_THREAD_NUM)
                max = n + FCFC_DATA_THREAD_NUM;
              for (int i = 0; i < nrcol; i++) {
                real *tmp = realloc(res[i], sizeof(real) * max);
                if (!tmp) {
                  P_ERR("failed to allocate memory for the data.\n");
                  FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
                }
                res[i] = tmp;
              }
            }
            for (int i = 0; i < nrcol; i++)
              memcpy(res[i] + n, pres + i * FCFC_DATA_THREAD_NUM,
                  sizeof(real) * FCFC_DATA_THREAD_NUM);
            n += FCFC_DATA_THREAD_NUM;
          }
          *pnum = 0;
        }
      }
    }
    nl = 0;
#endif

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

#ifdef OMP
  /* Record the rest of the private data. */
  for (int i = 0; i < nomp; i++) {
    if (pndata[i]) {
      /* Enlarge the memory for the data if necessary. */
      if (n + pndata[i] >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many objects in the file: `%s'.\n", fname);
          FCFC_ASCII_QUIT(FCFC_ERR_FILE);
        }
        max <<= 1;
        if (max < n + pndata[i]) max = n + pndata[i];
        for (int k = 0; k < nrcol; k++) {
          real *tmp = realloc(res[k], sizeof(real) * max);
          if (!tmp) {
            P_ERR("failed to allocate memory for the data.\n");
            FCFC_ASCII_QUIT(FCFC_ERR_MEMORY);
          }
          res[k] = tmp;
        }
      }
      for (int k = 0; k < nrcol; k++) {
        memcpy(res[k] + n, pdata + (k + i * nrcol) * FCFC_DATA_THREAD_NUM,
            sizeof(real) * pndata[i]);
      }
      n += pndata[i];
    }
  }
  /* Release memory for thread-private data structures. */
  CLEAN_ASCII_PTR_OMP;
#endif

  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    CLEAN_ASCII_PTR; return FCFC_ERR_FILE;
  }
  if (!n) {
    P_ERR("no object is read from file: `%s'\n", fname);
    CLEAN_ASCII_PTR; return FCFC_ERR_FILE;
  }

  if (verb) {
    printf("  Number of lines processed in total: %zu\n"
        "  Number of recorded objects: %zu\n", nline, n);
  }

  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);
  for (int i = 0; i < nrcol; i++) ast_destroy(ast_real[i]);
  free(ast_real);
  ast_destroy(ast_sel);
  ascii_arg_destroy(arg, nc); free(col); free(chunk);

  /* Release pre-allocated memory that is not necessary any more. */
#if     FCFC_SIMD  !=  FCFC_SIMD_NONE
  if (max != n + FCFC_NUM_REAL) {
    max = n + FCFC_NUM_REAL;
    for (int i = 0; i < nrcol; i++) {
      real *tmp = realloc(res[i], max * sizeof(real));
      if (!tmp) {
        P_ERR("failed to adjust memory for the data\n");
        return FCFC_ERR_MEMORY;
      }
      res[i] = tmp;
      memset(res[i] + n, 0, sizeof(real) * FCFC_NUM_REAL);
    }
  }
#else
  if (max != n) {
    for (int i = 0; i < nrcol; i++) {
      real *tmp = realloc(res[i], n * sizeof(real));
      if (tmp) res[i] = tmp;
    }
  }
#endif

  *num = n;
  *rout = res;
  return 0;
}

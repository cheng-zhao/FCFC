/*******************************************************************************
* read_fits.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#ifdef WITH_CFITSIO

#include "read_file.h"
#include "libast.h"
#include <fitsio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
#define FITS_ABORT {                                                    \
  P_ERR("cfitsio error: ");                                             \
  fits_report_error(stderr, status);                                    \
  status = 0;                                                           \
  if (fp) fits_close_file(fp, &status);                                 \
  return FCFC_ERR_FILE;                                                 \
}

#define CLEAN_FITS_PTR                                                  \
  if (fp) fits_close_file(fp, &status);                                 \
  fts_col_destroy(fcol);                                                \
  if (exp) {                                                            \
    for (int _i = 0; _i <= nrcol; _i++) if (exp[_i]) free(exp[_i]);     \
    free(exp);                                                          \
  }                                                                     \
  if (res) {                                                            \
    for (int _i = 0; _i < nrcol; _i++) if (res[_i]) free(res[_i]);      \
    free(res);                                                          \
  }                                                                     \
  if (ast_real) {                                                       \
    for (int _i = 0; _i < nrcol; _i++) ast_destroy(ast_real[_i]);       \
    free(ast_real);                                                     \
  }                                                                     \
  ast_destroy(ast_sel);

#ifdef OMP
  #define CLEAN_FITS_PTR_OMP                                            \
    if (nomp > 1) {                                                     \
      for (int _i = 0; _i < nomp - 1; _i++)                             \
        ast_destroy(ast_psel[_i]);                                      \
      for (int _i = 0; _i < (nomp - 1) * nrcol; _i++)                   \
        ast_destroy(ast_preal[_i]);                                     \
      free(ast_preal); free(ast_psel);                                  \
    }                                                                   \
    free(pdata); free(pndata);

  #define FCFC_FITS_QUIT(x) {                                           \
    printf(FMT_FAIL);                                                   \
    P_EXT("failed to read the FITS file\n");                            \
    exit(x);                                                            \
  }
#else
  #define CLEAN_FITS_PTR_OMP
#endif

/*============================================================================*\
                        Data structure for FITS columns
\*============================================================================*/

typedef struct {
  int ncol;             /* number of columns to be processed       */
  int *colnum;          /* column numbers in the FITS file         */
  int *dtype;           /* FITS data type of columns               */
  long *width;          /* widths of columns in bytes              */
  void **val;           /* pointer to the data read from fits file */
  char ***str;          /* pointer to strings read from fits file  */
} FTS_COL;


/*============================================================================*\
              Functions for manipulating the FITS column structure
\*============================================================================*/

/******************************************************************************
Function `fts_col_init`:
  Initialize the structure for storing FITS column information.
Return:
  Address of the structure on success; NULL on error.
******************************************************************************/
static inline FTS_COL *fts_col_init(void) {
  FTS_COL *fcol = calloc(1, sizeof(FTS_COL));
  if (!fcol) {
    P_ERR("failed to allocate memory for storing FITS column information\n");
    return NULL;
  }

  fcol->colnum = fcol->dtype = NULL;
  fcol->width = NULL;
  fcol->val = NULL;
  fcol->str = NULL;
  return fcol;
}

/******************************************************************************
Function `fts_col_insert`:
  Insert a column to the FITS column structure.
Arguments:
  * `fcol`:     structure for storing FITS column information;
  * `colnum`:   column number of the column to be appended.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int fts_col_insert(FTS_COL *fcol, const int colnum) {
  /* Check if there is enough allocated space. */
  if ((fcol->ncol & (fcol->ncol - 1)) == 0) {   /* ncol is 0 or power of 2 */
    if (INT_MAX / 2 < fcol->ncol) {
      P_ERR("too many columns required by the expressions\n");
      return FCFC_ERR_CFG;
    }
    int size = (fcol->ncol) ? fcol->ncol << 1 : 1;   /* double the size */
    int *tmp = realloc(fcol->colnum, size * sizeof(int));
    if (!tmp) {
      P_ERR("failed to allocate memory for recording FITS columns\n");
      return FCFC_ERR_MEMORY;
    }
    fcol->colnum = tmp;
  }

  fcol->colnum[fcol->ncol++] = colnum;
  return 0;
}

/******************************************************************************
Function `fts_col_destroy`:
  Deconstruct the structure for storing FITS column information.
Arguments:
  * `fcol`:     structure for storing FITS column information.
******************************************************************************/
static void fts_col_destroy(FTS_COL *fcol) {
  if (!fcol) return;
  if (fcol->colnum) free(fcol->colnum);
  if (fcol->dtype) free(fcol->dtype);
  if (fcol->width) free(fcol->width);
  if (fcol->val) {
    for (int i = 0; i < fcol->ncol; i++)
      if (fcol->val[i]) free(fcol->val[i]);
    free(fcol->val);
  }
  if (fcol->str) {
    for (int i = 0; i < fcol->ncol; i++) {
      if (fcol->str[i]) {
        free(fcol->str[i][0]);
        free(fcol->str[i]);
      }
    }
    free(fcol->str);
  }
  free(fcol);
}

/******************************************************************************
Function `fts_col_finalize`:
  Complete the structure for storing FITS column information.
Arguments:
  * `fp`:       pointer to the opened FITS table;
  * `fcol`:     structure for storing FITS column information;
  * `nrow`:     number of rows allocated for each column.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int fts_col_finalize(fitsfile *fp, FTS_COL *fcol, const long nrow) {
  if (!fcol->ncol) {
    P_ERR("no column to be read from the FITS file\n");
    return FCFC_ERR_CFG;
  }
  /* Reduce memory cost if applicable. */
  if ((fcol->ncol & (fcol->ncol - 1)) != 0) {
    int *tmp = realloc(fcol->colnum, fcol->ncol * sizeof(int));
    if (tmp) fcol->colnum = tmp;
  }

  /* Allocate memory for the rest column information. */
  if (!(fcol->dtype = calloc(fcol->ncol, sizeof(int))) ||
      !(fcol->width = calloc(fcol->ncol, sizeof(long))) ||
      !(fcol->val = malloc(fcol->ncol * sizeof(void *)))) {
    P_ERR("failed to allocate memory for storing FITS column information\n");
    return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < fcol->ncol; i++) fcol->val[i] = NULL;
  if (!(fcol->str = malloc(fcol->ncol * sizeof(char **)))) {
    P_ERR("failed to allocate memory for storing FITS column information\n");
    return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < fcol->ncol; i++) fcol->str[i] = NULL;

  /* Read column information. */
  for (int i = 0; i < fcol->ncol; i++) {
    int status = 0;
    long repeat = 0;
    if (fits_get_coltype(fp, fcol->colnum[i], fcol->dtype + i, &repeat,
        fcol->width + i, &status)) {
      P_ERR("cfitsio error: ");
      fits_report_error(stderr, status);
      return FCFC_ERR_FILE;
    }
    if ((fcol->dtype[i] == TSTRING && repeat != fcol->width[i]) ||
        (fcol->dtype[i] != TSTRING && repeat != 1)) {
      P_ERR("subcolumn of FITS table is not supported\n"
          "Please contact the author if you need this feature\n");
      return FCFC_ERR_FITS;
    }

    /* Check the column type and allocate memory for reading data. */
    switch (fcol->dtype[i]) {
      case TLOGICAL:
        if (!(fcol->val[i] = malloc(nrow * sizeof(char)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TBYTE:
        if (!(fcol->val[i] = malloc(nrow * sizeof(unsigned char)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TSHORT:
        if (!(fcol->val[i] = malloc(nrow * sizeof(short)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TINT:
        if (!(fcol->val[i] = malloc(nrow * sizeof(int)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TLONG:
        if (!(fcol->val[i] = malloc(nrow * sizeof(long)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TLONGLONG:
        if (!(fcol->val[i] = malloc(nrow * sizeof(long long)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TFLOAT:
        if (!(fcol->val[i] = malloc(nrow * sizeof(float)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TDOUBLE:
        if (!(fcol->val[i] = malloc(nrow * sizeof(double)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        break;
      case TSTRING:
        if (!(fcol->str[i] = malloc(nrow * sizeof(char *)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        for (long j = 0; j < nrow; j++) fcol->str[i][j] = NULL;
        if (!(fcol->str[i][0] = malloc(nrow * (fcol->width[i] + 1)))) {
          P_ERR("failed to allocate memory for FITS columns\n");
          return FCFC_ERR_MEMORY;
        }
        for (long j = 1; j < nrow; j++)
          fcol->str[i][j] = fcol->str[i][0] + j * (fcol->width[i] + 1);
        break;
      default:
        P_ERR("unsupported FITS column type: %d\n", fcol->dtype[i]);
        return FCFC_ERR_FITS;
    }
  }

  return 0;
}


/*============================================================================*\
                     Functions for processing FITS columns
\*============================================================================*/

/******************************************************************************
Function `exp_insert`:
  Insert a character to a string expression.
Arguments:
  * `str`:      address of the string;
  * `size`:     length of the string;
  * `pos`:      position to be inserted the character;
  * `c`:        the character to be inserted.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int exp_insert(char **str, size_t *size, const size_t pos,
    const char c) {
  /* Increase the length of the string if necessary. */
  if (pos + 1 >= *size) {
    if ((SIZE_MAX - 1) / 2 < *size) {
      P_ERR("unexpected length of the expression\n");
      return FCFC_ERR_UNKNOWN;
    }
    *size <<= 1;

    char *tmp = realloc(*str, *size * sizeof(char));
    if (!tmp) {
      P_ERR("failed to allocate memory for the expression\n");
      return FCFC_ERR_MEMORY;
    }
    *str = tmp;
  }

  (*str)[pos] = c;
  return 0;
}

/******************************************************************************
Function `print_exp_error`:
  Print the error message for the expression.
Arguments:
  * `fmtr`:     the formatter string;
  * `ptr`:      pointer to the position on error.
******************************************************************************/
static inline void print_exp_error(const char *exp, const char *ptr) {
  fprintf(stderr, "%s\n", exp);
  if (ptr > exp) {
    for (ptrdiff_t i = 0; i < ptr - exp; i++) fprintf(stderr, " ");
    fprintf(stderr, "^\n");
  }
}

/******************************************************************************
Function `get_fts_colnum`:
  Get the column index of a certain column in a FITS table.
Arguments:
  * `fp`:       pointer to the opened FITS file;
  * `name`:     name of the column;
  * `num`:      the output column number.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int get_fts_colnum(fitsfile *fp, char *name, int *num) {
  int cnum = 0, status = 0;
  char cname[FLEN_VALUE];
  memset(cname, 0, FLEN_VALUE);

  /* Get the column number and name. */
  if (fits_get_colname(fp, FCFC_FITS_CASESEN, name, cname, &cnum, &status)) {
    P_ERR("cfitsio error: ");
    fits_report_error(stderr, status);
    return FCFC_ERR_FITS;
  }
  /* Compare the input and retrieved names. */
  for (int i = 0; i < FLEN_VALUE; i++) {
    int diff = tolower(name[i]) - tolower(cname[i]);
    if (diff) {
      P_WRN("the FITS column name for %c%c%s%c is: '%s'\n",
          AST_VAR_FLAG, AST_VAR_START, name, AST_VAR_END, cname);
      break;
    }
    if (!name[i]) break;
  }
  *num = cnum;
  return 0;
}

/******************************************************************************
Function `get_fts_cols`:
  Get the information of fits columns required by an expression,
  and construct the corresponding libast-compatible expression.
Arguments:
  * `fp`:       pointer to the opened FITS file;
  * `exp`:      the input expression;
  * `fcol`:     structure for storing FITS column information;
  * `aexp`:     the output libast-compatible expression.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int get_fts_cols(fitsfile *fp, const char *exp, FTS_COL *fcol,
    char **aexp) {
  /* Allocate memory for the output expression. */
  size_t len = strlen(exp);
  char *res = malloc(len * sizeof(char));
  if (!res) {
    P_ERR("failed to allocate memory for parsing expressions\n");
    return FCFC_ERR_MEMORY;
  }

  size_t k = 0;
  bool col_name = false;
  for (const char *c = exp; *c; c++) {
    /* Copy the character directly if it is not a column name. */
    if (!col_name) {
      if (*c != AST_VAR_FLAG) {
        if (exp_insert(&res, &len, k++, *c)) {
          free(res); return FCFC_ERR_MEMORY;
        }
      }
      else if (*(++c) != AST_VAR_START) {
        P_ERR("name of a FITS column must be enclosed by %c%c:\n",
            AST_VAR_START, AST_VAR_END);
        print_exp_error(exp, c);
        free(res); return FCFC_ERR_CFG;
      }
      else {
        if (exp_insert(&res, &len, k++, AST_VAR_FLAG) ||
            exp_insert(&res, &len, k++, AST_VAR_START)) {
          free(res); return FCFC_ERR_MEMORY;
        }
        col_name = true;
      }
    }
    else {              /* col_name == true: parse column identifier */
      if (*c == AST_VAR_END) {
        P_ERR("empty column identifier:\n");
        print_exp_error(exp, c);
        free(res); return FCFC_ERR_CFG;
      }
      const char *p = c;
      size_t kk = k;
      while (*c && *c != AST_VAR_END) {
        /* Copy the column name to the output expression temporarily. */
        if (exp_insert(&res, &len, kk++, *c)) {
          free(res); return FCFC_ERR_MEMORY;
        }
        c++;
      }
      if (*c == '\0') {
        P_ERR("unbalanced column identifier:\n");
        print_exp_error(exp, p);
        free(res); return FCFC_ERR_CFG;
      }

      /* Look for the column number in the FITS table. */
      if (exp_insert(&res, &len, kk, '\0')) {
        free(res); return FCFC_ERR_MEMORY;
      }
      int colnum = 0;
      if (get_fts_colnum(fp, res + k, &colnum)) {
        free(res); return FCFC_ERR_FILE;
      }

      /* Check if the column number exists. */
      int idx = 0;
      for (; idx < fcol->ncol; idx++) {
        if (colnum == fcol->colnum[idx]) break;
      }
      if (idx == fcol->ncol) {          /* record this column */
        if (fts_col_insert(fcol, colnum)) {
          free(res); return FCFC_ERR_MEMORY;
        }
      }

      /* Set `idx + 1` as the new column identifier. */
      int cnt = snprintf(res + k, len - k, "%d", idx + 1);
      if (cnt < 0) {
        P_ERR("failed to construct the AST expression\n");
        return FCFC_ERR_MEMORY;
      }
      if (k + cnt + 2 >= len) {
        /* Increase the length of the new expression. */
        len = k + cnt + 2;
        char *tmp = realloc(res, len * sizeof(char));
        if (!tmp) {
          P_ERR("failed to allocate memory for the AST expression\n");
          free(res); return FCFC_ERR_MEMORY;
        }
        res = tmp;
        if (snprintf(res + k, len - k, "%d", idx + 1) < 0) {
          P_ERR("failed to construct the AST expression\n");
          return FCFC_ERR_MEMORY;
        }
      }
      k += cnt;
      res[k++] = AST_VAR_END;
      col_name = false;
    }
  }

  res[k++] = '\0';
  /* Reduce the memory cost if applicable. */
  if (k < len) {
    char *tmp = realloc(res, k * sizeof(char));
    if (tmp) res = tmp;
  }

  *aexp = res;
  return 0;
}


/*============================================================================*\
                   Functions for reading and processing data
\*============================================================================*/

/******************************************************************************
Function `load_fits_columns`:
  Read FITS columns.
Arguments:
  * `fp`:       pointer to the opened FITS file;
  * `fcol`:     structure for storing FITS column information;
  * `start`:    first row to be read;
  * `nrow`:     number of rows to be read.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int load_fits_columns(fitsfile *fp, FTS_COL *fcol, const long start,
    const long nrow) {
  for (int i = 0; i < fcol->ncol; i++) {
    int anynul = 0, status = 0;
    switch (fcol->dtype[i]) {
      case TLOGICAL:
        fits_read_col_log(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TBYTE:
        fits_read_col_byt(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TSHORT:
        fits_read_col_sht(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TINT:
        fits_read_col_int(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TLONG:
        fits_read_col_lng(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TLONGLONG:
        fits_read_col_lnglng(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TFLOAT:
        fits_read_col_flt(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TDOUBLE:
        fits_read_col_dbl(fp, fcol->colnum[i], start, 1, nrow, 0,
            fcol->val[i], &anynul, &status);
        break;
      case TSTRING:
        fits_read_col_str(fp, fcol->colnum[i], start, 1, nrow, NULL,
            fcol->str[i], &anynul, &status);
        break;
      default:
        P_ERR("unexpected FITS column type: %d\n", fcol->dtype[i]);
        return FCFC_ERR_FITS;
    }
    if (status) {
      P_ERR("cfitsio error: "); fits_report_error(stderr, status);
      return FCFC_ERR_FITS;
    }
    if (anynul) {
      P_ERR("undefined value found in column %d\n", fcol->colnum[i]);
      return FCFC_ERR_FITS;
    }
  }
  return 0;
}

/******************************************************************************
Function `fits_apply_sel`:
  Apply the selection criteria based on the AST and FITS columns.
Arguments:
  * `ast`:      abstract syntax tree for the selection criteria;
  * `fcol`:     structure for storing FITS column information;
  * `idx`:      index of the row to be processed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int fits_apply_sel(ast_t *ast, const FTS_COL *fcol, const long idx,
    bool *res) {
  for (long i = 0; i < ast->nvar; i++) {
    char *vc;
    bool vb;
    int vi;
    long vl;
    float vf;
    double vd;
    const long j = ast->vidx[i];
    const long k = j - 1;
    switch (fcol->dtype[k]) {
      case TFLOAT:
        vf = ((float *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vf, 0, AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TDOUBLE:
        vd = ((double *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vd, 0, AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TBYTE:
        vi = (int) ((unsigned char *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vi, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TSHORT:
        vi = (int) ((short *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vi, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TINT:
        vi = ((int *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vi, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TLONG:
        vl = ((long *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vl, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TLONGLONG:
        vl = (long) ((long long *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vl, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TLOGICAL:
        vc = (char *) fcol->val[k];
        if (vc[idx] == 1) vb = true;
        else if (vc[idx] == 0) vb = false;
        else {
          P_ERR("invalid boolean value in FITS column %d: %d\n",
              fcol->colnum[k], vc[idx]);
          return FCFC_ERR_FITS;
        }
        if (ast_set_var(ast, j, &vb, 0, AST_DTYPE_BOOL)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TSTRING:
        if (ast_set_var(ast, j, fcol->str[k][idx], fcol->width[k],
            AST_DTYPE_STRING)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      default:
        P_ERR("unexpected FITS column type: %d\n", fcol->dtype[i]);
        return FCFC_ERR_FITS;
    }
  }

  if (ast_eval(ast, res)) {
    P_AST_ERR(ast);
    return FCFC_ERR_AST;
  }
  return 0;
}

/******************************************************************************
Function `fits_read_real`:
  Read a floating-point number based on the AST and FITS columns.
Arguments:
  * `ast`:      abstract syntax tree for the number evaluation;
  * `fcol`:     structure for storing FITS column information;
  * `idx`:      index of the row to be processed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int fits_read_real(ast_t *ast, const FTS_COL *fcol, const long idx,
    real *res) {
  for (long i = 0; i < ast->nvar; i++) {
    int vi;
    long vl;
    float vf;
    double vd;
    const long j = ast->vidx[i];
    const long k = j - 1;
    switch (fcol->dtype[k]) {
      case TFLOAT:
        vf = ((float *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vf, 0, AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TDOUBLE:
        vd = ((double *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vd, 0, AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TBYTE:
        vi = (int) ((unsigned char *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vi, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TSHORT:
        vi = (int) ((short *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vi, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TINT:
        vi = ((int *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vi, 0, AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TLONG:
        vl = ((long *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vl, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case TLONGLONG:
        vl = (long) ((long long *) fcol->val[k])[idx];
        if (ast_set_var(ast, j, &vl, 0, AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      default:
        P_ERR("unexpected FITS column for numerical evaluation: %d\n",
            fcol->colnum[k]);
        return FCFC_ERR_FITS;
    }
  }

  if (ast_eval(ast, res)) {
    P_AST_ERR(ast);
    return FCFC_ERR_AST;
  }
  return 0;
}


/*============================================================================*\
                       Interface for reading a FITS file
\*============================================================================*/

/******************************************************************************
Function `read_fits_data`:
  Read a FITS file for the positions and weights.
Arguments:
  * `fname`:    filename of the input catalog;
  * `rcol_ids`: identifiers of columns for floating-point numbers;
  * `nrcol`:    number of columns for floating-point numbers;
  * `sel`:      data selection criteria;
  * `rout`:     address for storing the output floating-point columns;
  * `num`:      number of objects read in total;
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_fits_data(const char *fname, char *const *rcol_ids, const int nrcol,
    const char *sel, real ***rout, size_t *num, const int verb) {
  /* Validate arguments. */
  if (!fname) {
    P_ERR("the filename for FITS table is not specified\n");
    return FCFC_ERR_ARG;
  }
  if (!rcol_ids || !rout || !num) {
    P_ERR("FITS data reader not initialized\n");
    return FCFC_ERR_ARG;
  }
  if (nrcol < 3 || nrcol > 4) {
    P_ERR("unexpected number of floating-point columns to be read: %d\n",
        nrcol);
    return FCFC_ERR_ARG;
  }

  /* Open the file for reading. */
  int status = 0;
  fitsfile *fp = NULL;
  if (fits_open_data(&fp, fname, READONLY, &status)) FITS_ABORT;
  if (verb) printf("  Filename: %s\n", fname);

  /* Get the number of objects, and the optimal number of rows read at once. */
  long ntot = 0, nstep = 0;
  if (fits_get_num_rows(fp, &ntot, &status)) FITS_ABORT;
  if (fits_get_rowsize(fp, &nstep, &status)) FITS_ABORT;

  real **res = NULL;            /* arrays for the outputs */
  char **exp = NULL;            /* AST expressions        */
  ast_t **ast_real = NULL;      /* ASTs for the outputs   */
  ast_t *ast_sel = NULL;        /* AST for the selection  */
  FTS_COL *fcol = fts_col_init();
  if (!fcol) {
    CLEAN_FITS_PTR; return FCFC_ERR_MEMORY;
  }

  /* Read FITS columns and construct libast expressions. */
  if (!(exp = malloc((nrcol + 1) * sizeof(char *)))) {
    P_ERR("failed to allocate memory for the AST expressions\n");
    CLEAN_FITS_PTR; return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i <= nrcol; i++) exp[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
    if (get_fts_cols(fp, rcol_ids[i], fcol, exp + i)) {
      CLEAN_FITS_PTR; return FCFC_ERR_FITS;
    }
  }
  if (sel && *sel && get_fts_cols(fp, sel, fcol, exp + nrcol)) {
    CLEAN_FITS_PTR; return FCFC_ERR_FITS;
  }

  /* Construct ASTs. */
  if (!(ast_real = malloc(sizeof(ast_t *) * nrcol))) {
    P_ERR("failed to allocate memory for the Abstract Syntax Trees\n");
    CLEAN_FITS_PTR; return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < nrcol; i++) ast_real[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
    if (!(ast_real[i] = ast_init())) {
      P_AST_ERR(ast_real[i]); CLEAN_FITS_PTR; return FCFC_ERR_AST;
    }
    if (ast_build(ast_real[i], exp[i], FCFC_AST_REAL, true)) {
      P_AST_ERR(ast_real[i]); CLEAN_FITS_PTR; return FCFC_ERR_AST;
    }
    if (ast_real[i]->nvar == 0)
      P_WRN("the column expression indicates a constant: `%s'\n", rcol_ids[i]);
  }
  if (exp[nrcol]) {
    if (!(ast_sel = ast_init())) {
      P_AST_ERR(ast_sel); CLEAN_FITS_PTR; return FCFC_ERR_AST;
    }
    if (ast_build(ast_sel, exp[nrcol], AST_DTYPE_BOOL, true)) {
      P_AST_ERR(ast_sel); CLEAN_FITS_PTR; return FCFC_ERR_AST;
    }
    if (ast_sel->nvar == 0)
      P_WRN("the expression for data selection is a constant: `%s'\n", sel);
  }

  /* Check column information and allocate memory for file reading. */
  if (fts_col_finalize(fp, fcol, nstep)) {
    CLEAN_FITS_PTR; return FCFC_ERR_FITS;
  }

  /* Allocate memory for the outputs. */
  if (!(res = malloc(sizeof(real *) * nrcol))) {
    P_ERR("failed to allocate memory for the data\n");
    CLEAN_FITS_PTR; return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < nrcol; i++) res[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
#if     FCFC_SIMD  ==  FCFC_SIMD_NONE
    if (!(res[i] = malloc(ntot * sizeof(real))))
#else
    if (!(res[i] = malloc((ntot + FCFC_NUM_REAL) * sizeof(real))))
#endif
    {
      P_ERR("failed to allocate memory for the data\n");
      CLEAN_FITS_PTR; return FCFC_ERR_MEMORY;
    }
  }

  /* Dynamic allocations for OpenMP threads. */
#ifdef OMP
  const int nomp = omp_get_max_threads();
  /* Construct the ASTs for non-master threads. */
  ast_t **ast_preal = NULL, **ast_psel = NULL;
  if (nomp > 1) {
    if (!(ast_preal = malloc(sizeof(ast_t *) * (nomp - 1) * nrcol)) ||
        !(ast_psel = malloc(sizeof(ast_t *) * (nomp - 1)))) {
      P_ERR("failed to allocate memory for thread-private ASTs\n");
      FCFC_FITS_QUIT(FCFC_ERR_MEMORY);
    }
  }
  #pragma omp parallel for num_threads(nomp)
  for (int j = 0; j < nomp - 1; j++) {
    /* Positions and weights. */
    for (int i = 0; i < nrcol; i++) {
      int k = j * nrcol + i;
      if (!(ast_preal[k] = ast_init())) {
        P_AST_ERR(ast_preal[k]); FCFC_FITS_QUIT(FCFC_ERR_AST);
      }
      if (ast_build(ast_preal[k], exp[i], FCFC_AST_REAL, true)) {
        P_AST_ERR(ast_preal[k]); FCFC_FITS_QUIT(FCFC_ERR_AST);
      }
    }
    /* Selection. */
    if (exp[nrcol]) {
      if (!(ast_psel[j] = ast_init())) {
        P_AST_ERR(ast_psel[j]); FCFC_FITS_QUIT(FCFC_ERR_AST);
      }
      if (ast_build(ast_psel[j], exp[nrcol], AST_DTYPE_BOOL, true)) {
        P_AST_ERR(ast_psel[j]); FCFC_FITS_QUIT(FCFC_ERR_AST);
      }
    }
    else ast_psel[j] = NULL;
  }
  /* Construct the private data pool. */
  real *pdata = malloc(sizeof(real) * nomp * nrcol * FCFC_DATA_THREAD_NUM);
  size_t *pndata = calloc(nomp, sizeof(size_t));
  if (!pdata || !pndata) {
    P_ERR("failed to allocate memory for the thread-private data\n");
    FCFC_FITS_QUIT(FCFC_ERR_MEMORY);
  }
#endif

  /* Read the file by chunks. */
  long nread = 1;
  long nrest = ntot;
  size_t n = 0;
  while (nrest) {
    long nrow = (nstep < nrest) ? nstep : nrest;
    /* Read the necessary columns. */
    if (load_fits_columns(fp, fcol, nread, nrow)) {
      CLEAN_FITS_PTR; CLEAN_FITS_PTR_OMP; return FCFC_ERR_FITS;
    }

    /* Process the data row by row. */
#ifdef OMP
  #pragma omp parallel num_threads(nomp) firstprivate(ast_real,ast_sel)
    {
      /* Redirect pointers to the private pools. */
      const int tid = omp_get_thread_num();
      if (nomp > 1 && tid > 0) {
        ast_real = ast_preal + (tid - 1) * nrcol;
        ast_sel = ast_psel[tid - 1];
      }
      real *pres = pdata + tid * nrcol * FCFC_DATA_THREAD_NUM;
      size_t *pnum = pndata + tid;
      /* Process rows in parallel. */
  #pragma omp for
      for (long i = 0; i < nrow; i++) {
        /* Apply selection. */
        if (ast_sel) {
          bool keep = false;
          if (fits_apply_sel(ast_sel, fcol, i, &keep)) {
            FCFC_FITS_QUIT(FCFC_ERR_AST);
          }
          if (!keep) continue;
        }

        /* Record coordinates and weights in the private data pool. */
        for (int j = 0; j < nrcol; j++) {
          if (fits_read_real(ast_real[j], fcol, i,
              pres + j * FCFC_DATA_THREAD_NUM + *pnum)) {
            FCFC_FITS_QUIT(FCFC_ERR_AST);
          }
        }

        /* Record the private data and clear the pool if necessary. */
        if (++(*pnum) >= FCFC_DATA_THREAD_NUM) {
  #pragma omp critical
          {
            for (int k = 0; k < nrcol; k++) {
              memcpy(res[k] + n, pres + k * FCFC_DATA_THREAD_NUM,
                  sizeof(real) * FCFC_DATA_THREAD_NUM);
            }
            n += FCFC_DATA_THREAD_NUM;
          }
          *pnum = 0;
        }
      }
    }
#else
    for (long i = 0; i < nrow; i++) {
      /* Apply selection. */
      if (ast_sel) {
        bool keep = false;
        if (fits_apply_sel(ast_sel, fcol, i, &keep)) {
          CLEAN_FITS_PTR; return FCFC_ERR_FITS;
        }
        if (!keep) continue;
      }

      /* Read coordinates and weights. */
      for (int j = 0; j < nrcol; j++) {
        if (fits_read_real(ast_real[j], fcol, i, res[j] + n)) {
          CLEAN_FITS_PTR; return FCFC_ERR_FITS;
        }
      }
      n++;
    }
#endif

    nread += nrow;
    nrest -= nrow;
  }

#ifdef OMP
  /* Record the rest of the private data. */
  for (int i = 0; i < nomp; i++) {
    if (pndata[i]) {
      for (int k = 0; k < nrcol; k++) {
        memcpy(res[k] + n, pdata + (k + i * nrcol) * FCFC_DATA_THREAD_NUM,
            sizeof(real) * pndata[i]);
      }
      n += pndata[i];
    }
  }
  /* Release memory for thread-private data structures. */
  CLEAN_FITS_PTR_OMP;
#endif

  if (!n) {
    P_ERR("no object is read from file: `%s'\n", fname);
    CLEAN_FITS_PTR; return FCFC_ERR_FILE;
  }

  if (verb) {
    printf("  Number of rows processed in total: %zu\n"
        "  Number of recorded objects: %zu\n", ntot, n);
  }

  if (fits_close_file(fp, &status))
    P_WRN("failed to close file: `%s'\n", fname);
  for (int i = 0; i < nrcol; i++) ast_destroy(ast_real[i]);
  free(ast_real);
  ast_destroy(ast_sel);
  fts_col_destroy(fcol);
  for (int i = 0; i <= nrcol; i++) if (exp[i]) free(exp[i]);
  free(exp);

  /* Release pre-allocated memory that is not necessary any more. */
#if     FCFC_SIMD  ==  FCFC_SIMD_NONE
  if ((size_t) ntot != n) {
    for (int i = 0; i < nrcol; i++) {
      real *tmp = realloc(res[i], n * sizeof(real));
      if (tmp) res[i] = tmp;
    }
  }
#else
  if ((size_t) ntot != n + FCFC_NUM_REAL) {
    for (int i = 0; i < nrcol; i++) {
      real *tmp = realloc(res[i], (n + FCFC_NUM_REAL) * sizeof(real));
      if (tmp) res[i] = tmp;
      memset(res[i] + n, 0, sizeof(real) * FCFC_NUM_REAL);
    }
  }
#endif

  *num = n;
  *rout = res;
  return 0;
}

#endif

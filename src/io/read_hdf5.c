/*******************************************************************************
* read_hdf5.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#ifdef WITH_HDF5

#include "read_file.h"
#include "libast.h"
#include <hdf5.h>
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
#define CLEAN_HDF5_PTR                                                  \
  hdf5_col_destroy(hcol);                                               \
  H5Fclose(file_id);                                                    \
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
  #define CLEAN_HDF5_PTR_OMP                                            \
    if (nomp > 1) {                                                     \
      for (int _i = 0; _i < nomp - 1; _i++)                             \
        ast_destroy(ast_psel[_i]);                                      \
      for (int _i = 0; _i < (nomp - 1) * nrcol; _i++)                   \
        ast_destroy(ast_preal[_i]);                                     \
      free(ast_preal); free(ast_psel);                                  \
    }                                                                   \
    free(pdata); free(pndata);

  #define FCFC_HDF5_QUIT(x) {                                           \
    printf(FMT_FAIL);                                                   \
    P_EXT("failed to read the HDF5 file\n");                            \
    exit(x);                                                            \
  }
#else
  #define CLEAN_HDF5_PTR_OMP
#endif

/*============================================================================*\
                       Data structure for HDF5 "columns"
\*============================================================================*/

/* Data type identifiers. */
typedef enum {
  FCFC_HDF5_DTYPE_INT,
  FCFC_HDF5_DTYPE_LONG,
  FCFC_HDF5_DTYPE_FLOAT,
  FCFC_HDF5_DTYPE_DOUBLE,
  FCFC_HDF5_DTYPE_FIX_STR,
  FCFC_HDF5_DTYPE_VAR_STR
} FCFC_h5dtype_e;

/* Structure for storing the HDF5 column information. */
typedef struct {
  int ncol;                     /* number of columns to be processed    */
  char **dname;                 /* names of the relevant datasets       */
  int *idx;                     /* index of multi-dimensional data      */
  int *dim;                     /* major dimension of the data          */
  hid_t *dset_ids;              /* identifiers of the datasets          */
  hid_t *space_ids;             /* identifiers of the dataspaces        */
  FCFC_h5dtype_e *dtypes;       /* datatype classes of the columns      */
  size_t *dsizes;               /* sizes of datatypes in bytes          */
  hid_t *memtypes;              /* datatypes for storing the column     */
  hid_t memspace;               /* dataspaces for storing a column      */
  void **data;                  /* addresses for fixed-length data      */
  char ***vstr;                 /* addresses for variable-length string */
} HDF5_COL;


/*============================================================================*\
              Functions for manipulating the HDF5 column structure
\*============================================================================*/

/******************************************************************************
Function `hdf5_col_init`:
  Initialize the structure for storing HDF5 column information.
Arguments:
  * `nchunk`:   number of objects to be read at once.
Return:
  Address of the structure on success; NULL on error.
******************************************************************************/
static inline HDF5_COL *hdf5_col_init(const hsize_t nchunk) {
  HDF5_COL *hcol = calloc(1, sizeof(HDF5_COL));
  if (!hcol) {
    P_ERR("failed to allocate memory for storing HDF5 information\n");
    return NULL;
  }

  hcol->memspace = H5Screate_simple(1, &nchunk, NULL);
  if (hcol->memspace < 0) {
    P_ERR("failed to create dataspace for storing HDF5 data\n");
    free(hcol); return NULL;
  }

  hcol->dname = NULL;
  hcol->idx = hcol->dim = NULL;
  hcol->dset_ids = hcol->space_ids = hcol->memtypes = NULL;
  hcol->dtypes = NULL;
  hcol->dsizes = NULL;
  hcol->data = NULL;
  hcol->vstr = NULL;
  return hcol;
}

/******************************************************************************
Function `hdf5_col_insert`:
  Insert a column to the HDF5 column structure.
Arguments:
  * `hcol`:     structure for storing HDF5 column information;
  * `dname`:    name of the dataset;
  * `idx`:      index of multi-dimensional data.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int hdf5_col_insert(HDF5_COL *hcol, char *dname, const int idx) {
  /* Check if there is enough allocated space. */
  if ((hcol->ncol & (hcol->ncol - 1)) == 0) {   /* ncol is 0 or power of 2 */
    if (INT_MAX / 2 < hcol->ncol) {
      P_ERR("too many columns required by the expressions\n");
      return FCFC_ERR_CFG;
    }
    int size = (hcol->ncol) ? hcol->ncol << 1 : 1;   /* double the size */
    char **ntmp = realloc(hcol->dname, size * sizeof(char *));
    if (!ntmp) {
      P_ERR("failed to allocate memory for recording HDF5 information\n");
      return FCFC_ERR_MEMORY;
    }
    hcol->dname = ntmp;
    for (int i = hcol->ncol; i < size; i++) hcol->dname[i] = NULL;
    int *itmp = realloc(hcol->idx, size * sizeof(int));
    if (!itmp) {
      P_ERR("failed to allocate memory for recording HDF5 information\n");
      return FCFC_ERR_MEMORY;
    }
    hcol->idx = itmp;
  }

  hcol->dname[hcol->ncol] = dname;
  hcol->idx[hcol->ncol++] = idx;
  return 0;
}

/******************************************************************************
Function `hdf5_col_destroy`:
  Deconstruct the structure for storing HDF5 column information.
Arguments:
  * `hcol`:     structure for storing HDF5 column information.
******************************************************************************/
static void hdf5_col_destroy(HDF5_COL *hcol) {
  if (!hcol) return;
  if (hcol->dset_ids) {
    /* Close HDF5 datasets if applicable. */
    for (int i = 0; i < hcol->ncol; i++) {
      if (hcol->dset_ids[i] >= 0) {
        /* Check if the same dataset is closed earlier. */
        int j = 0;
        for (; j < i; j++) {
          if (!strcmp(hcol->dname[i], hcol->dname[j])) break;
        }
        if (j != i) continue;
        if (H5Dclose(hcol->dset_ids[i]) < 0)
          P_WRN("failed to close dataset: '%s'\n", hcol->dname[i]);
      }
    }
    free(hcol->dset_ids);
  }
  if (hcol->dname) {
    for (int i = 0; i < hcol->ncol; i++)
      if (hcol->dname[i]) free(hcol->dname[i]);
    free(hcol->dname);
  }
  if (hcol->idx) free(hcol->idx);
  if (hcol->dim) free(hcol->dim);
  if (hcol->space_ids) free(hcol->space_ids);
  if (hcol->dtypes) free(hcol->dtypes);
  if (hcol->dsizes) free(hcol->dsizes);
  if (hcol->memtypes) free(hcol->memtypes);
  if (hcol->data) {
    for (int i = 0; i < hcol->ncol; i++)
      if (hcol->data[i]) free(hcol->data[i]);
    free(hcol->data);
  }
  if (hcol->vstr) {
    for (int i = 0; i < hcol->ncol; i++)
      if (hcol->vstr[i]) free(hcol->vstr[i]);
    free(hcol->vstr);
  }
  free(hcol);
}

/******************************************************************************
Function `hdf5_col_finalize`:
  Complete the structure for storing HDF5 column information.
Arguments:
  * `file_id`:  identifier of the opened HDF5 file;
  * `hcol`:     structure for storing HDF5 column information;
  * `nchunk`:   number of objects to be read at once;
  * `ndata`:    number of data points in the HDF5 file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int hdf5_col_finalize(hid_t file_id, HDF5_COL *hcol,
    const hsize_t nchunk, hsize_t *ndata) {
  if (!hcol->ncol) {
    P_ERR("no data to be read from the HDF5 file\n");
    return FCFC_ERR_CFG;
  }
  /* Reduce memory cost if applicable. */
  if ((hcol->ncol & (hcol->ncol - 1)) != 0) {
    char **ntmp = realloc(hcol->dname, hcol->ncol * sizeof(char *));
    if (ntmp) hcol->dname = ntmp;
    int *itmp = realloc(hcol->idx, hcol->ncol * sizeof(int));
    if (itmp) hcol->idx = itmp;
  }

  /* Allocate memory for the rest column information. */
  if (!(hcol->dim = calloc(hcol->ncol, sizeof(int))) ||
      !(hcol->dset_ids = malloc(hcol->ncol * sizeof(hid_t)))) {
    P_ERR("failed to allocate memory for recording HDF5 file information\n");
    return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < hcol->ncol; i++) hcol->dset_ids[i] = -1;
  if (!(hcol->space_ids = malloc(hcol->ncol * sizeof(hid_t))) ||
      !(hcol->dtypes = calloc(hcol->ncol, sizeof(FCFC_h5dtype_e))) ||
      !(hcol->dsizes = calloc(hcol->ncol, sizeof(size_t))) ||
      !(hcol->memtypes = malloc(hcol->ncol * sizeof(hid_t))) ||
      !(hcol->data = malloc(hcol->ncol * sizeof(void *)))) {
    P_ERR("failed to allocate memory for recording HDF5 file information\n");
    return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < hcol->ncol; i++) hcol->data[i] = NULL;
  if (!(hcol->vstr = malloc(hcol->ncol * sizeof(char **)))) {
    P_ERR("failed to allocate memory for recording HDF5 file information\n");
    return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < hcol->ncol; i++) hcol->vstr[i] = NULL;

  /* Read column information. */
  hsize_t num = 0;              /* number of data points in the file */
  for (int i = 0; i < hcol->ncol; i++) {
    /* Check if the dataset is already opened. */
    bool open_dset = true;
    for (int j = 0; j < i; j++) {
      if (!strcmp(hcol->dname[i], hcol->dname[j])) {
        hcol->dset_ids[i] = hcol->dset_ids[j];
        hcol->space_ids[i] = hcol->space_ids[j];
        open_dset = false;
        break;
      }
    }
    /* Open the dataset, and get the dataspace. */
    if (open_dset) {
      hcol->dset_ids[i] = H5Dopen(file_id, hcol->dname[i], H5P_DEFAULT);
      if (hcol->dset_ids[i] < 0) {
        P_ERR("failed to open the HDF5 dataset: '%s'\n", hcol->dname[i]);
        return FCFC_ERR_HDF5;
      }
      hcol->space_ids[i] = H5Dget_space(hcol->dset_ids[i]);
      if (hcol->space_ids[i] < 0) {
        P_ERR("failed to get the dataspace of dataset: '%s'\n", hcol->dname[i]);
        return FCFC_ERR_HDF5;
      }
    }

    /* Read the dimensions of the dataset. */
    int ndim = H5Sget_simple_extent_ndims(hcol->space_ids[i]);
    if (ndim < 1 || ndim > 2) {
      P_ERR("unsupported number of HDF5 dataset dimension: %d\n", ndim);
      return FCFC_ERR_HDF5;
    }
    if (ndim == 1 && hcol->idx[i] != 0) {
      P_WRN("omitting the index %d of dataset '%s' as there is only"
          "one dimension\n", hcol->idx[i] + 1, hcol->dname[i]);
      hcol->idx[i] = 0;
    }
    /* Choose the longer dimension as the column direction. */
    hsize_t dims[2] = {1, 1};
    if (H5Sget_simple_extent_dims(hcol->space_ids[i], dims, NULL) < 0) {
      P_ERR("failed to get the dimensions of dataset: '%s'\n",
          hcol->dname[i]);
      return FCFC_ERR_HDF5;
    }
    hsize_t nrow, ncol;
    if (dims[0] >= dims[1]) {
      nrow = dims[0];
      ncol = dims[1];
      hcol->dim[i] = 0;
    }
    else {
      nrow = dims[1];
      ncol = dims[0];
      hcol->dim[i] = 1;
    }
    if (num == 0) num = nrow;
    else if (num != nrow) {
      P_ERR("ununique numbers of objects detected in different datasets: "
          "%zu vs %zu\n", (size_t) num, (size_t) nrow);
      return FCFC_ERR_HDF5;
    }
    /* Check the index of dataset. */
    if ((hsize_t) hcol->idx[i] >= ncol) {
      P_ERR("index (%d) exceeds the dimension (%zu) of dataset: '%s'\n",
          hcol->idx[i] + 1, (size_t) ncol, hcol->dname[i]);
      return FCFC_ERR_CFG;
    }

    /* Get the datatype and allocate memory for storing the data. */
    hid_t dtype_id = H5Dget_type(hcol->dset_ids[i]);
    if (dtype_id < 0) {
      P_ERR("failed to get the datatype of dataset: '%s'\n", hcol->dname[i]);
      return FCFC_ERR_HDF5;
    }
    H5T_class_t dtype = H5Tget_class(dtype_id);
    htri_t is_vstr;
    switch (dtype) {
      case H5T_INTEGER:
        if (!(hcol->dsizes[i] = H5Tget_size(dtype_id))) {
          P_ERR("failed to get the datatype size of dataset: '%s'\n",
              hcol->dname[i]);
          return FCFC_ERR_HDF5;
        }
        if (hcol->dsizes[i] > sizeof(int)) {
          hcol->dtypes[i] = FCFC_HDF5_DTYPE_LONG;
          hcol->memtypes[i] = H5T_NATIVE_LONG;
          if (!(hcol->data[i] = malloc(sizeof(long) * nchunk))) {
            P_ERR("failed to allocate memory for storing HDF5 data\n");
            return FCFC_ERR_MEMORY;
          }
        }
        else {
          hcol->dtypes[i] = FCFC_HDF5_DTYPE_INT;
          hcol->memtypes[i] = H5T_NATIVE_INT;
          if (!(hcol->data[i] = malloc(sizeof(int) * nchunk))) {
            P_ERR("failed to allocate memory for storing HDF5 data\n");
            return FCFC_ERR_MEMORY;
          }
        }
        break;
      case H5T_FLOAT:
        if (!(hcol->dsizes[i] = H5Tget_size(dtype_id))) {
          P_ERR("failed to get the datatype size of dataset: '%s'\n",
              hcol->dname[i]);
          return FCFC_ERR_HDF5;
        }
        if (hcol->dsizes[i] > sizeof(float)) {
          hcol->dtypes[i] = FCFC_HDF5_DTYPE_DOUBLE;
          hcol->memtypes[i] = H5T_NATIVE_DOUBLE;
          if (!(hcol->data[i] = malloc(sizeof(double) * nchunk))) {
            P_ERR("failed to allocate memory for storing HDF5 data\n");
            return FCFC_ERR_MEMORY;
          }
        }
        else {
          hcol->dtypes[i] = FCFC_HDF5_DTYPE_FLOAT;
          hcol->memtypes[i] = H5T_NATIVE_FLOAT;
          if (!(hcol->data[i] = malloc(sizeof(float) * nchunk))) {
            P_ERR("failed to allocate memory for storing HDF5 data\n");
            return FCFC_ERR_MEMORY;
          }
        }
        break;
      case H5T_STRING:
        /* Check if the length of string is fixed. */
        if ((is_vstr = H5Tis_variable_str(dtype_id)) < 0) {
          P_ERR("failed to check if the length of string is variable in "
              "dataset: '%s'\n", hcol->dname[i]);
          return FCFC_ERR_HDF5;
        }
        if ((hcol->memtypes[i] = H5Tcopy(H5T_C_S1)) < 0) {
          P_ERR("failed to set the datatype of string\n");
          return FCFC_ERR_HDF5;
        }
        if (is_vstr) {          /* the string length is variable */
          hcol->dtypes[i] = FCFC_HDF5_DTYPE_VAR_STR;
          if (H5Tset_size(hcol->memtypes[i], H5T_VARIABLE) < 0) {
            P_ERR("failed to set the datatype of variable-length string\n");
            return FCFC_ERR_HDF5;
          }
          if (!(hcol->vstr[i] = malloc(sizeof(char *) * nchunk))) {
            P_ERR("failed to allocate memory for storing HDF5 data\n");
            return FCFC_ERR_MEMORY;
          }
        }
        else {                  /* the string length is fixed */
          hcol->dtypes[i] = FCFC_HDF5_DTYPE_FIX_STR;
          /* Allocate memory also for the teriminator '\0'. */
          if (H5Tset_size(hcol->memtypes[i], hcol->dsizes[i] + 1) < 0) {
            P_ERR("failed to set the datatype for fixed-length string\n");
            return FCFC_ERR_HDF5;
          }
          if (!(hcol->data[i] = malloc((hcol->dsizes[i] + 1) * nchunk))) {
            P_ERR("failed to allocate memory for storing HDF5 data\n");
            return FCFC_ERR_MEMORY;
          }
        }
        break;
      default:
        P_ERR("unsupported datatype (%d) of dataset: '%s'\n", (int) dtype,
            hcol->dname[i]);
        return FCFC_ERR_HDF5;
    }
  }

  *ndata =  num;
  return 0;
}


/*============================================================================*\
                       Functions for processing HDF5 data
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
Function `get_hdf5_cols`:
  Get the information of HDF5 columns required by an expression,
  and construct the corresponding libast-compatible expression.
Arguments:
  * `exp`:      the input expression;
  * `hcol`:     structure for storing HDF5 column information;
  * `aexp`:     the output libast-compatible expression.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int get_hdf5_cols(const char *exp, HDF5_COL *hcol, char **aexp) {
  /* Allocate memory for the output expression. */
  size_t len = strlen(exp);
  char *res = malloc(len * sizeof(char));
  if (!res) {
    P_ERR("failed to allocate memory for parsing expressions\n");
    return FCFC_ERR_MEMORY;
  }

  size_t k = 0;
  bool dset_name = false;
  for (const char *c = exp; *c; c++) {
    /* Copy the character directly if it is not a dataset name. */
    if (!dset_name) {
      if (*c != AST_VAR_FLAG) {
        if (exp_insert(&res, &len, k++, *c)) {
          free(res); return FCFC_ERR_MEMORY;
        }
      }
      else if (*(++c) != AST_VAR_START) {
        P_ERR("name of a HDF5 dataset must be enclosed by %c%c:\n",
            AST_VAR_START, AST_VAR_END);
        print_exp_error(exp, c);
        free(res); return FCFC_ERR_CFG;
      }
      else {
        if (exp_insert(&res, &len, k++, AST_VAR_FLAG) ||
            exp_insert(&res, &len, k++, AST_VAR_START)) {
          free(res); return FCFC_ERR_MEMORY;
        }
        dset_name = true;
      }
    }
    else {              /* dset_name == true: parse dataset name */
      if (*c == AST_VAR_END || *c == FCFC_COL_IDX_START) {
        P_ERR("empty dataset name identifier:\n");
        print_exp_error(exp, c);
        free(res); return FCFC_ERR_CFG;
      }
      const char *nstart = c;   /* beginning of the dataset name */
      /* Find the end of the dataset name. */
      while (*c && *c != AST_VAR_END && *c != FCFC_COL_IDX_START) c++;
      if (*c == '\0') {
        P_ERR("unbalanced dataset name identifier:\n");
        print_exp_error(exp, nstart);
        free(res); return FCFC_ERR_CFG;
      }

      const char *nend = c;     /* end of the dataset name */
      int cidx = 0;
      if (*c == FCFC_COL_IDX_START) {  /* parse the dataset index */
        /* Find the end of the dataset index. */
        do c++; while (isdigit(*c));
        if (*c == '\0') {
          P_ERR("unbalanced column index identifier:\n");
          print_exp_error(exp, nend);
          free(res); return FCFC_ERR_CFG;
        }
        if (*c != FCFC_COL_IDX_END) {
          P_ERR("invalid column index:\n");
          print_exp_error(exp, nend + 1);
          free(res); return FCFC_ERR_CFG;
        }
        /* Convert the string to column index. */
        for (const char *p = nend + 1; p < c; p++) {
          int digit = *p - '0';
          if ((INT_MAX - digit) / 10 < cidx) {
            P_ERR("the column index is too large:\n");
            print_exp_error(exp, nend + 1);
            free(res); return FCFC_ERR_CFG;
          }
          cidx = cidx * 10 + digit;
        }
        if (!cidx) {
          P_ERR("invalid column index (starting from 1):\n");
          print_exp_error(exp, nend + 1);
          free(res); return FCFC_ERR_CFG;
        }
        cidx--;         /* internal index starting from 0. */
        c++;
      }

      /* Check if the dataset name and index exist already. */
      int idx = 0;
      for (; idx < hcol->ncol; idx++) {
        if (cidx == hcol->idx[idx] && !strncmp(nstart, hcol->dname[idx],
            nend - nstart) && !(hcol->dname[idx][nend - nstart])) break;
      }
      if (idx == hcol->ncol) {          /* record this dataset and index */
        char *dname = malloc((nend - nstart + 1) * sizeof(char));
        if (!dname) {
          P_ERR("failed to allocate memory for dataset name\n");
          free(res); return FCFC_ERR_MEMORY;
        }
        memcpy(dname, nstart, nend - nstart);
        dname[nend - nstart] = '\0';
        if (hdf5_col_insert(hcol, dname, cidx)) {
          free(res); free(dname); return FCFC_ERR_MEMORY;
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
      dset_name = false;

      if (*c != AST_VAR_END) {
        P_ERR("invalid dataset name identifier:\n");
        print_exp_error(exp, c);
        free(res); return FCFC_ERR_CFG;
      }
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
Function `load_hdf5_datasets`:
  Read HDF5 datasets.
Arguments:
  * `hcol`:     structure for storing HDF5 column information;
  * `start`:    first row to be read;
  * `nchunk`:   number of rows to be read at once.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int load_hdf5_datasets(HDF5_COL *hcol, const hsize_t start,
    const hsize_t nchunk) {
  hsize_t offset[2], count[2];
  const hsize_t stride[2] = {1, 1};
  const hsize_t block[2] = {1, 1};

  for (int i = 0; i < hcol->ncol; i++) {
    /* Set the hyperslab. */
    offset[hcol->dim[i]] = start;
    offset[1 - hcol->dim[i]] = hcol->idx[i];
    count[hcol->dim[i]] = nchunk;
    count[1 - hcol->dim[i]] = 1;
    if (H5Sselect_hyperslab(hcol->space_ids[i], H5S_SELECT_SET,
        offset, stride, count, block)) {
      P_ERR("failed to set the hyperslab for reading HDF5 file by chunks\n");
      return FCFC_ERR_HDF5;
    }

    /* Read data. */
    void *read = (hcol->dtypes[i] == FCFC_HDF5_DTYPE_VAR_STR) ?
        hcol->vstr[i] : hcol->data[i];
    if (H5Dread(hcol->dset_ids[i], hcol->memtypes[i], hcol->memspace,
        hcol->space_ids[i], H5P_DEFAULT, read) < 0) {
      P_ERR("failed to read dataset '%s'\n", hcol->dname[i]);
      return FCFC_ERR_HDF5;
    }
  }
  return 0;
}

/******************************************************************************
Function `hdf5_apply_sel`:
  Apply the selection criteria based on the AST and HDF5 columns.
Arguments:
  * `ast`:      abstract syntax tree for the selection criteria;
  * `hcol`:     structure for storing HDF5 column information;
  * `idx`:      index of the row to be processed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int hdf5_apply_sel(ast_t *ast, const HDF5_COL *hcol, const hsize_t idx,
    bool *res) {
  for (long i = 0; i < ast->nvar; i++) {
    const long j = ast->vidx[i];
    const long k = j - 1;
    switch (hcol->dtypes[k]) {
      case FCFC_HDF5_DTYPE_DOUBLE:
        if (ast_set_var(ast, j, ((double *) hcol->data[k]) + idx, 0,
            AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_FLOAT:
        if (ast_set_var(ast, j, ((float *) hcol->data[k]) + idx, 0,
            AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_LONG:
        if (ast_set_var(ast, j, ((long *) hcol->data[k]) + idx, 0,
            AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_INT:
        if (ast_set_var(ast, j, ((int *) hcol->data[k]) + idx, 0,
            AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_FIX_STR:
        if (ast_set_var(ast, j, ((char *) hcol->data[k]) +
            idx * (hcol->dsizes[k] + 1), hcol->dsizes[k], AST_DTYPE_STRING)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_VAR_STR:
        if (ast_set_var(ast, j, hcol->vstr[k][idx],
            strlen(hcol->vstr[k][idx]), AST_DTYPE_STRING)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
    }
  }

  if (ast_eval(ast, res)) {
    P_AST_ERR(ast);
    return FCFC_ERR_AST;
  }
  return 0;
}

/******************************************************************************
Function `hdf5_read_real`:
  Read a floating-point number based on the AST and HDF5 columns.
Arguments:
  * `ast`:      abstract syntax tree for the number evaluation;
  * `hcol`:     structure for storing HDF5 column information;
  * `idx`:      index of the row to be processed;
  * `res`:      the result.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int hdf5_read_real(ast_t *ast, const HDF5_COL *hcol, const hsize_t idx,
    real *res) {
  for (long i = 0; i < ast->nvar; i++) {
    const long j = ast->vidx[i];
    const long k = j - 1;
    switch(hcol->dtypes[k]) {
      case FCFC_HDF5_DTYPE_DOUBLE:
        if (ast_set_var(ast, j, ((double *) hcol->data[k]) + idx, 0,
            AST_DTYPE_DOUBLE)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_FLOAT:
        if (ast_set_var(ast, j, ((float *) hcol->data[k]) + idx, 0,
            AST_DTYPE_FLOAT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_LONG:
        if (ast_set_var(ast, j, ((long *) hcol->data[k]) + idx, 0,
            AST_DTYPE_LONG)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      case FCFC_HDF5_DTYPE_INT:
        if (ast_set_var(ast, j, ((int *) hcol->data[k]) + idx, 0,
            AST_DTYPE_INT)) {
          P_AST_ERR(ast); return FCFC_ERR_AST;
        }
        break;
      default:
        P_ERR("strings can not be used for numerical evaluation: '%s'\n",
            hcol->dname[k]);
        return FCFC_ERR_CFG;
    }
  }

  if (ast_eval(ast, res)) {
    P_AST_ERR(ast);
    return FCFC_ERR_AST;
  }
  return 0;
}


/*============================================================================*\
                       Interface for reading an HDF5 file
\*============================================================================*/

/******************************************************************************
Function `read_hdf5_data`:
  Read an HDF5 file for the positions and weights.
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
int read_hdf5_data(const char *fname, char *const *rcol_ids, const int nrcol,
    const char *sel, real ***rout, size_t *num, const int verb) {
  /* Validate arguments. */
  if (!fname) {
    P_ERR("the HDF5 filename is not specified\n");
    return FCFC_ERR_ARG;
  }
  if (!rcol_ids || !rout || !num) {
    P_ERR("HDF5 data reader not initialized\n");
    return FCFC_ERR_ARG;
  }
  if (nrcol < 3 || nrcol > 4) {
    P_ERR("unexpected number of floating-point columns to be read: %d\n",
        nrcol);
    return FCFC_ERR_ARG;
  }

  /* Suppress default HDF5 error messages if desired. */
  if (FCFC_SUP_HDF5_ERR_MSG) {
    if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0)
      P_WRN("unable to suppress default HDF5 error messages\n");
  }

  /* Open the file for reading. */
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    P_ERR("failed to open HDF5 file for reading: `%s'\n", fname);
    return FCFC_ERR_FILE;
  }
  if (verb) printf("  Filename: %s\n", fname);

  /* Initialize arrays. */
  real **res = NULL;            /* arrays for the outputs */
  char **exp = NULL;            /* AST expressions        */
  ast_t **ast_real = NULL;      /* ASTs for the outputs   */
  ast_t *ast_sel = NULL;        /* AST for the selection  */
  const hsize_t nchunk = FCFC_HDF5_CHUNK_SIZE;
  HDF5_COL *hcol = hdf5_col_init(nchunk);
  if (!hcol) {
    CLEAN_HDF5_PTR; return FCFC_ERR_MEMORY;
  }

  /* Parse expressions with HDF5 datasets and construct libast expressions. */
  if (!(exp = malloc((nrcol + 1) * sizeof(char *)))) {
    P_ERR("failed to allocate memory for the AST expressions\n");
    CLEAN_HDF5_PTR; return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i <= nrcol; i++) exp[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
    if (get_hdf5_cols(rcol_ids[i], hcol, exp + i)) {
      CLEAN_HDF5_PTR; return FCFC_ERR_HDF5;
    }
  }
  if (sel && *sel && get_hdf5_cols(sel, hcol, exp + nrcol)) {
    CLEAN_HDF5_PTR; return FCFC_ERR_HDF5;
  }

  /* Construct ASTs. */
  if (!(ast_real = malloc(sizeof(ast_t *) * nrcol))) {
    P_ERR("failed to allocate memory for the Abstract Syntax Trees\n");
    CLEAN_HDF5_PTR; return FCFC_ERR_MEMORY;
  }
  for (int i = 0; i < nrcol; i++) ast_real[i] = NULL;
  for (int i = 0; i < nrcol; i++) {
    if (!(ast_real[i] = ast_init())) {
      P_AST_ERR(ast_real[i]); CLEAN_HDF5_PTR; return FCFC_ERR_AST;
    }
    if (ast_build(ast_real[i], exp[i], FCFC_AST_REAL, true)) {
      P_AST_ERR(ast_real[i]); CLEAN_HDF5_PTR; return FCFC_ERR_AST;
    }
    if (ast_real[i]->nvar == 0)
      P_WRN("the column expression indicates a constant: `%s'\n", rcol_ids[i]);
  }
  if (exp[nrcol]) {
    if (!(ast_sel = ast_init())) {
      P_AST_ERR(ast_sel); CLEAN_HDF5_PTR; return FCFC_ERR_AST;
    }
    if (ast_build(ast_sel, exp[nrcol], AST_DTYPE_BOOL, true)) {
      P_AST_ERR(ast_sel); CLEAN_HDF5_PTR; return FCFC_ERR_AST;
    }
    if (ast_sel->nvar == 0)
      P_WRN("the expression for data selection is a constant: `%s'\n", sel);
  }

  /* Check dataset information and allocate memory for file reading. */
  hsize_t ntot;         /* number of objects in the HDF5 file */
  if (hdf5_col_finalize(file_id, hcol, nchunk, &ntot)) {
    CLEAN_HDF5_PTR; return FCFC_ERR_HDF5;
  }

  /* Allocate memory for the outputs. */
  if (!(res = malloc(sizeof(real *) * nrcol))) {
    P_ERR("failed to allocate memory for the data\n");
    CLEAN_HDF5_PTR; return FCFC_ERR_MEMORY;
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
      CLEAN_HDF5_PTR; return FCFC_ERR_MEMORY;
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
      FCFC_HDF5_QUIT(FCFC_ERR_MEMORY);
    }
  }
  #pragma omp parallel for num_threads(nomp)
  for (int j = 0; j < nomp - 1; j++) {
    /* Positions and weights. */
    for (int i = 0; i < nrcol; i++) {
      int k = j * nrcol + i;
      if (!(ast_preal[k] = ast_init())) {
        P_AST_ERR(ast_preal[k]); FCFC_HDF5_QUIT(FCFC_ERR_AST);
      }
      if (ast_build(ast_preal[k], exp[i], FCFC_AST_REAL, true)) {
        P_AST_ERR(ast_preal[k]); FCFC_HDF5_QUIT(FCFC_ERR_AST);
      }
    }
    /* Selection. */
    if (exp[nrcol]) {
      if (!(ast_psel[j] = ast_init())) {
        P_AST_ERR(ast_psel[j]); FCFC_HDF5_QUIT(FCFC_ERR_AST);
      }
      if (ast_build(ast_psel[j], exp[nrcol], AST_DTYPE_BOOL, true)) {
        P_AST_ERR(ast_psel[j]); FCFC_HDF5_QUIT(FCFC_ERR_AST);
      }
    }
    else ast_psel[j] = NULL;
  }
  /* Construct the private data pool. */
  real *pdata = malloc(sizeof(real) * nomp * nrcol * FCFC_DATA_THREAD_NUM);
  size_t *pndata = calloc(nomp, sizeof(size_t));
  if (!pdata || !pndata) {
    P_ERR("failed to allocate memory for the thread-private data\n");
    FCFC_HDF5_QUIT(FCFC_ERR_MEMORY);
  }
#endif

  /* Read the file by chunks. */
  hsize_t nread = 0;
  hsize_t nrest = ntot;
  size_t n = 0;
  while (nrest) {
    hsize_t nrow = nchunk;
    if (nchunk > nrest){
      nrow = nrest;
      /* Revise the dataspace in memory if there are not enough objects. */
      if ((hcol->memspace = H5Screate_simple(1, &nrow, NULL)) < 0) {
        P_ERR("failed to revise the dataspace for storing HDF5 data\n");
        CLEAN_HDF5_PTR; CLEAN_HDF5_PTR_OMP; return FCFC_ERR_HDF5;
      }
    }
    /* Read the necessary columns. */
    if (load_hdf5_datasets(hcol, nread, nrow)) {
      CLEAN_HDF5_PTR; CLEAN_HDF5_PTR_OMP; return FCFC_ERR_HDF5;
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
      for (hsize_t i = 0; i < nrow; i++) {
        /* Apply selection. */
        if (ast_sel) {
          bool keep = false;
          if (hdf5_apply_sel(ast_sel, hcol, i, &keep)) {
            FCFC_HDF5_QUIT(FCFC_ERR_AST);
          }
          if (!keep) continue;
        }

        /* Record coordinates and weights in the private data pool. */
        for (int j = 0; j < nrcol; j++) {
          if (hdf5_read_real(ast_real[j], hcol, i,
              pres + j * FCFC_DATA_THREAD_NUM + *pnum)) {
            FCFC_HDF5_QUIT(FCFC_ERR_AST);
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
    for (hsize_t i = 0; i < nrow; i++) {
      /* Apply selection. */
      if (ast_sel) {
        bool keep = false;
        if (hdf5_apply_sel(ast_sel, hcol, i, &keep)) {
          CLEAN_HDF5_PTR; return FCFC_ERR_HDF5;
        }
        if (!keep) continue;
      }

      /* Read coordinates and weights. */
      for (int j = 0; j < nrcol; j++) {
        if (hdf5_read_real(ast_real[j], hcol, i, res[j] + n)) {
          CLEAN_HDF5_PTR; return FCFC_ERR_HDF5;
        }
      }
      n++;
    }
#endif

    /* Release memory allocated for variable-length strings. */
    for (int i = 0; i < hcol->ncol; i++) {
      if (hcol->dtypes[i] == FCFC_HDF5_DTYPE_VAR_STR) {
        if (H5Dvlen_reclaim(hcol->memtypes[i], hcol->memspace, H5P_DEFAULT,
            hcol->vstr[i]) < 0) {
          P_WRN("failed to realease memory allocated for variable-length "
              "strings");
        }
      }
    }

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
  CLEAN_HDF5_PTR_OMP;
#endif

  if (!n) {
    P_ERR("no object is read from file: `%s'\n", fname);
    CLEAN_HDF5_PTR; return FCFC_ERR_FILE;
  }

  if (verb) {
    printf("  Number of rows processed in total: %zu\n"
        "  Number of recorded objects: %zu\n", (size_t) ntot, n);
  }

  hdf5_col_destroy(hcol);
  if (H5Fclose(file_id) < 0) P_WRN("failed to close file: `%s'\n", fname);
  for (int i = 0; i < nrcol; i++) ast_destroy(ast_real[i]);
  free(ast_real);
  ast_destroy(ast_sel);
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

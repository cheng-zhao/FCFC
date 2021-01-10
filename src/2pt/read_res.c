/*******************************************************************************
* 2pt/read_res.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_res.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

/*============================================================================*\
                         Function for memory allocation
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


/*============================================================================*\
                       Interface for reading pair counts
\*============================================================================*/

/******************************************************************************
Function `read_pair_count`:
  Read pair counts from a text file.
Arguments:
  * `fname`:    name of the input file for pair counts;
  * `nx`:       number of separation or s_perp bins;
  * `ny`:       number of mu or pi bins;
  * `res`:      array for storing the pair counts to be read from file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_pair_count(const char *fname, const int nx, const int ny,
    double *res) {
  const size_t ntot = nx * ny;
  if (nx <= 0 || ny <= 0 || !ntot) {
    P_ERR("invalid dimension for pair counts: (%d,%d)\n", nx, ny);
    return FCFC_ERR_ARG;
  }
  if (!res) {
    P_ERR("memory not allocated for storing pair counts\n");
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

  enum {
    PAIR_COUNT_FMT_NULL,        /* invalid pair count format   */
    PAIR_COUNT_FMT_ISO,         /* isotropic pair counts       */
    PAIR_COUNT_FMT_LIST,        /* 2-D pair counts as a list   */
    PAIR_COUNT_FMT_MATRIX       /* 2-D pair counts as a matrix */
  } fmt = PAIR_COUNT_FMT_NULL;

  int nheader = 0;
  size_t n, nread, nrest;
  n = nrest = 0;

  /* Check the header for the file format. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      nheader += 1;
      if (*p == FCFC_SAVE_COMMENT) {
        if (nheader == 1) {
          const char *line = "Created by " FCFC_CODE_NAME;
          if (strncmp(line, p + 2, strlen(line) + 1)) {
            P_ERR("invalid pair count file: `%s'\n", fname);
            fclose(fp); free(chunk);
            return FCFC_ERR_FILE;
          }
        }
        else if (nheader == 4) {
          int num = 0;
          if (sscanf(p + 2, "Rows: %d %*s bins", &num) == 1) {
            fmt = PAIR_COUNT_FMT_MATRIX;
            if (num != nx) {
              P_ERR("invalid number of rows in file: `%s'\n", fname);
              fclose(fp); free(chunk);
              return FCFC_ERR_FILE;
            }
          }
          else if (!strncmp(p + 2,
              "Columns: s_min(1) s_max(2) normalized_pair_count(3)", 51)) {
            fmt = PAIR_COUNT_FMT_ISO;
            if (ny != 1) {
              P_ERR("invalid second dimension for pair counts: %d\n", ny);
              return FCFC_ERR_ARG;
            }
          }
          else if (!strncmp(p + 2, "Columns: ", 9)) fmt = PAIR_COUNT_FMT_LIST;
        }
        else if (fmt == PAIR_COUNT_FMT_MATRIX && nheader == 6) {
          int num = 0;
          if (sscanf(p + 2, "Columns: %d %*s bins", &num) == 1) {
            if (num != ny) {
              P_ERR("invalid number of columns in file: `%s'\n", fname);
              fclose(fp); free(chunk);
              return FCFC_ERR_FILE;
            }
          }
          else fmt = PAIR_COUNT_FMT_NULL;
        }
        p = endl + 1;
        continue;
      }
      else if (isgraph(*p)) {
        if (fmt == PAIR_COUNT_FMT_NULL) {
          P_ERR("invalid format of pair count file: `%s'\n", fname);
          fclose(fp); free(chunk);
          return FCFC_ERR_FILE;
        }
        else break;
      }
      else {
        P_ERR("invalid starting character of line %d in file: `%s'\n",
            nheader, fname);
        fclose(fp); free(chunk);
        return FCFC_ERR_FILE;
      }
    }

    if (fmt != PAIR_COUNT_FMT_NULL) {           /* finish reading header */
      *endl = '\n';
      nrest = end - p;
      memmove(chunk, p, nrest);
      break;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        fclose(fp); free(chunk);
        return FCFC_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }
  nread = 0;

  if (fmt == PAIR_COUNT_FMT_MATRIX) {
    /* Start reading the pair counts by chunk. */
    do {
      char *p = chunk;
      char *end = p + nrest + nread;
      char *endl;
      if (nread && nread < csize - nrest) *end = '\n';

      /* Process lines in the chunk. */
      while ((endl = memchr(p, '\n', end - p))) {
        *endl = '\0';             /* replace '\n' by string terminator '\0' */
        while (isspace(*p)) ++p;  /* omit leading whitespaces */
        if (*p == '\0') {
          p = endl + 1;
          continue;
        }

        if (n == (size_t) nx) {
          P_ERR("too many rows of pair count file: `%s'\n", fname);
          fclose(fp); free(chunk);
          return FCFC_ERR_FILE;
        }

        /* Parse the line (there should be no header already). */
        for (int i = 0; i < ny; i++) {
          int cnt = 0;
          if (sscanf(p, "%lf%n", res + n + i * nx, &cnt) != 1 || !cnt) {
            P_ERR("failed to read line %zu of file: `%s'\n",
                nheader + n + 1, fname);
            fclose(fp); free(chunk);
            return FCFC_ERR_FILE;
          }
          p += cnt;
        }

        /* Continue with the next line. */
        n += 1;
        p = endl + 1;
      }

      /* The chunk cannot hold a full line. */
      if (p == chunk) {
        if (chunk_resize(&chunk, &csize)) {
          P_ERR("failed to allocate memory for reading the file by chunk\n");
          fclose(fp); free(chunk);
          return FCFC_ERR_MEMORY;
        }
        nrest += nread;
        continue;
      }

      /* Copy the remaining characters to the beginning of the chunk. */
      nrest = end - p;
      memmove(chunk, p, nrest);
    }
    while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp)));
    if (n != (size_t) nx) {
      P_ERR("unexpected number of records (%zu) of pair count file: `%s'\n",
          n, fname);
      fclose(fp); free(chunk);
      return FCFC_ERR_FILE;
    }
  }
  else {        /* fmt == PAIR_COUNT_FMT_LIST or fmt == PAIR_COUNT_FMT_ISO */
    const char *line_fmt = (fmt == PAIR_COUNT_FMT_ISO) ?
        "%*lf %*lf %lf" : "%*lf %*lf %*lf %*lf %lf";

    /* Start reading the pair counts by chunk. */
    do {
      char *p = chunk;
      char *end = p + nrest + nread;
      char *endl;
      if (nread && nread < csize - nrest) *end = '\n';

      /* Process lines in the chunk. */
      while ((endl = memchr(p, '\n', end - p))) {
        *endl = '\0';             /* replace '\n' by string terminator '\0' */
        while (isspace(*p)) ++p;  /* omit leading whitespaces */
        if (*p == '\0') {
          p = endl + 1;
          continue;
        }

        if (n == ntot) {
          P_ERR("too many rows of pair count file: `%s'\n", fname);
          fclose(fp); free(chunk);
          return FCFC_ERR_FILE;
        }

        /* Parse the line (there should be no header already). */
        if (sscanf(p, line_fmt, res + n) != 1) {
          P_ERR("failed to read line %zu of file: `%s'\n",
              nheader + n + 1, fname);
          return FCFC_ERR_FILE;
        }

        /* Continue with the next line. */
        n += 1;
        p = endl + 1;
      }

      /* The chunk cannot hold a full line. */
      if (p == chunk) {
        if (chunk_resize(&chunk, &csize)) {
          P_ERR("failed to allocate memory for reading the file by chunk\n");
          fclose(fp); free(chunk);
          return FCFC_ERR_MEMORY;
        }
        nrest += nread;
        continue;
      }

      /* Copy the remaining characters to the beginning of the chunk. */
      nrest = end - p;
      memmove(chunk, p, nrest);
    }
    while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp)));
    if (n != ntot) {
      P_ERR("unexpected number of records (%zu) of pair count file: `%s'\n",
          n, fname);
      fclose(fp); free(chunk);
      return FCFC_ERR_FILE;
    }
  }

  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    fclose(fp); free(chunk);
    return FCFC_ERR_FILE;
  }

  free(chunk);
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  return 0;
}

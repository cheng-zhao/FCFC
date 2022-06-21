/*******************************************************************************
* 2pt_box/read_res.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_res.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

/*============================================================================*\
                  Data structure for reading a pair count file
\*============================================================================*/

typedef struct {
  FILE *fp;                     /* pointer to the input file */
  char *chunk;                  /* chunk for storing file contents */
  size_t size;                  /* size of data stored in chunk */
  size_t capacity;              /* capacity of the chunk */
  const char *fname;            /* name of the input file */
} CNTFILE;


/*============================================================================*\
             Functions for manipulating the file reading structure
\*============================================================================*/

/******************************************************************************
Function `cfile_init`:
  Initialize the structure for reading pair counting files.
Arguments:
  * `fname`:    name of the input pair count file.
Return:
  Address of the structure on success; NULL on error.
******************************************************************************/
static CNTFILE *cfile_init(const char *fname) {
  CNTFILE *cfile = calloc(1, sizeof(CNTFILE));
  cfile->fname = fname;

  /* Open the file for reading. */
  if (!(cfile->fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    free(cfile); return NULL;
  }

  /* Allocate memory for the chunk. */
  cfile->capacity = FCFC_FILE_CHUNK;
  if (!(cfile->chunk = malloc(cfile->capacity * sizeof(char)))) {
    P_ERR("failed to allocate memory for reading the file by chunks\n");
    fclose(cfile->fp); free(cfile); return NULL;
  }

  return cfile;
}

/******************************************************************************
Function `cfile_destroy`:
  Deconstruct the structure for reading pair counting files.
Arguments:
  * `cfile`:    structure for reading pair count files.
******************************************************************************/
static void cfile_destroy(CNTFILE *cfile) {
  if (!cfile) return;
  if (cfile->fp) fclose(cfile->fp);
  if (cfile->chunk) free(cfile->chunk);
  free(cfile);
}

/******************************************************************************
Function `chunk_resize`:
  Enlarge the size of a chunk.
Arguments:
  * `cfile`:    structure for reading pair count files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int chunk_resize(CNTFILE *cfile) {
  if (FCFC_MAX_CHUNK / 2 < cfile->capacity) {
    P_ERR("cannot load the input file into memory: `%s'\n", cfile->fname);
    return FCFC_ERR_FILE;
  }
  cfile->capacity <<= 1;

  char *tmp = realloc(cfile->chunk, cfile->capacity * sizeof(char));
  if (!tmp) {
    P_ERR("failed to allocate memory for reading the file by chunks\n");
    return FCFC_ERR_MEMORY;
  }

  cfile->chunk = tmp;
  return 0;
}


/*============================================================================*\
                       Functions for reading pair counts
\*============================================================================*/

/******************************************************************************
Function `read_pair_ascii`:
  Read pair counts from a text file.
Arguments:
  * `cfile`:    structure for reading pair count files;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the pair counts or correlation function to be read.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_pair_ascii(CNTFILE *cfile, CF *cf, const int idx) {
  size_t nread = 0, nline = 1;
  do {
    char *p = cfile->chunk;
    char *end = p + cfile->size + nread;
    char *endl;

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      nline += 1;
      if (*p == FCFC_SAVE_COMMENT) {
        /* Determine the file format based on the 4th row. */
        if (nline == 4) {
          const char header_iso[] =
              "Columns: s_min(1) s_max(2) normalized_pair_count(3)";
          const char header_smu[] = "Columns: s_min(1) s_max(2) mu_min(3) "
              "mu_max(4) normalized_pair_count(5)";
          const char header_spi[] = "Columns: s_perp_min(1) s_perp_max(2) "
              "pi_min(3) pi_max(4) normalized_pair_count(5)";
          if ((cf->bintype == FCFC_BIN_ISO && strncmp(p + 2, header_iso,
                  strlen(header_iso))) ||
              (cf->bintype == FCFC_BIN_SMU && strncmp(p + 2, header_smu,
                  strlen(header_smu))) ||
              (cf->bintype == FCFC_BIN_SPI && strncmp(p + 2, header_spi,
                  strlen(header_spi)))) {
            P_ERR("unexpected separation bins in pair count file: `%s'\n",
                cfile->fname);
            return FCFC_ERR_FILE;
          }
        }
        p = endl + 1;
        continue;
      }
      else {
        if (nline != 5) {
          P_ERR("invalid number of header lines in pair count file: `%s'\n",
              cfile->fname);
          return FCFC_ERR_FILE;
        }
        break;
      }
    }

    /* The chunk cannot hold a full line. */
    if (p == cfile->chunk) {
      if (chunk_resize(cfile)) return FCFC_ERR_MEMORY;
      cfile->size += nread;
      continue;
    }

    if (nline != 5) {
      cfile->size = end - p;
      memmove(cfile->chunk, p, cfile->size);
    }
    else {
      *endl = '\n';
      cfile->size = end - p;
      memmove(cfile->chunk, p, cfile->size);
      break;
    }
  }
  while ((nread = fread(cfile->chunk + cfile->size, sizeof(char),
      cfile->capacity - cfile->size, cfile->fp)));

  /* Read pair counts. */
  const char *line_fmt = (cf->bintype == FCFC_BIN_ISO) ?
      "%*lf %*lf %lf" : "%*lf %*lf %*lf %*lf %lf";

  nread = nline = 0;
  do {
    char *p = cfile->chunk;
    char *end = p + cfile->size + nread;
    char *endl;

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;  /* omit leading whitespaces */
      if (*p == '\0') {
        p = endl + 1;
        continue;
      }

      if (nline == cf->ntot) {
        P_ERR("too many rows in pair count file: `%s'\n", cfile->fname);
        return FCFC_ERR_FILE;
      }

      /* Parse the line. */
      if (sscanf(p, line_fmt, cf->ncnt[idx] + nline) != 1) {
        P_ERR("failed to read line %zu of file: `%s'\n", nline + 5,
            cfile->fname);
        return FCFC_ERR_FILE;
      }

      /* Continue with the next line. */
      nline++;
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == cfile->chunk) {
      if (chunk_resize(cfile)) return FCFC_ERR_MEMORY;
      cfile->size += nread;
      continue;
    }

    /* Copy the remaining data to the beginning of the chunk. */
    cfile->size = end - p;
    memmove(cfile->chunk, p, cfile->size);
  }
  while ((nread = fread(cfile->chunk + cfile->size, sizeof(char),
      cfile->capacity - cfile->size, cfile->fp)));

  if (nline != cf->ntot) {
    P_ERR("unexpected number of bins (%zu) in pair count file: `%s'\n",
        nline, cfile->fname);
    return FCFC_ERR_FILE;
  }

  return 0;
}

/******************************************************************************
Function `read_pair_bin`:
  Read pair counts from a binary file.
Arguments:
  * `cfile`:    structure for reading pair count files;
  * `conf`:     structure for storing all configurations;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the pair counts or correlation function to be read.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_pair_bin(CNTFILE *cfile, CF *cf, const int idx) {
  uint64_t nblock = 0, nread = 0, size;
  bool single_prec = false, auto_cnt = false, with_wt = false;
  const uint64_t off = (cf->bintype == FCFC_BIN_SPI) ? 1 : 0;
  do {
    cfile->size += nread;

    /* Process blocks. */
    while (cfile->size >= sizeof(uint64_t)) {
      memcpy(&size, cfile->chunk, sizeof(uint64_t));
      if (cfile->size < 2 * sizeof(uint64_t) + size) break;

      char *p = cfile->chunk + sizeof(uint64_t);
      if (nblock == 0) {                /* pair count settings */
        if ((cf->bintype == FCFC_BIN_ISO && size != 7 * sizeof(int32_t)) ||
            (cf->bintype != FCFC_BIN_ISO && size != 8 * sizeof(int32_t))) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
        int32_t spec[8];
        memcpy(spec, p, size);
        if (spec[0] == 0) single_prec = true;
  #ifdef SINGLE_PREC
        if (spec[0] == 1) {
          P_WRN("pair count file was generated with double precision: `%s'\n",
              cfile->fname);
        }
  #else
        if (spec[0] == 0) {
          P_WRN("pair count file was generated with single precision: `%s'\n",
              cfile->fname);
        }
  #endif
  #ifdef WITH_MU_ONE
        if (spec[1] == 0 && cf->bintype == FCFC_BIN_SMU) {
          P_WRN("pair count file was generated without `-DWITH_MU_ONE`: `%s'\n",
              cfile->fname);
        }
  #else
        if (spec[1] == 1 && cf->bintype == FCFC_BIN_SMU) {
          P_WRN("pair count file was generated with `-DWITH_MU_ONE`: `%s'\n",
              cfile->fname);
        }
  #endif
        if (spec[2] != 1) {
          P_ERR("non-periodic pair counts in file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (spec[3] == FCFC_COUNT_AUTO) auto_cnt = true;
        if (spec[4] == FCFC_COUNT_WITH_WT) with_wt = true;
        if (spec[5] != (int32_t) cf->bintype) {
          P_ERR("unexpected binning scheme in pair count file: `%s'\n",
              cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (spec[6] != (int32_t) cf->ns) {
          P_ERR("unexpected number of separation bins (%d) in pair count file: "
              "`%s'\n", spec[6], cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (cf->bintype == FCFC_BIN_SMU && spec[7] != (int32_t) cf->nmu) {
          P_ERR("unexpected number of mu bins (%d) in pair count file: `%s'\n",
              spec[7], cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (cf->bintype == FCFC_BIN_SPI && spec[7] != (int32_t) cf->np) {
          P_ERR("unexpected number of pi bins (%d) in pair count file: `%s'\n",
              spec[7], cfile->fname);
        }
      }
      else if (nblock == 1) {           /* box sizes */
        if ((single_prec && size != 3 * sizeof(float)) ||
            (!single_prec && size != 3 * sizeof(double))) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (single_prec) {
          float b[3];
          memcpy(b, p, size);
          if (fabs(b[0] - cf->bsize[0]) > FLOAT_TOL ||
              fabs(b[1] - cf->bsize[1]) > FLOAT_TOL ||
              fabs(b[2] - cf->bsize[2]) > FLOAT_TOL) {
            P_ERR("unexpected box size (" OFMT_DBL "," OFMT_DBL "," OFMT_DBL
                ") in pair count file: `%s'\n", b[0], b[1], b[2], cfile->fname);
            return FCFC_ERR_FILE;
          }
        }
        else {
          double b[3];
          memcpy(b, p, size);
          if (fabs(b[0] - cf->bsize[0]) > REAL_TOL ||
              fabs(b[1] - cf->bsize[1]) > REAL_TOL ||
              fabs(b[2] - cf->bsize[2]) > REAL_TOL) {
            P_ERR("unexpected box size (" OFMT_DBL "," OFMT_DBL "," OFMT_DBL
                ") in pair count file: `%s'\n", b[0], b[1], b[2], cfile->fname);
            return FCFC_ERR_FILE;
          }
        }
      }
      else if (nblock == 2) {           /* separation bin edges */
        if ((single_prec && size != (cf->ns + 1) * sizeof(float)) ||
            (!single_prec && size != (cf->ns + 1) * sizeof(double))) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (single_prec) {
          for (int i = 0; i <= cf->ns; i++) {
            float sbin;
            memcpy(&sbin, p + i * sizeof(float), sizeof(float));
            if (fabs(sbin - cf->sbin_raw[i]) > FLOAT_TOL) {
              P_ERR("unexpected separation bin in pair count file: `%s'\n",
                  cfile->fname);
              return FCFC_ERR_FILE;
            }
          }
        }
        else {
          for (int i = 0; i <= cf->ns; i++) {
            double sbin;
            memcpy(&sbin, p + i * sizeof(double), sizeof(double));
            if (fabs(sbin - cf->sbin_raw[i]) > REAL_TOL) {
              P_ERR("unexpected separation bin in pair count file: `%s'\n",
                  cfile->fname);
              return FCFC_ERR_FILE;
            }
          }
        }
      }
      else if (nblock == 3 && cf->bintype == FCFC_BIN_SPI) {    /* pi bins */
        if ((single_prec && size != (cf->np + 1) * sizeof(float)) ||
            (!single_prec && size != (cf->np + 1) * sizeof(double))) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
        if (single_prec) {
          for (int i = 0; i <= cf->np; i++) {
            float pbin;
            memcpy(&pbin, p + i * sizeof(float), sizeof(float));
            if (fabs(pbin - cf->pbin_raw[i]) > FLOAT_TOL) {
              P_ERR("unexpected pi bin in pair count file: `%s'\n",
                  cfile->fname);
              return FCFC_ERR_FILE;
            }
          }
        }
        else {
          for (int i = 0; i <= cf->np; i++) {
            double pbin;
            memcpy(&pbin, p + i * sizeof(double), sizeof(double));
            if (fabs(pbin - cf->pbin_raw[i]) > REAL_TOL) {
              P_ERR("unexpected pi bin in pair count file: `%s'\n",
                  cfile->fname);
              return FCFC_ERR_FILE;
            }
          }
        }
      }
      else if (nblock == 3 + off) {     /* number of tracers */
        const size_t nsize = (with_wt) ? sizeof(double) : sizeof(uint64_t);
        const int num = (auto_cnt) ? 1 : 2;
        if (size != nsize * num) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
      }
      else if (nblock == 4 + off) {     /* normalization factor */
        if (size != sizeof(double)) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
      }
      else if (nblock == 5 + off) {     /* normalized pair counts */
        if (size != cf->ntot * sizeof(double)) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
        memcpy(cf->ncnt[idx], p, size);
      }
      else if (nblock == 6 + off) {     /* raw pair counts */
        const size_t nsize = (with_wt) ? sizeof(double) : sizeof(int64_t);
        if (size != cf->ntot * nsize) {
          P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
          return FCFC_ERR_FILE;
        }
      }
      else {
        P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
        return FCFC_ERR_FILE;
      }

      /* Check the size of block. */
      p += size;
      uint64_t size_check;
      memcpy(&size_check, p, sizeof(uint64_t));
      if (size != size_check) {
        P_ERR("invalid format of pair count file: `%s'\n", cfile->fname);
        return FCFC_ERR_FILE;
      }

      /* Continue with the next block. */
      p += sizeof(uint64_t);
      cfile->size -= 2 * sizeof(uint64_t) + size;
      memmove(cfile->chunk, p, cfile->size);
      nblock++;
    }

    /* Stop reading if all the blocks are processed. */
    if (nblock == 7 + off) break;

    /* Enlarge the size of the chunk if necessary. */
    if (cfile->size == cfile->capacity) {
      if (chunk_resize(cfile)) return FCFC_ERR_MEMORY;
      continue;
    }
  }
  while ((nread = fread(cfile->chunk + cfile->size, sizeof(char),
      cfile->capacity - cfile->size, cfile->fp)));

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
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the pair counts or correlation function to be read.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_pair_count(const char *fname, CF *cf, const int idx) {
  /* Validate function arguments. */
  if (!fname || !cf || idx < 0 || idx >= cf->npc) {
    P_ERR("pair counts not initialized\n");
    return FCFC_ERR_ARG;
  }

  /* Initialize file reading. */
  CNTFILE *cfile = cfile_init(fname);
  if (!cfile) return FCFC_ERR_FILE;

  /* Read the first header line. */
  size_t nread;
  int fmt = -1;
  while ((nread = fread(cfile->chunk + cfile->size, sizeof(char),
      cfile->capacity - cfile->size, cfile->fp))) {
    char *p = cfile->chunk;
    char *end = p + cfile->size + nread;
    char *endl = memchr(p, '\n', end - p);
    if (endl == NULL) {
      /* The chunk cannot hold a full line. */
      if (chunk_resize(cfile)) {
        P_ERR("failed to allocate memory for reading the file by chunk\n");
        cfile_destroy(cfile);
        return FCFC_ERR_MEMORY;
      }
      cfile->size += nread;
      continue;
    }

    /* Analyze the header. */
    int ver;
    if (sscanf(p + 2, "Created by " FCFC_CODE_NAME " %*s (%d) ; format = %d",
        &ver, &fmt) != 2) {
      P_ERR("unrecognized pair count file: `%s'\n", fname);
      cfile_destroy(cfile);
      return FCFC_ERR_FILE;
    }
    if (FCFC_VERNUM < ver) {
      P_WRN("the pair count file is generated by a future version of FCFC: "
          "`%s'\n", fname);
    }

    /* Copy the remaining data to the beginning of the chunk. */
    p = endl + 1;
    cfile->size = end - p;
    memmove(cfile->chunk, p, cfile->size);
    break;
  }

  switch (fmt) {
    case FCFC_OFMT_ASCII:
      if (read_pair_ascii(cfile, cf, idx)) {
        cfile_destroy(cfile);
        return FCFC_ERR_FILE;
      }
      break;
    case FCFC_OFMT_BIN:
      if (read_pair_bin(cfile, cf, idx)) {
        cfile_destroy(cfile);
        return FCFC_ERR_FILE;
      }
      break;
    default:
      P_ERR("unrecognized pair count file format: %d\n", fmt);
      cfile_destroy(cfile);
      return FCFC_ERR_FILE;
  }

  if (!feof(cfile->fp)) {
    P_ERR("unexpected end of file: `%s'\n", cfile->fname);
    cfile_destroy(cfile);
    return FCFC_ERR_FILE;
  }

  cfile_destroy(cfile);
  return 0;
}

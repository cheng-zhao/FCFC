/*******************************************************************************
* 2pt_box/save_res.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "write_file.h"
#include "save_res.h"
#include <stdlib.h>

/* Shortcuts for saving outputs. */
#define WRITE_LINE(...)                                                 \
  if (output_writeline(ofile, __VA_ARGS__)) return FCFC_ERR_FILE;

#define WRITE_DATA(src,size)                                            \
  if (output_write(ofile, &(size), sizeof(uint64_t)) ||                 \
      output_write(ofile, (src), (size)) ||                             \
      output_write(ofile, &(size), sizeof(uint64_t)))                   \
    return FCFC_ERR_FILE;


/*============================================================================*\
                   Functions for saving pair counts and 2PCFs
\*============================================================================*/

/******************************************************************************
Function `save_res_ascii`:
  Write the pair counts or correlation function to a text file.
Arguments:
  * `ofile`:    interface for writing to a file;
  * `conf`:     structure for storing all configurations;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the pair counts or correlation function to be saved;
  * `flag`:     identifier of the results to be saved.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int save_res_ascii(OFILE *ofile, const CONF *conf,
    const CF *cf, const int idx, const FCFC_otype_e flag) {
  /* Write the header for pair counts. */
  if (flag == FCFC_OUTPUT_PAIR_COUNT) {
    WRITE_LINE("%c Created by " FCFC_CODE_NAME " v" FCFC_VERSION
        " (%d) ; format = %d\n", FCFC_SAVE_COMMENT, FCFC_VERNUM, conf->ofmt);
    if (cf->pc_idx[0][idx] == cf->pc_idx[1][idx]) {     /* auto pairs */
      WRITE_LINE("%c Number of tracers: %zu (auto pair counts",
          FCFC_SAVE_COMMENT, cf->data[cf->pc_idx[0][idx]].n);
      if (cf->wt[idx]) {
        WRITE_LINE(", weighted number: " OFMT_DBL " ",
            cf->data[cf->pc_idx[0][idx]].wt);
      }
      WRITE_LINE(")\n");
    }
    else {                                              /* cross pairs */
      WRITE_LINE("%c Numbers of tracers: %zu and %zu",
          FCFC_SAVE_COMMENT, cf->data[cf->pc_idx[0][idx]].n,
          cf->data[cf->pc_idx[1][idx]].n);
      if (cf->wt[idx]) {
        WRITE_LINE(" (weighted numbers: " OFMT_DBL " and " OFMT_DBL " )",
            cf->data[cf->pc_idx[0][idx]].wt, cf->data[cf->pc_idx[1][idx]].wt);
      }
      WRITE_LINE("\n");
    }
    WRITE_LINE("%c Normalization (total number of pairs): " OFMT_DBL
        "\n", FCFC_SAVE_COMMENT, cf->norm[idx]);
  }

  /* Write row/column identifiers and the results. */
  if (flag == FCFC_OUTPUT_2PCF_INTEG) {
    if (conf->bintype == FCFC_BIN_SMU) {
      /* Header for 2PCF multipoles. */
      WRITE_LINE("%c Columns: s_cen(1) s_min(2) s_max(3)",
          FCFC_SAVE_COMMENT);
      for (int i = 0; i < cf->nl; i++) {
        WRITE_LINE(" xi_%d(%d)", cf->poles[i], 4 + i);
      }
      WRITE_LINE("\n");
      /* Results of 2PCF multipoles. */
      for (int i = 0; i < cf->ns; i++) {
        WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL,
            (cf->sbin_raw[i] + cf->sbin_raw[i + 1]) * 0.5,
            cf->sbin_raw[i], cf->sbin_raw[i + 1]);
        for (int j = 0; j < cf->nl; j++) {
          WRITE_LINE(" " OFMT_DBL, cf->mp[idx][i + j * cf->ns]);
        }
        WRITE_LINE("\n");
      }
    }
    else {
      /* Header for projected 2PCF. */
      WRITE_LINE("%c Columns: s_perp_cen(1) s_perp_min(2) s_perp_max(3)"
          " wp(4)\n", FCFC_SAVE_COMMENT);
      /* Results of projected 2PCF. */
      for (int i = 0; i < cf->ns; i++) {
        WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL "\n",
            (cf->sbin_raw[i] + cf->sbin_raw[i + 1]) * 0.5,
            cf->sbin_raw[i], cf->sbin_raw[i + 1], cf->wp[idx][i]);
      }
    }
  }
  else {        /* flag == FCFC_OUTPUT_PAIR_COUNT or FCFC_OUTPUT_2PCF_RAW */
    double **res = (flag == FCFC_OUTPUT_PAIR_COUNT) ? cf->ncnt : cf->cf;
    if (conf->bintype == FCFC_BIN_SMU) {
      /* Header for count(s,mu) or xi(s,mu) as a list. */
      WRITE_LINE("%c Columns: s_min(1) s_max(2) mu_min(3) mu_max(4) %s(5)\n",
          FCFC_SAVE_COMMENT, (flag == FCFC_OUTPUT_PAIR_COUNT) ?
          "normalized_pair_count" : "xi");
      /* Results of count(s,mu) or xi(s,mu) as a list. */
      for (int j = 0; j < cf->nmu; j++) {
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL " "
              OFMT_DBL "\n", cf->sbin_raw[i], cf->sbin_raw[i + 1],
              j / (double) cf->nmu, (j + 1) / (double) cf->nmu,
              res[idx][i + j * cf->ns]);
        }
      }
    }
    else if (conf->bintype == FCFC_BIN_SPI) {
      /* Header for count(s_perp,s_par) or xi(s_perp,s_par) as a list. */
      WRITE_LINE("%c Columns: s_perp_min(1) s_perp_max(2) pi_min(3) pi_max(4)"
          " %s(5)\n", FCFC_SAVE_COMMENT,
          (flag == FCFC_OUTPUT_PAIR_COUNT) ? "normalized_pair_count" : "xi");
      /* Results of count(s_perp,s_par) or xi(s_perp,s_par) as a list. */
      for (int j = 0; j < cf->np; j++) {
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL " "
              OFMT_DBL "\n", cf->sbin_raw[i], cf->sbin_raw[i + 1],
              cf->pbin_raw[j], cf->pbin_raw[j + 1], res[idx][i + j * cf->ns]);
        }
      }
    }
    else {      /* conf->bintype == FCFC_BIN_ISO */
      if (flag == FCFC_OUTPUT_PAIR_COUNT) {
        /* Header for isotropic pair counts. */
        WRITE_LINE("%c Columns: s_min(1) s_max(2) normalized_pair_count(3)\n",
            FCFC_SAVE_COMMENT);
        /* Results of isotropic pair counts. */
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL "\n",
              cf->sbin_raw[i], cf->sbin_raw[i + 1], cf->ncnt[idx][i]);
        }
      }
      else {    /* flag == FCFC_OUTPUT_2PCF_RAW */
        /* Header for isotropic 2PCF. */
        WRITE_LINE("%c Columns: s_cen(1) s_min(2) s_max(3) xi_iso(4)\n",
            FCFC_SAVE_COMMENT);
        /* Results of isotropic 2PCF. */
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL "\n",
              (cf->sbin_raw[i] + cf->sbin_raw[i + 1]) * 0.5,
              cf->sbin_raw[i], cf->sbin_raw[i + 1], cf->cf[idx][i]);
        }
      }
    }
  }

  return 0;
}

/******************************************************************************
Function `save_res_bin`:
  Write the pair counts to a binary file.
Arguments:
  * `ofile`:    interface for writing to a file;
  * `conf`:     structure for storing all configurations;
  * `cf`:       structure for correlation function evaluations;
  * `idx`:      index of the pair counts or correlation function to be saved;
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int save_res_bin(OFILE *ofile, const CONF *conf,
    const CF *cf, const int idx) {
  /* Write the header for pair counts. */
  WRITE_LINE("%c Created by " FCFC_CODE_NAME " v" FCFC_VERSION
      " (%d) ; format = %d\n", FCFC_SAVE_COMMENT, FCFC_VERNUM, conf->ofmt);

  /* Write the settings: precision, mu=1, box/lightcone, auto/cross pairs,
   * weights, bins. */
  uint64_t size = (cf->bintype == FCFC_BIN_ISO) ?
      7 * sizeof(int) : 8 * sizeof(int);
  int32_t spec[8];
#ifdef SINGLE_PREC
  spec[0] = 0;
#else
  spec[0] = 1;
#endif
#ifdef WITH_MU_ONE
  spec[1] = 1;
#else
  spec[1] = 0;
#endif
#ifdef FCFC_METRIC_PERIODIC
  spec[2] = 1;
#else
  spec[2] = 0;
#endif
  spec[3] = (cf->pc_idx[0][idx] == cf->pc_idx[1][idx]) ?
      FCFC_COUNT_AUTO : FCFC_COUNT_CROSS;
  spec[4] = (cf->wt[idx]) ? FCFC_COUNT_WITH_WT : FCFC_COUNT_NO_WT;
  spec[5] = cf->bintype;
  spec[6] = cf->ns;
  if (cf->bintype == FCFC_BIN_SMU) spec[7] = cf->nmu;
  else if (cf->bintype == FCFC_BIN_SPI) spec[7] = cf->np;

  WRITE_DATA(spec, size);

#ifdef FCFC_METRIC_PERIODIC
  /* Write the box sizes. */
  size = 3 * sizeof(real);
  WRITE_DATA(cf->bsize, size);
#endif

  /* Write the separation bin edges. */
  size = (cf->ns + 1) * sizeof(real);
  WRITE_DATA(cf->sbin_raw, size);

  if (cf->bintype == FCFC_BIN_SPI) {
    size = (cf->np + 1) * sizeof(real);
    WRITE_DATA(cf->pbin_raw, size);
  }

  /* Write the number of tracers and normalization factor. */
  if (cf->pc_idx[0][idx] == cf->pc_idx[1][idx]) {       /* auto pair counts */
    if (cf->wt[idx]) {
      size = sizeof(double);
      WRITE_DATA(&(cf->data[cf->pc_idx[0][idx]].wt), size);
    }
    else {
      size = sizeof(uint64_t);
      uint64_t num = cf->data[cf->pc_idx[0][idx]].n;
      WRITE_DATA(&num, size);
    }
  }
  else {                                                /* cross pair counts */
    if (cf->wt[idx]) {
      size = 2 * sizeof(double);
      double sumw[2] = {cf->data[cf->pc_idx[0][idx]].wt,
          cf->data[cf->pc_idx[1][idx]].wt};
      WRITE_DATA(sumw, size);
    }
    else {
      size = 2 * sizeof(uint64_t);
      uint64_t num[2] = {(uint64_t) cf->data[cf->pc_idx[0][idx]].n,
          (uint64_t) cf->data[cf->pc_idx[1][idx]].n};
      WRITE_DATA(num, size);
    }
  }
  size = sizeof(double);
  WRITE_DATA(cf->norm + idx, size);

  /* Write the normalized pair counts. */
  size = cf->ntot * sizeof(double);
  WRITE_DATA(cf->ncnt[idx], size);

  /* Write the raw pair counts. */
  if (cf->wt[idx]) {
    size = cf->ntot * sizeof(double);
    if (sizeof(double) == sizeof(COUNT)) {
      WRITE_DATA(cf->cnt[idx], size);
    }
    else {
      double *cnt = malloc(size);
      if (!cnt) {
        P_ERR("failed to allocate memory for saving pair counts\n");
        return FCFC_ERR_MEMORY;
      }
      for (size_t i = 0; i < cf->ntot; i++) cnt[i] = cf->cnt[idx][i].d;
      WRITE_DATA(cnt, size);
      free(cnt);
    }
  }
  else {
    size = cf->ntot * sizeof(int64_t);
    if (sizeof(int64_t) == sizeof(COUNT)) {
      WRITE_DATA(cf->cnt[idx], size);
    }
    else {
      int64_t *cnt = malloc(size);
      if (!cnt) {
        P_ERR("failed to allocate memory for saving pair counts\n");
        return FCFC_ERR_MEMORY;
      }
      for (size_t i = 0; i < cf->ntot; i++) cnt[i] = cf->cnt[idx][i].i;
      WRITE_DATA(cnt, size);
      free(cnt);
    }
  }

  return 0;
}

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
    const FCFC_otype_e flag) {
  /* Find the name of the output file. */
  const char *fname = NULL;
  switch (flag) {
    case FCFC_OUTPUT_PAIR_COUNT:
      if (idx < 0 || idx >= cf->npc) {
        P_ERR("invalid index of the pair counts to be saved: %d\n", idx);
        return FCFC_ERR_ARG;
      }
      fname = conf->pcout[idx];
      break;
    case FCFC_OUTPUT_2PCF_RAW:
      if (idx < 0 || idx >= cf->ncf) {
        P_ERR("invalid index of the correlation function to be saved: %d\n",
            idx);
        return FCFC_ERR_ARG;
      }
      fname = conf->cfout[idx];
      break;
    case FCFC_OUTPUT_2PCF_INTEG:
      if (idx < 0 || idx >= cf->ncf) {
        P_ERR("invalid index of the correlation function to be saved: %d\n",
            idx);
        return FCFC_ERR_ARG;
      }
      if (conf->bintype == FCFC_BIN_SMU) fname = conf->mpout[idx];
      else if (conf->bintype == FCFC_BIN_SPI) fname = conf->wpout[idx];
      else {
        P_ERR("integrated correlation function not available "
            "for binning scheme: %d\n", conf->bintype);
        return FCFC_ERR_ARG;
      }
      break;
    default:
      P_ERR("invalid identifier of the results to be saved: %d\n", flag);
      return FCFC_ERR_ARG;
  }

  /* Initialise the interface for writing files. */
  OFILE *ofile = output_init();
  if (!ofile) return FCFC_ERR_FILE;
  if (output_newfile(ofile, fname)) {
    output_destroy(ofile); return FCFC_ERR_FILE;
  }

  if (flag == FCFC_OUTPUT_PAIR_COUNT && conf->ofmt == FCFC_OFMT_BIN) {
    if (save_res_bin(ofile, conf, cf, idx)) {
      output_destroy(ofile); return FCFC_ERR_SAVE;
    }
  }
  else {
    if (save_res_ascii(ofile, conf, cf, idx, flag)) {
      output_destroy(ofile); return FCFC_ERR_SAVE;
    }
  }

  output_destroy(ofile);
  if (conf->verbose) printf("  Results written to file: `%s'\n", fname);
  return 0;
}

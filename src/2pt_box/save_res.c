/*******************************************************************************
* 2pt_box/save_res.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "write_file.h"
#include "save_res.h"

/* Shortcut for saving an output line. */
#define WRITE_LINE(...)                                                 \
  if (output_writeline(ofile, __VA_ARGS__)) {                           \
    output_destroy(ofile); return FCFC_ERR_FILE;                        \
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
    const fcfc_out_t flag) {
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
            "for binning scheme: %d", conf->bintype);
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

  /* Write the header for pair counts. */
  if (flag == FCFC_OUTPUT_PAIR_COUNT) {
    WRITE_LINE("%c Created by " FCFC_CODE_NAME "\n", FCFC_SAVE_COMMENT);
    if (cf->pc_idx[0][idx] == cf->pc_idx[1][idx]) {     /* auto pairs */
      WRITE_LINE("%c Number of tracers: %zu (auto pair counts)\n",
          FCFC_SAVE_COMMENT, cf->ndata[cf->pc_idx[0][idx]]);
    }
    else {                                              /* cross pairs */
      WRITE_LINE("%c Numbers of tracers: %zu and %zu\n",
          FCFC_SAVE_COMMENT, cf->ndata[cf->pc_idx[0][idx]],
          cf->ndata[cf->pc_idx[1][idx]]);
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
            (cf->sbin[i] + cf->sbin[i + 1]) * 0.5,
            cf->sbin[i], cf->sbin[i + 1]);
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
            (cf->sbin[i] + cf->sbin[i + 1]) * 0.5,
            cf->sbin[i], cf->sbin[i + 1], cf->wp[idx][i]);
      }
    }
  }
  else {        /* flag == FCFC_OUTPUT_PAIR_COUNT or FCFC_OUTPUT_2PCF_RAW */
    double **res = (flag == FCFC_OUTPUT_PAIR_COUNT) ? cf->ncnt : cf->cf;
    if (conf->bintype == FCFC_BIN_SMU) {
      if (conf->ostyle == FCFC_OUT_LIST) {
        /* Header for count(s,mu) or xi(s,mu) as a list. */
        WRITE_LINE("%c Columns: s_min(1) s_max(2) mu_min(3) mu_max(4) %s(5)\n",
            FCFC_SAVE_COMMENT, (flag == FCFC_OUTPUT_PAIR_COUNT) ?
            "normalized_pair_count" : "xi");
        /* Results of count(s,mu) or xi(s,mu) as a list. */
        for (int j = 0; j < cf->nmu; j++) {
          for (int i = 0; i < cf->ns; i++) {
            WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL " "
                OFMT_DBL "\n", cf->sbin[i], cf->sbin[i + 1],
                j / (double) cf->nmu, (j + 1) / (double) cf->nmu,
                res[idx][i + j * cf->ns]);
          }
        }
      }
      else {    /* conf->ostyle == FCFC_OUT_MATRIX */
        /* Header for count(s,mu) or xi(s,mu) as a matrix. */
        WRITE_LINE("%c Rows: %d s bins, edges:\n%c    ",
            FCFC_SAVE_COMMENT, cf->ns, FCFC_SAVE_COMMENT);
        for (int i = 0; i <= cf->ns; i++) {
          WRITE_LINE(" " OFMT_DBL, cf->sbin[i]);
        }
        WRITE_LINE("\n%c Columns: %d mu bins, edges:\n%c    ",
            FCFC_SAVE_COMMENT, cf->nmu, FCFC_SAVE_COMMENT);
        for (int i = 0; i <= cf->nmu; i++) {
          WRITE_LINE(" " OFMT_DBL, i / (double) cf->nmu);
        }
        WRITE_LINE("\n");
        /* Results of count(s,mu) or xi(s,mu) as a matrix. */
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL, res[idx][i]);
          for (int j = 1; j < cf->nmu; j++) {
            WRITE_LINE(" " OFMT_DBL, res[idx][i + j * cf->ns]);
          }
          WRITE_LINE("\n");
        }
      }
    }
    else if (conf->bintype == FCFC_BIN_SPI) {
      if (conf->ostyle == FCFC_OUT_LIST) {
        /* Header for count(s_perp,s_par) or xi(s_perp,s_par) as a list. */
        WRITE_LINE("%c Columns: s_perp_min(1) s_perp_max(2) pi_min(3) pi_max(4)"
            " %s(5)\n", FCFC_SAVE_COMMENT,
            (flag == FCFC_OUTPUT_PAIR_COUNT) ? "normalized_pair_count" : "xi");
        /* Results of count(s_perp,s_par) or xi(s_perp,s_par) as a list. */
        for (int j = 0; j < cf->np; j++) {
          for (int i = 0; i < cf->ns; i++) {
            WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL " "
                OFMT_DBL "\n", cf->sbin[i], cf->sbin[i + 1],
                cf->pbin[j], cf->pbin[j + 1], res[idx][i + j * cf->ns]);
          }
        }
      }
      else {    /* conf->ostyle == FCFC_OUT_MATRIX */
        /* Header for count(s_perp,s_par) or xi(s_perp,s_par) as a matrix. */
        WRITE_LINE("%c Rows: %d s_perp bins, edges:\n%c    ",
            FCFC_SAVE_COMMENT, cf->ns, FCFC_SAVE_COMMENT);
        for (int i = 0; i <= cf->ns; i++) {
          WRITE_LINE(" " OFMT_DBL, cf->sbin[i]);
        }
        WRITE_LINE("\n%c Columns: %d pi bins, edges:\n%c    ",
            FCFC_SAVE_COMMENT, cf->np, FCFC_SAVE_COMMENT);
        for (int i = 0; i <= cf->np; i++) {
          WRITE_LINE(" " OFMT_DBL, cf->pbin[i]);
        }
        WRITE_LINE("\n");
        /* Results of count(s_perp,s_par) or xi(s_perp,s_par) as a matrix. */
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL, res[idx][i]);
          for (int j = 1; j < cf->np; j++) {
            WRITE_LINE(" " OFMT_DBL, res[idx][i + j * cf->ns]);
          }
          WRITE_LINE("\n");
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
              cf->sbin[i], cf->sbin[i + 1], cf->ncnt[idx][i]);
        }
      }
      else {    /* flag == FCFC_OUTPUT_2PCF_RAW */
        /* Header for isotropic 2PCF. */
        WRITE_LINE("%c Columns: s_cen(1) s_min(2) s_max(3) xi_iso(4)\n",
            FCFC_SAVE_COMMENT);
        /* Results of isotropic 2PCF. */
        for (int i = 0; i < cf->ns; i++) {
          WRITE_LINE(OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL "\n",
              (cf->sbin[i] + cf->sbin[i + 1]) * 0.5,
              cf->sbin[i], cf->sbin[i + 1], cf->cf[idx][i]);
        }
      }
    }
  }

  output_destroy(ofile);
  if (conf->verbose) printf("  Results written to file: `%s'\n", fname);
  return 0;
}

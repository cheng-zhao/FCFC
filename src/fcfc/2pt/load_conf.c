/*******************************************************************************
* 2pt/load_conf.c: this file is part of the FCFC program.

* FCFC: Fast Correlation Function Calculator.

* Github repository:
        https://github.com/cheng-zhao/FCFC

* Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "libcfg.h"
#include "libast.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
/* Check existence of configuration parameters. */
#define CHECK_EXIST_PARAM(name, cfg, var)                       \
  if (!cfg_is_set((cfg), (var))) {                              \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return FCFC_ERR_CFG;                                        \
  }
#define CHECK_EXIST_ARRAY(name, cfg, var, num)                  \
  if (!(num = cfg_get_size((cfg), (var)))) {                    \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return FCFC_ERR_CFG;                                        \
  }

/* Check length of array. */
#define CHECK_ARRAY_LENGTH(name, cfg, var, fmt, num, nexp)      \
  if (num < (nexp)) {                                           \
    P_ERR("too few elements of " FMT_KEY(name) "\n");           \
    return FCFC_ERR_CFG;                                        \
  }                                                             \
  if (num > (nexp)) {                                           \
    P_WRN("omitting the following " FMT_KEY(name) ":");         \
    for (int i = (nexp); i < num; i++)                          \
      fprintf(stderr, " " fmt, (var)[i]);                       \
    fprintf(stderr, "\n");                                      \
  }

#define CHECK_STR_ARRAY_LENGTH(name, cfg, var, num, nexp)       \
  if (num < (nexp)) {                                           \
    P_ERR("too few elements of " FMT_KEY(name) "\n");           \
    return FCFC_ERR_CFG;                                        \
  }                                                             \
  if (num > (nexp)) {                                           \
    P_WRN("omitting the following " FMT_KEY(name) ":\n");       \
    for (int i = (nexp); i < num; i++)                          \
      fprintf(stderr, "  %s\n", (var)[i]);                      \
  }

/* Release memory for configuration parameters. */
#define FREE_ARRAY(x)           {if(x) free(x);}
#define FREE_STR_ARRAY(x)       {if(x) {if (*(x)) free(*(x)); free(x);}}

/* Print the warning and error messages. */
#define P_CFG_WRN(cfg)  cfg_pwarn(cfg, stderr, FMT_WARN);
#define P_CFG_ERR(cfg)  {                                       \
  cfg_perror(cfg, stderr, FMT_ERR);                             \
  cfg_destroy(cfg);                                             \
  return NULL;                                                  \
}


/*============================================================================*\
                    Functions called via command line flags
\*============================================================================*/

/******************************************************************************
Function `usage`:
  Print the usage of command line options.
******************************************************************************/
static void usage(void *args) {
  (void) args;
  printf(FCFC_LOGO "\nUsage: " FCFC_CODE_NAME " [OPTION]\n\
Compute the 2-point correlation functions of survey-like catalogs.\n\
  -h, --help\n\
        Display this message and exit\n\
  -V, --version\n\
        Display the version information\n\
  -t, --template\n\
        Print a template configuration file to the standard output and exit\n\
  -c, --conf            " FMT_KEY(CONFIG_FILE) "     String\n\
        Specify the configuration file (default: `%s')\n\
  -i, --input           " FMT_KEY(CATALOG) "         String array\n\
        Specify the input catalogs\n\
  -l, --label           " FMT_KEY(CATALOG_LABEL) "   Character array\n\
        Specify the labels of the input catalogs\n\
  -T, --type            " FMT_KEY(CATALOG_TYPE) "    Integer array\n\
        Type (format) of the input catalogs\n\
      --skip            " FMT_KEY(ASCII_SKIP) "      Long integer array\n\
        Numbers of lines to be skipped for the ASCII format input catalogs\n\
      --comment         " FMT_KEY(ASCII_COMMENT) "   Character array\n\
        Comment symbols for the ASCII format input catalogs\n\
  -f, --formatter       " FMT_KEY(ASCII_FORMATTER) " String array\n\
        Formatters for columns of ASCII format input catalogs\n\
  -x, --position        " FMT_KEY(POSITION) "        String array\n\
        Column indicator or expression for the 3-D positions of the inputs\n\
  -w, --weight          " FMT_KEY(WEIGHT) "          String array\n\
        Column indicator or expression for weights of the inputs\n\
  -s, --select          " FMT_KEY(SELECTION) "       String array\n\
        Expressions for sample selection criteria\n\
      --convert         " FMT_KEY(COORD_CONVERT) "   Boolean array\n\
        Indicate whether to apply coordinate conversion for the inputs\n\
  -d, --omega-m         " FMT_KEY(OMEGA_M) "         Double\n\
        Density parameter of matter at z = 0, for coordinate conversion\n\
      --omega-l         " FMT_KEY(OMEGA_LAMBDA) "    Double\n\
        Density parameter of Lambda at z = 0, for coordinate conversion\n\
      --eos-w           " FMT_KEY(DE_EOS_W) "        Double\n\
        Dark energy equation of state\n\
      --cmvdst-err      " FMT_KEY(CMVDST_ERR) "      Double\n\
        Specify the error allowed for comoving distance integration\n\
      --cmvdst-file     " FMT_KEY(Z_CMVDST_FILE) "   String\n\
        Specify the table for redshift to comoving distance conversion\n\
  -S, --data-struct     " FMT_KEY(DATA_STRUCT) "     Integer\n\
        Specify the data structure for pair counting\n\
  -B, --bin             " FMT_KEY(BINNING_SCHEME) "  Integer\n\
        Specify the binning scheme of the correlation functions\n\
  -p, --pair            " FMT_KEY(PAIR_COUNT) "      String array\n\
        Specify pairs to be counted or read, using the catalog labels\n\
  -P, --pair-output     " FMT_KEY(PAIR_COUNT_FILE) " String array\n\
        Specify the output files for pair counts\n\
  -e, --cf              " FMT_KEY(CF_ESTIMATOR) "    String array\n\
        Expressions for correlation function estimators based on pair counts\n\
  -E, --cf-output       " FMT_KEY(CF_OUTPUT_FILE) "  String array\n\
        Specify the output files for correlation functions\n\
  -m, --multipole       " FMT_KEY(MULTIPOLE) "       Integer array\n\
        Orders of Legendre multipoles of correlation functions to be evaluated\n\
  -M, --mp-output       " FMT_KEY(MULTIPOLE_FILE) "  String array\n\
        Specify the output files for correlation function multipoles\n\
  -u, --wp              " FMT_KEY(PROJECTED_CF) "    Boolean\n\
        Indicate whether to compute the projected correlation functions\n\
  -U, --wp-output       " FMT_KEY(PROJECTED_FILE) "  String array\n\
        Specify the output files for projected correlation functions\n\
      --s-file          " FMT_KEY(SEP_BIN_FILE) "    String\n\
        Specify the file defining edges of separation (or s_perp) bins\n\
      --s-min           " FMT_KEY(SEP_BIN_MIN) "     Double\n\
        Specify the lower limit of linear separation (or s_perp) bins\n\
      --s-max           " FMT_KEY(SEP_BIN_MAX) "     Double\n\
        Specify the upper limit of linear separation (or s_perp) bins\n\
      --s-step          " FMT_KEY(SEP_BIN_SIZE) "    Double\n\
        Specify the width of linear separation (or s_perp) bins\n\
      --mu-num          " FMT_KEY(MU_BIN_NUM) "      Integer\n\
        Specify the number of linear mu bins in the range [0,1"
#ifdef WITH_MU_ONE
      "]"
#else
      ")"
#endif
"\n\
      --pi-file         " FMT_KEY(PI_BIN_FILE) "     String\n\
        Specify the file defining edges of pi (a.k.a. s_para) bins\n\
      --pi-min          " FMT_KEY(PI_BIN_MIN) "      Double\n\
        Specify the lower limit of linear pi bins\n\
      --pi-max          " FMT_KEY(PI_BIN_MAX) "      Double\n\
        Specify the upper limit of linear pi bins\n\
      --pi-step         " FMT_KEY(PI_BIN_SIZE) "     Double\n\
        Specify the width of linear pi bins\n\
  -F, --out-format      " FMT_KEY(OUTPUT_FORMAT) "   Integer\n\
        Format of the output pair count files\n\
  -O, --overwrite       " FMT_KEY(OVERWRITE) "       Integer\n\
        Indicate whether to overwrite existing output files\n\
  -v, --verbose         " FMT_KEY(VERBOSE) "         Boolean\n\
        Indicate whether to display detailed standard outputs\n\
Consult the -t option for more information on the parameters\n\
Github repository: https://github.com/cheng-zhao/FCFC\n\
Licence: MIT\n",
    DEFAULT_CONF_FILE);
  exit(0);
}

/******************************************************************************
Function `version`:
  Print the version information.
******************************************************************************/
static void version(void *args) {
  (void) args;
  printf(FCFC_LOGO "\n\x1B[35C\x1B[33;1mv" FCFC_VERSION "\n\x1B[34C"
      FCFC_CODE_NAME "\x1B[0m\n"
      "\n- Parallelization schemes\n"
      "  * MPI: "
#ifdef MPI
      "enabled\n"
#else
      "disabled (enable with -DMPI)\n"
#endif
      "  * OpenMP: "
#ifdef OMP
      "enabled\n"
#else
      "disabled (enable with -DOMP)\n"
#endif
      "  * SIMD: "
#if             FCFC_SIMD  ==  FCFC_SIMD_NONE
      "disabled (enable with -DWITH_SIMD)"
#elif           FCFC_SIMD  ==  FCFC_SIMD_AVX
      "AVX"
  #ifdef        FCFC_SIMD_FMA
      " + FMA"
  #endif
#elif           FCFC_SIMD  ==  FCFC_SIMD_AVX2
      "AVX2"
  #ifdef        FCFC_SIMD_FMA
      " + FMA"
  #endif
#else        /* FCFC_SIMD  ==  FCFC_SIMD_AVX512 */
      "AVX-512F"
  #ifdef        FCFC_SIMD_AVX512DQ
      " + AVX-512DQ"
  #endif
#endif
      "\n- Compilation options\n"
      "  * Floating-point precision: "
#ifdef SINGLE_PREC
      "single (-DSINGLE_PREC enabled)\n"
#else
      "double (-DSINGLE_PREC disabled)\n"
#endif
      "  * (s,mu) pairs with mu = 1: "
#ifdef WITH_MU_ONE
      "included (-DWITH_MU_ONE enabled)\n"
#else
      "excluded (-DWITH_MU_ONE disabled)\n"
#endif
      "- External libraries\n"
      "  * CFITSIO: "
#ifdef WITH_CFITSIO
      "enabled\n"
#else
      "disabled (enable with -DWITH_CFITSIO)\n"
#endif
      "  * HDF5: "
#ifdef WITH_HDF5
      "enabled\n"
#else
      "disabled (enable with -DWITH_HDF5)\n"
#endif
      "\n\
- Copyright (c) 2020 -- 2022 Cheng ZHAO.\n\
- Github repository: https://github.com/cheng-zhao/FCFC\n\
- Licence: MIT\n");
  exit(0);
}

/******************************************************************************
Function `conf_template`:
  Print a template configuration file.
******************************************************************************/
static void conf_template(void *args) {
  (void) args;
  printf("# Configuration file for " FCFC_CODE_NAME " (default: `"
DEFAULT_CONF_FILE "').\n\
# Format: keyword = value # comment\n\
#     or: keyword = [element1, element2]\n\
#    see: https://github.com/cheng-zhao/libcfg for details.\n\
# Some of the entries allow expressions, see\n\
#         https://github.com/cheng-zhao/libast for details.\n\
# NOTE that command line options have priority over this file.\n\
# Unnecessary entries can be left unset.\n\
\n\
##########################################\n\
#  Specifications of the input catalogs  #\n\
##########################################\n\
\n\
CATALOG         = \n\
    # Filename of the input catalogs, string or string array.\n\
CATALOG_LABEL   = \n\
    # Label of the input catalogs, must be non-repetitive uppercase letters.\n\
    # Character, same dimension as `CATALOG`.\n\
    # If unset, catalogs are labelled in alphabetical order, i.e. [A,B,...].\n\
CATALOG_TYPE    = \n\
    # File format of the input catalogs (unset: %d).\n\
    # Integer, same dimension as `CATALOG`.\n\
    # Allowed values are:\n\
    # * %d: ASCII text file"
#ifdef WITH_CFITSIO
    ";\n    # * %d: FITS table"
#endif
#ifdef WITH_HDF5
    ";\n    # * %d: HDF5 file"
#endif
".\nASCII_SKIP      = \n\
    # Number of lines to be skipped for ASCII catalogs (unset: %ld).\n\
    # Long integer, same dimension as `CATALOG`.\n\
ASCII_COMMENT   = \n\
    # Character indicating comment lines for ASCII catalogs (unset: '%c%s.\n\
    # Character, same dimension as `CATALOG`.\n\
    # Empty character ('') means disabling comments.\n\
ASCII_FORMATTER = \n\
    # C99-style formatter for parsing lines of ASCII catalogs.\n\
    # String, same dimension as `DATA_CATALOG` (e.g. \"%%d %%ld %%f %%lf %%s\").\n\
    # If a column is suppressed by *, it is not counted for the column number.\n\
    # E.g., for \"%%d %%*s %%f\", the float number corresponds to column %c2.\n\
    # See https://en.cppreference.com/w/c/io/fscanf for details on the format.\n\
POSITION        = \n\
    # 3-D comoving coordinates, in the order of {RA,Dec,redshift} or {x,y,z}.\n\
    # String array, 3 times the length of `CATALOG`.\n\
    # They can be column indicator or expressions, e.g.,\n\
    #     \"(%c1 * %c%c10%c) %% 100\" / \"%c%cRA%c\" / \"%c%cgroup/dataset%c2%c%c\"\n\
    # Allowed values enclosed by %c%c%c:\n\
    # * long integer: column number of an ASCII file (starting from 1);\n\
    # * string: column name of a FITS file;\n\
    # * string%cinteger%c: dataset name and column index (starting from 1)\n\
    #                    of an HDF5 file.\n\
WEIGHT          = \n\
    # Weights for pair counts (unset: 1, i.e. no weight).\n\
    # Column indicator or expression, same dimension as `DATA_CATALOG`.\n\
SELECTION       = \n\
    # Selection criteria for the catalogs (unset: no selection).\n\
    # Logical expression, same dimension as `CATALOG` (e.g. \"%c3 > 0.5\").\n\
COORD_CONVERT   = \n\
    # Boolean option, same dimension as `CATALOG` (unset: %c).\n\
    # If a scalar value is given, it is applied to all catalogs.\n\
    # True (T) for converting the coordinates from {RA,Dec,redshift} to\n\
    # the comoving {x,y,z}, given the fiducial cosmology set below.\n\
\n\
##################################################\n\
#  Fiducial cosmology for coordinate conversion  #\n\
##################################################\n\
\n\
OMEGA_M         = \n\
    # Density parameter of matter at z = 0.\n\
    # Double-precision number.\n\
OMEGA_LAMBDA    = \n\
    # Density parameter of Lambda at z = 0 (unset: 1 - `OMEGA_M`).\n\
    # Double-precision number.\n\
DE_EOS_W        = \n\
    # Dark energy equation of state: w (unset: %g).\n\
    # Double-precision number.\n\
CMVDST_ERR      = \n\
    # Error for comoving distance evaluation (unset: %g).\n\
    # Double-precision number.\n\
Z_CMVDST_FILE   = \n\
    # Filename of a table for redshift to comoving distance conversion.\n\
    # It must be a text file with two columns: {redshift, comoving distance}.\n\
    # If this file is set, the cosmological parameters above are omitted.\n\
    # Lines starting with '%c' are omitted.\n\
\n\
################################################################\n\
#  Configurations for the 2-point correlation function (2PCF)  #\n\
################################################################\n\
\n\
DATA_STRUCT     = \n\
    # Data structure for evaluating pair counts, integer (unset: %d).\n\
    # Allowed values are:\n\
    # * %d: k-d tree;\n\
    # * %d: ball tree.\n\
BINNING_SCHEME  = \n\
    # Binning scheme of the 2PCFs, integer (unset: %d).\n\
    # Allowed values are:\n\
    # * %d: isotropic separation bins;\n\
    # * %d: (s, mu) bins (required by 2PCF multipoles);\n\
    # * %d: (s_perp, pi) bins (required by projected 2PCFs);\n\
PAIR_COUNT      = \n\
    # Identifiers of pairs to be counted or read, string or string array.\n\
    # Pairs are labelled by their source catalogs.\n\
    # E.g., \"DD\" denotes auto pairs from the catalog 'D',\n\
    # while \"DR\" denotes cross pairs from catalogs 'D' and 'R'.\n\
PAIR_COUNT_FILE = \n\
    # Name of the files for storing pair counts.\n\
    # String, same dimension as `PAIR_COUNT`.\n\
    # Depending on `OVERWRITE`, pair counts can be read from existing files.\n\
CF_ESTIMATOR    = \n\
    # Estimator of the 2PCFs to be evaluated, string or string array.\n\
    # It must be an expression with pair identifiers.\n\
CF_OUTPUT_FILE  = \n\
    # Name of the files for saving 2PCFs with the desired binning scheme.\n\
    # String, same dimension as `CF_ESTIMATOR`.\n\
MULTIPOLE       = \n\
    # Orders of Legendre multipoles to be evaluated, integer or integer array.\n\
MULTIPOLE_FILE  = \n\
    # Name of the files for saving 2PCF multipoles.\n\
    # String, same dimension as `CF_ESTIMATOR`.\n\
PROJECTED_CF    = \n\
    # Boolean option, indicate whether computing the projected 2PCFs (unset: %c).\n\
PROJECTED_FILE  = \n\
    # Name of the files for saving projected 2PCFs.\n\
    # String, same dimension as `CF_ESTIMATOR`.\n\
\n\
#############################\n\
#  Definitions of the bins  #\n\
#############################\n\
\n\
SEP_BIN_FILE    = \n\
    # Filename of the table defining edges of separation (or s_perp) bins.\n\
    # It mush be a text file with two columns, for the lower and upper limits\n\
    # of the distance bins, respectively.\n\
    # Lines starting with '%c' are omitted.\n\
SEP_BIN_MIN     = \n\
SEP_BIN_MAX     = \n\
SEP_BIN_SIZE    = \n\
    # Lower and upper limits, and width of linear separation (or s_perp) bins.\n\
    # Double-precision numbers. They are only used if `SEP_BIN_FILE` is unset.\n\
MU_BIN_NUM      = \n\
    # Number of linear mu bins in the range [0,1"
#ifdef WITH_MU_ONE
    "]"
#else
    ")"
#endif
", integer.\n\
PI_BIN_FILE     = \n\
    # Filename of the table defining edges of pi (a.k.a. s_para) bins.\n\
    # Lines starting with '%c' are omitted.\n\
PI_BIN_MIN      = \n\
PI_BIN_MAX      = \n\
PI_BIN_SIZE     = \n\
    # Lower and upper limits, and width of linear pi bins.\n\
    # Double-precision numbers. They are only used if `PI_BIN_FILE` is unset.\n\
\n\
####################\n\
#  Other settings  #\n\
####################\n\
\n\
OUTPUT_FORMAT   = \n\
    # Format of the output `PAIR_COUNT_FILE`, integer (unset: %d).\n\
    # Allowed values are:\n\
    # * %d: FCFC binary format;\n\
    # * %d: ASCII text format.\n\
OVERWRITE       = \n\
    # Flag indicating whether to overwrite existing files, integer (unset: %d).\n\
    # Allowed values are:\n\
    # * %d: quit the program when an output file exist;\n\
    # * %d: overwrite 2PCF files silently, but keep existing pair count files;\n\
    # * %d or larger: overwrite all files silently;\n\
    # * negative: notify for decisions, and the maximum allowed number of failed\n\
    #             trials are given by the absolute value of this number.\n\
VERBOSE         = \n\
    # Boolean option, indicate whether to show detailed outputs (unset: %c).\n",
      DEFAULT_FILE_TYPE, FCFC_FFMT_ASCII,
#ifdef WITH_CFITSIO
      FCFC_FFMT_FITS,
#endif
#ifdef WITH_HDF5
      FCFC_FFMT_HDF5,
#endif
      (long) DEFAULT_ASCII_SKIP,
      DEFAULT_ASCII_COMMENT ? DEFAULT_ASCII_COMMENT : '\'',
      DEFAULT_ASCII_COMMENT ? "')" : ")",
      AST_VAR_FLAG, AST_VAR_FLAG, AST_VAR_FLAG, AST_VAR_START, AST_VAR_END,
      AST_VAR_FLAG, AST_VAR_START, AST_VAR_END,
      AST_VAR_FLAG, AST_VAR_START, FCFC_COL_IDX_START, FCFC_COL_IDX_END,
      AST_VAR_END, AST_VAR_FLAG, AST_VAR_START, AST_VAR_END,
      FCFC_COL_IDX_START, FCFC_COL_IDX_END, AST_VAR_FLAG,
      DEFAULT_COORD_CNVT ? 'T' : 'F',
      (double) DEFAULT_DE_EOS_W, DEFAULT_CNVT_ERR, FCFC_READ_COMMENT,
      DEFAULT_STRUCT, FCFC_STRUCT_KDTREE, FCFC_STRUCT_BALLTREE,
      DEFAULT_BINNING, FCFC_BIN_ISO, FCFC_BIN_SMU, FCFC_BIN_SPI,
      DEFAULT_PROJECTED_CF ? 'T' : 'F', FCFC_READ_COMMENT, FCFC_READ_COMMENT,
      DEFAULT_OUTPUT_FORMAT, FCFC_OFMT_BIN, FCFC_OFMT_ASCII,
      DEFAULT_OVERWRITE, FCFC_OVERWRITE_NONE, FCFC_OVERWRITE_CFONLY,
      FCFC_OVERWRITE_ALL, DEFAULT_VERBOSE ? 'T' : 'F');
  exit(0);
}


/*============================================================================*\
                      Function for reading configurations
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing configurations.
Return:
  Address of the structure.
******************************************************************************/
static CONF *conf_init(void) {
  CONF *conf = calloc(1, sizeof *conf);
  if (!conf) return NULL;
  conf->fconf = conf->label = conf->comment = NULL;
  conf->fcnvt = conf->fsbin = conf->fpbin = NULL;
  conf->input = conf->fmtr = conf->pos = conf->wt = conf->sel = NULL;
  conf->pc = conf->pcout = conf->cf = conf->cfout = NULL;
  conf->mpout = conf->wpout = NULL;
  conf->ftype = conf->poles = NULL;
  conf->has_wt = conf->cnvt = NULL;
  conf->comp_pc = NULL;
  conf->skip = NULL;
  return conf;
}

/******************************************************************************
Function `conf_read`:
  Read configurations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  Interface of libcfg.
******************************************************************************/
static cfg_t *conf_read(CONF *conf, const int argc, char *const *argv) {
  if (!conf) {
    P_ERR("the structure for configurations is not initialised\n");
    return NULL;
  }
  cfg_t *cfg = cfg_init();
  if (!cfg) P_CFG_ERR(cfg);

  /* Functions to be called via command line flags. */
  const cfg_func_t funcs[] = {
    {'h', "help"        , usage         ,               NULL},
    {'V', "version"     , version       ,               NULL},
    {'t', "template"    , conf_template ,               NULL}
  };

  /* Configuration parameters. */
  const cfg_param_t params[] = {
    {'c', "conf"        , "CONFIG_FILE"    , CFG_DTYPE_STR , &conf->fconf   },
    {'i', "input"       , "CATALOG"        , CFG_ARRAY_STR , &conf->input   },
    {'l', "label"       , "CATALOG_LABEL"  , CFG_ARRAY_CHAR, &conf->label   },
    {'T', "type"        , "CATALOG_TYPE"   , CFG_ARRAY_INT , &conf->ftype   },
    { 0 , "skip"        , "ASCII_SKIP"     , CFG_ARRAY_LONG, &conf->skip    },
    { 0 , "comment"     , "ASCII_COMMENT"  , CFG_ARRAY_CHAR, &conf->comment },
    {'f', "formatter"   , "ASCII_FORMATTER", CFG_ARRAY_STR , &conf->fmtr    },
    {'x', "position"    , "POSITION"       , CFG_ARRAY_STR , &conf->pos     },
    {'w', "weight"      , "WEIGHT"         , CFG_ARRAY_STR , &conf->wt      },
    {'s', "select"      , "SELECTION"      , CFG_ARRAY_STR , &conf->sel     },
    { 0 , "convert"     , "COORD_CONVERT"  , CFG_ARRAY_BOOL, &conf->cnvt    },
    {'d', "omega-m"     , "OMEGA_M"        , CFG_DTYPE_DBL , &conf->omega_m },
    { 0 , "omega-l"     , "OMEGA_LAMBDA"   , CFG_DTYPE_DBL , &conf->omega_l },
    { 0 , "eos-w"       , "DE_EOS_W"       , CFG_DTYPE_DBL , &conf->dew     },
    { 0 , "cmvdst-err"  , "CMVDST_ERR"     , CFG_DTYPE_DBL , &conf->ecnvt   },
    { 0 , "cmvdst-file" , "Z_CMVDST_FILE"  , CFG_DTYPE_STR , &conf->fcnvt   },
    {'S', "data-struct" , "DATA_STRUCT"    , CFG_DTYPE_INT , &conf->dstruct },
    {'B', "bin"         , "BINNING_SCHEME" , CFG_DTYPE_INT , &conf->bintype },
    {'p', "pair"        , "PAIR_COUNT"     , CFG_ARRAY_STR , &conf->pc      },
    {'P', "pair-output" , "PAIR_COUNT_FILE", CFG_ARRAY_STR , &conf->pcout   },
    {'e', "cf"          , "CF_ESTIMATOR"   , CFG_ARRAY_STR , &conf->cf      },
    {'E', "cf-output"   , "CF_OUTPUT_FILE" , CFG_ARRAY_STR , &conf->cfout   },
    {'m', "multipole"   , "MULTIPOLE"      , CFG_ARRAY_INT , &conf->poles   },
    {'M', "mp-output"   , "MULTIPOLE_FILE" , CFG_ARRAY_STR , &conf->mpout   },
    {'u', "wp"          , "PROJECTED_CF"   , CFG_DTYPE_BOOL, &conf->wp      },
    {'U', "wp-output"   , "PROJECTED_FILE" , CFG_ARRAY_STR , &conf->wpout   },
    { 0 , "s-file"      , "SEP_BIN_FILE"   , CFG_DTYPE_STR , &conf->fsbin   },
    { 0 , "s-min"       , "SEP_BIN_MIN"    , CFG_DTYPE_DBL , &conf->smin    },
    { 0 , "s-max"       , "SEP_BIN_MAX"    , CFG_DTYPE_DBL , &conf->smax    },
    { 0 , "s-step"      , "SEP_BIN_SIZE"   , CFG_DTYPE_DBL , &conf->ds      },
    { 0 , "mu-num"      , "MU_BIN_NUM"     , CFG_DTYPE_INT , &conf->nmu     },
    { 0 , "pi-file"     , "PI_BIN_FILE"    , CFG_DTYPE_STR , &conf->fpbin   },
    { 0 , "pi-min"      , "PI_BIN_MIN"     , CFG_DTYPE_DBL , &conf->pmin    },
    { 0 , "pi-max"      , "PI_BIN_MAX"     , CFG_DTYPE_DBL , &conf->pmax    },
    { 0 , "pi-step"     , "PI_BIN_SIZE"    , CFG_DTYPE_DBL , &conf->dpi     },
    {'F', "out-format"  , "OUTPUT_FORMAT"  , CFG_DTYPE_INT , &conf->ofmt    },
    {'O', "overwrite"   , "OVERWRITE"      , CFG_DTYPE_INT , &conf->ovwrite },
    {'v', "verbose"     , "VERBOSE"        , CFG_DTYPE_BOOL, &conf->verbose }
  };

  /* Register functions and parameters. */
  if (cfg_set_funcs(cfg, funcs, sizeof(funcs) / sizeof(funcs[0])))
      P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);
  if (cfg_set_params(cfg, params, sizeof(params) / sizeof(params[0])))
      P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read configurations from command line options. */
  int optidx;
  if (cfg_read_opts(cfg, argc, argv, FCFC_PRIOR_CMD, &optidx))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read parameters from configuration file. */
  if (!cfg_is_set(cfg, &conf->fconf)) conf->fconf = DEFAULT_CONF_FILE;
  if (access(conf->fconf, R_OK))
    P_WRN("cannot access the configuration file: `%s'\n", conf->fconf);
  else if (cfg_read_file(cfg, conf->fconf, FCFC_PRIOR_FILE)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  return cfg;
}


/*============================================================================*\
                      Functions for parameter verification
\*============================================================================*/

/******************************************************************************
Function `check_input`:
  Check whether an input file can be read.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_input(const char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR("the input " FMT_KEY(%s) " is not set\n", key);
    return FCFC_ERR_CFG;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot access " FMT_KEY(%s) ": `%s'\n", key, fname);
    return FCFC_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `check_output`:
  Check whether an output file can be written.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file;
  * `ovwrite`:  option for overwriting exisiting files;
  * `force`:    flag for overwriting files without notification.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_output(char *fname, const char *key, int ovwrite,
    const int force) {
  if (!fname || *fname == '\0') {
    P_ERR("the output file " FMT_KEY(%s) " is not set\n", key);
    return FCFC_ERR_CFG;
  }

  /* Check if the file exists. */
  if (!access(fname, F_OK)) {
    if (ovwrite < 0) {                          /* ask for decision */
      P_WRN("the output file " FMT_KEY(%s) " exists: `%s'\n", key, fname);
      char confirm = 0;
      for (int i = 0; i != ovwrite; i--) {
        fprintf(stderr, "Are you going to overwrite it? (y/n): ");
        if (scanf("%c", &confirm) != 1) continue;
        int c;
        while((c = getchar()) != '\n' && c != EOF) continue;
        if (confirm == 'n' || confirm == 'N') {
          ovwrite = force - 1;
          break;
        }
        else if (confirm == 'y' || confirm == 'Y') {
          ovwrite = force;
          break;
        }
      }
      if (confirm != 'y' && confirm != 'Y' &&
          confirm != 'n' && confirm != 'N') {
        P_ERR("too many failed inputs\n");
        return FCFC_ERR_FILE;
      }
    }

    if (ovwrite <= FCFC_OVERWRITE_NONE) {       /* not overwriting */
      P_ERR("abort to avoid overwriting " FMT_KEY(%s) ": `%s'\n", key, fname);
      return FCFC_ERR_FILE;
    }
    else if (ovwrite >= force) {                /* force overwriting */
      P_WRN(FMT_KEY(%s) " will be overwritten: `%s'\n", key, fname);
    }
    else {                                      /* this is an input file */
      if (access(fname, R_OK)) {
        P_ERR("cannot access " FMT_KEY(%s) ": `%s'\n", key, fname);
        return FCFC_ERR_FILE;
      }
      return FCFC_ERR_SAVE;     /* indicate that the file will be read */
    }

    /* Check file permission for overwriting. */
    if (access(fname, W_OK)) {
      P_ERR("cannot write to file: `%s'\n", fname);
      return FCFC_ERR_FILE;
    }
  }
  /* Check the path permission. */
  else {
    char *end;
    if ((end = strrchr(fname, FCFC_PATH_SEP)) != NULL) {
      *end = '\0';
      if (access(fname, X_OK)) {
        P_ERR("cannot access the directory `%s'\n", fname);
        return FCFC_ERR_FILE;
      }
      *end = FCFC_PATH_SEP;
    }
  }
  return 0;
}

/******************************************************************************
Function `check_cosmo`:
  Verify cosmological parameters for coordinate conversion.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_cosmo(const cfg_t *cfg, CONF *conf) {
  /* Check OMEGA_M. */
  CHECK_EXIST_PARAM(OMEGA_M, cfg, &conf->omega_m);
  if (conf->omega_m <= 0 || conf->omega_m > 1) {
    P_ERR(FMT_KEY(OMEGA_M) " must be > 0 and <= 1\n");
    return FCFC_ERR_CFG;
  }

  /* Check OMEGA_LAMBDA */
  if (!cfg_is_set(cfg, &conf->omega_l)) {
    conf->omega_l = 1 - conf->omega_m;
    conf->omega_k = 0;
  }
  else if (conf->omega_l < 0) {
    P_ERR(FMT_KEY(OMEGA_LAMBDA) " must be >= 0\n");
    return FCFC_ERR_CFG;
  }
  else conf->omega_k = 1 - conf->omega_m - conf->omega_l;

  /* Check DE_EOS_W. */
  if (!cfg_is_set(cfg, &conf->dew)) conf->dew = DEFAULT_DE_EOS_W;
  else if (conf->dew > -1 / (double) 3) {
    P_ERR(FMT_KEY(DE_EOS_W) " must be <= -1/3\n");
    return FCFC_ERR_CFG;
  }

  /* Finally, make sure that H^2 (z) > 0. */
  double w3 = conf->dew * 3;
  double widx = w3 + 1;
  if (conf->omega_k * pow(conf->omega_l * (-widx), widx / w3) <=
      conf->omega_l * w3 * pow(conf->omega_m, widx / w3)) {
    P_ERR("negative H^2 given the cosmological parameters\n");
    return FCFC_ERR_CFG;
  }
  return 0;
}

/******************************************************************************
Function `conf_verify`:
  Verify configuration parameters.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int conf_verify(const cfg_t *cfg, CONF *conf) {
  int e, num;

  /* CATALOG */
  CHECK_EXIST_ARRAY(CATALOG, cfg, &conf->input, conf->ninput);
  if (conf->ninput > 26) {
    P_ERR("there can be at most 26 catalogs set via " FMT_KEY(CATALOG) "\n");
    return FCFC_ERR_CFG;
  }
  for (int i = 0; i < conf->ninput; i++) {
    if ((e = check_input(conf->input[i], "CATALOG"))) return e;
  }

  /* CATALOG_LABEL */
  num = cfg_get_size(cfg, &conf->label);
  if (!num) {
    char *tmp = realloc(conf->label, conf->ninput * sizeof *tmp);
    if (!tmp) {
      P_ERR("failed to allocate memory for " FMT_KEY(CATALOG_LABEL) "\n");
      return FCFC_ERR_MEMORY;
    }
    conf->label = tmp;
    for (int i = 0; i < conf->ninput; i++) conf->label[i] = 'A' + i;
  }
  else {
    CHECK_ARRAY_LENGTH(CATALOG_LABEL, cfg, conf->label, "%c", num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->label[i] < 'A' || conf->label[i] > 'Z') {
        P_ERR("invalid " FMT_KEY(CATALOG_LABEL) ": %c\n", conf->label[i]);
        return FCFC_ERR_CFG;
      }
    }
    /* Check duplicates. */
    for (int i = 0; i < conf->ninput - 1; i++) {
      for (int j = i + 1; j < conf->ninput; j++) {
        if (conf->label[i] == conf->label[j]) {
          P_ERR("duplicate " FMT_KEY(CATALOG_LABEL) ": %c\n", conf->label[i]);
          return FCFC_ERR_CFG;
        }
      }
    }
  }

  /* CATALOG_TYPE */
  if ((num = cfg_get_size(cfg, &conf->ftype))) {
    CHECK_ARRAY_LENGTH(CATALOG_TYPE, cfg, conf->ftype, "%d", num, conf->ninput);
  }
  conf->ascii = false;
  for (int i = 0; i < conf->ninput; i++) {
    int type = conf->ftype ? conf->ftype[i] : DEFAULT_FILE_TYPE;
    switch (type) {
      case FCFC_FFMT_ASCII:
        conf->ascii = true;
        break;
      case FCFC_FFMT_FITS:
#ifdef WITH_CFITSIO
        break;
#else
        P_ERR("FITS format is not enabled\n"
            "Please re-compile the code with option -DWITH_CFITSIO\n");
        return FCFC_ERR_CFG;
#endif
      case FCFC_FFMT_HDF5:
#ifdef WITH_HDF5
        break;
#else
        P_ERR("HDF5 format is not enabled\n"
            "Please re-compile the code with opetion -DWITH_HDF5\n");
        return FCFC_ERR_CFG;
#endif
      default:
        P_ERR("invalid " FMT_KEY(CATALOG_TYPE) ": %d\n", type);
        return FCFC_ERR_CFG;
    }
  }

  /* Check format settings for ASCII catalogues. */
  if (conf->ascii) {
    /* ASCII_SKIP */
    if ((num = cfg_get_size(cfg, &conf->skip))) {
      CHECK_ARRAY_LENGTH(ASCII_SKIP, cfg, conf->skip, "%ld", num, conf->ninput);
      for (int i = 0; i < conf->ninput; i++) {
        int type = conf->ftype ? conf->ftype[i] : DEFAULT_FILE_TYPE;
        if (type == FCFC_FFMT_ASCII && conf->skip[i] < 0) {
          P_ERR(FMT_KEY(ASCII_SKIP) " cannot be negative\n");
          return FCFC_ERR_CFG;
        }
      }
    }

    /* ASCII_COMMENT */
    if ((num = cfg_get_size(cfg, &conf->comment))) {
      CHECK_ARRAY_LENGTH(ASCII_COMMENT, cfg, conf->comment, "%c",
          num, conf->ninput);
      for (int i = 0; i < conf->ninput; i++) {
        int type = conf->ftype ? conf->ftype[i] : DEFAULT_FILE_TYPE;
        if (type == FCFC_FFMT_ASCII) {
          if (conf->comment[i] && !isgraph(conf->comment[i])) {
            P_ERR("invalid " FMT_KEY(ASCII_COMMENT) ": '%c' (ASCII code: %d)\n",
                conf->comment[i], conf->comment[i]);
            return FCFC_ERR_CFG;
          }
        }
      }
    }

    /* ASCII_FORMATTER */
    CHECK_EXIST_ARRAY(ASCII_FORMATTER, cfg, &conf->fmtr, num);
    CHECK_STR_ARRAY_LENGTH(ASCII_FORMATTER, cfg, conf->fmtr, num, conf->ninput);
  }

  /* POSITION */
  CHECK_EXIST_ARRAY(POSITION, cfg, &conf->pos, num);
  CHECK_STR_ARRAY_LENGTH(POSITION, cfg, conf->pos, num, conf->ninput * 3);

  /* WEIGHT */
  if (!(conf->has_wt = malloc(conf->ninput * sizeof(bool)))) {
    P_ERR("failed to allocate memory for " FMT_KEY(WEIGHT) "\n");
    return FCFC_ERR_MEMORY;
  }
  if ((num = cfg_get_size(cfg, &conf->wt))) {
    CHECK_STR_ARRAY_LENGTH(WEIGHT, cfg, conf->wt, num, conf->ninput);
    for (int i = 0; i < conf->ninput; i++) {
      /* Disable weighting if the weight is 1 or empty. */
      if ((conf->wt[i][0] == '1' && conf->wt[i][1] == '\0') ||
          (((conf->wt[i][0] == '\'' && conf->wt[i][1] == '\'') ||
          (conf->wt[i][0] == '"' && conf->wt[i][1] == '"')) &&
          conf->wt[i][2] == '\0')) conf->has_wt[i] = false;
      else conf->has_wt[i] = true;
    }
  }
  else {
    for (int i = 0; i < conf->ninput; i++) conf->has_wt[i] = false;
  }

  /* SELECTION */
  if ((num = cfg_get_size(cfg, &conf->sel))) {
    CHECK_STR_ARRAY_LENGTH(SELECTION, cfg, conf->sel, num, conf->ninput);
  }

  /* COORD_CONVERT */
  if ((num = cfg_get_size(cfg, &conf->cnvt))) {
    if (num == 1 && conf->ninput > 1) {
      bool *tmp = realloc(conf->cnvt, conf->ninput * sizeof(bool));
      if (!tmp) {
        P_ERR("failed to allocate memory for " FMT_KEY(COORD_CONVERT) "\n");
        return FCFC_ERR_MEMORY;
      }
      conf->cnvt = tmp;
      for (int i = 1; i < conf->ninput; i++) conf->cnvt[i] = conf->cnvt[0];
    }
    else if (num < conf->ninput) {
      P_ERR("too few elements of " FMT_KEY(COORD_CONVERT) "\n");
      return FCFC_ERR_CFG;
    }
    if (num > conf->ninput) {
      P_WRN("omitting the following " FMT_KEY(COORD_CONVERT) ":");
      for (int i = conf->ninput; i < num; i++)
        fprintf(stderr, " %c", conf->cnvt[i] ? 'T' : 'F');
      fprintf(stderr, "\n");
    }
  }

  /* Check the fiducial cosmology. */
  conf->has_cnvt = false;
  if (!conf->cnvt && DEFAULT_COORD_CNVT == true) conf->has_cnvt = true;
  else if (conf->cnvt) {
    for (int i = 0; i < conf->ninput; i++) {
      if (conf->cnvt[i]) {
        conf->has_cnvt = true;
        break;
      }
    }
  }
  if (conf->has_cnvt) {
    /* Check Z_CMVDST_FILE. */
    if (cfg_is_set(cfg, &conf->fcnvt)) {
      if ((e = check_input(conf->fcnvt, "Z_CMVDST_FILE"))) return e;
    }
    else {
      if ((e = check_cosmo(cfg, conf))) return e;
      /* Check CMVDST_ERR. */
      if (!cfg_is_set(cfg, &conf->ecnvt)) conf->ecnvt = DEFAULT_CNVT_ERR;
      if (conf->ecnvt < DBL_EPSILON) {
        P_ERR(FMT_KEY(CMVDST_ERR) " is smaller than the machine epsilon.\n");
        return FCFC_ERR_CFG;
      }
    }
  }

  /* OVERWRITE */
  if (!cfg_is_set(cfg, &conf->ovwrite)) conf->ovwrite = DEFAULT_OVERWRITE;

  /* DATA_STRUCT */
  if (!cfg_is_set(cfg, &conf->dstruct)) conf->dstruct = DEFAULT_STRUCT;
  switch (conf->dstruct) {
    case FCFC_STRUCT_KDTREE:
    case FCFC_STRUCT_BALLTREE:
      break;
    default:
      P_ERR("invalid " FMT_KEY(DATA_STRUCT) ": %d\n", conf->dstruct);
      return FCFC_ERR_CFG;
  }

  /* BINNING_SCHEME */
  if (!cfg_is_set(cfg, &conf->bintype)) conf->bintype = DEFAULT_BINNING;
  switch (conf->bintype) {
    case FCFC_BIN_ISO:
    case FCFC_BIN_SMU:
    case FCFC_BIN_SPI:
      break;
    default:
      P_ERR("invalid " FMT_KEY(BINNING_SCHEME) ": %d\n", conf->bintype);
      return FCFC_ERR_CFG;
  }

  /* PAIR_COUNT */
  CHECK_EXIST_ARRAY(PAIR_COUNT, cfg, &conf->pc, conf->npc);
  /* Simple validation. */
  for (int i = 0; i < conf->npc; i++) {
    char *s = conf->pc[i];
    if (s[0] < 'A' || s[0] > 'Z' || s[1] < 'A' || s[1] > 'Z' || s[2]) {
      P_ERR("invalid " FMT_KEY(PAIR_COUNT) ": %s\n", s);
      return FCFC_ERR_CFG;
    }
  }
  /* Check duplicates. */
  for (int i = 0; i < conf->npc - 1; i++) {
    for (int j = i + 1; j < conf->npc; j++) {
      if (conf->pc[i][0] == conf->pc[j][0] &&
          conf->pc[i][1] == conf->pc[j][1]) {
        P_ERR("duplicate " FMT_KEY(PAIR_COUNT) ": %s\n", conf->pc[i]);
        return FCFC_ERR_CFG;
      }
    }
  }

  if (!(conf->comp_pc = malloc(conf->npc * sizeof(bool)))) {
    P_ERR("failed to allocate memory for checking pair counts\n");
    return FCFC_ERR_MEMORY;
  }

  /* PAIR_COUNT_FILE */
  CHECK_EXIST_ARRAY(PAIR_COUNT_FILE, cfg, &conf->pcout, num);
  CHECK_STR_ARRAY_LENGTH(PAIR_COUNT_FILE, cfg, conf->pcout, num, conf->npc);
  for (int i = 0; i < conf->npc; i++) {
    e = check_output(conf->pcout[i], "PAIR_COUNT_FILE", conf->ovwrite,
        FCFC_OVERWRITE_ALL);
    if (!e) conf->comp_pc[i] = true;
    else if (e == FCFC_ERR_SAVE) conf->comp_pc[i] = false;
    else return e;

    /* Check if the labels exist if evaluating pair counts. */
    if (conf->comp_pc[i]) {
      int label_found = 0;
      for (int j = 0; j < conf->ninput; j++) {
        if (conf->pc[i][0] == conf->label[j]) label_found += 1;
        if (conf->pc[i][1] == conf->label[j]) label_found += 1;
      }
      if (label_found != 2) {
        P_ERR("catalog label not found for " FMT_KEY(PAIR_COUNT) ": %s\n",
            conf->pc[i]);
        return FCFC_ERR_CFG;
      }
    }
  }

  /* CF_ESTIMATOR */
  if ((conf->ncf = cfg_get_size(cfg, &conf->cf))) {
    for (int i = 0; i < conf->ncf; i++) {
      if (!conf->cf[i] || !(*conf->cf[i])) {
        P_ERR("unexpected empty " FMT_KEY(CF_ESTIMATOR) "\n");
        return FCFC_ERR_CFG;
      }
    }
    /* CF_OUTPUT_FILE */
    CHECK_EXIST_ARRAY(CF_OUTPUT_FILE, cfg, &conf->cfout, num);
    CHECK_STR_ARRAY_LENGTH(CF_OUTPUT_FILE, cfg, conf->cfout, num, conf->ncf);
    for (int i = 0; i < conf->ncf; i++) {
      if ((e = check_output(conf->cfout[i], "CF_OUTPUT_FILE", conf->ovwrite,
          FCFC_OVERWRITE_CFONLY))) return e;
    }

    if (conf->bintype == FCFC_BIN_SMU) {
      /* MULTIPOLE */
      if ((conf->npole = cfg_get_size(cfg, &conf->poles))) {
        /* Sort multipoles and remove duplicates. */
        if (conf->npole > 1) {
          /* 5-line insertion sort from https://doi.org/10.1145/3812.315108 */
          for (int i = 1; i < conf->npole; i++) {
            int tmp = conf->poles[i];
            for (num = i; num > 0 && conf->poles[num - 1] > tmp; num--)
              conf->poles[num] = conf->poles[num - 1];
            conf->poles[num] = tmp;
          }
          /* Remove duplicates from the sorted array. */
          num = 0;
          for (int i = 1; i < conf->npole; i++) {
            if (conf->poles[i] != conf->poles[num]) {
              num++;
              conf->poles[num] = conf->poles[i];
            }
          }
          conf->npole = num + 1;
        }
        if (conf->poles[0] < 0 || conf->poles[conf->npole - 1] > FCFC_MAX_ELL) {
          P_ERR(FMT_KEY(MULTIPOLE) " must be between 0 and %d\n", FCFC_MAX_ELL);
          return FCFC_ERR_CFG;
        }

        /* MULTIPOLE_FILE */
        CHECK_EXIST_ARRAY(MULTIPOLE_FILE, cfg, &conf->mpout, num);
        CHECK_STR_ARRAY_LENGTH(MULTIPOLE_FILE, cfg, conf->mpout,
            num, conf->ncf);
        for (int i = 0; i < conf->ncf; i++) {
          if ((e = check_output(conf->mpout[i], "MULTIPOLE_FILE",
              conf->ovwrite, FCFC_OVERWRITE_CFONLY))) return e;
        }
      }
    }
    else if (conf->bintype == FCFC_BIN_SPI) {
      /* PROJECTED_CF */
      if (!cfg_is_set(cfg, &conf->wp)) conf->wp = DEFAULT_PROJECTED_CF;
      if (conf->wp) {
        /* PROJECTED_FILE */
        CHECK_EXIST_ARRAY(PROJECTED_FILE, cfg, &conf->wpout, num);
        CHECK_STR_ARRAY_LENGTH(PROJECTED_FILE, cfg, conf->wpout,
            num, conf->ncf);
        for (int i = 0; i < conf->ncf; i++) {
          if ((e = check_output(conf->wpout[i], "PROJECTED_FILE",
              conf->ovwrite, FCFC_OVERWRITE_CFONLY))) return e;
        }
      }
    }
  }

  /* SEP_BIN_FILE */
  if (cfg_is_set(cfg, &conf->fsbin)) {
    if ((e = check_input(conf->fsbin, "SEP_BIN_FILE"))) return e;
  }
  else {
    /* SEP_BIN_MIN, SEP_BIN_MAX, SEP_BIN_SIZE */
    CHECK_EXIST_PARAM(SEP_BIN_MIN, cfg, &conf->smin);
    CHECK_EXIST_PARAM(SEP_BIN_MAX, cfg, &conf->smax);
    CHECK_EXIST_PARAM(SEP_BIN_SIZE, cfg, &conf->ds);
    if (!isfinite(conf->ds) || conf->ds <= 0) {
      P_ERR(FMT_KEY(SEP_BIN_SIZE) " must be finite and positive\n");
      return FCFC_ERR_CFG;
    }
    if (conf->smin + conf->ds > conf->smax + DOUBLE_TOL) {
      P_ERR(FMT_KEY(SEP_BIN_MIN) " + " FMT_KEY(SEP_BIN_SIZE)
          " cannot be larger than " FMT_KEY(SEP_BIN_MAX) "\n");
      return FCFC_ERR_CFG;
    }
    double smax = conf->smin;
    conf->nsbin = 0;
    while (smax < conf->smax - DOUBLE_TOL) {
      smax += conf->ds;
      if (++conf->nsbin > FCFC_MAX_NSBIN) {
        P_ERR("too many separations bins given " FMT_KEY(SEP_BIN_MIN) ", "
            FMT_KEY(SEP_BIN_MAX) ", and " FMT_KEY(SEP_BIN_SIZE) "\n");
        return FCFC_ERR_CFG;
      }
    }
    if (smax > conf->smax + DOUBLE_TOL) {
      P_WRN("reduce " FMT_KEY(SEP_BIN_MAX) " to " OFMT_DBL " given "
          FMT_KEY(SEP_BIN_MIN) " and " FMT_KEY(SEP_BIN_SIZE) "\n", smax);
    }
    conf->smax = smax;
  }

  if (conf->npole) {
    /* MU_BIN_NUM */
    CHECK_EXIST_PARAM(MU_BIN_NUM, cfg, &conf->nmu);
    if (conf->nmu <= 1) {
      P_ERR(FMT_KEY(MU_BIN_NUM) " must be larger than 1\n");
      return FCFC_ERR_CFG;
    }
    if (conf->nmu > FCFC_MAX_NMU) {
      P_ERR(FMT_KEY(MU_BIN_NUM) " cannot be larger than %d\n", FCFC_MAX_NMU);
      return FCFC_ERR_CFG;
    }
  }
  else if (conf->bintype == FCFC_BIN_ISO) conf->nmu = 1;

  if (conf->bintype == FCFC_BIN_SPI) {
    /* PI_BIN_FILE */
    if (cfg_is_set(cfg, &conf->fpbin)) {
      if ((e = check_input(conf->fpbin, "PI_BIN_FILE"))) return e;
    }
    else {
      /* PI_BIN_MIN, PI_BIN_MAX, PI_BIN_SIZE */
      CHECK_EXIST_PARAM(PI_BIN_MIN, cfg, &conf->pmin);
      CHECK_EXIST_PARAM(PI_BIN_MAX, cfg, &conf->pmax);
      CHECK_EXIST_PARAM(PI_BIN_SIZE, cfg, &conf->dpi);
      if (!isfinite(conf->dpi) || conf->dpi <= 0) {
        P_ERR(FMT_KEY(PI_BIN_SIZE) " must be finite and positive\n");
        return FCFC_ERR_CFG;
      }
      if (conf->pmin + conf->dpi > conf->pmax + DOUBLE_TOL) {
        P_ERR(FMT_KEY(PI_BIN_MIN) " + " FMT_KEY(PI_BIN_SIZE)
            " cannot be larger than " FMT_KEY(PI_BIN_MAX) "\n");
        return FCFC_ERR_CFG;
      }
      double pmax = conf->pmin;
      conf->npbin = 0;
      while (pmax < conf->pmax - DOUBLE_TOL) {
        pmax += conf->dpi;
        if (++conf->npbin > FCFC_MAX_NSBIN) {
          P_ERR("too many pi bins given " FMT_KEY(PI_BIN_NUM) ", "
              FMT_KEY(PI_BIN_MAX) ", and " FMT_KEY(PI_BIN_SIZE) "\n");
          return FCFC_ERR_CFG;
        }
      }
      if (pmax > conf->pmax + DOUBLE_TOL) {
        P_WRN("reduce " FMT_KEY(PI_BIN_MAX) " to " OFMT_DBL " given "
            FMT_KEY(PI_BIN_MIN) " and " FMT_KEY(PI_BIN_SIZE) "\n", pmax);
      }
      conf->pmax = pmax;
    }
  }

  /* OUTPUT_FORMAT */
  if (!cfg_is_set(cfg, &conf->ofmt)) conf->ofmt = DEFAULT_OUTPUT_FORMAT;
  switch (conf->ofmt) {
    case FCFC_OFMT_BIN:
    case FCFC_OFMT_ASCII:
      break;
    default:
      P_ERR("invalid " FMT_KEY(OUTPUT_STYLE) ": %d\n", conf->ofmt);
      return FCFC_ERR_CFG;
  }

  /* VERBOSE */
  if (!cfg_is_set(cfg, &conf->verbose)) conf->verbose = DEFAULT_VERBOSE;

  return 0;
}


/*============================================================================*\
                      Function for printing configurations
\*============================================================================*/

/******************************************************************************
Function `conf_print`:
  Print configuration parameters.
Arguments:
  * `conf`:     structure for storing configurations;
  * `ntask`:    number of MPI tasks;
  * `nthread`:  number of OpenMP threads.
******************************************************************************/
static void conf_print(const CONF *conf
#ifdef MPI
    , const int ntask
#endif
#ifdef OMP
    , const int nthread
#endif
    ) {
  /* Configuration file */
  printf("\n  CONFIG_FILE     = %s", conf->fconf);

  /* Input catalogs. */
  printf("\n  CATALOG         = %s", conf->input[0]);
  for (int i = 1; i < conf->ninput; i++)
    printf("\n                    %s", conf->input[i]);
  printf("\n  CATALOG_LABEL   = '%c'", conf->label[0]);
  for (int i = 1; i < conf->ninput; i++) printf(" , '%c'", conf->label[i]);

  const char *ftype[] = {"ASCII", "FITS", "HDF5"};
  const int ntype = sizeof(ftype) / sizeof(ftype[0]);
  if (!conf->ftype) {
    printf("\n  CATALOG_TYPE    = %d (%s)",
        DEFAULT_FILE_TYPE,
        DEFAULT_FILE_TYPE < ntype ? ftype[DEFAULT_FILE_TYPE] : "unknown");
  }
  else {
    printf("\n  CATALOG_TYPE    = %d (%s)",
        conf->ftype[0],
        conf->ftype[0] < ntype ? ftype[conf->ftype[0]] : "unknown");
    for (int i = 1; i < conf->ninput; i++) {
      printf("\n                    %d (%s)",
          conf->ftype[i],
          conf->ftype[i] < ntype ? ftype[conf->ftype[i]] : "unknown");
    }
  }

  if (conf->ascii) {
    if (!conf->skip) {
      printf("\n  ASCII_SKIP      = %ld", (long) DEFAULT_ASCII_SKIP);
    }
    else {
      printf("\n  ASCII_SKIP      = %ld", conf->skip[0]);
      for (int i = 1; i < conf->ninput; i++) printf(" , %ld", conf->skip[i]);
    }

    if (!conf->comment) {
      if (DEFAULT_ASCII_COMMENT == 0) printf("\n  ASCII_COMMENT   = ''");
      else printf("\n  ASCII_COMMENT   = '%c'", DEFAULT_ASCII_COMMENT);
    }
    else {
      int type = conf->ftype ? conf->ftype[0] : DEFAULT_FILE_TYPE;
      if (type != FCFC_FFMT_ASCII || conf->comment[0] == 0)
        printf("\n  ASCII_COMMENT   = ''");
      else printf("\n  ASCII_COMMENT   = '%c'", conf->comment[0]);
      for (int i = 1; i < conf->ninput; i++) {
        type = conf->ftype ? conf->ftype[i] : DEFAULT_FILE_TYPE;
        if (type != FCFC_FFMT_ASCII || conf->comment[i] == 0) printf(" , ''");
        else printf(" , '%c'", conf->comment[i]);
      }
    }

    printf("\n  ASCII_FORMATTER = %s", conf->fmtr[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf("\n                    %s", conf->fmtr[i]);
  }

  printf("\n  POSITION        = %s , %s , %s",
      conf->pos[0], conf->pos[1], conf->pos[2]);
  for (int i = 1; i < conf->ninput; i++) {
    printf("\n                    %s , %s , %s",
        conf->pos[i * 3], conf->pos[i * 3 + 1], conf->pos[i * 3 + 2]);
  }

  if (conf->wt) {
    printf("\n  WEIGHT          = %s", conf->wt[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf("\n                    %s", conf->wt[i]);
  }

  if (conf->sel) {
    printf("\n  SELECTION       = %s", conf->sel[0]);
    for (int i = 1; i < conf->ninput; i++)
      printf("\n                    %s", conf->sel[i]);
  }

  if (conf->cnvt) {
    printf("\n  COORD_CONVERT   = %c", conf->cnvt[0] ? 'T' : 'F');
    for (int i = 1; i < conf->ninput; i++)
      printf(" , %c", conf->cnvt[i] ? 'T' : 'F');
  }
  else printf("\n  COORD_CONVERT   = %c", DEFAULT_COORD_CNVT ? 'T' : 'F');

  /* Fiducial cosmology. */
  if (conf->has_cnvt) {
    if (conf->fcnvt) printf("\n  Z_CMVDST_FILE   = %s", conf->fcnvt);
    else {
      printf("\n  OMEGA_M         = " OFMT_DBL, conf->omega_m);
      printf("\n  OMEGA_LAMBDA    = " OFMT_DBL, conf->omega_l);
      if (conf->dew != -1)
        printf("\n  DE_EOS_W        = " OFMT_DBL, conf->dew);
      printf("\n  CMVDST_ERR      = " OFMT_DBL, conf->ecnvt);
    }
  }

  /* 2PCF configurations. */
  const char *tname[] = {"k-d tree", "ball tree"};
  const int ntname = sizeof(tname) / sizeof(tname[0]);
  printf("\n  DATA_STRUCT     = %d (%s)", conf->dstruct,
      conf->dstruct < ntname ? tname[conf->dstruct] : "unknown");

  const char *bname[] = {"s", "s & mu", "s_perp & pi"};
  const int nbname = sizeof(bname) / sizeof(bname[0]);
  printf("\n  BINNING_SCHEME  = %d (%s)", conf->bintype,
      conf->bintype < nbname ? bname[conf->bintype] : "unknown");
  printf("\n  PAIR_COUNT      = %s", conf->pc[0]);
  for (int i = 1; i < conf->npc; i++) printf(" , %s", conf->pc[i]);

  printf("\n  PAIR_COUNT_FILE = <%c> %s",
      conf->comp_pc[0] ? 'W' : 'R', conf->pcout[0]);
  for (int i = 1; i < conf->npc; i++) {
    printf("\n                    <%c> %s",
        conf->comp_pc[i] ? 'W' : 'R', conf->pcout[i]);
  }

  if (conf->ncf) {
    printf("\n  CF_ESTIMATOR    = %s", conf->cf[0]);
    for (int i = 1; i < conf->ncf; i++)
      printf("\n                    %s", conf->cf[i]);
    printf("\n  CF_OUTPUT_FILE  = %s", conf->cfout[0]);
    for (int i = 1; i < conf->ncf; i++)
      printf("\n                    %s", conf->cfout[i]);

    if (conf->bintype == FCFC_BIN_SMU && conf->npole) {
      printf("\n  MULTIPOLE       = %d", conf->poles[0]);
      for (int i = 1; i < conf->npole; i++) printf(" , %d", conf->poles[i]);
      printf("\n  MULTIPOLE_FILE  = %s", conf->mpout[0]);
      for (int i = 1; i < conf->ncf; i++)
        printf("\n                    %s", conf->mpout[i]);
    }

    if (conf->bintype == FCFC_BIN_SPI) {
      printf("\n  PROJECTED_CF    = %c", conf->wp ? 'T' : 'F');
      if (conf->wp) {
        printf("\n  PROJECTED_FILE  = %s", conf->wpout[0]);
        for (int i = 1; i < conf->ncf; i++)
          printf("\n                    %s", conf->wpout[i]);
      }
    }
  }

  /* Bin definitions. */
  if (conf->fsbin) printf("\n  SEP_BIN_FILE    = %s", conf->fsbin);
  else {
    printf("\n  SEP_BIN_MIN     = " OFMT_DBL, conf->smin);
    printf("\n  SEP_BIN_MAX     = " OFMT_DBL, conf->smax);
    printf("\n  SEP_BIN_SIZE    = " OFMT_DBL " (%d bins)",
        conf->ds, conf->nsbin);
  }
  if (conf->bintype == FCFC_BIN_SMU)
    printf("\n  MU_BIN_NUM      = %d", conf->nmu);

  if (conf->bintype == FCFC_BIN_SPI) {
    if (conf->fpbin) printf("\n  PI_BIN_FILE     = %s", conf->fpbin);
    else {
      printf("\n  PI_BIN_MIN      = " OFMT_DBL, conf->pmin);
      printf("\n  PI_BIN_MAX      = " OFMT_DBL, conf->pmax);
      printf("\n  PI_BIN_SIZE     = " OFMT_DBL " (%d bins)",
          conf->dpi, conf->npbin);
    }
  }

  /* Others. */
  const char *sname[] = {"binary", "ASCII"};
  const int nsname = sizeof(sname) / sizeof(sname[0]);
  printf("\n  OUTPUT_FORMAT   = %d (%s)", conf->ofmt,
      conf->ofmt < nsname ? sname[conf->ofmt] : "unknown");
  printf("\n  OVERWRITE       = %d", conf->ovwrite);

#ifdef MPI
  printf("\n  MPI_NUM_TASKS   = %d", ntask);
#endif
#ifdef OMP
  printf("\n  OMP_NUM_THREADS = %d", nthread);
#endif
  printf("\n");
}


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments;
  * `ntask`:    number of MPI tasks;
  * `nthread`:  number of OpenMP threads.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv
#ifdef MPI
    , const int ntask
#endif
#ifdef OMP
    , const int nthread
#endif
    ) {
  CONF *conf = conf_init();
  if (!conf) return NULL;

  cfg_t *cfg = conf_read(conf, argc, argv);
  if (!cfg) {
    conf_destroy(conf);
    return NULL;
  }

  printf("Loading configurations ...");
  fflush(stdout);

  if (conf_verify(cfg, conf)) {
    if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
    conf_destroy(conf);
    cfg_destroy(cfg);
    return NULL;
  }

  if (conf->verbose)
    conf_print(conf
#ifdef MPI
        , ntask
#endif
#ifdef OMP
        , nthread
#endif
        );

  if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
  cfg_destroy(cfg);

  printf(FMT_DONE);
#ifdef MPI
  fflush(stdout);
#endif
  return conf;
}

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf) {
  if (!conf) return;
  FREE_STR_ARRAY(conf->input);
  FREE_ARRAY(conf->label);
  FREE_ARRAY(conf->ftype);
  FREE_ARRAY(conf->skip);
  FREE_ARRAY(conf->comment);
  FREE_STR_ARRAY(conf->fmtr);
  FREE_STR_ARRAY(conf->pos);
  FREE_STR_ARRAY(conf->wt);
  FREE_STR_ARRAY(conf->sel);
  FREE_ARRAY(conf->has_wt);
  FREE_ARRAY(conf->cnvt);
  FREE_ARRAY(conf->fcnvt);
  FREE_STR_ARRAY(conf->pc);
  FREE_ARRAY(conf->comp_pc);
  FREE_STR_ARRAY(conf->pcout);
  FREE_STR_ARRAY(conf->cf);
  FREE_STR_ARRAY(conf->cfout);
  FREE_ARRAY(conf->poles);
  FREE_STR_ARRAY(conf->mpout);
  FREE_STR_ARRAY(conf->wpout);
  FREE_ARRAY(conf->fsbin);
  FREE_ARRAY(conf->fpbin);
  free(conf);
}

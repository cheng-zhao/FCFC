/*******************************************************************************
* libcfg.h: this file is part of the libcfg library.

* libcfg: C library for parsing command line option and configuration files.

* Github repository:
        https://github.com/cheng-zhao/libcfg

* Copyright (c) 2019 Cheng Zhao <zhaocheng03@gmail.com>
 
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

#ifndef _LIBCFG_H_
#define _LIBCFG_H_

#include <stdio.h>
#include <stdbool.h>

/*============================================================================*\
                           Definitions for data types
\*============================================================================*/

typedef enum {
  CFG_DTYPE_NULL,
  CFG_DTYPE_BOOL,
  CFG_DTYPE_CHAR,
  CFG_DTYPE_INT,
  CFG_DTYPE_LONG,
  CFG_DTYPE_FLT,
  CFG_DTYPE_DBL,
  CFG_DTYPE_STR,
  CFG_ARRAY_BOOL,
  CFG_ARRAY_CHAR,
  CFG_ARRAY_INT,
  CFG_ARRAY_LONG,
  CFG_ARRAY_FLT,
  CFG_ARRAY_DBL,
  CFG_ARRAY_STR
} cfg_dtype_t;

#define CFG_DTYPE_INVALID(x)    ((x) < CFG_DTYPE_BOOL || (x) > CFG_ARRAY_STR)
#define CFG_DTYPE_IS_ARRAY(x)   ((x) >= CFG_ARRAY_BOOL && (x) <= CFG_ARRAY_STR)

/*============================================================================*\
                         Definitions for string lengths
\*============================================================================*/
#define CFG_MAX_NAME_LEN        128
#define CFG_MAX_LOPT_LEN        128
#define CFG_MAX_FILENAME_LEN    1024

/*============================================================================*\
                          Definitions for the formats
\*============================================================================*/
#define CFG_SYM_EQUAL           '='
#define CFG_SYM_ARRAY_START     '['
#define CFG_SYM_ARRAY_END       ']'
#define CFG_SYM_ARRAY_SEP       ','
#define CFG_SYM_COMMENT         '#'
#define CFG_SYM_NEWLINE         '\\'

#define CFG_CMD_FLAG            '-'
#define CFG_CMD_ASSIGN          '='


/*============================================================================*\
                         Definition of data structures
\*============================================================================*/

/* Main entry for all configuration parameters and command line functions. */
typedef struct {
  int npar;             /* number of verified configuration parameters  */
  int nfunc;            /* number of verified command line functions    */
  void *params;         /* data structure for storing parameters        */
  void *funcs;          /* data structure for storing function pointers */
  void *error;          /* data structure for storing error messages    */
} cfg_t;

/* Interface for registering configuration parameters. */
typedef struct {
  int opt;                      /* short command line option            */
  char *lopt;                   /* long command line option             */
  char *name;                   /* name of the parameter                */
  cfg_dtype_t dtype;            /* data type of the parameter           */
  void *var;                    /* variable for the retrieved value     */
} cfg_param_t;

/* Interface for registering command line functions. */
typedef struct {
  int opt;                      /* short command line option            */
  char *lopt;                   /* long command line option             */
  void (*func) (void *);        /* pointer to the function              */
  void *args;                   /* pointer to the arguments             */
} cfg_func_t;


/*============================================================================*\
                            Definition of functions
\*============================================================================*/

/******************************************************************************
Function `cfg_init`:
  Initialise the entry for all parameters and command line functions.
Return:
  The address of the structure.
******************************************************************************/
cfg_t *cfg_init(void);

/******************************************************************************
Function `cfg_set_params`:
  Verify and register configuration parameters.
Arguments:
  * `cfg`:      entry for all configuration parameters;
  * `param`:    stucture for the input configuration parameters;
  * `npar`:     number of input configuration parameters.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cfg_set_params(cfg_t *cfg, const cfg_param_t *param, const int npar);

/******************************************************************************
Function `cfg_set_funcs`:
  Verify and register command line functions.
Arguments:
  * `cfg`:      entry for all command line functions;
  * `func`:     stucture for the input functions;
  * `nfunc`:    number of input functions.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cfg_set_funcs(cfg_t *cfg, const cfg_func_t *func, const int nfunc);

/******************************************************************************
Function `cfg_read_opts`:
  Parse command line options.
Arguments:
  * `cfg`:      entry for the configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments;
  * `prior`:    priority of values set via command line options;
  * `optidx`:   position of the first unparsed argument.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cfg_read_opts(cfg_t *cfg, const int argc, char *const *argv,
    const int prior, int *optidx);

/******************************************************************************
Function `cfg_read_file`:
  Read configuration parameters from a file.
Arguments:
  * `cfg`:      entry for the configurations;
  * `fname`:    name of the input file;
  * `prior`:    priority of values read from this file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int cfg_read_file(cfg_t *cfg, const char *fname, const int prior);

/******************************************************************************
Function `cfg_is_set`:
  Check if a variable is set via the command line or files.
Arguments:
  * `cfg`:      entry of all configurations;
  * `var`:      address of the variable.
Return:
  True if the variable is set; false otherwise.
******************************************************************************/
bool cfg_is_set(const cfg_t *cfg, const void *var);

/******************************************************************************
Function `cfg_get_size`:
  Return the number of elements for the parsed array.
Arguments:
  * `cfg`:      entry of all configurations;
  * `var`:      address of the variable.
Return:
  The number of array elements on success; 0 on error.
******************************************************************************/
int cfg_get_size(const cfg_t *cfg, const void *var);

/******************************************************************************
Function `cfg_destroy`:
  Release memory allocated for the configuration parameters.
Arguments:
  * `cfg`:      pointer to the entry of all configurations.
******************************************************************************/
void cfg_destroy(cfg_t *cfg);

/******************************************************************************
Function `cfg_perror`:
  Print the error message if there is an error.
Arguments:
  * `cfg`:      entry of all configurations;
  * `fp`:       output file stream to write to;
  * `msg`:      string to be printed before the error message.
******************************************************************************/
void cfg_perror(const cfg_t *cfg, FILE *fp, const char *msg);

/******************************************************************************
Function `cfg_pwarn`:
  Print the warning messages if there is any, and clean the warnings.
Arguments:
  * `cfg`:      entry of all configurations;
  * `fp`:       output file stream to write to;
  * `msg`:      string to be printed before the error message.
******************************************************************************/
void cfg_pwarn(cfg_t *cfg, FILE *fp, const char *msg);

#endif

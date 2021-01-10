# Configuration parameters for `FCFC_2PT_BOX`

Template configuration file: [`fcfc_2pt_box.conf`](../etc/fcfc_2pt_box.conf).

## Table of Contents

-   [Specifications of the input catalogues](#specifications-of-the-input-catalogs)
-   [Configurations for the 2-point correlation function](#configurations-for-the-2-point-correlation-function)
-   [Definitions of the bins](#definitions-of-the-bins)
-   [Other settings](#other-settings)

## Specifications of the input catalogues

### `CATALOG`

Filename of the input catalogues. They can be either strings or string arrays.

*Examples*

```nginx
CATALOG = /path/to/input_data.dat
CATALOG = [ data_catalog.dat, \
            rand_catalog.dat ]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CATALOG_LABEL`

Label of the input catalogs. They must be non-repetitive uppercase letters, and are used for indicating the sources of the pairs to be counted. For instance, the two sources of a pair "`DR`" are the catalogues labelled by "`D`" and "`R`", respectively.

If this parameter is not set, the input catalogues are labelled in alphabetical order, i.e. [`A`, `B`, ...].

*Examples*

```nginx
CATALOG_LABEL = A
CATALOG_LABEL = [D,R,S]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CATALOG_TYPE`

Format of the input catalogs. They must be integers or integer arrays, depending on the dimension of  `CATALOG`. The allowed values are:

-   `0`: for ASCII format text files (default);
-   `1`: for FITS tables (*coming soon*).

In particular, FITS tables are supported via the CFITSIO library, so the compilation flag `-DWITH_CFITSIO` has to be enabled for reading FITS files.

*Examples*

```nginx
CATALOG_TYPE = 0
CATALOG_TYPE = [0,1]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `ASCII_SKIP`

Number of lines to be skipped for ASCII format input files. They can be non-negative long integers or long integer arrays, depending on the dimension of `CATALOG`.

*Examples*

```nginx
ASCII_SKIP = 10
ASCII_SKIP = [0,5]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `ASCII_COMMENT`

Indicator of comment lines for ASCII format input files. They can be characters or character arrays, depending on the dimension of `CATALOG`. If the first non-whitespace character of a line is the specified character, then the whole line of the input catalogue is omitted. If empty characters (`''`) are supplied, then no line is treated as comments.

*Examples*

```nginx
ASCII_COMMENT = '#'
ASCII_COMMENT = ['', '!']
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `ASCII_FORMATTER`

C99-style formatter specifying the format of columns of ASCII format input files. They must be strings or string arrays, depending on the dimension of  `CATALOG`. Note however that only formatters with the following argument types are supported (see  [cppreference.com](https://en.cppreference.com/w/c/io/fscanf)  for details):

-   `int *`
-   `long *`
-   `float *`
-   `double *`
-   `char *`

*Examples*

```nginx
ASCII_FORMATTER = "%d %ld %f %lf %s"  # for int, long, float, double, and string types
ASCII_FORMATTER = "%*d,%10s,%[123]"
        # Column separators are ',';
        # The first column is treated as an integer, but is omitted;
        # The second column is a string with 10 characters;
        # The third column is a string composed of characters '1', '2', and '3'.
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `POSITION`

3-D coordinates, in the order of {*x*,  *y*,  *z*} or {RA, Dec, redshift}, where RA and Dec must be in degrees. The dimension of this parameter must be at least 3 times that of `CATALOG`.

The strings must be column indicators, or expressions, which are parsed using the [libast](https://github.com/cheng-zhao/libast) library. Columns are indicated by <span><code>${&bull;}</code></span>, where <span><code>&bull;</code></span> must be a number for ASCII format files, and a string for FITS tables. For instance, the 3rd column of an ASCII file can be indicated by `$3`, and the "RA" column of a FITS table can be indicated by `${RA}`. Note that if there are more than one digits for the ASCII column numbers, the number must be enclosed by braces, e.g, `${10}`. And if an ASCII column is omitted via the formatter (e.g. `%*lf`), it is not counted for the column number.

Moreover, expressions are supported for pre-processing the columns, with some basic arithmetic operators, and mathematical functions (see  [libast](https://github.com/cheng-zhao/libast)  for details).

*Examples*

```nginx
POSITION = [${RA}*180/3.1415927, ${DEC}, ${Z}, ${RA}, ${DEC}, ${Z}]
POSITION = [($1+1000)%1000, ($2+1000)%1000, ($3+1000)%1000]
POSITION = [$1, $2, ($3+$6*(1+0.6)/(100*sqrt(0.31*(1+0.6)**3+0.69))+1000)%1000]
      # The last expression implies real to redshift space conversion
      # with periodic boundary conditions given the box size 1000 Mpc/h,
      # in a FlatLCDM cosmology with Omega_m = 0.31, at redshift 0.6.
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SELECTION`

Selection criteria for the input catalogues. They can be column indicators or expressions. Numerical, bitwise, and logical expressions are all supported. Only objects with columns fulfilling the conditions are kept.

*Examples*

```nginx
SELECTION = isfinite($1) && $3 == "YES" && log($4) < 1.0
SELECTION = [($5 & 1) != 0, $3 >= 0.5]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `BOX_SIZE`

Side length of the box that catalogues are placed in. It must be a double-precision floating-point number. All distances are evaluated following periodic boundary conditions with this box size.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Configurations for the 2-point correlation function

### `BINNING_SCHEME`

Binning scheme of the 2PCF to be evaluated. It must be an integer, and the allowed values are:

-   `0`: isotropic separation bins for *&xi;*(*s*)
-   `1`: anisotropic separation bins for *&xi;*(*s*, *&mu;*) (required by 2PCF multipoles);
-   `2`: 2-D separation bins for *&xi;*(*s*<sub>perp</sub>, *&pi;*) (required by the projected 2PCF).

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PAIR_COUNT`

Identifiers of pair counts, which are consist of `CATALOG_LABEL`. The labels indicate the source catalogues of the pairs. This parameter must be two-character strings or string arrays.

*Examples*

```nginx
PAIR_COUNT = DD          # auto pair counts of catalog 'D'
PAIR_COUNT = [DD,DR,RR]  # "DR" denotes cross pairs from catalogs 'D' and 'R'
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PAIR_COUNT_FILE`

Name of the files for pair counts. It can be a string or string array, depending on the dimension of `PAIR_COUNT`. If a specified file exists, then the pair counts are read from this file; otherwise the pair counts are evaluated and saved to the file.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CF_ESTIMATOR`

Estimators of the 2PCFs to be evaluated. It can be expressions consist of `PAIR_COUNT`. In particular, `@@` denotes the analytical random-random pair counts.

*Examples*

```nginx
CF_ESTIMATOR = DD / @@ - 1            # natural estimator
CF_ESTIMATOR = (DD - 2*DS + SS) / @@  # Landy-Szalay estimator for reconstructed catalogues
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CF_OUTPUT_FILE`

Name of the output files for 2PCFs. It can be a string or string array, depending on the dimension of `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `MULTIPOLE`

Legendre multipoles of the power spectra to be evaluated. It must be a non-negative integer, or integer arrays. The current maximum supported *&ell;* is `6`. Note that these multipoles are evaluated for all 2PCFs defined by `CF_ESTIMATOR`.

*Examples*

```nginx
MULTIPOLE = [0,2,4]  # for monopole. quadrupole, and hexadecapole
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `MULTIPOLE_FILE`

Name of the output files for 2PCF multipoles. It can be a string or string array, depending on the dimension of `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PROJECTED_CF`

Indicate whether the projected 2PCF is required. It must be a boolean value, and is applied to all 2PCFs defined by `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PROJECTED_FILE`

Name of the output files for the projected 2PCFs. It can be a string or string array, depending on the dimension of `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Definitions of the bins

*Note*: all bins defined in this section are left-closed and right-open, i.e., the lower boundary is included in the bin, while the upper boundary is not.

### `SEP_BIN_FILE`

Name of the file defining separation (either *s* or *s*<sub>perp</sub> bins. The file should be an ASCII table, with the first two columns being the lower and upper boundary of each separation bin, respectively.

If this parameter is not set, then linear separation bins defined by `SEP_BIN_MIN`, `SEP_BIN_MAX`, and `SEP_BIN_SIZE` are used.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SEP_BIN_MIN`

Lower boundary of separation (either *s* or *s*<sub>perp</sub> bins. It is the minimum separation of interest. It must be a double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SEP_BIN_MAX`

Upper boundary of separation (either *s* or *s*<sub>perp</sub> bins. Only separations below this value are recorded. It must be a double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SEP_BIN_SIZE`

The width of all separation bins. It must be a positive double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `MU_BIN_NUM`

Number of *&mu;* bins in the range [0, 1). It must be a positive integer.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_FILE`

Name of the file defining *&pi;* bins. The file should be an ASCII table, with the first two columns being the lower and upper boundary of each *&pi;* bin, respectively.

If this parameter is not set, then linear separation bins defined by `PI_BIN_MIN`, `PI_BIN_MAX`, and `PI_BIN_SIZE` are used.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_MIN`

Lower boundary of *&pi;* bins. It is the minimum separation of interest. It must be a double-precision floating-point number.


<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_MAX`

Upper boundary of *&pi;* bins. Only separations below this value are recorded. It must be a double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_SIZE`

The width of all *&pi;* bins. It must be a positive double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SQ_DIST_PREC`

Precision when distributing (squared) distances into separation bins. It must be an integer, and digits of squared distances after 10<sup>`SQ_DIST_PREC`</sup> are truncated. For instance, the squared distance 1234.5678 is truncated to 1230 if `SQ_DIST_PREC` = 1, or 1234.5 of `SQ_DIST_PREC` = &minus;1. If `SQ_DIST_PREC` is outside the range [&minus;4,4], then no truncation is applied.

The truncation is for fast separation bin search. It does not necessarily result in the loss of accuracy for pair counts. For instance, if the edges of the distance bins are all integers &mdash; which is true for most practical cases &mdash; then the results are still exact if `SQ_DIST_PREC` = 0.

Note that only *&pi;* is truncated with out being squared.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Other settings

### `OUTPUT_STYLE`

Define the format of the outputs. It must be an integer, and the allowed values are:

-   `0`: flatten 2-D results and save them as 1-D lists;
-   `1`: save 2-D results as matrices.

More details of the output file formats can be found in the header of the files.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `OVERWRITE`

An integer value indicating whether to overwrite existing files. Allowed values are

-   `0`: quit the program when an output file exist;
-   postive: force overwriting output files whenever possible;
-   negative: notify at most this number (absolute value) of times, for asking whether overwriting existing files.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `VERBOSE`

A boolean value indicating whether to show detailed standard outputs.

<sub>[\[TOC\]](#table-of-contents)</sub>


# Configuration parameters for `FCFC_2PT_BOX`

Template configuration file: [`fcfc_2pt_box.conf`](../etc/fcfc_2pt_box.conf).

## Table of Contents

-   [Specifications of the input catalogues](#specifications-of-the-input-catalogs)
-   [Configurations for the 2-point correlation function](#configurations-for-the-2-point-correlation-function)
-   [Definitions of the bins](#definitions-of-the-bins)
-   [Other settings](#other-settings)

## Specifications of the input catalogues

### `CATALOG` (`-i` / `--input`)

Filename of the input catalogues. They can be either strings or string arrays.

*Examples*

```nginx
CATALOG = /path/to/input_data.dat
CATALOG = [ data_catalog.dat, \
            rand_catalog.dat ]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CATALOG_LABEL` (`-l` / `--label`)

Label of the input catalogs. They must be non-repetitive uppercase letters, and are used for indicating the sources of the pairs to be counted. For instance, the two sources of a pair "`DR`" are the catalogues labelled by "`D`" and "`R`", respectively.

If this parameter is not set, the input catalogues are labelled in alphabetical order, i.e. [`A`, `B`, ...].

*Examples*

```nginx
CATALOG_LABEL = A
CATALOG_LABEL = [D,R,S]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CATALOG_TYPE` (`-T` / `--type`)

Format of the input catalogs. They must be integers or integer arrays, depending on the dimension of  `CATALOG`. The allowed values are:

-   `0` (default): for ASCII format text files;
-   `1`: for FITS tables;
-   `2`: for HDF5 files.

In particular, FITS and HDF5 formats are supported through external libraries (see [README.md](../README.md#compilation) for details).

*Examples*

```nginx
CATALOG_TYPE = 0
CATALOG_TYPE = [0,1]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `ASCII_SKIP` (`--skip`)

Number of lines to be skipped for ASCII format input files. They can be non-negative long integers or long integer arrays, depending on the dimension of `CATALOG`.

*Examples*

```nginx
ASCII_SKIP = 10
ASCII_SKIP = [0,5]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `ASCII_COMMENT` (`--comment`)

Indicator of comment lines for ASCII format input files. They can be characters or character arrays, depending on the dimension of `CATALOG`. If the first non-whitespace character of a line is the specified character, then the whole line of the input catalogue is omitted. If empty characters (`''`) are supplied, then no line is treated as comments.

*Examples*

```nginx
ASCII_COMMENT = '#'
ASCII_COMMENT = ['', '!']
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `ASCII_FORMATTER` (`-f` / `--formatter`)

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

### `POSITION` (`-x` / `--position`)

3-D coordinates, in the order of {*x*,  *y*,  *z*} or {RA, Dec, redshift}, where RA and Dec must be in degrees. The dimension of this parameter must be at least 3 times that of `CATALOG`.

The strings must be column indicators, or expressions, which are parsed using the [libast](https://github.com/cheng-zhao/libast) library. The syntaxes of column indicators for different file formats are summarised below:

| File format | Column indicator                              | Description                                                                                                                                           | Example                        |
|-------------|-----------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------|
| ASCII       | &dollar;*number* or &dollar;{*number*}        | *number* indicates the index of a column (starting from 1, omitted columns are not counted), and must be enclosed by braces if it is not less than 10 | `$3`, `${12}`                  |
| FITS        | &dollar;{*name*}                              | *name* indicates the name of a column (column with multiple values is not supported so far)                                                           | `${RA}`                        |
| HDF5        | &dollar;{*name*} or &dollar;{*name*(*index*)} | *name* indicates the full name of a dataset, and *index* indicates the index of the shorter dimension (starting from 1) if the the dataset is in 2D   | `${/data/x}`, `${position(1)}` |

Moreover, expressions are supported for pre-processing the columns, with some basic arithmetic operators, and mathematical functions (see  [libast](https://github.com/cheng-zhao/libast)  for details).

*Examples*

```nginx
POSITION = [(${X}+1000)%1000, (${Y}+1000)%1000, (${Z}+1000)%1000]
POSITION = [$1, $2, ($3+$6*(1+0.6)/(100*sqrt(0.31*(1+0.6)**3+0.69))+1000)%1000]
      # The second expression implies real to redshift space conversion
      # with periodic boundary conditions given the box size 1000 Mpc/h,
      # in a FlatLCDM cosmology with Omega_m = 0.31, at redshift 0.6.
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `WEIGHT` (`-w` / `--weight`)

Weights for pair counts. They can be column indicators or expressions, depending on the dimension of `CATALOG`. If this parameter is not set, no weights are applied to any of the catalogues.

*Examples*

```nginx
WEIGHT = ${WEIGHT}
WEIGHT = [1, $4 * $5]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SELECTION` (`-s` / `--select`)

Selection criteria for the input catalogues. They can be column indicators or expressions. Numerical, bitwise, and logical expressions are all supported. Only objects with columns fulfilling the conditions are kept.

*Examples*

```nginx
SELECTION = isfinite($1) && $3 == "YES" && log($4) < 1.0
SELECTION = [($5 & 1) != 0, $3 >= 0.5]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `BOX_SIZE` (`-b` / `--box`)

Side lengths of the box that catalogues are placed in. It can be a double-precision floating-point number for cubes, or 3-element double numbers for arbitrary cuboids. All distances are evaluated following periodic boundary conditions with the corresponding side lengths.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Configurations for the 2-point correlation function

### `DATA_STRUCT` (`-S` / `--data-struct`)

Data structure used for pair count evaluations. It must be an integer, and the allowed values are:

-   `0` (default): *k*-d tree;
-   `1`: ball tree.

In most cases the choice of data structure does not significantly affect the performance.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `BINNING_SCHEME` (`-B` / `--bin`)

Binning scheme of the 2PCF to be evaluated. It must be an integer, and the allowed values are:

-   `0` (default): isotropic separation bins for *&xi;*(*s*)
-   `1`: anisotropic separation bins for *&xi;*(*s*, *&mu;*) (required by 2PCF multipoles);
-   `2`: 2-D separation bins for *&xi;*(*s*<sub>perp</sub>, *&pi;*) (required by the projected 2PCF).

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PAIR_COUNT` (`-p` / `--pair`)

Identifiers of pair counts, which are consist of `CATALOG_LABEL`. The labels indicate the source catalogues of the pairs. This parameter must be two-character strings or string arrays.

*Examples*

```nginx
PAIR_COUNT = DD          # auto pair counts of catalog 'D'
PAIR_COUNT = [DD,DR,RR]  # "DR" denotes cross pairs from catalogs 'D' and 'R'
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PAIR_COUNT_FILE` (`-P` / `--pair-output`)

Name of the files for pair counts. It can be a string or string array, depending on the dimension of `PAIR_COUNT`. If a specified file exists, then the pair counts are read from this file; otherwise the pair counts are evaluated and saved to the file.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CF_ESTIMATOR` (`-e` / `--cf`)

Estimators of the 2PCFs to be evaluated. It can be expressions consist of `PAIR_COUNT`. In particular, `@@` denotes the analytical random-random pair counts.

*Examples*

```nginx
CF_ESTIMATOR = DD / @@ - 1            # natural estimator
CF_ESTIMATOR = (DD - 2*DS + SS) / @@  # Landy-Szalay estimator for reconstructed catalogues
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `CF_OUTPUT_FILE` (`-E` / `--cf-output`)

Name of the output files for 2PCFs. It can be a string or string array, depending on the dimension of `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `MULTIPOLE` (`-m` / `--multipole`)

Legendre multipoles of the power spectra to be evaluated. It must be a non-negative integer, or integer arrays. The current maximum supported *&ell;* is `6`. Note that these multipoles are evaluated for all 2PCFs defined by `CF_ESTIMATOR`.

*Examples*

```nginx
MULTIPOLE = [0,2,4]  # for monopole. quadrupole, and hexadecapole
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### `MULTIPOLE_FILE` (`-M` / `--mp-output`)

Name of the output files for 2PCF multipoles. It can be a string or string array, depending on the dimension of `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PROJECTED_CF` (`-u` / `--wp`)

Indicate whether the projected 2PCF is required. It must be a boolean value, and is applied to all 2PCFs defined by `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PROJECTED_FILE` (`-U` / `--wp-output`)

Name of the output files for the projected 2PCFs. It can be a string or string array, depending on the dimension of `CF_ESTIMATOR`.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Definitions of the bins

*Note*: all bins defined in this section are left-closed and right-open (unless otherwise stated), i.e., the lower boundary is included in the bin, while the upper boundary is not.

### `SEP_BIN_FILE` (`--s-file`)

Name of the file defining separation (either *s* or *s*<sub>perp</sub> bins. The file should be an ASCII table, with the first two columns being the lower and upper boundary of each separation bin, respectively.

If this parameter is not set, then linear separation bins defined by `SEP_BIN_MIN`, `SEP_BIN_MAX`, and `SEP_BIN_SIZE` are used.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SEP_BIN_MIN` (`--s-min`)

Lower boundary of separation (either *s* or *s*<sub>perp</sub> bins. It is the minimum separation of interest. It must be a double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SEP_BIN_MAX` (`--s-max`)

Upper boundary of separation (either *s* or *s*<sub>perp</sub> bins. Only separations below this value are recorded. It must be a double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `SEP_BIN_SIZE` (`--s-step`)

The width of all separation bins. It must be a positive double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `MU_BIN_NUM` (`--mu-num`)

Number of *&mu;* bins in the range [0, 1) or [0, 1] (if the program is compiled with `WITH_MU_ONE = T`, see [README.md](../README.md#compilation) for details). It must be a positive integer.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_FILE` (`--pi-file`)

Name of the file defining *&pi;* bins. The file should be an ASCII table, with the first two columns being the lower and upper boundary of each *&pi;* bin, respectively.

If this parameter is not set, then linear separation bins defined by `PI_BIN_MIN`, `PI_BIN_MAX`, and `PI_BIN_SIZE` are used.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_MIN` (`--pi-min`)

Lower boundary of *&pi;* bins. It is the minimum separation of interest. It must be a double-precision floating-point number.


<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_MAX` (`--pi-max`)

Upper boundary of *&pi;* bins. Only separations below this value are recorded. It must be a double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `PI_BIN_SIZE` (`--pi-step`)

The width of all *&pi;* bins. It must be a positive double-precision floating-point number.

<sub>[\[TOC\]](#table-of-contents)</sub>


## Other settings

### `OUTPUT_FORMAT` (`-F` / `--out-format`)

Format of the output pair count files (`PAIR_COUNT_FILE`). It must be an integer, and the allowed values are:

-   `0` (default): FCFC binary format;
-   `1`: ASCII text format.

The binary format stores all metadata of the pair counts, and can be read with the script [`read_pair_count.py`](../scripts/read_pair_count.py).

<sub>[\[TOC\]](#table-of-contents)</sub>

### `OVERWRITE` (`-O` / `--overwrite`)

An integer value indicating whether to overwrite existing files. Allowed values are

-   `0` (default): quit the program when an output file exist;
-   `1`: overwrite 2PCF files silently, but keep (read from) existing pair count files;
-   `2` or larger: overwrite all output files silently whenever possible;
-   negative: notify at most this number (absolute value) of times, for asking whether overwriting existing files.

<sub>[\[TOC\]](#table-of-contents)</sub>

### `VERBOSE` (`-v` / `--verbose`)

A boolean value indicating whether to show detailed standard outputs.

<sub>[\[TOC\]](#table-of-contents)</sub>


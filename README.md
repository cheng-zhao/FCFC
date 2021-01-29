# Fast Correlation Function Calculator (FCFC)

![GitHub](https://img.shields.io/github/license/cheng-zhao/FCFC.svg)
![Codacy grade](https://img.shields.io/codacy/grade/4a85732bb6264027aefac7f002550cdd.svg)

<img src="doc/logo/FCFC_logo.svg" align="right" />

## Table of Contents

-   [Introduction](#introduction)
-   [Compilation](#compilation)
-   [Components and configurations](#components-and-configurations)
-   [Acknowledgements](#acknowledgements)
-   [References](#references)

## Introduction

**F**ast **C**orrelation **F**unction **C**alculator (FCFC) is a C toolkit for computing cosmological correlation functions. It is designed in the hope of being (both time and space) efficient, portable, and user-friendly.

So far the following products are supported:
-   Isotropic 2-point correlation function (2PCF): *&xi;*(*s*);
-   Anisotropic 2PCF: *&xi;*(*s*, *&mu;*);
-   2-D 2PCF: *&xi;*(*s*<sub>perp</sub>, *s*<sub>para</sub>), also known as *&xi;*(*s*<sub>perp</sub>, *&pi;*);
-   2PCF multipoles: *&xi;*<sub>*&ell;*</sub>(*s*);
-   Projected 2PCF: *w*<sub>*p*</sub>(*s*<sub>perp</sub>).

This program is compliant with the ISO C99 and IEEE POSIX.1-2008 standards, and no external library is mandatory. Parallelisation can be enabled with [OpenMP](https://www.openmp.org). Thus it is compatible with most of the modern C compilers and operating systems. It is written by Cheng Zhao (&#36213;&#25104;), and is distributed under the [MIT license](LICENSE.txt).

If you use this program in research that results in publications, please cite the following paper:

> Zhao et al. 2020, [arXiv:2007.08997](https://ui.adsabs.harvard.edu/abs/2020arXiv200708997Z/abstract)

<sub>[\[TOC\]](#table-of-contents)</sub>

## Compilation

The following command should compile the code for most of the time:

```bash
make
```

To compile only a certain component of the program, the name of the component can be supplied via

```bash
make [COMPONENT_NAME]
```

<sub>[\[TOC\]](#table-of-contents)</sub>

## Components and configurations

FCFC comes along with several components for different tasks. They are served as separate executables, and have to be supplied the corresponding configurations, either via command line options or a text file with configuration parameters.

The list of available command line options can be consulted using the `-h` or `--help` flags, and a template configuration file can be printed via the `-t` or `--template` flags.

An introduction of the components and the corresponding configuration parameters are listed below:

| Component    | Description                                                     | Configuration parameters               |
|:------------:|-----------------------------------------------------------------|:--------------------------------------:|
| FCFC_2PT     | Compute 2PCF for survey-like data                               | [FCFC_2PT.md](doc/FCFC_2PT.md)         |
| FCFC_2PT_BOX | Compute 2PCF for periodic simulation boxes<sup>[*](#tab1)</sup> | [FCFC_2PT_BOX.md](doc/FCFC_2PT_BOX.md) |

<span id="tab1">*: treat the 3<sup>rd</sup> dimension (*z*-direction) as the line of sight</span>

<sub>[\[TOC\]](#table-of-contents)</sub>

## Acknowledgements

This program benefits from the following open-source projects:
-   [Fast Cubic Spline Interpolation](https://doi.org/10.5281/zenodo.3611922) (see also [arXiv:2001.09253](https://arxiv.org/abs/2001.09253))
-   [https://github.com/andralex/MedianOfNinthers](https://github.com/andralex/MedianOfNinthers) (see also [this paper](http://dx.doi.org/10.4230/LIPIcs.SEA.2017.24))
-   [https://github.com/swenson/sort](https://github.com/swenson/sort)

<sub>[\[TOC\]](#table-of-contents)</sub>

## References

TBE

<sub>[\[TOC\]](#table-of-contents)</sub>

## TODO

-   FITS and HDF5 file formats
-   Jackknife covariance estimation
-   3-point correlation functions
-   Approximate correlation function calculators
-   More data structures for potentially faster correlation function evaluation

<sub>[\[TOC\]](#table-of-contents)</sub>

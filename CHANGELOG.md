# FCFC Change Log

## [1.1.0] &ndash; 2025-06-11

-   Fixed an early-return bug affecting tiny catalogs.
-   Corrected separation bins of the saved projected correlation data ([PR #4](pull/4)).
-   Revised the normalization factor for weighted pair counts ([PR #3](pull/3))
-   Added a `Makefile` for benchmark code compilation.
-   Fixed a macro typo reported in [Issue #2](issues/2).
-   Updated the FCFC paper reference in the README.

## [1.0.1] &ndash; 2022-06-26

-   Refactorized the implementation of MPI and OpenMP parallelisms.
-   Optimized MPI communications with custom structs, and fixed the deadlock due to multiple broadcasts.

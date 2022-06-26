# Set the compiler and default compilation flags.
# * CC is used when WITH_MPI != T;
# * MPICC is used when WITH_MPI == T.
CC := gcc
MPICC := mpicc
CFLAGS = -O3 -march=native -flto

# Set `WITH_MPI := T` to enable MPI.
WITH_MPI := 

# Set `WITH_OMP := T` to enable OpenMP.
# The compiler flag for OpenMP can be set via `OMP_FLAG`.
WITH_OMP := T
OMP_FLAG := -fopenmp

# Set `WITH_SIMD := T` to enable SIMD.
# Do not forget to set the CPU type (e.g. `-march=native`) in `CFLAGS`.
WITH_SIMD := T

# Set `SINGLE_PREC := T` to enable single-precision.
SINGLE_PREC := 

# Set `WITH_MU_ONE := T` to include (s,mu) pairs with mu = 1.
WITH_MU_ONE := 

# Set `WITH_CFITSIO := T` to enable FITS file format.
# The paths for cfitsio header (fitsio.h) and library (libcfitsio.{a,so,dylib})
# files can be set via `CFITSIO_INC_DIR` and `CFITSIO_LIB_DIR`, respectively.
WITH_CFITSIO := 
CFITSIO_INC_DIR := 
CFITSIO_LIB_DIR := 

# Set `WITH_HDF5 := T` to enable HDF5 file format.
# The paths for libhdf5 header (hdf5.h) and library (libhdf5.{a,so,dylib})
# files can be set via `HDF5_INC_DIR` and `HDF5_LIB_DIR`, respectively.
WITH_HDF5 := 
HDF5_INC_DIR := 
HDF5_LIB_DIR := 


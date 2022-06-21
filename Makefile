# Load the custom compilation options.
include options.mk

# Set flag for the POSIX standard.
CFLAGS += -std=c99 -D_POSIX_C_SOURCE=200809 -Wall

# Set flags according to the custom options.
ifeq ($(strip $(WITH_MPI)), T)
  C_CMPLR = $(MPICC)
  CFLAGS += -DMPI
else
  C_CMPLR = $(CC)
endif

ifeq ($(strip $(WITH_OMP)), T)
  CFLAGS += -DOMP
  ifneq ($(OMP_FLAG),)
    CFLAGS += $(OMP_FLAG)
  else
    CFLAGS += -fopenmp
  endif
endif

ifeq ($(strip $(WITH_SIMD)), T)
  CFLAGS += -DWITH_SIMD
endif

ifeq ($(strip $(SINGLE_PREC)), T)
  CFLAGS += -DSINGLE_PREC
endif

ifeq ($(strip $(WITH_MU_ONE)), T)
  CFLAGS += -DWITH_MU_ONE
endif

ifeq ($(strip $(WITH_CFITSIO)), T)
  CFLAGS += -DWITH_CFITSIO
  LIBS += -lcfitsio
  ifneq ($(CFITSIO_INC_DIR),)
    INCL += -I$(CFITSIO_INC_DIR)
  endif
  ifneq ($(CFITSIO_LIB_DIR),)
    LIBS += -L$(CFITSIO_LIB_DIR)
  endif
endif

ifeq ($(strip $(WITH_HDF5)), T)
  CFLAGS += -DWITH_HDF5
  LIBS += -lhdf5
  ifneq ($(HDF5_INC_DIR),)
    INCL += -I$(HDF5_INC_DIR)
  endif
  ifneq ($(HDF5_LIB_DIR),)
    LIBS += -L$(HDF5_LIB_DIR)
  endif
endif

# Set the components to be compiled
TARGETS = FCFC_2PT_BOX FCFC_2PT

# Set mandatory options (do not edit below this line).
LIBS += -lm
INCL += -Isrc/fcfc/common -Isrc/io -Isrc/lib -Isrc/math -Isrc/tree
SRCS = $(wildcard src/fcfc/common/*.c src/*/*.c)


all: $(TARGETS)

FCFC_2PT: src/fcfc/2pt
	$(C_CMPLR) $(CFLAGS) -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

FCFC_2PT_BOX: src/fcfc/2pt_box
	$(C_CMPLR) $(CFLAGS) -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

clean:
	rm $(TARGETS)

CC = gcc
CFLAGS = -std=c99 -D_POSIX_C_SOURCE=200809 -Wall -O3 -march=native
LIBS = -lm
INCL = -Icommon
SRCS = $(wildcard common/*.c)
CFLAGS += -DBENCHMARK_TIMING

TARGETS = BENCHMARK_struct_serial BENCHMARK_struct_avx BENCHMARK_struct_single_serial BENCHMARK_struct_single_avx BENCHMARK_hist_serial BENCHMARK_hist_avx BENCHMARK_hist_avx2 BENCHMARK_hist_single_serial BENCHMARK_hist_single_avx BENCHMARK_hist_single_avx2

# Benchmarks with AVX512
#TARGETS += BENCHMARK_struct_avx512 BENCHMARK_struct_single_avx512 BENCHMARK_hist_avx512 BENCHMARK_hist_single_avx512

# Benchmarks for OpenMP strategy
#TARGETS += BENCHMARK_struct_omp_dynamic BENCHMARK_struct_omp_static


.PHONY: all
all: $(TARGETS)

BENCHMARK_struct_serial: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TREE_PREALLOC -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_struct_avx: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TREE_PREALLOC -DWITH_AVX -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_struct_avx512: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TREE_PREALLOC -DWITH_AVX512 -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_struct_single_serial: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TREE_PREALLOC -DSINGLE_PREC -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_struct_single_avx: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TREE_PREALLOC -DSINGLE_PREC -DWITH_AVX -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_struct_single_avx512: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TREE_PREALLOC -DSINGLE_PREC -DWITH_AVX512 -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<


BENCHMARK_hist_serial: histogram
	$(CC) $(CFLAGS) -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_avx: histogram
	$(CC) $(CFLAGS) -DWITH_AVX -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_avx2: histogram
	$(CC) $(CFLAGS) -DWITH_AVX2 -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_avx512: histogram
	$(CC) $(CFLAGS) -DWITH_AVX512 -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_single_serial: histogram
	$(CC) $(CFLAGS) -DSINGLE_PREC -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_single_avx: histogram
	$(CC) $(CFLAGS) -DSINGLE_PREC -DWITH_AVX -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_single_avx2: histogram
	$(CC) $(CFLAGS) -DSINGLE_PREC -DWITH_AVX2 -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_hist_single_avx512: histogram
	$(CC) $(CFLAGS) -DSINGLE_PREC -DWITH_AVX512 -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<


BENCHMARK_struct_omp_dynamic: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TIMING -DBENCHMARK_OMP_DYNAMIC -fopenmp -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

BENCHMARK_struct_omp_static: struct
	$(CC) $(CFLAGS) -DBENCHMARK_TIMING -DBENCHMARK_OMP_STATIC -fopenmp -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<


.PHONY: clean
clean:
	rm -f $(TARGETS)

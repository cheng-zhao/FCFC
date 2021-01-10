CC = gcc
CFLAGS = -std=c99 -O3 -march=native -Wall
LIBS = -lm
INCL = -Isrc -Iio -Ilib -Imath -Itree
SRCS = $(wildcard src/*.c io/*.c lib/*.c math/*.c tree/*.c)
TARGETS = FCFC_2PT FCFC_2PT_BOX

# Settings for OpenMP (comment the following line to disable OpenMP)
LIBS += -DOMP -fopenmp

all: $(TARGETS)

FCFC_2PT: src/2pt
	$(CC) $(CFLAGS) -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

FCFC_2PT_BOX: src/2pt_box
	$(CC) $(CFLAGS) -o $@ $(SRCS) $(wildcard $</*.c) $(LIBS) $(INCL) -I$<

clean:
	rm $(TARGETS)

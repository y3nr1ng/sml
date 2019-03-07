#---------------------------------------------------------------------------
#	For Windows 7 (without OpenMP, best performance)
#
# CC     := gcc
# CFLAGS := -O2 -Wall -fopenmp -static -I/home/thhsieh/win64/include
# TIFF   := -L/home/thhsieh/win64/lib -ltiff.dll -ljbig -ljpeg
# LAPACK := -llapack -lblas -lgfortran

#---------------------------------------------------------------------------
#	For MacOS
#
# CC     := gcc
# CFLAGS := -O2 -fopenmp -Wall -I/sw/include
# TIFF   := -L/sw/lib -ltiff
# LAPACK := -llapack -lblas

#---------------------------------------------------------------------------
#	For Linux
#
# MKLDIR := /opt/intel/mkl
# CC     := icc
# CFLAGS := -O2 -openmp -g -Wall -wd981
# TIFF   := -ltiff
# LAPACK := -L$(MKLDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5

MKLDIR := /opt/intel/mkl
CC     := gcc
CFLAGS := -O2 -fopenmp -g -Wall
TIFF   := -ltiff
LAPACK := -L$(MKLDIR)/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lgomp

#---------------------------------------------------------------------------
#	For debug purpose.
#
CFLAGS += -DDEBUG
MD5LIB := -lcrypto

#---------------------------------------------------------------------------
#	Code building rules.
#
BIN := pix
OBJ := main.o frame.o spot.o spotfit.o output.o frameIO.o nlinfit.o \
       reginal.o matrix.o solver.o debug.o
LIB := $(TIFF) $(LAPACK) $(MD5LIB) -lz -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $<

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) -o $(BIN) $(CFLAGS) $(OBJ) $(LIB)

clean:
	rm -f $(OBJ) $(BIN)

main.o:     main.c     pix.h
frame.o:    frame.c    pix.h
frameIO.o:  frameIO.c  pix.h
spot.o:     spot.c     pix.h
spotfit.o:  spotfit.c  pix.h
output.o:   output.c   pix.h
nlinfit.o:  nlinfit.c  pix.h
reginal.o:  reginal.c  pix.h
matrix.o:   matrix.c   pix.h
nlinfit.o:  nlinfit.c  pix.h
debug.o:    debug.c    pix.h

#---------------------------------------------------------------------------
#	For Windows 7 (without MPI and OpenMP, best performance)
#
# CC     := gcc
# MPICC  := gcc
# CFLAGS := -O2 -DUSE_OMP -Wall -fopenmp -static -I/home/thhsieh/win64/include
# TIFF   := -L/home/thhsieh/win64/lib -ltiff.dll -ljbig -ljpeg
# LAPACK := -llapack -lblas -lgfortran

#---------------------------------------------------------------------------
#	For MacOS without MPI
#
# CC     := gcc
# MPICC  := gcc
# CFLAGS := -O2 -DUSE_OMP -fopenmp -Wall -I/sw/include
# TIFF   := -L/sw/lib -ltiff
# LAPACK := -llapack -lblas

#---------------------------------------------------------------------------
#	For Linux without MPI
#
# CC     := icc
# MPICC  := icc
# CFLAGS := -O2 -DUSE_OMP -openmp -g -Wall -wd981
# TIFF   := -ltiff
# LAPACK := -L/opt/intel/Compiler/11.1/069/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5

CC     := gcc
MPICC  := gcc
CFLAGS := -O2 -DUSE_OMP -fopenmp -g -Wall
TIFF   := -ltiff -ljpeg

#---------------------------------------------------------------------------
#	For Linux with MPI
#
# CC     := icc
# MPICC  := mpiicc
# CFLAGS := -O3 -DUSE_MPI -openmp -Wall -wd981
# TIFF   := -ltiff
# LAPACK := -L/opt/intel/Compiler/11.1/069/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5

BIN := pix
SRC := main.c frame.c spot.c spotfit.c output.c \
       frameIO.c matIO.c nlinfit.c reginal.c matrix.c solver.c dsysv.c
OBJ := $(SRC:.c=.o)
LIB := $(TIFF) $(LAPACK) -lz -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $<

all: $(BIN)

$(BIN): $(OBJ)
	$(MPICC) -o $(BIN) $(CFLAGS) $(OBJ) $(LIB)

clean:
	rm -f $(OBJ) $(BIN)

main.o:     main.c     pix.h
	$(MPICC) -c $(CFLAGS) $<
frame.o:    frame.c    pix.h
frameIO.o:  frameIO.c  pix.h
spot.o:     spot.c     pix.h
spotfit.o:  spotfit.c  pix.h
output.o:   output.c   pix.h
nlinfit.o:  nlinfit.c  pix.h
matIO.o:    matIO.c    pix.h
reginal.o:  reginal.c  pix.h
matrix.o:   matrix.c   pix.h
nlinfit.o:  nlinfit.c  pix.h
dsysv.o:    dsysv.c    pix.h

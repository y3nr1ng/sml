CC     := gcc
CFLAGS := -O2 -fopenmp -Wall
# CC     := icc
# CFLAGS := -O3 -openmp -Wall -wd981

BIN := pclst
SRC := main.c spot.c pclst.c xcor.c
OBJ := $(SRC:.c=.o)
LIB := -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $<

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) -o $(BIN) $(CFLAGS) $(OBJ) $(LIB)

clean:
	rm -f $(OBJ) $(BIN)

main.o:     main.c     pclst.h
spot.o:     spot.c     pclst.h
pclst.o:    pclst.c    pclst.h
xcor.o:     xcor.c     pclst.h

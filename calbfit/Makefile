CC     := gcc
CFLAGS := -O2 -Wall

BIN := calbfit
SRC := main.c calb.c func.c nlinfit.c tool.c dsysv.c
OBJ := $(SRC:.c=.o)
LIB := -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $<

all: $(OBJ)
	$(CC) -o $(BIN) $(CFLAGS) $(OBJ) $(LIB)

main.o:    main.c    calb.h
calb.o:    calb.c    calb.h
func.o:    func.c    calb.h
nlinfit.o: nlinfit.c calb.h
tool.o:    tool.c    calb.h
dsysv.o:   dsysv.c   calb.h


clean:
	rm -f $(BIN) $(OBJ)

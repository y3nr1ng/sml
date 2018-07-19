ARCH   := LINUX
CC     := gcc
CFLAGS := -O2 -Wall

ifeq ($(ARCH),APPLE)
    JPEG_INC   := -I/sw/include
    JPEG_LIB   := -L/sw/lib -ljpeg
else
    JPEG_INC   := 
    JPEG_LIB   := -ljpeg
endif

SRC := main.c aimg.c jpeg.c output.c
OBJ := $(SRC:.c=.o)
EXE := aimg
INC := $(JPEG_INC)
LIB := $(JPEG_LIB) -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $<

all: $(EXE)

.PHONY: clean GLFONT

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(INC) -o $(EXE) $(OBJ) $(LIB)

main.o:   main.c   aimg.h
aimg.o:   aimg.c   aimg.h
jpeg.o:   jpeg.c   aimg.h
output.o: output.c aimg.h

clean:
	rm -f $(OBJ) $(EXE)

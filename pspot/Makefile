ARCH   := LINUX
CC     := gcc
CFLAGS := -O2 -Wall

ifeq ($(ARCH),APPLE)
    OPENGL_LIB := -framework GLUT -framework OpenGL
    OPENGL_INC :=
    JPEG_INC   := -I/sw/include
    JPEG_LIB   := -L/sw/lib -ljpeg
else
ifeq ($(ARCH),CYGWIN)
#    OPENGL_LIB := -lglut32 -lglu32 -lopengl32
    OPENGL_LIB := /usr/bin/glut32.dll -lglu32 -lopengl32
    OPENGL_INC := -I/usr/include/opengl
    JPEG_INC   := 
    JPEG_LIB   := -ljpeg
else
    OPENGL_LIB := -lglut -lGLU
    OPENGL_INC :=
    JPEG_INC   := 
    JPEG_LIB   := -ljpeg
endif
endif

SRC := main.c pspot.c xcor.c opengl.c jpeg.c
OBJ := $(SRC:.c=.o) glfont/glfont.o
EXE := pspot
INC := $(JPEG_INC) $(OPENGL_INC)
LIB := $(OPENGL_LIB) $(JPEG_LIB) -lm

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $<

all: GLFONT $(EXE)

.PHONY: clean GLFONT

GLFONT:
	ARCH=$(ARCH) make -C glfont

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(INC) -o $(EXE) $(OBJ) $(LIB)

main.o:   main.c   pspot.h
pspot.o:  pspot.c  pspot.h
xcor.o:   xcor.c   pspot.h
opengl.o: opengl.c pspot.h
jpeg.o:   jpeg.c   pspot.h

clean:
	make -C glfont clean
	rm -f $(OBJ) $(EXE)

CC     := gcc
CFLAGS := -O2 -Wall

ifeq ($(ARCH),APPLE)
    OPENGL_LIB := -framework GLUT -framework OpenGL
    OPENGL_INC :=
else
ifeq ($(ARCH),CYGWIN)
    OPENGL_LIB := -lglut32 -lopengl32
    OPENGL_INC := -I/usr/include/opengl
else
    OPENGL_LIB := -lglut -lGLU
    OPENGL_INC :=
endif
endif

SRC := glfont.c
OBJ := $(SRC:.c=.o)

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPENGL_INC) $(INC) $<

all: $(OBJ)

glfont.o: glfont.c glfont.h

clean:
	rm -f $(OBJ) $(EXE)

# Makefile for libfastell

SRC = fastell.f
OBJ = fastell.o

%.o: %.f
	/opt/homebrew/bin/gfortran -O5 -c -o $@ $^

all: libfastell.a

libfastell.a: $(OBJ)
	/opt/homebrew/bin/gcc-ar-15 r $@ $^

clean:
	rm -f *.o *.a


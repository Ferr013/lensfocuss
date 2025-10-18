COMPILER = /opt/homebrew/bin/gcc-15
SOURCE_LIBS = -Ilib/
OSX_OPT = -Llib/ -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL lib/libraylib.a
OSX_OUT = -o "bin/lensfocus"
EXTRA = -std=c99 -O3 -Wno-deprecated-non-prototype
CFILES = src/main.c

QFITSFILES = src/qfits/*.c
QFITS_OUT = -o "bin/qfits"

FASTELL_INC = -Ilib/fastell/
FASTELL_LIB = -Llib/fastell/ -lfastell -lgfortran

.PHONY = all
all: lensfocus run

lensfocus:
	$(COMPILER) $(QFITSFILES) $(CFILES) $(SOURCE_LIBS) $(FASTELL_INC) $(OSX_OUT) $(OSX_OPT) $(EXTRA) $(FASTELL_LIB)

qfits:
	$(COMPILER) $(QFITSFILES) $(SOURCE_LIBS) $(QFITS_OUT)

run:
	./bin/lensfocus

clean:
	rm -f bin/lensfocus
	rm -f bin/qfits
	rm -f src/*.o
	rm -f src/*.dSYM
	rm -f lib/*.o
	rm -f lib/*.dSYM

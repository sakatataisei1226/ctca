.PHONY: all install clean

PREFIX=..
CC=mpiicc
FC=mpiifort
CFLAGS=-I.
FFLAGS=-I.
LDFLAGS=-L.
LIBS=-lctca
F_LIBS=-lctca_f -lctca_mod

FILES = ctca ctca_f ctca_mod

OBJECTS = $(FILES:=.o)
TARGETS = $(FILES:%=lib%.a)

install: all
	/usr/bin/mkdir -p ${PREFIX}/lib
	/usr/bin/mkdir -p ${PREFIX}/include
	/bin/cp libctca.a libctca_f.a libctca_mod.a ${PREFIX}/lib
	/bin/cp ctca.h ctca.mod ${PREFIX}/include

all: ${TARGETS} ${OBJECTS}

lib%.a: %.o
	ar rc $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $^
%.o: %.f90
	$(FC) $(FFLAGS) -c $^

clean:
	/bin/rm -f -r *.o *.mod *.a ../include ../lib

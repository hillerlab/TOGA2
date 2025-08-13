# Copyright 2016 Peter Schwede

SILENCE=@

CC = cc
LD = cc

LIB = criterion
LIBLOC = ${HOME}/usr/local/lib/
INC = ${LIBLOC}../include

CFLAGS = -c -g -O0 -Wextra -Wall -DDEBUG -I../src/ -I${INC} --std=c11

LDLIBS = -L${LIBLOC} -l${LIB} -lm -lc

srcfiles := $(filter-out ../src/Cesar.c, $(wildcard ../src/*.c) $(wildcard *.c))
objects := $(srcfiles:%.c=%.o)


default: all


.PHONY: all
all: AllTests
	AllTests -j1
	#AllTests


.PHONY: valgrind
valgrind: AllTests
	valgrind -q --track-origins=yes --leak-check=yes AllTests


.PHONY: AllTests
AllTests: ${objects} ${srcfiles}
	${CC} -o $@ ${objects} ${LDLIBS}


.PHONY: clean
clean:
	${SILENCE}rm -f ${tests} *.o AllTests cesar_matrix.log

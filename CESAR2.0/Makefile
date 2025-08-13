CC=gcc
LDLIBS=-lm -lc
LDFLAGS=${LDLIBS}
FLAGS_DEBUG=-g -Wall -DDEBUG -DLOGLEVEL=0 -Wextra -Wno-unused-function
VERSION=-DVERSION="\"$(shell git log -n1 --pretty=format:%h%d) added cd_acc/do\""

#CFLAGS=-O0 -std=c11 ${VERSION} ${FLAGS_DEBUG}
#CFLAGS=-O2 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5 -std=c11 ${VERSION} ${FLAGS_DEBUG}
CFLAGS=-O3 -std=c11 ${VERSION}

srcfiles=$(wildcard src/*.c)
objects=$(srcfiles:%.c=%.o)


# compile (fast) production version
default: cesar
.PHONY: cesar
cesar: CESAR
	mv CESAR cesar

CESAR: ${objects}
	${CC} ${CFLAGS} -o $@ $^ ${LDLIBS}


# compile debugging version
.PHONY: debug
debug: ${objects}
	${CC} ${FLAGS_DEBUG} -o $@ $^ ${LDLIBS}


# dev setup
workspace:
	virtualenv venv


# dev doc
.PHONY: doc
doc:
	mkdir -p doc/
	doxygen doxygen.conf


.PHONY: test
test: ${objects} cesar
	tests.sh ./cesar


.PHONY: test/valgrind
test/valgrind: ${objects} cesar
	test/run.sh cesar
	make -C test valgrind

# cleaners
.PHONY: clean
clean: cleanCode cleanTest
	rm -rf doc
	mkdir doc

cleanCode:
	rm -fv *.o
	rm -fv src/*.o
	rm -fv src/*.gch
	rm -f cesar
	rm -f debug

cleanTest:
	make -i -C test clean

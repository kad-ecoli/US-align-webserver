CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=-static# -lm

PROGRAM=cif2pdb USalign pdb2fasta

all: ${PROGRAM}

cif2pdb: cif2pdb.cpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2fasta: pdb2fasta.cpp pstream.h basic_fun.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

USalign: USalign.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

USalign.zip: USalign
	zip -9r USalign.zip USalign readme.txt

clean:
	rm -f ${PROGRAM}

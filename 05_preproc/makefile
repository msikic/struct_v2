
CC=gcc
CXX=g++
CFLAGS=-O -Wall
CXXFLAGS=-O -Wall
CPPFLAGS= -Istride
LDFLAGS=  -Lstride -lstride -llapack

OBJECTS=main.o readpdb.o calcstride.o geomdescr.o linefit.o line.o vec.o mlexception.o


preproc:  ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} ${LDFLAGS} -o ../preproc

clean:
	rm -f ${OBJECTS}

linefittest: linefittest.cc linefit.cc vec.cc mlexception.cc\
 linefit.h vec.h mlexception.h
	gcc -g linefittest.cc linefit.cc vec.cc mlexception.cc\
 -llapack -o linefittest

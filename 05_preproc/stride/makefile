
LIBNAME=stride
LIBFILE=lib${LIBNAME}.a
LIBTEST=${LIBNAME}test

LIBOBJS=stride.o stride_cmd.o stride_pdb.o stride_chain.o\
 stride_util.o stride_phipsi.o stride_geom.o stride_map.o stride_hbond.o\
 stride_helix.o stride_sheet.o stride_sheet_mpat.o stride_sheet_fasn.o
OBJECTS=${LIBTEST}.o ${LIBOBJS}

CFLAGS=-O -Wall
CPPFLAGS=
LDFLAGS=-L. -l${LIBNAME} -lm

${LIBTEST}: makefile ${LIBTEST}.o ${LIBFILE}
	${CC} ${CFLAGS} ${LIBTEST}.o ${LDFLAGS} -o $@	

${LIBFILE}: makefile ${LIBOBJS}
	ar crs $@ ${LIBOBJS}

clean:
	rm -f ${OBJECTS} ${LIBTEST}

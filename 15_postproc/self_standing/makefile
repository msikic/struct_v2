
CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= gcc   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = pdbtfm

OPTM =  -O

VPATH = ./

PLATFORM = $(shell uname)



LIBS =  -lm   -llapack -lblas -lm

ifeq (${PLATFORM},Linux)
	PD =  -DLINUX
endif

DEFINES =  $(PD)  
SOURCES.F = 
SOURCES.f = 
SOURCES.c =   pdb_find_tfm.c utils.c


AUX = pdb.h utils.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)



$(NAME): $(OBJECTS)
	$(FCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.F.o:
	$(FCOMPL)   $(FFLAGS) -c   $(DEFINES)  $<  $(OPTM) $(INCLUDE)
.f.o:
	$(FCOMPL)  $(FFLAGS) -c   $(DEFINES)  $<  $(OPTM) $(INCLUDE)
.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 



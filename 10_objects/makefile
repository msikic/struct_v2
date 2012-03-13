CFLAGS= -m64   # options for C compiler
FFLAGS =
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall -fopenmp
CPLUSCOMPL = g++

HERE = $(shell pwd)

PREP_NAME = ../preproc
TEST_NAME = ../urchin_test
NAME = ../struct
OPTM =  -g


VPATH = ../02_parsers:../06_main:\
        ../04_geometric_match:../05_preproc:\
	../01_data_structures:../03_output:../07_utils:\
	../09_surface_fitting:../15_postproc


PLATFORM = $(shell uname)


LIBS =  -lm -lc  -L/usr/lib  -llapack  
INCLUDE = -I ../00_include



DEFINES =   # -DDMALLOC # -DDEBUG
SOURCES.F = 
SOURCES.f = 
PREP_SOURCES.c = struct_aa_code.c  struct_geometry.c\
		struct_geometric_descriptor.c struct_linalg.c  \
	        struct_preproc_out.c \
		struct_pp_ivana.c struct_read_pdb.c \
		struct_utils.c 

TEST_SOURCES.c = struct_urchin_test.c struct_quaternion.c\
		struct_utils.c

SOURCES.c =  struct.c struct_aa_code.c  struct_complement_match.c \
	     struct_descr_out.c	struct_geometry.c \
	     struct_linalg.c  struct_map.c struct_map2rotation.c \
	     struct_postproc.c struct_quaternion.c struct_read_pdb.c  \
	     struct_read_cmd_file.c struct_representation.c \
	     struct_recursive_output.c struct_scoring_fn.c \
	     struct_table.c  struct_table_input.c struct_input.c struct_utils.c	\
             hungarian.c 



AUX =  makefile  taxtrie.h  README \
       struct_comparison.h  struct_geometry.h  struct.h  struct_pdb.h  \
       struct_tree.h  struct_utils.h hungarian.h

PREP_SOURCES = $(PREP_SOURCES.c)  $(PREP_SOURCES.F) $(PREP_SOURCES.f)
PREP_OBJECTS =  $(PREP_SOURCES.F:.F=.o)  $(PREP_SOURCES.f:.f=.o)\
		 $(PREP_SOURCES.c:.c=.o)

TEST_SOURCES = $(TEST_SOURCES.c)  $(TEST_SOURCES.F) $(TEST_SOURCES.f)
TEST_OBJECTS = $(TEST_SOURCES.F:.F=.o)  $(TEST_SOURCES.f:.f=.o)\
		 $(TEST_SOURCES.c:.c=.o)

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

$(PREP_NAME): $(PREP_OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(PREP_OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)

$(TEST_NAME): $(TEST_OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(TEST_OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)



.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm  $(OBJECTS) 

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 



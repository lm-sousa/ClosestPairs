#######################################################################
## Makefile
#######################################################################

# program name
EXEC = gen_points
EXECDIR = aux

# Source
SRC = gen_points.f90

OBJ = $(SRC:.f90=.o)


#################################################
##
## READ COMPILATION FLAGS FROM
## makefile.cfg
##
#################################################
include ../makefile.cfg

# Compiler (GCC)
FC = gfortran
# flags moved to ../makefile.cfg

.PHONY: all clean

all: $(EXEC);
	ln -sf $(EXECDIR)/$(EXEC) ../$(EXEC)

$(EXEC) : $(OBJ) ;
	$(FC) $(FFLAGS) $^ -o $@

$(filter-out $(MODOBJ),%.o) : $(filter-out $(MODSRC),%.f90) $(MODOBJ) ;
	@echo "compiling object $<"
	$(FC) $(FFLAGS) -c $<

clean : ;
	rm -f *.o *.mod $(EXEC)
	rm -f ../$(EXEC)

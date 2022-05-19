#######################################################################
## Makefile
#######################################################################

# program name
EXEC = closest_mpi
EXECDIR = parallel

# Source
SRC = closest_mpi.f90

OBJ = $(SRC:.f90=.o)
MODDIR = ../mod
MODSRC = $(MODDIR)/closest_mod.f90 $(MODDIR)/sort_mod.f90
MODOBJ = closest_mod.o sort_mod.o


#################################################
##
## READ COMPILATION FLAGS FROM
## makefile.cfg
##
#################################################
include ../makefile.cfg

# Compiler (GCC)
FC = mpif90
# flags moved to ../makefile.cfg
FFLAGS += -fallow-argument-mismatch

.PHONY: all clean

all: $(EXEC);
	ln -sf $(EXECDIR)/$(EXEC) ../$(EXEC)

$(EXEC) : $(MODOBJ) $(OBJ) ;
	$(FC) $(FFLAGS) $^ -o $@

# module compilation
%mod.o : $(MODDIR)/%mod.f90 ;
	@echo "compiling module $<"
	$(FC) $(FFLAGS) -c $<

$(filter-out $(MODOBJ),%.o) : $(filter-out $(MODSRC),%.f90) $(MODOBJ) ;
	@echo "compiling object $<"
	$(FC) $(FFLAGS) -c $<

clean : ;
	rm -f *.o *.mod $(EXEC)
	rm -f ../$(EXEC)

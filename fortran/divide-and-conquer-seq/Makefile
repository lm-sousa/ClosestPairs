#######################################################################
## Makefile
#######################################################################

# program name
EXEC = closest_dc
EXECDIR = divide-and-conquer-seq

# Source
SRC = closest_dc.f90

# Other sources
# SRC +=

OBJ = $(SRC:.f90=.o)
MODDIR = ../mod
MODSRC = $(MODDIR)/closest_mod.f90 $(MODDIR)/sort_mod.f90
# MODOBJ = $(MODSRC:.f90=.o)
MODOBJ = closest_mod.o sort_mod.o


#################################################
##
## READ COMPILATION FLAGS FROM
## makefile.cfg
##
#################################################
include ../makefile.cfg

# Compiler (GCC)
FC = gfortran

.PHONY: all clean

all: $(EXEC) ;
	ln -sf $(EXECDIR)/$(EXEC) ../$(EXEC)

$(EXEC) : $(MODOBJ) $(OBJ) ;
	$(FC) $(FFLAGS) $^ -o $@

# module compilation
closest_mod.o : $(MODDIR)/closest_mod.f90;
	@echo "compiling module $<"
	$(FC) $(FFLAGS) -c $<

sort_mod.o : $(MODDIR)/sort_mod.f90;
	@echo "compiling module $<"
	$(FC) $(FFLAGS) -c $<

$(filter-out $(MODOBJ),%.o) : $(filter-out $(MODSRC),%.f90) $(MODOBJ) ;
	@echo "compiling object $<"
	$(FC) $(FFLAGS) -c $<

clean : ;
	rm -f *.o *.mod $(EXEC)
	rm -f ../$(EXEC)

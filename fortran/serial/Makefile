#######################################################################
## Makefile
#######################################################################

# program name
EXEC = closest_serial
EXECDIR = serial

# Source
SRC = closest_serial.f90

OBJ = $(SRC:.f90=.o)
MODDIR = ../mod
MODSRC = $(MODDIR)/closest_mod.f90
MODOBJ = closest_mod.o


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

$(EXEC) : $(MODOBJ) $(OBJ) ;
	$(FC) $(FFLAGS) $^ -o $@

# module compilation
closest_mod.o : $(MODDIR)/closest_mod.f90;
	@echo "compiling module $<"
	$(FC) $(FFLAGS) -c $<

$(filter-out $(MODOBJ),%.o) : $(filter-out $(MODSRC),%.f90) $(MODOBJ) ;
	@echo "compiling object $<"
	$(FC) $(FFLAGS) -c $<

clean : ;
	rm -f *.o *.mod $(EXEC) $(EXEC2) $(EXECMPI)
	rm -f ../$(EXEC)

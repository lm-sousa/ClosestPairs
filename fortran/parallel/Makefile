#######################################################################
## Makefile
#######################################################################

# program name
EXEC = closest_mpi
EXECDIR = parallel

# Source
SRC = closest_mpi.f90

# Other sources
# SRC +=

OBJ = $(SRC:.f90=.o)
MODDIR = ../mod
MODSRC = $(MODDIR)/closest_mod.f90 $(MODDIR)/sort_mod.f90 $(MODDIR)/mpi_mod.f90
# MODOBJ = $(MODSRC:.f90=.o)
MODOBJ = closest_mod.o sort_mod.o mpi_mod.o


# Compiler (GCC)
FC = mpif90
# move flags to general Makefile
FFLAGS = -O0

## Debug
FFLAGS += -g -fbacktrace
FFLAGS += -Wall -Wextra -Wconversion 
# FFLAGS += -pedantic
FFLAGS += -ffpe-trap=zero,overflow,underflow
FFLAGS += -fbounds-check
FFLAGS += -fallow-argument-mismatch 
# FFLAGS += -Wargument-mismatch -Werror

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

#######################################################################
## Makefile
#######################################################################

# program name
EXEC = gen_points
EXECDIR = aux

# Source
SRC = gen_points.f90

# Other sources
# SRC +=

OBJ = $(SRC:.f90=.o)
# MODDIR = ../mod
# MODSRC = $(MODDIR)/closest_mod.f90
# # MODOBJ = $(MODSRC:.f90=.o)
# MODOBJ = closest_mod.o


#################################################
##
## READ COMPILATION FLAGS FROM
## makefile.cfg
##
#################################################
include ../makefile.cfg

# Compiler (GCC)
# FC = gfortran
# flags moved to ../makefile.cfg
# FFLAGS += -O0
# ## Debug
# FFLAGS += -g -fbacktrace
# FFLAGS += -Wall -Wextra -Wconversion -pedantic
# FFLAGS += -ffpe-trap=zero,overflow,underflow
# FFLAGS += -fbounds-check

.PHONY: all clean

all: $(EXEC);
	ln -sf $(EXECDIR)/$(EXEC) ../$(EXEC)

$(EXEC) : $(OBJ) ;
	$(FC) $(FFLAGS) $^ -o $@

# module compilation
# closest_mod.o : $(MODDIR)/closest_mod.f90;
# 	@echo "compiling module $<"
# 	$(FC) $(FFLAGS) -c $<
# 
$(filter-out $(MODOBJ),%.o) : $(filter-out $(MODSRC),%.f90) $(MODOBJ) ;
	@echo "compiling object $<"
	$(FC) $(FFLAGS) -c $<

clean : ;
	rm -f *.o *.mod $(EXEC)
	# rm -f *__genmod.f90 *_genmod.mod
	rm -f ../$(EXEC)

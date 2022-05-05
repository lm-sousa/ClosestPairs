#######################################################################
## Makefile
#######################################################################

# program name
EXEC = closest_serial
EXEC2 = gen_points
EXECDIR = .

# Source
SRC = closest_serial.f90
SRC2 = gen_points.f90

# Other sources
# SRC +=

OBJ = $(SRC:.f90=.o)
MODSRC = closest_mod.f90
MODOBJ = $(MODSRC:.f90=.o)

# Compiler (GCC)
FC = gfortran
FFLAGS += -O0

## Debug
FFLAGS += -g -fbacktrace
FFLAGS += -Wall -Wextra -Wconversion -pedantic
FFLAGS += -ffpe-trap=zero,overflow,underflow
FFLAGS += -fbounds-check

.PHONY: all clean

all: $(EXEC) $(EXEC2)

# $(EXEC) : $(SRC) 
# 	$(FC) $(FFLAGS) $(SRC) -o $(EXEC)
# $(EXEC) : $(EXEC).o $(FC) $(FFLAGS) $^ -o $@
# 
$(EXEC) : $(MODOBJ) $(OBJ) ;
	$(FC) $(FFLAGS) $^ -o $@
# %.o : %.f90 ;
# 	$(FC) $(FFLAGS) -c $<
#
$(EXEC2) : $(SRC2) 
	$(FC) $(FFLAGS) $(SRC2) -o $(EXEC2)

# module compilation
closest_mod.o : closest_mod.f90;
	@echo "compiling module $<"
	$(FC) $(FFLAGS) -c $<

$(filter-out $(MODOBJ),%.o) : $(filter-out $(MODSRC),%.f90) $(MODOBJ) ;
	@echo "compiling object $<"
	$(FC) $(FFLAGS) -c $<

clean : ;
	rm -f *.o *.mod $(EXEC) $(EXEC2)
	# rm -f *__genmod.f90 *_genmod.mod
	# rm -f ../$(EXEC)

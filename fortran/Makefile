########################################################################
## Makefile
########################################################################

### SOURCES DIRECTORIES
SRCDIR1 = serial
SRCDIR2 = parallel
SRCDIR3 = aux

.PHONY : all clean

all: closest_serial closest_mpi gen_points ; 

closest_serial : ;
	rm -f closest_serial
	@cd $(SRCDIR1) ; make

closest_mpi : ;
	rm -f closest_mpi
	@cd $(SRCDIR2) ; make

gen_points : ;
	rm -f gen_points
	@cd $(SRCDIR3) ; make

clean_srcs : ;
	@cd $(SRCDIR1) ; make clean 
	@cd ../
	@cd $(SRCDIR2) ; make clean 
	@cd ../
	@cd $(SRCDIR3) ; make clean 

clean : clean_srcs ;

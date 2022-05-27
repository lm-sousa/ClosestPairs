########################################################################
## Makefile
########################################################################

### SOURCES DIRECTORIES
SRCDIR1 = serial
SRCDIR2 = parallel
SRCDIR3 = aux
SRCDIR4 = divide-and-conquer-seq
SRCDIR5 = divide-and-conquer-mpi

.PHONY : all clean

all: closest_bf closest_mpi gen_points closest_dc closest_dc_mpi; 

closest_bf : ;
	rm -f closest_bf
	@cd $(SRCDIR1) ; make

closest_mpi : ;
	rm -f closest_mpi
	@cd $(SRCDIR2) ; make

gen_points : ;
	rm -f gen_points
	@cd $(SRCDIR3) ; make

closest_dc : ;
	rm -f closest_dc
	@cd $(SRCDIR4) ; make

closest_dc_mpi : ;
	rm -f closest_dc_mpi
	@cd $(SRCDIR5) ; make

clean_srcs : ;
	@cd $(SRCDIR1) ; make clean 
	@cd ../
	@cd $(SRCDIR2) ; make clean 
	@cd ../
	@cd $(SRCDIR3) ; make clean 
	@cd ../
	@cd $(SRCDIR4) ; make clean 
	@cd ../
	@cd $(SRCDIR5) ; make clean 

clean : clean_srcs ;

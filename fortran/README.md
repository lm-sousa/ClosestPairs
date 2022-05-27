# ClosestPairs Fortran implementation

## Directory structure:

 - `aux/` : auxiliary programs (e.g., point generator);
 - `serial/` : serial (brute force) version of ClosestPairs;
 - `parallel/` : parallel version of ClosestPairs (requires `mpif90`);
 - `divide-and-conquer-seq` : sequential Divide-and-Conquer approach;
 - `divide-and-conquer-mpi` : parallel version of Divide-and-Conquer approach;

## Instructions:

Compilation (requires `gfortran` and an MPI wrapper `mpif90`):

 - `make` :  compiles all tools.
 - `make clean`: removes compiled objects and programs.
 
One can also compile separate tools by:

 - `make gen_points`
 - `make closest_bf`
 - `make closest_mpi`
 - `make closest_dc`
 - `make closest_dc_mpi`
 
Or heading into a specific directory and running `make`.
 
Compilation flags can be defined in makefile.cfg (flags for parallel code 
defined infividually).
 
## Checklist:

- [x] Serial (brute-force) code.
- [x] Random point generator.
- [x] Parallel (MPI) version.
- [x] Divide and conquer method (sequential)
- [x] Divide and conquer method (parallel)
- [ ] Parallel quicksort/mergesort

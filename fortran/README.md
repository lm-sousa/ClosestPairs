# ClosestPairs Fortran implementation

## Directory structure:

 - `aux/` : auxiliary programs (e.g., point generator);
 - `serial/` : serial (brute force) version of ClosestPairs;
 - `parallel/` : parallel version of ClosestPairs (requires `mpif90`).

## Instructions:

Compilation (requires `gfortran` and an MPI wrapper `mpif90`):

 - `make` :  compiles all tools.
 - `make clean`: removes compiled objects and programs.
 
One can also compile separate tools by:

 - `make closest_serial`
 - `make closest_mpi`
 
Or heading into `serial` or `parallel` directories and running `make`.
 
Compilation flags can be defined in base Makefile.
 
## Checklist:

- [x] Serial (brute-force) code.
- [x] Random point generator.
- [ ] Parallel (MPI) version.

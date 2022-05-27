#!/bin/sh
make clean && make
echo ""
echo "------------------BruteForce----------------------"
 ./closest_bf points.dat out.dat 
echo ""
echo "------------------BruteForceMPI-------------------"
 mpirun -n 4 ./closest_mpi points.dat out.dat
echo ""
echo "----------------------DaC-------------------------"
 ./closest_dc points.dat out.dat 
echo ""
echo "--------------------DaC-MPI-----------------------"
 mpirun -n 4 ./closest_dc_mpi points.dat out.dat 

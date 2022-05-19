!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO OBTAIN THE TWO CLOSEST POINTS FROM A LIST OF POINTS (2D)
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! REQUIREMENTS !!!
!! - [x] Brute force algorithm
!!   - [x] Allow for 1, 2 and 3D points to be specified
!!   - [x] Read points from a file specified as an argument
!! - [x] external code generates points randomly
!! - [x] Coordinates are DOUBLE.
!! - [x] Resulting points written on a separate file.
!!!!!!!!!!!!!!!!!!!

PROGRAM closest_serial

! USE iso_fortran_env
USE closest_mod  ! contains variable declarations, subroutines

IMPLICIT NONE

! variable declaration: all in closest_mod

! start measuring cputime
CALL measure_cpu_time(cput)

! read command line arguments
CALL parseArguments(infile, outfile)

! read from file 
np = 0
ndim = 2
CALL readPoints(infile, np, ndim, domain, all_pts, verb=0)

! calculate closest
ALLOCATE (p1(ndim))
ALLOCATE (p2(ndim))
p1 = 0.d0
p2 = 0.d0
mindist = 0.d0

cputlabel = 'getClosest'
CALL measure_cpu_time(cput1, msg=cputlabel)
CALL getClosest(np, ndim, all_pts, p1, p2, mindist, verb=0)
CALL measure_cpu_time(cput1, msg=cputlabel)

DEALLOCATE (all_pts)

! write results (stdout - optional)
CALL writeResultsToStd(mindist, ndim, p1, p2)

! write results (file)
CALL writeResultsToFile(outfile, ndim, p1, p2)

! output cputime
CALL measure_cpu_time(cput)

DEALLOCATE (p1)
DEALLOCATE (p2)

END PROGRAM closest_serial

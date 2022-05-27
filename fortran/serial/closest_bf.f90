!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO OBTAIN THE TWO CLOSEST POINTS FROM A LIST OF POINTS (2D)
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM closest_bf

! USE iso_fortran_env
USE closest_mod  ! contains variable declarations, subroutines

IMPLICIT NONE

! variable declaration: all in closest_mod

verb = 0  ! verbosity level (global to all subroutines)

! start measuring cputime
CALL measure_cpu_time(cput)

! read command line arguments
CALL parseArguments(infile, outfile)

! read from file 
np = 0
ndim = 2
CALL readPoints(infile, np, ndim, domain, all_pts)

! calculate closest
ALLOCATE (p1(ndim))
ALLOCATE (p2(ndim))
p1 = 0.d0
p2 = 0.d0
mindist = 0.d0

cputlabel = 'getClosest'
CALL measure_cpu_time(cput1, msg=cputlabel)
CALL getClosest(np, all_pts, p1, p2, mindist)
CALL measure_cpu_time(cput1, msg=cputlabel)

DEALLOCATE (all_pts)

! write results (stdout - optional)
CALL writeResultsToStd(mindist, p1, p2)

! write results (file)
CALL writeResultsToFile(outfile, p1, p2)

! output cputime
CALL measure_cpu_time(cput)

DEALLOCATE (p1)
DEALLOCATE (p2)

END PROGRAM closest_bf

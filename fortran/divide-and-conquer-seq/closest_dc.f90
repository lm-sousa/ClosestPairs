!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO OBTAIN THE TWO CLOSEST POINTS FROM A LIST OF POINTS (2D)
!!  DIVIDE AND CONQUER VERSION
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM closest_dc

USE closest_mod  ! contains variable declarations, subroutines
USE sort_mod  ! contains sort subroutines

IMPLICIT NONE

! variable declaration: see closest_mod
INTEGER, ALLOCATABLE, DIMENSION(:) :: x1_isort, x2_isort  
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:, :) :: x1_sort, x2_sort, dummy

verb = 0

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

! sort in x1
ALLOCATE (dummy(np, ndim))
ALLOCATE (x1_isort(np))
CALL quickSortMod (np, all_pts(:, 1), dummy, x1_isort)
ALLOCATE (x1_sort(np, ndim))
x1_sort = all_pts(x1_isort, :)
DEALLOCATE (x1_isort)
DEALLOCATE (all_pts)

! sort in x2
ALLOCATE (x2_isort(np))
CALL quickSortMod (np, x1_sort(:, 2), dummy, x2_isort)
ALLOCATE (x2_sort(np, ndim))
x2_sort = x1_sort(x2_isort, :)
DEALLOCATE (x2_isort)
DEALLOCATE (dummy)

! check sort
! DO i = 1, np-1
!    IF (x1_sort(i, 1).GT.x1_sort(i+1, 1)) THEN
!        WRITE(*,*) "NOT SORTED IN X1", i, x1_sort(i+1, 1)
!    END IF
! END DO
! 
! DO i = 1, np-1
!    IF (x2_sort(i, 2).GT.x2_sort(i+1, 2)) THEN
!        WRITE(*,*) "NOT SORTED IN X2", i
!    END IF
! END DO

! READ(*,*)


cputlabel = 'getClosestDC'
CALL measure_cpu_time(cput1, msg=cputlabel)
CALL closestPairDC (np, x1_sort, x2_sort, p1, p2, mindist, 1, 1)
CALL measure_cpu_time(cput1, msg=cputlabel)


! write results (stdout - optional)
CALL writeResultsToStd(mindist, p1, p2)

! write results (file)
CALL writeResultsToFile(outfile, p1, p2)

! output cputime
CALL measure_cpu_time(cput)

DEALLOCATE (p1)
DEALLOCATE (p2)

END PROGRAM closest_dc

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
! instrumentation
cput = 0
cpufile = 'cputime_dc.csv'

CALL CPU_TIME (cput(1))
! read command line arguments
CALL parseArguments(infile, outfile)
CALL CPU_TIME (cput(2))

! read from file 
np = 0
ndim = 2
CALL readPoints(infile, np, ndim, domain, all_pts)
CALL CPU_TIME (cput(3))

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

CALL CPU_TIME (cput(4))

CALL closestPairDC (np, x1_sort, x2_sort, p1, p2, mindist, 1, 1)
CALL CPU_TIME (cput(5))

! write results (stdout - optional)
CALL writeResultsToStd(mindist, p1, p2)
CALL CPU_TIME (cput(6))

! write results (file)
CALL writeResultsToFile(outfile, p1, p2)
CALL CPU_TIME (cput(7))

DEALLOCATE (p1)
DEALLOCATE (p2)

CALL CPU_TIME (cput(8))

! write to cputime.txt
cput = cput - cput(1)
cput = cput * 1000.  ! miliseconds

! stdout
WRITE(*,*) 
WRITE(*,*) "Time output (ms):"
WRITE(*, '(" Argument parsing: ", F12.3)') cput(2)
WRITE(*, '(" Read points:      ", F12.3)') cput(3) - cput(2)
WRITE(*, '(" Sort (x1, x2):    ", F12.3)') cput(4) - cput(3)
WRITE(*, '(" Closest Pair DaC: ", F12.3)') cput(5) - cput(4)
WRITE(*, '(" Write to stdout:  ", F12.3)') cput(6) - cput(5)
WRITE(*, '(" Write to file:    ", F12.3)') cput(7) - cput(6)
WRITE(*, '(" ** Total time:    ", F12.3, " **")') cput(8)

INQUIRE (FILE=cpufile, EXIST=fexist)
IF (.NOT.fexist) THEN
   OPEN(8, FILE=cpufile, STATUS='NEW', ACTION='WRITE')
   WRITE(8, *) "parseargs, readpts, sort, getclosest, writeStd, writeFile, total"
ELSE
   OPEN(8, FILE=cpufile, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
END IF

! WRITE MESSAGE
WRITE(8, 1080) cput(2), cput(3) - cput(2), cput(4) - cput(3), &
   cput(5) - cput(4), cput(6) - cput(5), cput(7) - cput(6), cput(8)
1080 FORMAT(F12.3, 6(', ', F12.3))

CLOSE(8)

END PROGRAM closest_dc

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

! instrumentation
cput = 0
cpufile = 'cputime_bf.csv'

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

CALL getClosest(np, all_pts, p1, p2, mindist)
CALL CPU_TIME (cput(4))

DEALLOCATE (all_pts)

! write results (stdout - optional)
CALL writeResultsToStd(mindist, p1, p2)
CALL CPU_TIME (cput(5))

! write results (file)
CALL writeResultsToFile(outfile, p1, p2)
CALL CPU_TIME (cput(6))

! output cputime

DEALLOCATE (p1)
DEALLOCATE (p2)
CALL CPU_TIME (cput(7))

! write to cputime.txt
cput = cput - cput(1)
cput = cput * 1000.  ! miliseconds

! stdout
WRITE(*,*) 
WRITE(*,*) "Time output (ms):"
WRITE(*, '(" Argument parsing: ", F12.3)') cput(2)
WRITE(*, '(" Read points:      ", F12.3)') cput(3) - cput(2)
WRITE(*, '(" getClosest:       ", F12.3)') cput(4) - cput(3)
WRITE(*, '(" Write to stdout:  ", F12.3)') cput(5) - cput(4)
WRITE(*, '(" Write to file:    ", F12.3)') cput(6) - cput(5)
WRITE(*, '(" ** Total time:    ", F12.3, " **")') cput(7)

! to file
INQUIRE (FILE=cpufile, EXIST=fexist)
IF (.NOT.fexist) THEN
   OPEN(8, FILE=cpufile, STATUS='NEW', ACTION='WRITE')
   WRITE(8, *) "parseargs, readpts, getclosest, writeStd, writeFile, total"
ELSE
   OPEN(8, FILE=cpufile, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
END IF

! WRITE MESSAGE
WRITE(8, 1080) cput(2), cput(3) - cput(2), cput(4) - cput(3), &
   cput(5) - cput(4), cput(6) - cput(5), cput(7)
1080 FORMAT(F12.3, 5(', ', F12.3))

CLOSE(8)

END PROGRAM closest_bf

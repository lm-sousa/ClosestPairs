!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO GENERATE LIST OF POINTS USING PSEUDO-RANDOM GENERATOR
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! REQUIREMENTS !!!
!! 1. Inputs: Number of dimensions, domain, number of points
!! 2. Outputs: same as above + list of points
!! 3. Coordinates must be DBLE. (output as unformatted?)
!!!!!!!!!!!!!!!!!!!

PROGRAM gen_points

USE iso_Fortran_env

IMPLICIT NONE

! variable declaration
INTEGER :: i
INTEGER :: np, ndim
INTEGER :: ioerr
CHARACTER(len=10) :: outfile = "points.dat"
!! REAL :: xmin, ymin, zmin, xmax, ymax, zmax
REAL, ALLOCATABLE, DIMENSION(:,:) :: domain
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: points

! variable initilisation
np = 0
ndim = 2  ! default
ioerr= 0

! get inputs
WRITE(*,*) "Number of dimensions (1, 2 or 3): "
READ(*,*) ndim
WRITE(*,*) "Number of points: "
READ(*,*) np

! domain definition
WRITE(*, *) "Domain definition: "
ALLOCATE (domain(ndim, 2))
! depending on ndim, check inputs
IF ((ndim.GE.1).AND.(ndim.LE.3)) THEN
   WRITE(*,*) "xmin, xmax"
   READ(*,*) domain(1, 1), domain(1, 2) 
   IF (ndim.GE.2) THEN
      WRITE(*,*) "ymin, ymax"
      READ(*,*) domain(2, 1), domain(2, 2)
   END IF
   IF (ndim.EQ.3) THEN
       WRITE(*,*) "zmin, zmax"
       READ(*,*) domain(3, 1), domain(3, 2)
   END IF
ELSE IF (ndim.GT.3) THEN
   WRITE(*,*) "Max. 3 dimensions. Exiting..."
   STOP
ELSE
   STOP
END IF

ALLOCATE (points(np, ndim))
points = 0.d0

WRITE(*, *) "Generating with pseudo-random routine, uniform distribution..."
CALL RANDOM_NUMBER(points)
! re-scale
WRITE(*, *) "Re-scaling points to domain..."
DO i=1,ndim
   points(:, i) = domain(i, 1) + (domain(i, 2) - domain(i, 1)) * points(:, i) 
END DO
! stdout: SELECT CASE (ndim)
!    CASE (1)
!       WRITE(*,101) points
!    CASE (2)
!       WRITE(*,102) points
!    CASE (3)
!       WRITE(*,103) points
! END SELECT stdout
! ! WRITE(*,100) points
! 101 FORMAT (F7.3)
! 102 FORMAT (F7.3, ', ', F7.3)
! 103 FORMAT (2(F7.3, ', '), F7.3)

! write to file, NOTE: .dat or unformatted?
WRITE(*,*) "Writing points to file..."
OPEN(UNIT=8, FILE=outfile, STATUS='REPLACE', ACTION='WRITE')
! header
WRITE(8,110) np
WRITE(8,111) ndim
WRITE(8,112) domain(1, :)
IF (ndim.GE.2) THEN
   WRITE(8,113) domain(2, :)
   IF (ndim.EQ.3) THEN
       WRITE(8,114) domain(3, :)
   END IF
END IF
110 FORMAT (I10, T20, '# np, number of points')
111 FORMAT (I4, T20, '# ndim, number of dimensions')
112 FORMAT (F7.2, ', ', F7.2, T20, '# xmin, xmax')
113 FORMAT (F7.2, ', ', F7.2, T20, '# ymin, ymax')
114 FORMAT (F7.2, ', ', F7.2, T20, '# zmin, zmax')

! points
DO i = 1, np
   fileout: SELECT CASE (ndim)
      CASE (1)
         WRITE(8,201) points(i, :)
      CASE (2)
         WRITE(8,202) points(i, :)
      CASE (3)
         WRITE(8,203) points(i, :)
   END SELECT fileout
END DO
201 FORMAT (E23.15E3)
202 FORMAT (E23.15E3, ', ', E23.15E3)
203 FORMAT (2(E23.15E3, ', '), E23.15E3)

STOP
END PROGRAM gen_points

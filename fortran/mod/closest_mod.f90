!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  MODULE WITH SERIAL SUBROUTINES FOR ClosestPoints
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE closest_mod
!
! Purpose:
!    To declare data for closest_serial and closest_mpi, and subroutines
!    for closest point calculation and I/O

USE iso_fortran_env  ! module for compiler-independent precision declaration

IMPLICIT NONE
SAVE

REAL(KIND=REAL64), PARAMETER :: LARGE = huge(LARGE)
   ! largest double precision value

INTEGER :: i, j, k
INTEGER :: np  ! number of points 
INTEGER :: ndim  ! number of dimensions

REAL :: cput, cput1 = -1.  ! cputime
REAL, ALLOCATABLE, DIMENSION(:,:) :: domain  ! point domain limits
REAL(KIND=REAL64) :: mindist  ! minimum distance
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: p1, p2  ! closest points
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: all_pts  ! list of points

! String
CHARACTER(len=80) :: infile, outfile  ! filenames for points I/O
CHARACTER(len=80) :: cputlabel  ! label for cputime
      

CONTAINS
   SUBROUTINE getClosest (npts, dims, points, pt1, pt2, min_dist, verb)
      ! Get the closest pair of points, calculate distances between all points.

      IMPLICIT NONE

      ! inputs
      INTEGER, INTENT(IN) :: npts
         ! number of points
      INTEGER, INTENT(IN) :: dims
         ! number of dimensions
      INTEGER, INTENT(IN), OPTIONAL :: verb
         ! verbosity level
      REAL(KIND=REAL64), DIMENSION(npts, dims), INTENT(IN) :: points
         ! array of points

      INTEGER, DIMENSION(2) :: closest_pts = 0
         ! indices (from array) of pair of closest points
      REAL(KIND=REAL64) :: dist
         ! distance between points

      REAL(KIND=REAL64), DIMENSION(dims), INTENT(OUT) :: pt1, pt2
         ! coordinates of points 1 and 2 of closest pair
      REAL(KIND=REAL64), INTENT(OUT) :: min_dist
         ! minimum distance (= distance between closest pair)

      pt1 = 0.d0
      pt2 = 0.d0
      min_dist = LARGE
      dist = 0.d0

      IF (verb.EQ.1) THEN
         WRITE(*,*) "Points: "
         DO i = 1, npts
            SELECT CASE (dims)
               CASE (1)
                  WRITE(*,101) points(i,:)
               CASE (2)
                  WRITE(*,102) points(i,:)
               CASE (3)
                  WRITE(*,103) points(i,:)
            END SELECT
         END DO
      END IF
      101 FORMAT (ES23.15E3)
      102 FORMAT (ES23.15E3, ', ', ES23.15E3)
      103 FORMAT (2(ES23.15E3, ', '), ES23.15E3)

      WRITE(*, *) "Calculating closest pair..."
      DO i = 1, npts
         DO j = i+1, npts
            dist = 0.d0
            DO k = 1, ndim
               dist = dist + (points(i,k) - points(j,k))**2
            END DO
            dist = SQRT(dist)
            dist = dist
            IF (dist.LT.min_dist) THEN
               min_dist = dist
               closest_pts(1) = i
               closest_pts(2) = j
            END IF
         END DO
      END DO

      IF (verb.EQ.1) THEN
         WRITE(*,*) "Distances: "
         WRITE(*,201) dist
         201 FORMAT (10(1X, F10.3))
      END IF

      pt1 = points(closest_pts(1), :)
      pt2 = points(closest_pts(2), :)

   END SUBROUTINE

   SUBROUTINE getClosestOld (npts, dims, points, pt1, pt2, min_dist, verb)
      ! Same as @see getClosest: get the closest pair of points, calculate 
      ! distances between all points.
      ! Store distances in array (this version is used for debugging).

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: npts
         ! number of points.
      INTEGER, INTENT(IN) :: dims
         ! number of dimensions.
      INTEGER, INTENT(IN) :: verb
         ! verbosity level (0, 1).
      REAL(KIND=REAL64), DIMENSION(npts, dims), INTENT(IN) :: points

      INTEGER, DIMENSION(2) :: closest_pts = 0
         ! indices (from array) of pair of closest points
      REAL(KIND=REAL64), DIMENSION(npts, npts) :: dist
         ! array with distances between points

      REAL(KIND=REAL64), DIMENSION(dims), INTENT(OUT) :: pt1, pt2
         ! coordinates of points 1 and 2 of closest pair
      REAL(KIND=REAL64), INTENT(OUT) :: min_dist
         ! minimum distance (= distance between closest pair)

      pt1 = 0.d0
      pt2 = 0.d0
      min_dist = LARGE
      dist = 0.d0

      IF (verb.EQ.1) THEN
         WRITE(*,*) "Points: "
         DO i = 1, npts
            SELECT CASE (dims)
               CASE (1)
                  WRITE(*,101) points(i,:)
               CASE (2)
                  WRITE(*,102) points(i,:)
               CASE (3)
                  WRITE(*,103) points(i,:)
            END SELECT
         END DO
      END IF
      101 FORMAT (ES23.15E3)
      102 FORMAT (ES23.15E3, ', ', ES23.15E3)
      103 FORMAT (2(ES23.15E3, ', '), ES23.15E3)

      ! Store distances on dist array (not needed, useful for debugging)
      DO i = 1, npts
         DO j = i+1, npts
            dist(i,j) = 0.d0
            DO k = 1, ndim
               dist(i,j) = dist(i,j) + (points(i,k) - points(j,k))**2
            END DO
            dist(i,j) = SQRT(dist(i,j))
            dist(j,i) = dist(i,j)
         END DO
      END DO

      IF (verb.EQ.1) THEN
         WRITE(*,*) "Distances: "
         WRITE(*,201) dist
         201 FORMAT (10(1X, F10.3))
      END IF

      closest_pts = MINLOC(dist)

      min_dist = dist(closest_pts(1), closest_pts(2))
      pt1 = points(closest_pts(1), :)
      pt2 = points(closest_pts(2), :)

   END SUBROUTINE

   SUBROUTINE parseArguments (arg1, arg2)
      ! Parse arguments used in the program's execution.

      INTEGER :: nargs, ioerr
         ! number of arguments, error flag

      CHARACTER(len=80), INTENT(OUT) :: arg1, arg2
         ! arguments 1 and 2

      nargs = 0
      ioerr = 0

      ! read command line arguments
      nargs = COMMAND_ARGUMENT_COUNT()
      IF (nargs.ne.2) THEN
          
         IF (nargs.eq.0) THEN
            WRITE(*,*) 'Must supply input and output filenames.'
         ELSE IF (nargs.eq.1) THEN
            WRITE(*,*) 'Missing either input or output file.'
         ELSE IF (nargs.gt.2) THEN
            WRITE(*,*) 'Too many arguments...'
         END IF
      
         WRITE(*,*) 'Execution: ./closest_brute <input_points.dat> '&
                    '<output_file.dat>'
         STOP 'Exiting...'
      END IF
      
      CALL GET_COMMAND_ARGUMENT(NUMBER=1, VALUE=arg1, STATUS=ioerr)
      IF (ioerr.NE.0) WRITE(*,*) 'Error on input file arg.'
      CALL GET_COMMAND_ARGUMENT(NUMBER=2, VALUE=arg2, STATUS=ioerr)
      IF (ioerr.NE.0) WRITE(*,*) 'Error on output file arg.'
      
      ! trim filenames
      arg1 = TRIM(arg1)
      arg2 = TRIM(arg2)

      WRITE(*, *) "Input file (points):        ", TRIM(arg1)
      WRITE(*, *) "Output file (closest pair): ", TRIM(arg2)

   END SUBROUTINE parseArguments

   SUBROUTINE readPoints (filename, npts, dims, dom, points, verb)
      ! Read points in input file, filename.

      IMPLICIT NONE

      INTEGER, INTENT(IN), OPTIONAL :: verb
         ! verbosity level
      CHARACTER(len=*), INTENT(IN) :: filename
         ! name of input file with list of points.

      INTEGER, INTENT(OUT) :: npts, dims
         ! number of points and dimensions
      REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: dom
         ! domain limits (min, max) for each dimension.
      REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: points
         !! array of points

      REAL :: cputime = -1.

      cputlabel = 'readPoints'
      CALL measure_cpu_time(cputime, cputlabel)
      
      npts = 0
      dims = 2

      WRITE(*,*) 'Reading points from ', TRIM(filename) ! try with formats
      
      ! file follows a specific convention
      OPEN(8, FILE=TRIM(filename), STATUS='OLD', ACTION='READ')
      READ(8,*) npts
      READ(8,*) dims

      ALLOCATE(dom(dims, 2))
      DO i=1, dims
         READ(8, *) dom(i, :)
      END DO
      
      ALLOCATE (points(npts, dims))
      points = 0.d0
      ! read points
      DO i = 1, npts
         READ(8,*) points(i,:)
      END DO
      CLOSE(8)
      
      WRITE(*,*) "Number of points: ", npts
      WRITE(*,*) "Number of dimensions: ", dims
      IF (verb.EQ.1) THEN
          WRITE(*,*) "Points: "
          DO i = 1, npts
             SELECT CASE (dims)
                CASE (1)
                   WRITE(*,101) points(i,:)
                CASE (2)
                   WRITE(*,102) points(i,:)
                CASE (3)
                   WRITE(*,103) points(i,:)
             END SELECT
          END DO
      END IF
      101 FORMAT (ES23.15E3)
      102 FORMAT (ES23.15E3, ', ', ES23.15E3)
      103 FORMAT (2(ES23.15E3, ', '), ES23.15E3)

      CALL measure_cpu_time(cputime, cputlabel)

   END SUBROUTINE

   SUBROUTINE writeResultsToStd (min_dist, dims, pt1, pt2)
      ! Write results (distance and closest pair) to stdout

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: dims
         ! number of dimensions.
      REAL(KIND=REAL64), INTENT(IN) :: min_dist
         ! distance of closest pair
      REAL(KIND=REAL64), DIMENSION(dims), INTENT(IN) :: pt1, pt2
         ! first and second point of closest pair

      REAL :: cputime = -1.

      cputlabel = 'writeToStd'
      CALL measure_cpu_time(cputime, cputlabel)

      WRITE(*,100) min_dist
      100 FORMAT (" Closest points' distance: ", F7.3)
      SELECT CASE (dims)
         CASE (1)
            WRITE(*,101) "p1:", pt1
            WRITE(*,101) "p2:", pt2
         CASE (2)
            WRITE(*,102) "p1:", pt1
            WRITE(*,102) "p2:", pt2
         CASE (3)
            WRITE(*,103) "p1:", pt1
            WRITE(*,103) "p2:", pt2
      END SELECT
      101 FORMAT (1X, A, 1X, ES23.15E3)
      102 FORMAT (1X, A, 1X, ES23.15E3, ', ', ES23.15E3)
      103 FORMAT (1X, A, 1X, 2(ES23.15E3, ', '), ES23.15E3)

      CALL measure_cpu_time(cputime, cputlabel)

   END SUBROUTINE writeResultsToStd

   SUBROUTINE writeResultsToFile (filename, dims, pt1, pt2)
      ! Write results (distance and closest pair) to output file

      IMPLICIT NONE

      CHARACTER(len=*) :: filename
         ! name of output file
      INTEGER, INTENT(IN) :: dims
         ! number of dimensions.
      REAL(KIND=REAL64), DIMENSION(dims), INTENT(IN) :: pt1, pt2
         ! first and second point of closest pair

      REAL :: cputime = -1.0

      cputlabel = 'writeToFile'
      CALL measure_cpu_time(cputime, cputlabel)

      filename = TRIM(filename)

      OPEN(9, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
         SELECT CASE (dims)
            CASE (1)
               WRITE(9,101) pt1
               WRITE(9,101) pt2
            CASE (2)
               WRITE(9,102) pt1
               WRITE(9,102) pt2
            CASE (3)
               WRITE(9,103) pt1
               WRITE(9,103) pt2
         END SELECT
      CLOSE(9)
      101 FORMAT (E23.15E3)
      102 FORMAT (E23.15E3, ', ', E23.15E3)
      103 FORMAT (2(E23.15E3, ', '), E23.15E3)

      CALL measure_cpu_time(cputime, cputlabel)

   END SUBROUTINE writeResultsToFile

   SUBROUTINE measure_cpu_time (cputime0, msg)
      IMPLICIT NONE

      REAL, INTENT(INOUT) :: cputime0
      REAL :: cputime1, diff

      LOGICAL :: fexist
      CHARACTER(len=80), OPTIONAL :: msg

      IF (0>cputime0) THEN
         CALL CPU_TIME ( cputime0 )
      ELSE
         CALL CPU_TIME ( cputime1 )
         diff = cputime1 - cputime0
         IF (PRESENT(msg)) THEN
            WRITE(*,1000) TRIM(msg), diff
         ELSE
            WRITE(*, 1003) diff
         END IF

         INQUIRE (FILE='cputime.txt', EXIST=fexist)
         IF (.NOT.fexist) THEN
            OPEN(8, FILE='cputime.txt', STATUS='NEW', ACTION='WRITE')
         ELSE
            OPEN(8, FILE='cputime.txt', STATUS='OLD', ACTION='WRITE')
         END IF
         IF (PRESENT(msg)) THEN
            WRITE(8, 1001) TRIM(msg), diff
         ELSE
            WRITE(8, 1002) diff
         END IF
         CLOSE(8)
      ENDIF
      1000 FORMAT (" ** ", A, ": Took ", F5.2, " seconds.")
      1001 FORMAT (A, ", ", F5.2)
      1002 FORMAT ("main, ", F5.2)
      1003 FORMAT ("* Took ", F5.2, " seconds.")

      RETURN
   END SUBROUTINE measure_cpu_time

        
END MODULE closest_mod

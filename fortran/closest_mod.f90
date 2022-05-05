!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  MODULE WITH SERIAL SUBROUTINES FOR ClosestPoints
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! CONTAINS !!!
!! 1. Variable declarations and initializations
!! 2. Subroutine for distance calculation and closest pair determination
!!!!!!!!!!!!!!!!!!!

MODULE closest_mod

USE iso_fortran_env

IMPLICIT NONE
SAVE

! variable declaration
INTEGER :: i, j, k
INTEGER :: np, ndim
REAL(KIND=REAL64) :: mindist
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: p1, p2
REAL, ALLOCATABLE, DIMENSION(:,:) :: domain
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: all_pts

CHARACTER(len=80) :: infile, outfile  ! filenames for points I/O


CONTAINS
    SUBROUTINE getClosest (npts, dims, points, pt1, pt2, min_dist, verb)

       IMPLICIT NONE

       ! inputs
       INTEGER, INTENT(IN) :: npts, dims, verb
       REAL(KIND=REAL64), DIMENSION(np, ndim), INTENT(IN) :: points

       ! intermediate
       INTEGER, DIMENSION(2) :: closest_pts = 0
       REAL(KIND=REAL64), DIMENSION(np, np) :: dist

       ! outputs
       REAL(KIND=REAL64), DIMENSION(dims), INTENT(OUT) :: pt1, pt2
       REAL(KIND=REAL64), INTENT(OUT) :: min_dist

       pt1 = 0.d0
       pt2 = 0.d0
       min_dist = 0.d0
       dist = 9999.d0

       ! Store distances on dist array (not needed, useful for debugging)
       ! NOTE: Could repeatedly evaluate minimum, at each i iteration
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

        INTEGER :: nargs, ioerr

        CHARACTER(len=80), INTENT(OUT) :: arg1, arg2

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
        
           WRITE(*,*) 'Execution: ./closest_brute <input_points.dat> <output_file.dat>'
           STOP 'Exiting...'
        END IF
        
        CALL GET_COMMAND_ARGUMENT(NUMBER=1, VALUE=arg1, STATUS=ioerr)
        IF (ioerr.NE.0) WRITE(*,*) 'Error on input file arg.'
        CALL GET_COMMAND_ARGUMENT(NUMBER=2, VALUE=arg2, STATUS=ioerr)
        IF (ioerr.NE.0) WRITE(*,*) 'Error on output file arg.'
        
        ! trim filenames
        arg1 = TRIM(arg1)
        arg2 = TRIM(arg2)

    END SUBROUTINE parseArguments

    SUBROUTINE readPoints (filename, npts, dims, dom, points, verb)

        IMPLICIT NONE

        ! input
        INTEGER, INTENT(IN), OPTIONAL :: verb
        CHARACTER(len=*), INTENT(IN) :: filename

        ! output
        INTEGER, INTENT(OUT) :: npts, dims
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: dom
        REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: points
        

        npts = 0
        dims = 2
        ! filename = TRIM(filename)
        ! read from file 
        WRITE(*,*) 'Reading points from ', TRIM(filename) ! try with formats
        
        OPEN(8, FILE=TRIM(filename), STATUS='OLD', ACTION='READ')
        READ(8,*) npts  ! number of points
        READ(8,*) dims  ! number of dimensions

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
    END SUBROUTINE

    !!! write results to stdout
    SUBROUTINE writeResultsToStd (min_dist, dims, pt1, pt2)

        IMPLICIT NONE

        !inputs
        INTEGER, INTENT(IN) :: dims
        REAL(KIND=REAL64), INTENT(IN) :: min_dist
        REAL(KIND=REAL64), DIMENSION(dims), INTENT(IN) :: pt1, pt2

        WRITE(*,100) min_dist
        100 FORMAT ("Closest points' distance: ", F7.3)
        SELECT CASE (dims)
           CASE (1)
              WRITE(*,101) p1
              WRITE(*,101) p2
           CASE (2)
              WRITE(*,102) p1
              WRITE(*,102) p2
           CASE (3)
              WRITE(*,103) p1
              WRITE(*,103) p2
        END SELECT
        101 FORMAT (ES23.15E3)
        102 FORMAT (ES23.15E3, ', ', ES23.15E3)
        103 FORMAT (2(ES23.15E3, ', '), ES23.15E3)

    END SUBROUTINE writeResultsToStd

    !!! write results to file
    SUBROUTINE writeResultsToFile (filename, dims, pt1, pt2)

       IMPLICIT NONE

       ! inputs
       CHARACTER(len=*) :: filename
       INTEGER, INTENT(IN) :: dims
       REAL(KIND=REAL64), DIMENSION(dims), INTENT(IN) :: pt1, pt2

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

    END SUBROUTINE writeResultsToFile
        
END MODULE closest_mod

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
INTEGER :: ndim  ! number of dimensions, global used for all subroutines
INTEGER :: verb ! verbosity level

REAL :: cput, cput1 = -1.  ! cputime
REAL, ALLOCATABLE, DIMENSION(:,:) :: domain  ! point domain limits
REAL(KIND=REAL64) :: mindist  ! minimum distance
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: p1, p2  ! closest points
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: all_pts  ! list of points

! String
CHARACTER(len=80) :: infile, outfile  ! filenames for points I/O
CHARACTER(len=80) :: cputlabel  ! label for cputime
      

CONTAINS
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
      ! arg1 = TRIM(arg1)
      ! arg2 = TRIM(arg2)

      WRITE(*, *) "Input file (points):        ", TRIM(arg1)
      WRITE(*, *) "Output file (closest pair): ", TRIM(arg2)

   END SUBROUTINE parseArguments

   SUBROUTINE readPoints (filename, npts, dims, dom, points)
      ! Read points in input file, filename.

      IMPLICIT NONE

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

   SUBROUTINE writeResultsToStd (min_dist, pt1, pt2)
      ! Write results (distance and closest pair) to stdout

      IMPLICIT NONE

      ! INTEGER, INTENT(IN) :: dims
         ! number of dimensions.
      REAL(KIND=REAL64), INTENT(IN) :: min_dist
         ! distance of closest pair
      REAL(KIND=REAL64), DIMENSION(ndim), INTENT(IN) :: pt1, pt2
         ! first and second point of closest pair

      REAL :: cputime = -1.

      cputlabel = 'writeToStd'
      CALL measure_cpu_time(cputime, cputlabel)

      WRITE(*,100) min_dist
      100 FORMAT (" Closest points' distance: ", ES23.15E3)
      SELECT CASE (ndim)
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

   SUBROUTINE writeResultsToFile (filename, pt1, pt2)
      ! Write results (distance and closest pair) to output file

      IMPLICIT NONE

      CHARACTER(len=*) :: filename
         ! name of output file
      ! INTEGER, INTENT(IN) :: dims
         ! number of dimensions.
      REAL(KIND=REAL64), DIMENSION(ndim), INTENT(IN) :: pt1, pt2
         ! first and second point of closest pair

      REAL :: cputime = -1.0

      cputlabel = 'writeToFile'
      CALL measure_cpu_time(cputime, cputlabel)

      filename = TRIM(filename)

      OPEN(9, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
         SELECT CASE (ndim)
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

   SUBROUTINE getClosest (npts, points, pt1, pt2, min_dist)
      ! Get the closest pair of points, calculate distances between all points.

      IMPLICIT NONE

      ! inputs
      INTEGER, INTENT(IN) :: npts
         ! number of points
      REAL(KIND=REAL64), DIMENSION(npts, ndim), INTENT(IN) :: points
         ! array of points

      INTEGER, DIMENSION(2) :: closest_pts = 0
         ! indices (from array) of pair of closest points
      REAL(KIND=REAL64) :: dist
         ! distance between points
      ! REAL(KIND=REAL64) :: ldist
         ! distance between points (single direction)

      REAL(KIND=REAL64), DIMENSION(ndim), INTENT(OUT) :: pt1, pt2
         ! coordinates of points 1 and 2 of closest pair
      REAL(KIND=REAL64), INTENT(OUT) :: min_dist
         ! minimum distance (= distance between closest pair)

      pt1 = 0.d0
      pt2 = 0.d0
      min_dist = LARGE

      dist = 0.d0
      ! ldist = 0.d0

      IF (verb.EQ.1) THEN
         WRITE(*,*) "Points: "
         DO i = 1, npts
            SELECT CASE (ndim)
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

      IF (verb.EQ.1) WRITE(*, *) "Calculating closest pair..."
      DO i = 1, npts
         DO j = i+1, npts
            dist = calcDist (ndim, points(i,:), points(j,:))
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

   RECURSIVE SUBROUTINE closestPairDC (npts, px1, px2, pt1, pt2, &
           min_dist, level, sublevel)
      ! Get the closest pair of points through a divide and conquer approach

      IMPLICIT NONE

      ! inputs
      INTEGER, INTENT(IN) :: npts
         ! number of points
      INTEGER, INTENT(IN) :: level, sublevel
         ! divide and conquer level
      REAL(KIND=REAL64), DIMENSION(npts, ndim), INTENT(IN) :: px1, px2
         ! array of points, sorted in x1 and x2

      INTEGER :: imid, il, ir
         ! index of mid-point
      ! NOTE: replace with arrays?
      REAL(KIND=REAL64) :: ldist, rdist, sdist, L = 0.d0
         ! distance between points (global)
      REAL(KIND=REAL64), DIMENSION(ndim) :: pr, qr, pl, ql, pm, qm
         ! p, q closest points from left and right domains
      REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:, :) :: px2_r, px2_l

      REAL(KIND=REAL64), DIMENSION(ndim), INTENT(OUT) :: pt1, pt2
         ! coordinates of points 1 and 2 of closest pair
      REAL(KIND=REAL64), INTENT(OUT) :: min_dist
         ! minimum distance (= distance between closest pair)

      IF (npts.LE.3) THEN
         CALL getClosest(npts, px1, pt1, pt2, min_dist)
         RETURN
      END IF

      ! partition P at midpoint in x1
      imid = npts / 2 
      ! middle point (x1 from imid, or between imid and imid+1?)
      ! L = points(imid, 1)
      L = (px1(imid, 1) + px1(imid+1, 1)) / 2.

      ! divide points in px2 (x2-sorted array)
      ALLOCATE(px2_l(imid, ndim))
      ALLOCATE(px2_r(npts-imid, ndim))
      il = 0
      ir = 0

      DO i = 1, npts
         IF ((px2(i, 1).LT.L).AND.(il.LT.imid)) THEN
            il = il + 1
            px2_l(il, :) = px2(i, :)
         ELSE
            ir = ir + 1
            px2_r(ir, :) = px2(i, :)
         END IF
      END DO

      ! left side
      CALL closestPairDC(imid, px1(1:imid,:), px2_l, pl, ql, &
          ldist, level+1, sublevel*2-1)
      ! right side
      CALL closestPairDC(npts - imid, px1(imid+1:,:), px2_r, pr, qr, &
          rdist, level+1, sublevel*2)

      ! re-calculate line to avoid modifications from lower levels
      ! L = points(imid, 1)
      L = (px1(imid, 1) + px1(imid+1, 1)) / 2.
      ! delta
      min_dist = MIN(ldist, rdist)
      ! pm, qm (strip)
      CALL closestPairStrip (npts, px2, L, min_dist, pm, qm, &
                             sdist)

      ! return closest pair
      IF (sdist.LT.min_dist) THEN
          pt1 = pm
          pt2 = qm
          min_dist = sdist
      ELSE IF (rdist.LT.ldist) THEN
          pt1 = pr 
          pt2 = qr
      ELSE
          pt1 = pl 
          pt2 = ql
      END IF

   END SUBROUTINE closestPairDC

   SUBROUTINE partitionLR (points, npoints, nl, nr, mid, left, &
                           right)

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: npoints, nl, nr
      REAL(KIND=REAL64), INTENT(IN) :: mid
      REAL(KIND=REAL64), DIMENSION(npoints, ndim), INTENT(IN) :: points

      INTEGER :: il, ir

      REAL(KIND=REAL64), DIMENSION(nl, ndim), INTENT(OUT) :: left
      REAL(KIND=REAL64), DIMENSION(nr, ndim), INTENT(OUT) :: right

      il = 0
      ir = 0
      DO i = 1, npoints
         IF ((points(i, 1).LT.mid).AND.(il.LT.nl)) THEN
            il = il + 1
            left(il, :) = points(i, :)
         ELSE
            ir = ir + 1
            right(ir, :) = points(i, :)
         END IF
      END DO
   END SUBROUTINE

   SUBROUTINE closestPairStrip (npoints, points, strip_loc, dlt,&
           ps, qs, strip_min)

      IMPLICIT NONE

      ! inputs
      INTEGER, INTENT(IN) :: npoints
         ! number of points
      REAL(KIND=REAL64), DIMENSION(npoints, ndim), INTENT(IN) :: points
         ! array of points
      REAL(KIND=REAL64), INTENT(IN) :: strip_loc, dlt
         ! strip location in x1
      

      INTEGER :: np_strip
         ! number of points in opposite and same side of strip
      REAL(KIND=REAL64) :: dist
         ! distance between points (global)
      REAL(KIND=REAL64), DIMENSION(npoints, ndim) :: strip
         ! points in strip

      REAL(KIND=REAL64), DIMENSION(ndim), INTENT(OUT) :: ps, qs
         ! coordinates of points 1 and 2 of closest pair
      REAL(KIND=REAL64), INTENT(OUT) :: strip_min
         ! minimum distance (= distance between closest pair)

      strip_min = LARGE
      strip(:, :) = LARGE
      ps = 0.d0
      qs = 0.d0
      np_strip = 0
      ! get strip points NOTE: Subroutine?
      DO i = 1, npoints
         IF (ABS(points(i, 1) - strip_loc).LT.dlt) THEN
            np_strip = np_strip + 1
            strip(np_strip, :) = points(i, :)
         END IF
      END DO

      CALL getClosestStrip(np_strip, strip(:np_strip, :), strip_loc, dlt, &
         ps, qs, strip_min)

   END SUBROUTINE closestPairStrip

   SUBROUTINE getClosestStrip (strip_n, strip_pts, strip_x1, max_dlt,&
           p1s, p2s, mindist_strip)

      IMPLICIT NONE

      ! inputs
      INTEGER, INTENT(IN) :: strip_n
         ! number of points
      REAL(KIND=REAL64), DIMENSION(strip_n, ndim), INTENT(IN) :: strip_pts
         ! array of points
      REAL(KIND=REAL64), INTENT(IN) :: strip_x1, max_dlt
         ! strip location in x1
      

      INTEGER :: np_opp, np_same = 0
         ! number of points in opposite and same side of strip
      REAL(KIND=REAL64) :: dist, dx1, dx2
         ! distance between points (global)
      REAL(KIND=REAL64), DIMENSION(ndim), INTENT(OUT) :: p1s, p2s
         ! coordinates of points 1 and 2 of closest pair
      REAL(KIND=REAL64), INTENT(OUT) :: mindist_strip
         ! minimum distance (= distance between closest pair)

      dist = LARGE
      mindist_strip = LARGE
      p1s = 0.d0
      p2s = 0.d0

      IF (strip_n.LT.2) THEN ! no need to calculate
         RETURN
      ELSE IF (strip_n.EQ.2) THEN ! needed?
         mindist_strip = calcDist (ndim, strip_pts(1, :), strip_pts(2, :))
         p1s = strip_pts(1, :)
         p1s = strip_pts(2, :)
      ELSE
         ! compare all points (for debug, or 3D - not needed for 2D)
         IF (ndim.NE.2) THEN
         CALL getClosest(strip_n, strip_pts, p1s, p2s, mindist_strip)
         ELSE
            ! NOTE: CLOSEST STRIP 2D
            DO i=1, strip_n - 1
               np_opp = 0  ! opposite side of L
               np_same = 0  ! same side of L
               dx1 = 0.d0
               dx2 = 0.d0
               dist = LARGE
               DO j=i+1, strip_n
                  ! exit, skip conditions (max. comparisons)
                  IF ((np_opp.EQ.4).AND.(np_same.LT.3)) CONTINUE
                  IF ((np_opp.LT.4).AND.(np_same.EQ.3)) CONTINUE
                  IF ((np_opp.EQ.4).AND.(np_same.EQ.3)) EXIT

                  dx2 = ABS(strip_pts(j, 2) - strip_pts(i, 2)) 
                  IF (dx2.GT.max_dlt) THEN
                     EXIT
                  ELSE 
                     dx1 = ABS(strip_pts(i, 1) - strip_pts(j, 1))
                     ! if larger than delta, skip
                     IF (dx1.LT.max_dlt) THEN
                         ! calculate distances
                         dist = SQRT(dx1*dx1 + dx2*dx2)
                         ! if on opposite side
                         IF (((strip_pts(i, 1) - strip_x1)*&
                            (strip_x1 - strip_pts(j, 1))).GT.0) THEN
                            np_opp = np_opp + 1 ! increase counter of opposite points
                         ! if on same side
                         ELSE
                             np_same = np_same + 1
                         END IF
                     END IF

                     IF (dist.LT.mindist_strip) THEN
                        mindist_strip = dist
                        p1s = strip_pts(i, :)
                        p2s = strip_pts(j, :)
                     END IF
                  END IF
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE getClosestStrip


   REAL(KIND=REAL64) FUNCTION calcDist (nd, point1, point2)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nd  ! number of dimensions
      REAL(KIND=REAL64), DIMENSION(nd), INTENT(IN) :: point1, point2

      REAL(KIND=REAL64) :: diff

      calcDist = 0.d0
      DO k = 1, nd
         diff = point1(k) - point2(k)
         calcDist = calcDist + diff * diff
      END DO
      calcDist = SQRT(calcDist)

   END FUNCTION calcDist


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
      1000 FORMAT (" ** ", A, ": Took ", F9.6, " seconds.")
      1001 FORMAT (A, ", ", F9.6)
      1002 FORMAT ("main, ", F6.6)
      1003 FORMAT ("* Took ", F9.6, " seconds.")

      RETURN
   END SUBROUTINE measure_cpu_time

END MODULE closest_mod

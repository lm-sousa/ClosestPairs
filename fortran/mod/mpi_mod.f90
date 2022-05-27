!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  MODULE WITH SERIAL SUBROUTINES FOR ClosestPoints
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mpi_mod
!
! Purpose:
!    To declare data for closest_serial and closest_mpi, and subroutines
!    for closest point calculation and I/O

USE iso_fortran_env  ! module for compiler-independent precision declaration
USE closest_mod

IMPLICIT NONE
SAVE

INCLUDE 'mpif.h'  ! subroutines for MPI

! NOTE: mpi variable declarations (separate module?)
INTEGER :: rank                       ! current rank
INTEGER :: nranks                     ! number of ranks
INTEGER :: source                     ! var for sources
INTEGER :: tgt                        ! var for sources
INTEGER :: msgsize                    ! message sizes
INTEGER :: tag                        ! communication tags
INTEGER :: ierr                       ! error flag
INTEGER :: mpi_stat(MPI_STATUS_SIZE)  ! status array
INTEGER :: root 
   ! root rank

CONTAINS

   RECURSIVE SUBROUTINE closestPairDC_MPI (npts, px1, px2, pt1, pt2, &
           min_dist, level, sublevel, available_ranks)
      ! Get the closest pair of points through a divide and conquer approach
   
      ! USE closest_mod
   
      IMPLICIT NONE
   
      ! inputs
      INTEGER, INTENT(IN) :: npts
         ! number of points
      ! INTEGER, INTENT(IN) :: dims
         ! number of dimensions
      INTEGER, INTENT(IN) :: level, sublevel, available_ranks
         ! divide and conquer level
      REAL(KIND=REAL64), DIMENSION(npts, ndim), INTENT(IN) :: px1, px2
         ! array of points, sorted in x1 and x2
   
      INTEGER :: imid, il, ir
         ! index of mid-point
      INTEGER :: ldivs, rdivs
   
      ! MPI
      INTEGER :: tgt_rank
      INTEGER :: commtag
      INTEGER :: src
      INTEGER :: msglen
   
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

      il = 0
      ir = 0
      ldivs = 0
      rdivs = 0
      imid = npts/2  ! first division
   
      ! exec brute force (with px1 or px2)
      IF (npts.LE.3) THEN
         CALL getClosest(npts, px1, pt1, pt2, min_dist)
         RETURN
      END IF
   
      ! if available_ranks > 1: divide again
      IF (available_ranks.GT.1) THEN
         ldivs = 2**INT(LOG(REAL(available_ranks))/LOG(REAL(2)))
         rdivs = available_ranks - ldivs
         IF (rdivs.GT.0) THEN
            imid = INT(REAL(ldivs) / REAL(available_ranks) * npts)  ! index for proportional division
         ELSE
            ldivs = ldivs/2
            rdivs = ldivs ! next division
         END IF
      END IF
   
      ! middle point (x1 from imid, or between imid and imid+1?)
      L = (px1(imid, 1) + px1(imid+1, 1)) / 2.
      ! L = points(imid, 1)

      ! NOTE: subroutine?
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

      ! if rdivs > 1, we need to send to another rank
      IF (rdivs.GT.0) THEN
         tgt_rank = rank + ldivs
         commtag = tgt_rank
         CALL MPI_Send (rank, 1, MPI_INTEGER, tgt_rank, commtag, &
            MPI_COMM_WORLD, ierr)
         ! Send rdivs (as available_ranks), npts, px1, px2
         commtag = tgt_rank + 100
         CALL MPI_Send (rdivs, 1, MPI_INTEGER, tgt_rank, commtag, &
            MPI_COMM_WORLD, ierr)
         commtag = tgt_rank + 200
         CALL MPI_Send (ir, 1, MPI_INTEGER, tgt_rank, commtag, &
            MPI_COMM_WORLD, ierr)
         msglen = ir * ndim
         commtag = tgt_rank + 300
         CALL MPI_Send (px1(imid+1:, :), msglen, MPI_DOUBLE_PRECISION, tgt_rank,&
            commtag, MPI_COMM_WORLD, ierr)
         commtag = tgt_rank + 400
         CALL MPI_Send (px2_r, msglen, MPI_DOUBLE_PRECISION, tgt_rank,&
            commtag, MPI_COMM_WORLD, ierr)
         commtag = tgt_rank + 500
         CALL MPI_Send (level+1, 1, MPI_INTEGER, tgt_rank, commtag, &
            MPI_COMM_WORLD, ierr)
         commtag = tgt_rank + 600
         CALL MPI_Send (sublevel*2, 1, MPI_INTEGER, tgt_rank, commtag, &
            MPI_COMM_WORLD, ierr)
      END IF

      ! left side
      CALL closestPairDC_MPI (imid, px1(1:imid, :), px2_l, pl, ql, &
          ldist, level+1, sublevel*2-1, ldivs)
   
      IF (rdivs.GT.0) THEN
         ! Receive pr, qr, rdist from tgt rank
         src = tgt_rank
         commtag = tgt_rank + 1000
         CALL MPI_Recv (rdist, 1, MPI_DOUBLE_PRECISION, src, commtag, &
            MPI_COMM_WORLD, mpi_stat, ierr)
         commtag = tgt_rank + 1100
         CALL MPI_Recv (pr, ndim, MPI_DOUBLE_PRECISION, src, commtag, &
            MPI_COMM_WORLD, mpi_stat, ierr)
         commtag = tgt_rank + 1200
         CALL MPI_Recv (qr, ndim, MPI_DOUBLE_PRECISION, src, commtag, &
            MPI_COMM_WORLD, mpi_stat, ierr)
      ELSE
         ! else execute here
         CALL closestPairDC_MPI (ir, px1(imid+1:, :), px2_r, pr, qr, &
             rdist, level+1, sublevel*2, 1)
      END IF
   
      ! re-obtain L, to avoid modifications from lower levels
      L = (px1(imid, 1) + px1(imid+1, 1)) / 2.
      ! L = points(imid, 1)

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
   
      ! if host != myrank: send results to host (needed?)
   
   END SUBROUTINE closestPairDC_MPI
        
END MODULE mpi_mod

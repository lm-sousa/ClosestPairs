!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO OBTAIN THE TWO CLOSEST POINTS FROM A LIST OF POINTS (2D)
!!
!!  version: 1.0
!!  date: 2022-05
!!  encoding: utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!
!!! IMPROVEMENTS !!!
!! - improve consistency of tags
!! - quickSort for large array sizes
!! - non-blocking communication operations? 
!! - do not store distances on getClosest
!!!!!!!!!!!!!!!!!!!!

PROGRAM closest_mpi

USE sort_mod  ! sorting subroutines
USE mpi_mod  ! mpi vars and subroutines (imports closest_mod, mpi.h)

IMPLICIT NONE

! variable declaration: general in closest_mod

INTEGER :: rank_np  
   ! no of points in current rank
INTEGER :: closest_rank, closest_rank_strip  
   ! rank with minimum distance, strip dist
INTEGER :: extra_pts  
   ! modulo of np/nranks
INTEGER :: i1
   ! extra index
INTEGER, ALLOCATABLE, DIMENSION(:) :: x1_isort, x2_isort
   ! x1, x2-ordered index lists
INTEGER, ALLOCATABLE, DIMENSION(:) :: np_ranks, np_strips
   ! no points for each rank
INTEGER, ALLOCATABLE, DIMENSION(:) :: scatter_counts, scatter_displace
   ! additional arrays for scatter operation
REAL(KIND=REAL64) :: delta  
   ! extent for strip search
REAL(KIND=REAL64) :: strip_mindist
   ! minimum distance within strip
REAL(KIND=REAL64) :: strip_location
   ! location (x1) of strip
REAL(KIND=REAL64), DIMENSION(2) :: reduce_mindist, reduce_strip_mindist
REAL(KIND=REAL64), DIMENSION(2) :: strip_mindist_out, mindist_out
   ! outputs for MPI_MINLOC from MPI_Reduce
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: strips
   ! location of vertical strips
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: p1_strip, p2_strip
   ! results from  strip getClosest
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: dummy
   ! dummy array (for sorted list output)
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: x1_sort, x2_sort
   ! sorted arrays
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: all_strips
   ! points for all strips
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: rank_pts, strip_pts
   ! points for each rank, strip

verb = 0
root = 0  ! root rank is 0

closest_rank = 0
mindist_out = 0.d0
reduce_mindist = 0.d0 ! mindist for reduce operation with MPI_MINLOC (dimension 2,
                   ! value in dim 1, rank in dim 2 (DBLE)

strip_mindist = LARGE ! large value
strip_location = 0.d0

CALL MPI_Init (ierr)
CALL MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
CALL MPI_Comm_size (MPI_COMM_WORLD, nranks, ierr)

! rank 0 operations
IF (rank.EQ.root) CALL measure_cpu_time(cput)

! read command line arguments (all ranks, so that any rank can output results)
 CALL parseArguments (infile, outfile)

IF (rank.EQ.root) THEN
   
   ! read from file 
   np = 0
   ndim = 2
   CALL readPoints (infile, np, ndim, domain, all_pts)

END IF

! broadcast ndim
tag = root
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
CALL MPI_Bcast (ndim, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

IF (rank.EQ.root) THEN
   
   ! SEQUENTIAL SORT
   ! NOTE: Explore parallel sort (mergesort, quicksort)
   ! sort in x1
   ALLOCATE (dummy(np))
   ALLOCATE (x1_isort(np))
   WRITE(*, *) "Sorting points in x1..."
   CALL quickSortMod (np, all_pts(:, 1), dummy, x1_isort)
   ALLOCATE (x1_sort(np, ndim))
   x1_sort = all_pts(x1_isort, :)
   DEALLOCATE (x1_isort)
   DEALLOCATE (all_pts)
   
   ! sort in x2
   ALLOCATE (x2_isort(np))
   WRITE(*, *) "Sorting points in x2..."
   CALL quickSortMod (np, x1_sort(:, 2), dummy, x2_isort)
   ALLOCATE (x2_sort(np, ndim))
   x2_sort = x1_sort(x2_isort, :)
   DEALLOCATE (x2_isort)
   DEALLOCATE (dummy)

   ! calculate points to send to each process (balanced approach)
   WRITE(*, *) "Calculating number of points for each rank..."
   ! NOTE: at least 2 points should be in each rank
   ALLOCATE (np_ranks(nranks))
   np_ranks(:) = INT(np/nranks)
   extra_pts = MODULO(np, nranks)
   DO i = 2, extra_pts + 1
      np_ranks(i) = np_ranks(i) + 1       
   END DO

END IF

! Broadcast np (or scatter)
ALLOCATE (scatter_counts(nranks))
ALLOCATE (scatter_displace(nranks))
scatter_counts = 1
scatter_displace = [ (i, i=0, nranks-1) ]

CALL MPI_Scatterv (np_ranks, scatter_counts, scatter_displace, MPI_INTEGER, rank_np, 1, &
   MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

ALLOCATE (rank_pts(rank_np, ndim))
! ScatterV
IF (rank==root) THEN 
   scatter_counts = np_ranks
   scatter_displace = 0
   DO i = 1, nranks-1
      scatter_displace(i+1) = scatter_displace(i) + np_ranks(i)
   END DO
END IF
CALL MPI_Scatterv (x1_sort(:, 1), scatter_counts, scatter_displace, MPI_DOUBLE_PRECISION, &
   rank_pts(:, 1), rank_np, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
CALL MPI_Scatterv (x1_sort(:, 2), scatter_counts, scatter_displace, MPI_DOUBLE_PRECISION, &
   rank_pts(:, 2), rank_np, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
! IF METHOD == divide and conquer (NOTE: Need to cross x2 with x1)
! CALL MPI_Scatterv (x2_sort(:, 1), scatter_counts, scatter_displace, MPI_DOUBLE_PRECISION, &
!    rank_pts(:, 1), rank_np, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
! CALL MPI_Scatterv (x2_sort(:, 2), scatter_counts, scatter_displace, MPI_DOUBLE_PRECISION, &
!    rank_pts(:, 2), rank_np, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
! DO i = 1, rank_np
!    WRITE(*,'("R", I2, " npr=", 2(1X, F9.5))') rank, rank_pts(i, :)
! END DO

! calculate closest pair for each rank
ALLOCATE (p1(ndim))
ALLOCATE (p2(ndim))
p1 = 0.d0
p2 = 0.d0
mindist = 0.d0
! CALL measure_cpu_time(cput1, msg=cputlabel)
CALL getClosest(rank_np, rank_pts, p1, p2, mindist)
! CALL measure_cpu_time(cput1, msg=cputlabel)
! 
! Reduce and obtain delta and rank with closest pair (with MIN_loc)
! delta is closest distance
reduce_mindist(1) = mindist
reduce_mindist(2) = DBLE(rank)
CALL MPI_Allreduce (reduce_mindist, mindist_out, 1, MPI_2DOUBLE_PRECISION, &
                    MPI_MINLOC, MPI_COMM_WORLD, ierr)

delta = mindist_out(1)
! round to nearest integer
closest_rank = NINT(mindist_out(2))

IF (rank.EQ.root) THEN
   WRITE(*, *) "Location of delta: RANK ", closest_rank, " , dist = ", delta
END IF

DEALLOCATE (rank_pts)

! root checks for points within +-delta of vertical strips
i1 = 0  ! auxiliary index
IF (rank.EQ.root) THEN
   ALLOCATE (strips(nranks))
   strips = LARGE
   ! Get strip locations
   DO i = 2, nranks
      i1 = i1 + np_ranks(i-1)
      ! NOTE: mid-way between last point of rank i-1, and first of i
      strips(i) = (x1_sort(i1, 1) + x1_sort(i1+1, 1)) / 2.d0
      ! strips(i-1) = all_pts(i1, 1)
   END DO

   ! get points within +-delta
   ALLOCATE (np_strips(nranks))
   ALLOCATE (all_strips(np, ndim))
   all_strips = LARGE
   WRITE(*, *) "Calculating np for each strip..."
   np_strips = 0
   i1 = 0
   DO j = 2, nranks
      DO i = 1, np
         IF (ABS(x2_sort(i, 1) - strips(j)).LT.delta) THEN
            i1 = i1 + 1
            np_strips(j) = np_strips(j) + 1
            all_strips(i1, :) = x2_sort(i, :)
         END IF
      END DO
      ! IF np_strips(j) < 2 eliminate points (no need for execution)
      IF (np_strips(j).LT.2) THEN
         i1 = i1 - np_strips(j)
         np_strips(j) = 0
         all_strips(i1:, :) = LARGE ! remove strip points
      END IF
   END DO
END IF

! Allocate strip points (re-use rank_pts, rank_np)
! ScatterV of np_strips
scatter_counts = 0  ! do not operate on
scatter_counts(2:) = 1
scatter_displace = [ (i, i=0, nranks-1) ]

rank_np = 0
CALL MPI_Scatterv (np_strips, scatter_counts, scatter_displace, MPI_INTEGER, rank_np, 1, &
   MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

! Scatter strip locations
strip_location = LARGE
CALL MPI_Scatterv (strips, scatter_counts, scatter_displace, MPI_DOUBLE_PRECISION, &
   strip_location, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)

ALLOCATE (rank_pts(rank_np, ndim))
! Scatter of strip points
IF (rank==root) THEN 
   scatter_counts = np_strips
   scatter_displace = 0
   DO i = 1, nranks-1
      scatter_displace(i+1) = scatter_displace(i) + np_strips(i)
   END DO
END IF

CALL MPI_Scatterv (all_strips(:i1, 1), scatter_counts, scatter_displace, &
   MPI_DOUBLE_PRECISION, rank_pts(:, 1), rank_np, MPI_DOUBLE_PRECISION, &
   root, MPI_COMM_WORLD, ierr)
CALL MPI_Scatterv (all_strips(:i1, 2), scatter_counts, scatter_displace, &
   MPI_DOUBLE_PRECISION, rank_pts(:, 2), rank_np, MPI_DOUBLE_PRECISION, &
   root, MPI_COMM_WORLD, ierr)

! calculate closest pair in strip for each rank
ALLOCATE (p1_strip(ndim))
ALLOCATE (p2_strip(ndim))
p1_strip = 0.d0
p2_strip = 0.d0
strip_mindist = LARGE

! 2D strip calculation NOTE: Replace for subroutine
IF (rank_np.GE.2) THEN
   CALL getClosestStrip(rank_np, rank_pts, strip_location, delta, &
      p1_strip, p2_strip, strip_mindist)
END IF

! reduce location and value
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
reduce_strip_mindist(1) = strip_mindist
reduce_strip_mindist(2) = DBLE(rank)
CALL MPI_Allreduce (reduce_strip_mindist, strip_mindist_out, 1, &
                 MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD, ierr)

! separate mindist for strip
mindist = strip_mindist_out(1)
closest_rank_strip = INT(strip_mindist_out(2))

IF (mindist.LT.delta) THEN
   closest_rank = closest_rank_strip
   p1 = p1_strip
   p2 = p2_strip
ELSE
   mindist = delta
END IF

DEALLOCATE (p1_strip)
DEALLOCATE (p2_strip)

! output from rank with closest pair
IF (rank.EQ.closest_rank) THEN
   ! write results (stdout - optional)
   CALL writeResultsToStd(mindist, p1, p2)
   ! write results (file)
   CALL writeResultsToFile(outfile, p1, p2)
END IF

! finish measuring cputime
IF (rank.EQ.root) CALL measure_cpu_time(cput)
 
! deallocate
IF (rank_np.GE.2) DEALLOCATE (rank_pts)
IF (rank.EQ.root) THEN
   DEALLOCATE (x1_sort)
   DEALLOCATE (x2_sort)
END IF
DEALLOCATE (p1)
DEALLOCATE (p2)

CALL MPI_Finalize(ierr)

END PROGRAM closest_mpi

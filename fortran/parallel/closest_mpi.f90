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

! USE MPI
USE closest_mod  ! contains variable declarations, subroutines
USE sort_mod  ! sorting subroutines

IMPLICIT NONE

INCLUDE 'mpif.h'  ! subroutines for MPI

! variable declaration: general in closest_mod

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
INTEGER :: rank_np  
   ! no of points in current rank
INTEGER :: rank_strip_np  
   ! no of points within strip in current rank
INTEGER :: closest_rank, strip_closest_rank  
   ! rank with minimum distance, strip dist
INTEGER :: extra_pts  
   ! modulo of np/nranks
INTEGER :: ipoint  
   ! index of some point
INTEGER :: i1, i2  
   ! start and end index for partitioning
INTEGER :: np_opp, np_same 
   ! np to opposite and same side of strip
INTEGER, ALLOCATABLE, DIMENSION(:) :: x1_isort, x2_isort  
   ! x1, x2-ordered index lists
INTEGER, ALLOCATABLE, DIMENSION(:) :: np_ranks, np_strips 
   ! no points for each rank
REAL(KIND=REAL64) :: delta  
   ! extent for strip search
REAL(KIND=REAL64) :: strip_mdist
   ! minimum distance within strip
REAL(KIND=REAL64) :: dx1, dx2
   ! distances in x and y
REAL(KIND=REAL64) :: dist
   ! distance between two points 
REAL(KIND=REAL64) :: strip_x1
   ! location (x1) of strip
REAL(KIND=REAL64), DIMENSION(2) :: mdist_reduce, smdist_reduce, mdist_out
REAL(KIND=REAL64), DIMENSION(2) :: smdist_out
   ! outputs for MPI_MINLOC from MPI_Reduce
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: strips
   ! location of vertical strips
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: dummy
   ! dummy array (for sorted list output)
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: rank_pts, strip_pts  
   ! points for each rank, strip

root = 0  ! root rank is 0

closest_rank = 0
mdist_out = 0.d0
mdist_reduce = 0.d0 ! mindist for reduce operation with MPI_MINLOC (dimension 2,
                   ! value in dim 1, rank in dim 2 (DBLE)

strip_mdist = LARGE ! large value
strip_x1 = 0.d0

CALL MPI_Init (ierr)
CALL MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
CALL MPI_Comm_size (MPI_COMM_WORLD, nranks, ierr)

! WRITE(*,*) "Current rank: ", rank
! WRITE(*,*) "Number of ranks: ", nranks

! rank 0 operations
IF (rank.EQ.root) THEN
   
   ! start measuring cputime
   CALL measure_cpu_time(cput)

   ! read command line arguments
   CALL parseArguments (infile, outfile)
   
   ! read from file 
   np = 0
   ndim = 2
   CALL readPoints (infile, np, ndim, domain, all_pts, verb=0)
   
   !! Divide in vertical strips after sorting, for a balanced distribution
   ! sort in x1 (get sorted indices), serial
   ! NOTE: Explore parallel sort (heapsort, quicksort)
   ALLOCATE (x1_isort(np))
   ALLOCATE (dummy(np))
   WRITE(*, *) "Sorting points in x1..."
   ! NOTE: high-order method for large arrays, prefer quickSort
   ! CALL selectionSort (all_pts(:, 1), np, dummy, x1_isort)
   CALL quickSortMod(np, all_pts(:, 1), dummy, x1_isort)
   ! WRITE(*, *) "Pre-sort: ", all_pts(:, 1)
   ! WRITE(*, *) "indices: ", x1_isort
   ! order points
   all_pts = all_pts(x1_isort, :)
   ! WRITE(*, *) "Post-sort: ", all_pts(:, 1)
   DEALLOCATE(dummy)
   DEALLOCATE(x1_isort)

   ! calculate points to send to each process (balanced approach)
   WRITE(*, *) "Calculating number of points for each rank..."
   ALLOCATE (np_ranks(nranks))
   np_ranks(:) = INT(np/nranks)
   extra_pts = MODULO(np, nranks)
   DO i = 2, extra_pts + 1
      np_ranks(i) = np_ranks(i) + 1       
   END DO
   ! WRITE(*, *) "Rank points: ", np_ranks(:)

   rank_np = np_ranks(1)
   ! send number of points to each rank
   DO i = 2, nranks
      tag = i - 1
      CALL MPI_Send (np_ranks(i), 1, MPI_INTEGER, i-1, tag,&
                     MPI_COMM_WORLD, ierr)
   END DO

ELSE

   ! receive number of points
   source = root
   tag = rank
   CALL MPI_Recv (rank_np, 1, MPI_INTEGER, source, tag, &
                  MPI_COMM_WORLD, mpi_stat, ierr)

END IF

WRITE(*, *) "RANK: ", rank, ". np: ", rank_np

! broadcast ndim
! WRITE(*, *) "Broadcasting ndim..."
CALL MPI_Bcast (ndim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
! NOTE: Barrier needed?
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
! WRITE(*, *) "RANK = ", rank, "ndim = ", ndim

msgsize = rank_np * ndim  ! msgsize for communication of points
ALLOCATE (rank_pts(rank_np, ndim))

IF (rank.EQ.root) THEN

   i1 = 1
   i2 = np_ranks(1)
   rank_pts = all_pts(i1:i2, :)
   ! WRITE(*, *) rank_pts
   ! get points and send to other processes, calculate strip loc (for root)
   ! number of vertical strips = nranks - 1
   ALLOCATE (strips(nranks - 1))
   ! WRITE(*, *) "Calculating strip locations..."
   DO i = 2, nranks
      ! NOTE: strip location is mid-way between last point of rank i-1, and 
      !       first of i
      i1 = i2 + 1
      ! WRITE(*, *) "Index1: ", i1
      ! WRITE(*, *) "Index2: ", i2
      strips(i-1) = (all_pts(i2, 1) + all_pts(i1, 1)) / 2.d0
      i2 = i2 + np_ranks(i)
      tag = 100 + i - 1
      msgsize = np_ranks(i) * ndim
      ! WRITE(*, *) "Sending points to rank ", i - 1, ". Tag = ", tag
      CALL MPI_Send (all_pts(i1:i2, :), msgsize, &
                     MPI_DOUBLE_PRECISION, i-1, tag,&
                     MPI_COMM_WORLD, ierr)
   END DO

   ! WRITE(*, *) "Strip locations: ", strips
   DEALLOCATE (np_ranks)

ELSE

   source = root
   tag = rank + 100
   msgsize = rank_np * ndim
   ! WRITE(*, *) "RANK: ", rank, ". Receiving subd points from source. Tag = ", &
   !             tag
   CALL MPI_Recv (rank_pts, msgsize, MPI_DOUBLE_PRECISION, source, tag, &
                  MPI_COMM_WORLD, mpi_stat, ierr)
END IF

! WRITE(*, *) "RANK: ", rank, ". Points (x1): ", rank_pts(:, 1)
! WRITE(*, *) "RANK: ", rank, ". Points (x2): ", rank_pts(:, 2)

! NOTE: Barrier needed here?
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
! CALL MPI_Bcast (strips, nranks-1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)

! calculate closest pair for each rank
ALLOCATE (p1(ndim))
ALLOCATE (p2(ndim))
p1 = 0.d0
p2 = 0.d0
mindist = 0.d0
! WRITE(*, *) "RANK ", rank, ": Calculating closest point (of ", rank_np, ")..."
WRITE(cputlabel, "('RANK ', I3, ': getClosest')") rank ! pass rank to cputlabel
CALL measure_cpu_time(cput1, msg=cputlabel)
CALL getClosest(rank_np, ndim, rank_pts, p1, p2, mindist, verb=0)

WRITE(*, '("    RANK ", I3, ": mindist = ", ES11.4)') rank, mindist
WRITE(*, '("    RANK ", I3, ": p1 = ", 2(ES23.15E3, 1X))') rank, p1
WRITE(*, '("    RANK ", I3, ": p2 = ", 2(ES23.15E3, 1X))') rank, p2

CALL measure_cpu_time(cput1, msg=cputlabel)

! wait for all to finish
! CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
! Reduce and obtain delta and rank with closest pair (with MIN_loc)
! delta is closest distance
! WRITE(*, *) "RANK ", rank, ": reducing..."
mdist_reduce(1) = mindist
! NOTE: Avoid DBLE, define double type
mdist_reduce(2) = DBLE(rank)
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
! CALL MPI_Reduce (mdist_reduce, mdist_out, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC,&
!                  root, MPI_COMM_WORLD, ierr)
CALL MPI_Allreduce (mdist_reduce, mdist_out, 1, MPI_2DOUBLE_PRECISION, &
                    MPI_MINLOC, MPI_COMM_WORLD, ierr)

delta = mdist_out(1)
! round to nearest integer
closest_rank = NINT(mdist_out(2))

! IF (rank.EQ.root) THEN
!    WRITE(*, *) "Location of mindist: RANK ", closest_rank, " , dist = ", delta
! END IF

msgsize = ndim
tgt = root
IF (closest_rank.NE.root) THEN
   IF (rank.EQ.closest_rank) THEN
      ! WRITE(*, *) "Sending p1, ", p1
      ! WRITE(*, *) "from RANK ", rank, " , to ", root
      CALL MPI_Send (p1, msgsize, MPI_DOUBLE_PRECISION, tgt, 10000,&
                     MPI_COMM_WORLD, ierr)
      ! WRITE(*, *) "Sending p2, from RANK ", rank, " , to ", root
      CALL MPI_Send (p2, msgsize, MPI_DOUBLE_PRECISION, tgt, 10001,&
                     MPI_COMM_WORLD, ierr)
   ELSE IF ((rank.EQ.root).AND.(closest_rank.NE.root)) THEN
      msgsize = ndim
      source = closest_rank
      ! get p1 and p2
      ! WRITE(*, *) "Receiving p1, from RANK ", closest_rank
      CALL MPI_Recv (p1, msgsize, MPI_DOUBLE_PRECISION, source, 10000, &
                     MPI_COMM_WORLD, mpi_stat, ierr)
      ! WRITE(*, *) "P1: ", p1
      ! WRITE(*, *) "Receiving p2, from RANK ", closest_rank
      CALL MPI_Recv (p2, msgsize, MPI_DOUBLE_PRECISION, source, 10001, &
                     MPI_COMM_WORLD, mpi_stat, ierr)
      ! WRITE(*, *) "P2: ", p2
   END IF
END IF


! update tag and target
tag = 200
tgt = 0

! root checks for points within +-delta of vertical strips
! NOTE: improve approach, this is practically sequential...
IF (rank.EQ.root) THEN
   ! get number of points within +-delta
   ALLOCATE (np_strips(nranks-1))
   WRITE(*, *) "Calculating np for each strip..."
   np_strips = 0
   DO j = 1, nranks-1
      DO i = 1, np
         IF (all_pts(i, 1).GT.(strips(j) + delta)) THEN
            EXIT
         ELSE IF (all_pts(i, 1).GE.(strips(j) - delta)) THEN
            np_strips(j) = np_strips(j) + 1
         END IF
      END DO

      ! WRITE(*, *) "  NP Strips = ", np_strips

      ! send number of points
      tgt = j
      ! tgt = tgt + 1
      tag = tag + 1
      ! WRITE(*, *) "Sending np to RANK ", tgt, ", TAG ", tag
      CALL MPI_Send (np_strips(j), 1, MPI_INTEGER, tgt, tag,&
                     MPI_COMM_WORLD, ierr)
      IF (np_strips(j).GT.1) THEN
         tag = tag + 100
         ! WRITE(*, *) "Sending strips to RANK ", tgt, ", TAG ", tag
         CALL MPI_Send (strips(j), 1, MPI_DOUBLE_PRECISION, tgt, tag,&
                        MPI_COMM_WORLD, ierr)
         ALLOCATE (strip_pts(np_strips(j), ndim))
         strip_pts = 0.d0
         ipoint = 0
         ! get points within +-delta of each strip
         WRITE(*, *) "Gathering points for each strip..."
         DO i = 1, np
            IF (all_pts(i, 1).GT.(strips(j) + delta)) THEN
               EXIT
            ELSE IF (all_pts(i, 1).GE.(strips(j) - delta)) THEN
               ipoint = ipoint + 1
               strip_pts(ipoint, :) = all_pts(i, :)
            END IF
         END DO

         ! 4. Send points to processes 
         tag = tag + 100
         msgsize = np_strips(j) * ndim
         ! WRITE(*, *) "Send points to RANK ", tgt, " TAG ", tag
         CALL MPI_Send (strip_pts, msgsize, MPI_DOUBLE_PRECISION, tgt, tag,&
                        MPI_COMM_WORLD, ierr)
         DEALLOCATE (strip_pts)
         tag = tag - 200

      END IF
   END DO
   DEALLOCATE (np_strips)
ELSE

   tag = 200 + rank
   source = root
   rank_strip_np = 0
   ! receive np
   ! WRITE(*, *) "RANK ", rank, " Recv np strip . TAG ", tag
   CALL MPI_Recv (rank_strip_np, 1, MPI_INTEGER, source, tag, &
                  MPI_COMM_WORLD, mpi_stat, ierr)
   ! WRITE(*, *) "    RANK ", rank, ": Strip NP ", rank_strip_np

   ! process only if rank_strip_np > 1
   IF (rank_strip_np.GT.1) THEN
      tag = tag + 100
      ! receive strip location
      ! WRITE(*, *) "RANK ", rank, " Recv strip location. TAG ", tag
      CALL MPI_Recv (strip_x1, 1, MPI_DOUBLE_PRECISION, source, tag, &
                     MPI_COMM_WORLD, mpi_stat, ierr)

      ! receive strip_pts
      ALLOCATE (strip_pts(rank_strip_np, ndim))
      msgsize = rank_np * ndim
      tag = tag + 100
      ! WRITE(*, *) "RANK ", rank, " Recv strip_pts . TAG ", tag
      CALL MPI_Recv (strip_pts, msgsize, MPI_DOUBLE_PRECISION, source, tag, &
                     MPI_COMM_WORLD, mpi_stat, ierr)

      ! sort with x2 
      ALLOCATE (x2_isort(rank_strip_np))
      ALLOCATE (dummy(rank_strip_np))
      ! WRITE(*, *) "Sorting according to x2..."
      ! WRITE(*, *) "OG strip_pts", strip_pts(:, 2)
      ! CALL selectionSort (strip_pts(:, 2), rank_strip_np, dummy, x2_isort)
      CALL quickSortMod(rank_strip_np, strip_pts(:, 2), dummy, x2_isort)
      ! apply sort
      strip_pts = strip_pts(x2_isort, :)
      ! WRITE(*, *) "Sorted strip_pts", strip_pts(:, 2)
      DEALLOCATE (dummy) 
      DEALLOCATE (x2_isort) 

      ! Calculate distances (max. 4 in opposite side, 3 in same side)
      WRITE(*, *) "Calculating distances within strip..."
      cput1 = -1.
      WRITE(cputlabel, "('RANK ', I3, ': getClosestStrip')") rank
      CALL measure_cpu_time(cput1, msg=cputlabel)
      ! strip_mdist = delta
      ! strip_mdist = LARGE
      DO i=1, rank_strip_np - 1
         np_opp = 0
         np_same = 0
         dx1 = 0.d0
         dx2 = 0.d0
         DO j=i+1, rank_strip_np
            ! WRITE(*, *) "np_opp = ", np_opp
            ! WRITE(*, *) "np_same = ", np_same

            ! exit, skip conditions (max. comparisons)
            IF ((np_opp.EQ.4).AND.(np_same.LT.3)) CONTINUE
            IF ((np_opp.LT.4).AND.(np_same.EQ.3)) CONTINUE
            IF ((np_opp.EQ.4).AND.(np_same.EQ.3)) EXIT

            dx2 = ABS(strip_pts(j, 2) - strip_pts(i, 2)) 
            IF (dx2.GT.delta) THEN
               ! WRITE(*, *) "Eliminated due to dx2: ", j
               EXIT
            ELSE 
               dx1 = ABS(strip_pts(i, 2) - strip_pts(j, 2))
               ! if larger than delta, skip
               IF (dx1.LE.delta) THEN
                   ! WRITE(*, *) "Point with lower distance..."
                   ! calculate distances
                   dist = SQRT(dx1**2 + dx2**2)
                   ! if on opposite side
                   IF (((strip_pts(i, 1) - strip_x1)*&
                      (strip_x1 - strip_pts(j, 1))).GT.0) THEN
                      np_opp = np_opp + 1 ! increase counter of opposite points
                   ! if on same side
                   ELSE
                       np_same = np_same + 1
                   END IF
               END IF
               ! replace ClosestPair (for strip)
               IF (dist.LT.strip_mdist) THEN
                  strip_mdist = dist
                  p1 = strip_pts(i, :)
                  p2 = strip_pts(j, :)
               END IF
            END IF
         END DO
      END DO
      CALL measure_cpu_time(cput1, msg=cputlabel)
   END IF
END IF

! 7. reduce location and value
! WRITE(*, *) "Barrier..."
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
smdist_reduce(1) = strip_mdist
smdist_reduce(2) = DBLE(rank)
! WRITE(*, *) "Reducing..."
! CALL MPI_Reduce (smdist_reduce, smdist_out, 1, &
!                  MPI_2DOUBLE_PRECISION, MPI_MINLOC, root, MPI_COMM_WORLD, 
!                  ierr)
CALL MPI_Allreduce (smdist_reduce, smdist_out, 1, &
                 MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD, ierr)

mindist = smdist_out(1)
strip_closest_rank = INT(smdist_out(2))

! NOTE: condition mindist.LT.delta works as all but the rank with the actual 
! minimum still have mindist = delta. If Allreduce had been used, this would 
! not be necessary
! IF (mindist.LT.delta) THEN

IF (rank.EQ.strip_closest_rank) THEN
   tag = 500
   msgsize = ndim
   CALL MPI_Send (p1, msgsize, MPI_DOUBLE_PRECISION, tgt, tag,&
                  MPI_COMM_WORLD, ierr)
   tag = 600
   CALL MPI_Send (p2, msgsize, MPI_DOUBLE_PRECISION, tgt, tag,&
                  MPI_COMM_WORLD, ierr)
ELSE IF (rank.EQ.root) THEN
   tag = 500
   msgsize = ndim
   source = strip_closest_rank
   CALL MPI_Recv (p1, msgsize, MPI_DOUBLE_PRECISION, source, tag, &
                  MPI_COMM_WORLD, mpi_stat, ierr)
   tag = 600
   CALL MPI_Recv (p2, msgsize, MPI_DOUBLE_PRECISION, source, tag, &
                  MPI_COMM_WORLD, mpi_stat, ierr)
END IF
! END IF

IF (rank.EQ.root) THEN
   ! write results (stdout - optional)
   CALL writeResultsToStd(delta, ndim, p1, p2)
   ! write results (file)
   CALL writeResultsToFile(outfile, ndim, p1, p2)
   ! finish measuring cputime
   CALL measure_cpu_time(cput)
END IF

! deallocate
IF (rank.EQ.root) THEN
   DEALLOCATE (all_pts)
END IF
DEALLOCATE (p1)
DEALLOCATE (p2)

CALL MPI_Finalize(ierr)

END PROGRAM closest_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO OBTAIN THE TWO CLOSEST POINTS FROM A LIST OF POINTS (2D)
!!
!!  version: 1.0
!!  date: 2022-05
!!  encoding: utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM closest_dc_mpi

! USE closest_mod  ! contains variable declarations, subroutines
USE mpi_mod  ! mpi variables and routines
USE sort_mod  ! sorting subroutines (sequential)

IMPLICIT NONE


! variable declaration: general in closest_mod

INTEGER :: lvl, sublvl, free_ranks
   ! level, sublevel of DaC, available ranks
INTEGER, ALLOCATABLE, DIMENSION(:) :: x1_isort, x2_isort  
   ! x1, x2-ordered index lists
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:) :: dummy
   ! dummy array (for sorted list output)
REAL(KIND=REAL64), ALLOCATABLE, DIMENSION(:,:) :: x1_sort, x2_sort
   ! points for each rank, strip

verb = 0  ! verbosity level
root = 0  ! root rank is 0

! init MPI, get ranks
CALL MPI_Init (ierr)
CALL MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
CALL MPI_Comm_size (MPI_COMM_WORLD, nranks, ierr)

! rank 0 operations
IF (rank.EQ.root) THEN
   
   ! start measuring cputime
   CALL measure_cpu_time(cput)

   ! read command line arguments
   CALL parseArguments (infile, outfile)
   
   ! read from file 
   np = 0
   ndim = 2
   CALL readPoints (infile, np, ndim, domain, all_pts)

END IF

tag = root
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
CALL MPI_Bcast (ndim, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

ALLOCATE (p1(ndim))
ALLOCATE (p2(ndim))

IF (rank.EQ.root) THEN
   ! SEQUENTIAL SORT
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
   CALL quickSortMod (np, x1_sort(:, 2), dummy, x2_isort)
   ALLOCATE (x2_sort(np, ndim))
   x2_sort = x1_sort(x2_isort, :)
   DEALLOCATE (x2_isort)
   DEALLOCATE (dummy)

   lvl = 1
   sublvl = 1
   free_ranks = nranks

   ! Execute closest
   cput1 = -1.
   cputlabel = 'closestPairDC_MPI'
   CALL measure_cpu_time(cput1, cputlabel)
   CALL closestPairDC_MPI (np, x1_sort, x2_sort, p1, p2, &
       mindist, lvl, sublvl, free_ranks)
   CALL measure_cpu_time(cput1, cputlabel)

   DEALLOCATE (x1_sort)
   DEALLOCATE (x2_sort)

ELSE
   tag = rank
   CALL MPI_Recv (source, 1, MPI_INTEGER, MPI_ANY_SOURCE, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)
   tag = rank+100
   CALL MPI_Recv (free_ranks, 1, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)
   tag = rank+200
   CALL MPI_Recv (np, 1, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)
   ALLOCATE (x1_sort(np, ndim))
   tag = rank+300
   msgsize = np*ndim
   CALL MPI_Recv (x1_sort, msgsize, MPI_DOUBLE_PRECISION, source, &
      tag, MPI_COMM_WORLD, mpi_stat, ierr)
   ALLOCATE (x2_sort(np, ndim))
   tag = rank+400
   CALL MPI_Recv (x2_sort, msgsize, MPI_DOUBLE_PRECISION, source, &
      tag, MPI_COMM_WORLD, mpi_stat, ierr)
   ! receive lvl, sublvl
   tag = rank+500
   CALL MPI_Recv (lvl, 1, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)
   tag = rank+600
   CALL MPI_Recv (sublvl, 1, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)

   ! Execute closest
   CALL closestPairDC_MPI (np, x1_sort, x2_sort, p1, p2, &
       mindist, lvl, sublvl, free_ranks)

   DEALLOCATE (x1_sort)
   DEALLOCATE (x2_sort)

   ! Send results to host rank
   tgt = source
   tag = rank + 1000
   WRITE(*, *) "Sending results to...", tgt, tag
   CALL MPI_Send (mindist, 1, MPI_DOUBLE_PRECISION, tgt, tag, &
      MPI_COMM_WORLD, ierr)
   tag = rank + 1100
   CALL MPI_Send (p1, ndim, MPI_DOUBLE_PRECISION, tgt, tag, &
      MPI_COMM_WORLD, ierr)
   tag = rank + 1200
   CALL MPI_Send (p2, ndim, MPI_DOUBLE_PRECISION, tgt, tag, &
      MPI_COMM_WORLD, ierr)
END IF

IF (rank.EQ.root) THEN
   ! write results (stdout - optional)
   CALL writeResultsToStd(mindist, p1, p2)
   ! write results (file)
   CALL writeResultsToFile(outfile, p1, p2)
   ! finish measuring cputime
   CALL measure_cpu_time(cput)
END IF

DEALLOCATE (p1)
DEALLOCATE (p2)

CALL MPI_Finalize(ierr)

END PROGRAM closest_dc_mpi

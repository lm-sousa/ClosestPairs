!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  PROGRAM TO OBTAIN THE TWO CLOSEST POINTS FROM A LIST OF POINTS (2D)
!!
!!  version: 1.0
!!  date: 2022-05
!!  encoding: utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM closest_dc_mpi

USE mpi_mod  ! mpi variables and routines
USE sort_mod  ! sorting subroutines (sequential)

IMPLICIT NONE


! variable declaration: general in closest_mod (loaded by mpi_mod)

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
! instrumentation
cput = 0
cpufile = 'cputime_dc_mpi.csv'
commfile = 'commtime_dc_mpi.csv'

CALL CPU_TIME(cput(1))
CALL MPI_Init (ierr)
CALL MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
CALL MPI_Comm_size (MPI_COMM_WORLD, nranks, ierr)
CALL CPU_TIME(cput(2))

! rank 0 operations
IF (rank.EQ.root) THEN
   ! read command line arguments
   CALL parseArguments (infile, outfile)
   ! read from file 
   np = 0
   ndim = 2
   CALL readPoints (infile, np, ndim, domain, all_pts)
END IF
   
tag = root
CALL MPI_Barrier (MPI_COMM_WORLD, ierr)
CALL CPU_TIME(cput(3))
CALL MPI_Bcast (ndim, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
CALL CPU_TIME(cput(4))

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

   CALL CPU_TIME(cput(5))

   lvl = 1
   sublvl = 1
   free_ranks = nranks

   ! Execute closest
   CALL closestPairDC_MPI (np, x1_sort, x2_sort, p1, p2, &
       mindist, lvl, sublvl, free_ranks)
   CALL CPU_TIME(cput(6))

   DEALLOCATE (x1_sort)
   DEALLOCATE (x2_sort)

ELSE
   ! SOURCE: MPI_ANY_SOURCE
   CALL CPU_TIME(cput(11))
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
   ! receive lvl, sublvl (optional)
   tag = rank+500
   CALL MPI_Recv (lvl, 1, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)
   tag = rank+600
   CALL MPI_Recv (sublvl, 1, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, mpi_stat, ierr)

   CALL CPU_TIME(cput(12))
   cput(13) = cput(13) + cput(12) - cput(11)

   ! Execute closest
   CALL closestPairDC_MPI (np, x1_sort, x2_sort, p1, p2, &
       mindist, lvl, sublvl, free_ranks)


   DEALLOCATE (x1_sort)
   DEALLOCATE (x2_sort)

   ! Send results to host rank
   CALL CPU_TIME(cput(11))
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
   CALL CPU_TIME(cput(12))
   cput(14) = cput(14) + cput(12) - cput(11)
END IF

CALL CPU_TIME(cput(7))
IF (rank.EQ.root) THEN
   ! write results (stdout - optional)
   CALL writeResultsToStd(mindist, p1, p2)
END IF
CALL CPU_TIME(cput(8))
IF (rank.EQ.root) THEN
   ! write results (file)
   CALL writeResultsToFile(outfile, p1, p2)
END IF
CALL CPU_TIME(cput(9))

DEALLOCATE (p1)
DEALLOCATE (p2)
CALL CPU_TIME(cput(10))

cput_red = cput
DO i = 2, 9
   cput(i) = cput_red(i) - cput_red(i-1)
END DO
CALL MPI_Reduce (cput, cput_red, MEASUREMENTS, MPI_REAL, &
   MPI_MAX, root, MPI_COMM_WORLD, ierr)

commtime = cput(13) + cput(14) + cput(2)
CALL MPI_Reduce (commtime, commtime_red, 1, MPI_REAL, &
   MPI_MAX, root, MPI_COMM_WORLD, ierr)

IF (rank.EQ.root) THEN
   ! write to cputime.txt
   cput = cput * 1000.  ! miliseconds
   cput_red = cput_red * 1000.  ! miliseconds
   commtime_red = commtime_red * 1000.
   
   ! stdout
   WRITE(*,*) 
   WRITE(*,*) "Time output (ms):"
   WRITE(*, '(" MPI Init:         ", F12.3)') cput_red(2)
   WRITE(*, '(" Arg + Read points:", F12.3)') cput_red(3)
   WRITE(*, '(" MPI Bcast:        ", F12.3)') cput_red(4)
   WRITE(*, '(" Sort:             ", F12.3)') cput(5)
   WRITE(*, '(" getClosest:       ", F12.3)') cput(6)
   WRITE(*, '(" MPI Send:         ", F12.3)') cput_red(13)
   WRITE(*, '(" MPI Recv:         ", F12.3)') cput_red(14)
   WRITE(*, '(" Strip calculation:", F12.3)') cput_red(15)
   WRITE(*, '(" Write to stdout:  ", F12.3)') cput_red(8)
   WRITE(*, '(" Write to file:    ", F12.3)') cput_red(9)
   WRITE(*, '(" ** Total MPI time:", F12.3, " **")') commtime_red
   WRITE(*, '(" ** Total time:    ", F12.3, " **")') cput_red(10)
   
   ! to file
   INQUIRE (FILE=cpufile, EXIST=fexist)
   IF (.NOT.fexist) THEN
      OPEN(8, FILE=cpufile, STATUS='NEW', ACTION='WRITE')
      WRITE(8, *) "np, nranks, args+readpts, sort, getclosest, strip, writeStd, writeFile, total"
      OPEN(9, FILE=commfile, STATUS='NEW', ACTION='WRITE')
      WRITE(9, *) "np, nranks, init, bcast, send, recv, total"
   ELSE
      OPEN(8, FILE=cpufile, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
      OPEN(9, FILE=commfile, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
   END IF
   
   ! WRITE MESSAGE
   WRITE(8, 1080) np, nranks, &
    cput_red(3), &
    cput_red(5), &
    cput_red(6), &
    cput_red(15), &
    cput_red(8), &
    cput_red(9), &
    cput_red(10)
   CLOSE(8)
   WRITE(9, 1090) np, nranks, &
      cput_red(2), &
      cput_red(4), &
      cput_red(13),&
      cput_red(14), &
      commtime_red
   CLOSE(9)
   1080 FORMAT(I7, ', ', I2, 7(', ', F11.3))
   1090 FORMAT(I7, ', ', I2, 5(', ', F11.3))
   
END IF

CALL MPI_Finalize(ierr)

END PROGRAM closest_dc_mpi

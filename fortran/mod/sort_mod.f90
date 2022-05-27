!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  MODULE WITH SORTING SUBROUTINES FOR ClosestPoints
!!
!!  date: 2022-05
!!  encoding utf-8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE sort_mod

USE iso_fortran_env  ! module for compiler-independent precision declaration
USE closest_mod  ! global variables

IMPLICIT NONE

CONTAINS

   SUBROUTINE selectionSort (array, asize, sorted, isorted)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: asize
      REAL(KIND=REAL64), DIMENSION(asize), INTENT(IN) :: array

      REAL(KIND=REAL64), DIMENSION(asize), INTENT(OUT) :: sorted
      INTEGER, DIMENSION(asize), INTENT(OUT) :: isorted

      INTEGER :: iptr, itemp  ! index, pointer to smallest value
      REAL(KIND=REAL64) :: temp  ! temporary variable for swapping
      REAL :: cputime = -1.0

      cputlabel = 'selectionSort'
      CALL measure_cpu_time(cputime, cputlabel)

      temp = 0.d0
      sorted = array
      isorted = [ (i, i=1,asize) ]

      ! Sort the data
      DO i = 1, asize - 1
         iptr = i
         DO j = i+1, asize
            IF (sorted(j).LT.sorted(iptr)) THEN
                iptr = j
            END IF
         END DO

         ! swap
         IF ( i.NE.iptr ) THEN
            temp = sorted(i)
            sorted(i) = sorted(iptr)
            sorted(iptr) = temp
            ! indices
            itemp = isorted(i)
            isorted(i) = isorted(iptr)
            isorted(iptr) = itemp
         END IF
      END DO

      CALL measure_cpu_time(cputime, cputlabel)

   END SUBROUTINE selectionSort

   SUBROUTINE quickSortMod (n, arr, sorted, isorted)
   
      ! Quicksort algorithm 
      ! From Numerical Recipes (1992) and (2007) - F77 and C++
      ! Implementation combines parts from both versions
      ! Mod: array is not replaced. Sorted array stored in 'sorted', 
      ! indices in isorted
   
      INTEGER, INTENT(IN) :: n
         ! number of points in array (1D) arr
   
      REAL(KIND=REAL64), DIMENSION(n), INTENT(IN) :: arr
         ! array to sort, will get replaced
   
      INTEGER, DIMENSION(n), INTENT(OUT) :: isorted
         ! array of sorted indices
      REAL(KIND=REAL64), DIMENSION(n), INTENT(OUT) :: sorted
         ! sorted arr
   
      INTEGER, PARAMETER :: M = 7
         ! size of subarray sorted by straight insertion
      INTEGER, PARAMETER :: NSTACK = 64
         ! required auxiliary storage
   
      ! INTEGER :: i, j, k
      INTEGER :: ir, jstack, l, b
      INTEGER, DIMENSION(NSTACK) :: istack
      REAL :: cputime
      REAL(KIND=REAL64) :: a

      cputime = -1.
      cputlabel = 'quicksort'
      CALL measure_cpu_time(cputime, cputlabel)
   
      jstack = 0
      l = 1
      ir = n
      sorted = arr
      isorted = [ (i, i=1,n) ]
   
      outer: DO 
         ! insertion sort when subarray is small enough
         IF ((ir - l).LT.M) THEN
            DO j = l+1, ir
               a = sorted(j)
               b = isorted(j)
               DO i = j-1, l, -1
                  IF (sorted(i).LE.a) EXIT
                  sorted(i+1) = sorted(i)
                  isorted(i+1) = isorted(i)
               END DO
               ! i = l - 1
               sorted(i+1) = a
               isorted(i+1) = b
            END DO
            IF (jstack.EQ.0) EXIT outer
            ir = istack(jstack)
            l = istack(jstack-1)
            jstack = jstack - 2
         ELSE
            k = (l + ir) / 2
            CALL swap(sorted(k), sorted(l+1))
            CALL swap_int(isorted(k), isorted(l+1))
            IF (sorted(l).GT.sorted(ir)) THEN
               CALL swap(sorted(l), sorted(ir))
               CALL swap_int(isorted(l), isorted(ir))
            END IF
            IF (sorted(l+1).GT.sorted(ir)) THEN
               CALL swap(sorted(l+1), sorted(ir))
               CALL swap_int(isorted(l+1), isorted(ir))
            END IF
            IF (sorted(l).GT.sorted(l+1)) THEN
               CALL swap(sorted(l), sorted(l+1))
               CALL swap_int(isorted(l), isorted(l+1))
            END IF
            i = l + 1
            j = ir
            a = sorted(l+1)
            b = isorted(l+1)
            DO
               DO
                  i = i + 1
                  IF (sorted(i).GE.a) EXIT
               END DO
               DO
                  j = j - 1
                  IF (sorted(j).LE.a) EXIT
               END DO
               IF (j.LT.i) EXIT
               CALL swap(sorted(i), sorted(j))
               CALL swap_int(isorted(i), isorted(j))
            END DO
            sorted(l+1) = sorted(j)
            sorted(j) = a
            isorted(l+1) = isorted(j)
            isorted(j) = b
            jstack = jstack + 2
            ! Push pointers to larger subsorteday on stack, process smaller 
            ! subarray immediately
            IF (jstack.GT.NSTACK) WRITE(*, *) "NSTACK too small in sort"
            IF ((ir - i + 1).GE.(j-1)) THEN
               istack(jstack) = ir
               istack(jstack - 1) = i
               ir = j - 1
            ELSE
               istack(jstack) = j - 1
               istack(jstack - 1) = l
               l = i
            END IF
         END IF
      END DO outer

      CALL measure_cpu_time(cputime, cputlabel)
   
   END SUBROUTINE quickSortMod

   SUBROUTINE swap (val1, val2)

      ! swap two values (double precision)

      IMPLICIT NONE

      REAL(KIND=REAL64), INTENT(INOUT) :: val1
      REAL(KIND=REAL64), INTENT(INOUT) :: val2
      REAL(KIND=REAL64) :: t

      t = val2
      val2 = val1
      val1 = t

   END SUBROUTINE swap

   SUBROUTINE swap_int (val1, val2)

      ! swap two values (integers)

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: val1
      INTEGER, INTENT(INOUT) :: val2
      INTEGER :: t

      t = val2
      val2 = val1
      val1 = t

   END SUBROUTINE swap_int

END MODULE sort_mod

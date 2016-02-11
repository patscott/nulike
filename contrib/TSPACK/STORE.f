      DOUBLE PRECISION FUNCTION STORE (X)
      USE omp_lib
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C On input:
C
C       X = Value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       STORE = Value of X after it has been stored and
C               possibly truncated or rounded to the single
C               precision word length.
C
C Modules required by STORE:  None
C
C Modified by Pat Scott Feb 01 2016 to make threadsafe
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
!$omp critical
      Y = X
      STORE = Y
!$omp end critical
      RETURN
      END

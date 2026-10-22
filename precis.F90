!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE DEFINING THE FLOATING NUMBERS PRECISION : DOUBLE OR TRIPLE !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE precis
  IMPLICIT NONE

  !! SINGLE PRECISION 32-BIT !!
!  INTEGER,PARAMETER :: prec = SELECTED_REAL_KIND(p=6,r=37)

  !! DOUBLE PRECISION 64-BIT !!
  ! INTEGER,PARAMETER :: prec = SELECTED_REAL_KIND(p=15,r=307)

  !! QUADRUPLE PRECISION 128-BIT !!
 INTEGER,PARAMETER :: prec = SELECTED_REAL_KIND(p=33,r=4931)

END MODULE precis
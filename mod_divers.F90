!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE CONTAINING DIFFERENT SECONDARY FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mod_divers
  USE mod_var
  IMPLICIT NONE

    REAL(prec), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230

CONTAINS

    !-----------------------------------------!
    ! SKIP kn LINES IN THE READING OF FILE id !
    !-----------------------------------------!
    SUBROUTINE jump(kn,id)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: kn,id
        INTEGER :: ij
        CHARACTER(LEN=1) :: ab

        DO ij=1,kn
        READ(id,*) ab
        print *,ab
        END DO

    END SUBROUTINE jump

    SUBROUTINE unique_sort(sizef,size,final,val)
    implicit none
    integer :: i = 0, min_val, max_val, size, sizef
    integer, dimension(size) :: val, unique
    integer, dimension(:), allocatable :: final

    print *,val
    min_val = minval(val,val>0)-1
    max_val = maxval(val)
    do while (min_val<max_val)
        i = i+1
        min_val = minval(val, mask=val>min_val)
        unique(i) = min_val
    enddo
    IF(ALLOCATED(final)) DEALLOCATE(final)
    allocate(final(i), source=unique(1:i))   !<-- Or, just use unique(1:i) 
    print *, final
    print *,i
    sizef = i

    END SUBROUTINE unique_sort

    FUNCTION det(v1,v2)
        IMPLICIT NONE
        REAL(prec),DIMENSION(2), INTENT(IN) ::v1,v2
        REAL(prec) :: det
        det = v1(1)*v2(2)-v1(2)*v2(1)
    END FUNCTION det

    FUNCTION normalized(vect)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: vect
        REAL(prec), DIMENSION(2) :: normalized

        normalized = vect/norm2(vect)
    END FUNCTION

    FUNCTION norm2(vect)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: vect
        REAL(prec) :: norm2
        norm2 = sqrt(vect(1)**2+vect(2)**2)
    END FUNCTION

    Function CheckIntersect(P1,P2, Q1,Q2)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2) :: P1,P2, Q1,Q2 !(Ax,Bx,Ay,By)
        REAL(prec) :: t,u
        REAL(prec) :: a,b,c,d, k1,k2, det
        LOGICAL  :: CheckIntersect

        CheckIntersect = .FALSE.

        a = P2(1) - P1(1); b = Q1(1) - Q2(1)
        c = P2(2) - P1(1); d = Q1(2) - Q2(2)

        k1= Q1(1) - P1(1); k2= Q1(2) - P1(2)
        det = a*d -b*c
        IF(abs(det) .LT. eps0) THEN 
            CheckIntersect = .FALSE.
            RETURN
        END IF

        t = (k1*d - k2*b)/det;       
        u = (k2*a - k1*c)/det
        print *,t,u
        print *,( eps0 < t .AND. t < 1), (eps0 < u .AND. u < 1)
        IF(( eps0 < t .AND. t < 1) .AND. (eps0 < u .AND. u < 1) ) CheckIntersect = .TRUE.

    END Function CheckIntersect

END MODULE mod_divers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE CONTAINING DIFFERENT SECONDARY FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mod_divers
  USE mod_var
  IMPLICIT NONE

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

END MODULE mod_divers
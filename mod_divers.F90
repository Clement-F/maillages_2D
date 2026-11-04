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

    FUNCTION det(v1,v2,v3)
        IMPLICIT NONE

        REAL(prec),DIMENSION(:), INTENT(IN) ::v1,v2
        REAL(prec),DIMENSION(:), OPTIONAL, INTENT(IN) :: v3
        REAL(prec) :: det

        det = v1(1)*v2(2)-v1(2)*v2(1)
        IF(present(v3)) &
            & det = v1(1)*v2(2)*v3(3) + v1(2)*v2(3)*v3(1) + v1(3)*v2(1)*v3(2) - v1(3)*v2(2)*v3(1) - v2(3)*v3(2)*v1(1) - v3(3)*v1(2)*v2(1)
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
        ! print *,t,u
        ! print *,( eps0 < t .AND. t < 1), (eps0 < u .AND. u < 1)
        IF(( eps0 < t .AND. t < 1) .AND. (eps0 < u .AND. u < 1) ) CheckIntersect = .TRUE.

    END Function CheckIntersect

    FUNCTION Signed_Surface(A,B,C)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C
        REAL(prec) :: Signed_Surface

        Signed_Surface = 0.5_prec*(det(B-A,C-A))
    END FUNCTION Signed_Surface

    FUNCTION Coordonnee_barycentrique(A,B,C,P)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C,P
        REAL(prec), DIMENSION(3) :: Coordonnee_barycentrique

        Coordonnee_barycentrique(1) = Signed_Surface(B,C,P)/Signed_Surface(A,B,C)
        Coordonnee_barycentrique(2) = Signed_Surface(C,A,P)/Signed_Surface(A,B,C)
        Coordonnee_barycentrique(3) = Signed_Surface(A,B,P)/Signed_Surface(A,B,C)

    END FUNCTION Coordonnee_barycentrique

    FUNCTION Is_inside_triangle(A,B,C,P)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C,P
        LOGICAL :: Is_inside_triangle

        Is_inside_triangle = .FALSE.
        IF(ALL(Coordonnee_barycentrique(A,B,C,P) >= 0.0_prec)) Is_inside_triangle = .TRUE.
        IF(ALL(Coordonnee_barycentrique(A,B,C,P) <= 0.0_prec)) STOP "Point undefined"

    END FUNCTION Is_inside_triangle


    SUBROUTINE Calc_center(Msh)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER :: iTri
        REAL(prec), DIMENSION(2) :: A, B, C
        
        IF(ALLOCATED(Msh%Center)) DEALLOCATE(Msh%Center)
        ALLOCATE(Msh%Center(Msh%Nb_Cell,2))

        ! Circumscribed circle Msh%Center
        DO iTri = 1, Msh%Nb_Cell
            ! A = Msh%Coord(Msh%Sommets(iTri,1),:)
            ! B = Msh%Coord(Msh%Sommets(iTri,2),:)
            ! C = Msh%Coord(Msh%Sommets(iTri,3),:)
            ! Msh%Center(iTri,:) = Circumcenter(A, B, C)

            ! IF(.NOT. Is_inside_triangle(A, B, C, Msh%Center(iTri,:)) .AND. Msh%CIsOnBdr(iTri) ) THEN
            !     print *, "Circumcenter is outside the triangle for cell ", iTri
            !     Msh%Center(iTri,1) = Sum(Msh%Coord(Msh%Sommets(iTri,1:3),1)) / 3.0_prec
            !     Msh%Center(iTri,2) = Sum(Msh%Coord(Msh%Sommets(iTri,1:3),2)) / 3.0_prec
            ! END IF

            Msh%Center(iTri,1) = Sum(Msh%Coord(Msh%Sommets(iTri,1:3),1)) / 3.0_prec
            Msh%Center(iTri,2) = Sum(Msh%Coord(Msh%Sommets(iTri,1:3),2)) / 3.0_prec
        END DO

    END SUBROUTINE Calc_center

    FUNCTION Circumcenter(A, B, C)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A, B, C
        REAL(prec), DIMENSION(2) :: circumcenter
        REAL(prec) :: D

        D = 2.0_prec * (A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))

        circumcenter(1) = ((A(1)**2 + A(2)**2)*(B(2)-C(2)) + (B(1)**2 + B(2)**2)*(C(2)-A(2)) + (C(1)**2 + C(2)**2)*(A(2)-B(2))) / D
        circumcenter(2) = ((A(1)**2 + A(2)**2)*(C(1)-B(1)) + (B(1)**2 + B(2)**2)*(A(1)-C(1)) + (C(1)**2 + C(2)**2)*(B(1)-A(1))) / D

    END FUNCTION Circumcenter

    SUBROUTINE Permute_Vertex_onCell(iVer,jVer,nCell, Msh)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: iVer,jVer,nCell
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER :: loc_i,loc_j

        print *,"permute ",iVer,jVer,'in cell',nCell

        loc_i = FINDLOC(Msh%Sommets(nCell,:),iVer,DIM=1)
        loc_j = FINDLOC(Msh%Sommets(nCell,:),jVer,DIM=1)
        Msh%Sommets(nCell,loc_i) = jVer
        MSh%Sommets(nCell,loc_j) = iVer
        ! print *,Msh%Sommets(nCell,:)

    END SUBROUTINE

    SUBROUTINE Insert_Vertex_onCell(iVer,jVer,nCell,Msh,bef_aft)   ! Insert vertex i behind Vertex j in Sommets(nCell,:)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: iVer,jVer,nCell
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER, DIMENSION(Msh%NbSom(nCell)) :: temp_Sommets
        INTEGER :: loc_j
        Character(len=3), INTENT(IN), OPTIONAL :: bef_aft

        ! print *,"insert ",iVer,jVer,'in cell',nCell
        ! write(*,'("insert ",i3,1x,i3,1x,"in cell",i3)') iVer,jVer,nCell
        ! print *,Msh%Sommets(nCell,1:Msh%NbSom(nCell))
        temp_Sommets = Msh%Sommets(nCell,1:Msh%NbSom(nCell))
        loc_j = FINDLOC(Msh%Sommets(nCell,:),jVer,DIM=1)

        ! print *,loc_j
        IF(bef_aft == 'bef') THEN
            loc_j = loc_j - 1
        END IF
        ! print *,loc_j
        Msh%Sommets(nCell,1:loc_j) = temp_Sommets(1:loc_j)
        Msh%Sommets(nCell,loc_j+1) = iVer
        Msh%Sommets(nCell,loc_j+2:Msh%NbSom(nCell)+1) = temp_Sommets(loc_j+1:Msh%NbSom(nCell))
        Msh%NbSom(nCell) = Msh%NbSom(nCell) + 1

    END SUBROUTINE Insert_Vertex_onCell

    SUBROUTINE IsDelauney(Msh)
        IMPLICIT NONE
        TYPE(MESH) :: Msh
        LOGICAL :: IsDelaun
        INTEGER :: nCell, iCell, iVer,ipVer, jVer,jpVer, tVer
        INTEGER :: i,j

        IsDelaun = .TRUE.

        DO nCell=1,Msh%Nb_Cell
            DO i=1,3
            iCell = Msh%CNeigh(nCell,i)
            iVer  = Msh%Sommets(nCell,mod(i  ,3)+1)
            ipVer = Msh%Sommets(nCell,mod(i+1,3)+1)
            IF(iCell .NE. 0) THEN
                DO j =1,3
                    jVer  = Msh%Sommets(iCell,mod(j  ,3)+1)
                    jpVer = Msh%Sommets(iCell,mod(j+1,3)+1)
                    IF((jVer == iVer .AND. jpVer == ipVer) .OR.  (jVer == ipVer .AND. jpVer == iVer) ) THEN

                        tVer = Msh%Sommets(iCell,mod(j+2,3)+1)
                        IsDelaun = .NOT. IsInCircum(nCell,tVer,Msh)
                        IF(.NOT. IsDelaun) THEN 
                            Print *,"This Mesh is Delauney : ",IsDelaun
                            write(*,fmt='(" Vertex ",i3," Is inside the circum circle of ",i3)') tVer,nCell
                        RETURN
                        END IF
                    END IF

                END DO
            END IF
            END DO
        END DO

        Print *,"This Mesh is Delauney : ",IsDelaun

        Contains 

    END SUBROUTINE IsDelauney

    FUNCTION IsInCircum(nCell,tVer,Msh)
        IMPLICIT NONE
        INTEGER :: nCell, tVer
        TYPE(Mesh) :: Msh
        REAL(prec), DIMENSION(2) :: c1,c2,c3, P
        REAL(prec), DIMENSION(3) :: v1,v2,v3
        LOGICAL :: IsInCircum

        c1 = Msh%Coord(Msh%Sommets(nCell,1),:)
        c2 = Msh%Coord(Msh%Sommets(nCell,2),:)
        c3 = Msh%Coord(Msh%Sommets(nCell,3),:)
        P  = Msh%Coord(tVer,:)

        v1(1:2) = c1-P; v1(3) = (norm2(c1-P))**2
        v2(1:2) = c2-P; v2(3) = (norm2(c2-P))**2
        v3(1:2) = c3-P; v3(3) = (norm2(c3-P))**2

        IF(det(v1,v2,v3) >=0) THEN; IsInCircum = .TRUE.
        ELSE ; IsInCircum = .FALSE.
        END IF

    END FUNCTION IsInCircum

    FUNCTION hash_find(iVer1, iVer2, Msh) RESULT(j_hsh)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: iVer1, iVer2
        TYPE(MESH), INTENT(IN) :: Msh
        INTEGER :: j_hsh, key

        key = iVer1 + iVer2
        ! print *,key
        j_hsh = Msh%Head(key)

        DO WHILE (j_hsh /= 0)
            ! print *,"looping through hash list, j_hsh:", j_hsh
            IF ((Msh%LstObj(j_hsh,1) == iVer1 .AND. Msh%LstObj(j_hsh,2) == iVer2) .OR. &
                (Msh%LstObj(j_hsh,1) == iVer2 .AND. Msh%LstObj(j_hsh,2) == iVer1)) THEN
                RETURN
            END IF
            j_hsh = Msh%LstObj(j_hsh,5)
        END DO

        j_hsh = 0

    END FUNCTION hash_find

    SUBROUTINE hash_add(iVer1, iVer2, iTri, Msh, i_hsh)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: iVer1, iVer2, iTri, i_hsh
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER :: key


        ! i_hsh is a initial guess. it should be 0 by default, 
        ! but if you have already run hash_find, you can bypass the find to imput the element at i_hsh

        ! if(i_hsh==0)  i_hsh = hash_find(hsh,iVer1,iVer2); ! check if the element isn't in the hash_list already

        IF (i_hsh == 0) THEN

            IF(Msh%NbrObj >= Msh%Nb_Edg) THEN
                print *, "Error: Exceeded maximum number of edges."
                STOP
            END IF
            key = iVer1 + iVer2
            Msh%NbrObj = Msh%NbrObj + 1
            Msh%LstObj(Msh%NbrObj,1) = iVer1
            Msh%LstObj(Msh%NbrObj,2) = iVer2
            Msh%LstObj(Msh%NbrObj,3) = iTri
            Msh%LstObj(Msh%NbrObj,5) = Msh%Head(key)
            Msh%Head(key) = Msh%NbrObj

            Msh%Edg(Msh%NbrObj,1) = iVer1    
            Msh%Edg(Msh%NbrObj,2) = iVer2
            Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
            Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1
            Msh%EdgNeigh(Msh%NbrObj,1) = iTri

        ELSE
            Msh%LstObj(i_hsh,4) = iTri
            Msh%EdgNeigh(i_hsh,2) = iTri
        END IF

    END SUBROUTINE hash_add


    FUNCTION IsConvex(nCell, Msh)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nCell
        TYPE(MESH), INTENT(IN) :: Msh
        REAL(prec),DIMENSION(2) :: Coord_i,Coord_ip, Coord_j,Coord_jp, Coord_k
        REAL(prec),DIMENSION(2) :: i_To_j, j_To_k
        LOGICAL :: IsConvex
        INTEGER :: iVer,ipVer,jVer,jpVer,kVer
        INTEGER :: i,j

        IsConvex = .TRUE.

        DO i=1,Msh%NbSom(nCell)

        iVer = Msh%Sommets(nCell,mod(i  ,Msh%NbSom(nCell))+1)
        ipVer= Msh%Sommets(nCell,mod(i+1,Msh%NbSom(nCell))+1)
        Coord_i = Msh%Coord(iVer,:)
        Coord_ip= Msh%Coord(ipVer,:)

        DO j=i,Msh%NbSom(nCell)

        jVer = Msh%Sommets(nCell,mod(j  ,Msh%NbSom(nCell))+1)
        jpVer= Msh%Sommets(nCell,mod(j+1,Msh%NbSom(nCell))+1)
        Coord_j = Msh%Coord(jVer,:)
        Coord_jp= Msh%Coord(jpVer,:)

        IF(CheckIntersect(Coord_i,Coord_ip, Coord_j,Coord_jp)) THEN
            IsConvex = .FALSE.
            RETURN
        END IF

        END DO
        END DO

        DO i=1,Msh%NbSom(nCell)

            iVer = Msh%Sommets(nCell,mod(i  ,Msh%NbSom(nCell))+1)
            jVer = Msh%Sommets(nCell,mod(i+1,Msh%NbSom(nCell))+1)
            kVer = Msh%Sommets(nCell,mod(i+2,Msh%NbSom(nCell))+1)
            Coord_i = Msh%Coord(iVer,:)
            Coord_j = Msh%Coord(jVer,:)
            Coord_k = Msh%Coord(kVer,:)
            i_To_j = Coord_j-Coord_i
            j_To_k = Coord_k-Coord_j

            IF(det(i_To_j,j_To_k) < -eps0) THEN
                IsConvex = .FALSE.
                RETURN
            END IF

        END DO

    END FUNCTION IsConvex

END MODULE mod_divers
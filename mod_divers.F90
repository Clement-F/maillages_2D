!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE CONTAINING DIFFERENT SECONDARY FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mod_divers
  USE mod_var
  IMPLICIT NONE

    REAL(prec), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230

CONTAINS

    SUBROUTINE jump(kn,id)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: kn,id
        INTEGER :: ij
        CHARACTER(LEN=1) :: ab

        DO ij=1,kn
        READ(id,*) ab
        ! print *,ab
        END DO

    END SUBROUTINE jump

    FUNCTION det(v1,v2,v3) ! Deteminant de 2/3 vecteurs de taille 2/3
        IMPLICIT NONE

        REAL(prec),DIMENSION(:), INTENT(IN) ::v1,v2
        REAL(prec),DIMENSION(:), OPTIONAL, INTENT(IN) :: v3
        REAL(prec) :: det

        det = v1(1)*v2(2)-v1(2)*v2(1)
        IF(present(v3)) &
            & det = v1(1)*v2(2)*v3(3) + v1(2)*v2(3)*v3(1) + v1(3)*v2(1)*v3(2) - v1(3)*v2(2)*v3(1) - v2(3)*v3(2)*v1(1) - v3(3)*v1(2)*v2(1)
    END FUNCTION det

    FUNCTION normalized(vect) ! normalise un vecteur de taille 2
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: vect
        REAL(prec), DIMENSION(2) :: normalized

        normalized = vect/norm2(vect)
    END FUNCTION

    FUNCTION norm2(vect)    ! donne la norme d'un vecteur de taille 2
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: vect
        REAL(prec) :: norm2
        norm2 = sqrt(vect(1)**2+vect(2)**2)
    END FUNCTION

    Function CheckIntersect(P1,P2, Q1,Q2)       ! Check if the segment (P1,P2) intersects with the segment (Q1,Q2)
        IMPLICIT NONE
        REAL(prec), DIMENSION(2) :: P1,P2, Q1,Q2 !(Ax,Bx,Ay,By)
        REAL(prec) :: t,u
        REAL(prec) :: a,b,c,d, k1,k2, det
        LOGICAL  :: CheckIntersect

        CheckIntersect = .FALSE.

        a = P2(1) - P1(1); b = Q1(1) - Q2(1)
        c = P2(2) - P1(2); d = Q1(2) - Q2(2)
        k1= Q1(1) - P1(1); k2= Q1(2) - P1(2)

        det = a*d -b*c

        IF(abs(det) .LT. eps0) THEN 
            CheckIntersect = .FALSE.
            RETURN
        END IF

        t = (k1*d - k2*b)/det;       
        u = (k2*a - k1*c)/det

        IF(( eps0 < t .AND. t < 1) .AND. (eps0 < u .AND. u < 1) ) CheckIntersect = .TRUE.

    END Function CheckIntersect

    FUNCTION Signed_Surface(A,B,C)              ! Calcul la surface signée du triangle ABC
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C
        REAL(prec) :: Signed_Surface

        Signed_Surface = 0.5_prec*(det(B-A,C-A))
    END FUNCTION Signed_Surface

    FUNCTION Coordonnee_barycentrique(A,B,C,P)  ! Calcul les coordonnées barycentriques du point P dans le triangle ABC
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C,P
        REAL(prec), DIMENSION(3) :: Coordonnee_barycentrique

        Coordonnee_barycentrique(1) = Signed_Surface(B,C,P)/Signed_Surface(A,B,C)
        Coordonnee_barycentrique(2) = Signed_Surface(C,A,P)/Signed_Surface(A,B,C)
        Coordonnee_barycentrique(3) = Signed_Surface(A,B,P)/Signed_Surface(A,B,C)

    END FUNCTION Coordonnee_barycentrique

    FUNCTION Is_inside_triangle(A,B,C,P)        ! Check if the point P is inside the triangle ABC
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C,P
        LOGICAL :: Is_inside_triangle

        Is_inside_triangle = .FALSE.
        IF(ALL(Coordonnee_barycentrique(A,B,C,P) >= 0.0_prec)) Is_inside_triangle = .TRUE.
        IF(ALL(Coordonnee_barycentrique(A,B,C,P) <= 0.0_prec)) STOP "Point undefined"

    END FUNCTION Is_inside_triangle

    SUBROUTINE Calc_center(Msh)                 ! Calcul les centres des mailles du maillage (circumcenter ou barycentre)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell) :: iTri
        INTEGER :: i
        REAL(prec), DIMENSION(2) :: A, B, C
        
        ! Circumscribed circle Msh%Center
        DO i = 1, Msh%Nb_Cell
            iTri = Msh%Cells(i)

            ! A = Msh%Coord(iTri%Sommets(1),:)
            ! B = Msh%Coord(iTri%Sommets(2),:)
            ! C = Msh%Coord(iTri%Sommets(3),:)
            ! iTri%Center(:) = Circumcenter(A, B, C)

            ! IF(.NOT. Is_inside_triangle(A, B, C, Msh%Center(iTri,:)) .AND. Msh%CIsOnBdr(iTri) ) THEN
            !     print *, "Circumcenter is outside the triangle for cell ", iTri
            !     Msh%Center(iTri,1) = Sum(Msh%Coord(Msh%Sommets(iTri,1:3),1)) / 3.0_prec
            !     Msh%Center(iTri,2) = Sum(Msh%Coord(Msh%Sommets(iTri,1:3),2)) / 3.0_prec
            ! END IF

            iTri%Center(1) = Sum(Msh%Coord(iTri%Sommets(1:3),1)) / 3.0_prec
            iTri%Center(2) = Sum(Msh%Coord(iTri%Sommets(1:3),2)) / 3.0_prec

            Msh%Cells(iTri%index) = iTri
        END DO

    END SUBROUTINE Calc_center

    FUNCTION Circumcenter(A, B, C)              ! Calcul le circumcenter du triangle ABC
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A, B, C
        REAL(prec), DIMENSION(2) :: circumcenter
        REAL(prec) :: D

        D = 2.0_prec * (A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))

        circumcenter(1) = ((A(1)**2 + A(2)**2)*(B(2)-C(2)) + (B(1)**2 + B(2)**2)*(C(2)-A(2)) + (C(1)**2 + C(2)**2)*(A(2)-B(2))) / D
        circumcenter(2) = ((A(1)**2 + A(2)**2)*(C(1)-B(1)) + (B(1)**2 + B(2)**2)*(A(1)-C(1)) + (C(1)**2 + C(2)**2)*(B(1)-A(1))) / D

    END FUNCTION Circumcenter

    SUBROUTINE Permute_Vertex_onCell(iVer,jVer, nCell)  ! Permute vertex i and j in Sommets(nCell,:)
        IMPLICIT NONE
        TYPE(Cell) :: nCell
        INTEGER,INTENT(IN) :: iVer,jVer
        INTEGER :: loc_i,loc_j

        ! print *,"permute ",iVer,jVer,'in cell',nCell

        loc_i = FINDLOC(nCell%Sommets(:),iVer,DIM=1)
        loc_j = FINDLOC(nCell%Sommets(:),jVer,DIM=1)
        nCell%Sommets(loc_i) = jVer
        nCell%Sommets(loc_j) = iVer
        ! print *,Msh%Sommets(nCell,:)

    END SUBROUTINE

    SUBROUTINE Insert_Vertex_onCell(iVer,jVer,nCell,bef_aft)   ! Insert vertex i behind/after Vertex j in Sommets(nCell,:)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: iVer,jVer
        TYPE(Cell) :: nCell
        INTEGER, DIMENSION(nCell%NbSom) :: temp_Sommets
        INTEGER :: loc_j
        Character(len=3), INTENT(IN), OPTIONAL :: bef_aft

        temp_Sommets = nCell%Sommets(1:nCell%NbSom)
        loc_j = FINDLOC(nCell%Sommets(:),jVer,DIM=1)

        IF(bef_aft == 'bef') THEN
            loc_j = loc_j - 1
        END IF

        nCell%Sommets(1:loc_j) = temp_Sommets(1:loc_j)
        nCell%Sommets(loc_j+1) = iVer
        nCell%Sommets(loc_j+2:nCell%NbSom+1) = temp_Sommets(loc_j+1:nCell%NbSom)

        nCell%NbSom = nCell%NbSom + 1

    END SUBROUTINE Insert_Vertex_onCell

    SUBROUTINE IsDelauney(Msh)                  ! Check if the mesh is Delauney, if not print the first vertex inside the circumcircle of a triangle
        IMPLICIT NONE
        TYPE(MESH) :: Msh
        TYPE(Cell) :: nCell
        LOGICAL :: IsDelaun
        INTEGER :: iCell, iVer,ipVer, jVer,jpVer, tVer
        INTEGER :: i,j,n

        IsDelaun = .TRUE.

        DO n=1,Msh%Nb_Cell
            nCell = Msh%Cells(n)
            DO i=1,3
            iCell = nCell%CNeigh(i)
            iVer  = nCell%Sommets(mod(i  ,3)+1)
            ipVer = nCell%Sommets(mod(i+1,3)+1)
            IF(iCell .NE. 0) THEN
                DO j =1,3
                    jVer  = nCell%Sommets(mod(j  ,3)+1)
                    jpVer = nCell%Sommets(mod(j+1,3)+1)
                    IF((jVer == iVer .AND. jpVer == ipVer) .OR.  (jVer == ipVer .AND. jpVer == iVer) ) THEN

                        tVer = nCell%Sommets(mod(j+2,3)+1)
                        IsDelaun = .NOT. IsInCircum(nCell,tVer,Msh)
                        IF(.NOT. IsDelaun) THEN 
                            Print *,"This Mesh is Delauney : ",IsDelaun
                            ! write(*,fmt='(" Vertex ",i3," Is inside the circum circle of ",i3)') tVer,nCell
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

    FUNCTION IsInCircum(nCell,tVer,Msh)         ! Check if the vertex tVer is inside the circumcircle of the triangle nCell
        IMPLICIT NONE
        INTEGER,    INTENT(IN) :: tVer
        TYPE(Mesh), INTENT(IN) :: Msh
        TYPE(Cell), INTENT(IN) :: nCell
        REAL(prec), DIMENSION(2) :: c1,c2,c3, P
        REAL(prec), DIMENSION(3) :: v1,v2,v3
        LOGICAL :: IsInCircum

        c1 = Msh%Coord(nCell%Sommets(1),:)
        c2 = Msh%Coord(nCell%Sommets(2),:)
        c3 = Msh%Coord(nCell%Sommets(3),:)
        P  = Msh%Coord(tVer,:)

        v1(1:2) = c1-P; v1(3) = (norm2(c1-P))**2
        v2(1:2) = c2-P; v2(3) = (norm2(c2-P))**2
        v3(1:2) = c3-P; v3(3) = (norm2(c3-P))**2

        IF(det(v1,v2,v3) >=0) THEN; IsInCircum = .TRUE.
        ELSE ; IsInCircum = .FALSE.
        END IF

    END FUNCTION IsInCircum

    FUNCTION hash_find(iVer1, iVer2, Msh) RESULT(j_hsh)     ! Give the index of the edge (iVer1,iVer2) in the hash table, or 0 if not found
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

    SUBROUTINE hash_add(iVer1, iVer2, iTri, Msh, i_hsh)     ! Add the edge (iVer1,iVer2) with triangle iTri to the hash table, or update the second triangle assigned to the edge if it already exists
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: iVer1, iVer2, iTri, i_hsh
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER :: key


        ! i_hsh is a initial guess. it should be 0 by default, 
        ! but if you have already run hash_find, you can bypass the find to imput the element at i_hsh

        ! if(i_hsh==0)  i_hsh = hash_find(hsh,iVer1,iVer2); ! check if the element isn't in the hash_list already

        IF (i_hsh == 0) THEN

            IF(Msh%NbrObj >= Msh%Nb_Edg) THEN
                print *, "Error: Exceeded maximum number of edges. Nb edges :",Msh%NbrObj 
                print *, "Processed ",iTri," Numbers of Cells out of ", Msh%Nb_Cell
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

    FUNCTION IsConvex(nCell, Msh)               ! Check if the cell nCell is convex
        IMPLICIT NONE
        TYPE(Cell), INTENT(IN) :: nCell
        TYPE(MESH), INTENT(IN) :: Msh
        REAL(prec),DIMENSION(2) :: Coord_i,Coord_ip, Coord_j,Coord_jp, Coord_k
        REAL(prec),DIMENSION(2) :: i_To_j, j_To_k
        LOGICAL :: IsConvex
        INTEGER :: iVer,ipVer,jVer,jpVer,kVer
        INTEGER :: i,j

        IsConvex = .TRUE.

        DO i=1,nCell%NbSom

        iVer = nCell%Sommets(mod(i  ,nCell%NbSom)+1)
        ipVer= nCell%Sommets(mod(i+1,nCell%NbSom)+1)
        Coord_i = Msh%Coord(iVer,:)
        Coord_ip= Msh%Coord(ipVer,:)

        DO j=i,nCell%NbSom

        jVer = nCell%Sommets(mod(j  ,nCell%NbSom)+1)
        jpVer= nCell%Sommets(mod(j+1,nCell%NbSom)+1)
        Coord_j = Msh%Coord(jVer,:)
        Coord_jp= Msh%Coord(jpVer,:)

        IF(CheckIntersect(Coord_i,Coord_ip, Coord_j,Coord_jp)) THEN
            IsConvex = .FALSE.
            RETURN
        END IF

        END DO
        END DO

        DO i=1,nCell%NbSom

            iVer = nCell%Sommets(mod(i  ,nCell%NbSom)+1)
            jVer = nCell%Sommets(mod(i+1,nCell%NbSom)+1)
            kVer = nCell%Sommets(mod(i+2,nCell%NbSom)+1)
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

    SUBROUTINE Fuse_cells(iCell,jCell,Msh)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(CELL), INTENT(INOUT) :: iCell,jCell
        INTEGER :: index_edg, Common_Edg, size_list
        INTEGER :: Ver_1, Ver_2, Ver_temp, index_temp
        INTEGER, DIMENSION(iCell%NbSom + jCell%NbSom -2) :: Sommets_temp

        print *,'-------------------------'
        ! Choix delete
        print *, "FUSING CELLS ", iCell%index, " WITH ", jCell%index
        IF(iCell%index == jCell%index) STOP "FUSE CELLS ERROR : SAME CELLS"
        ! jCell%index = max(iCell%index,jCell%index); To_Rec_ind= min(iCell%index,jCell%index)
        ! jCell = Msh%Cells(jCell%index); iCell = Msh%Cells(To_Rec_ind)

        ! CNeigh(jEdg) = Cell_Nghb ; Cell2Edg(jEdg) = j_hsh
        Common_Edg = iCell%Cell2Edg(FINDLOC(iCell%CNeigh(:), jCell%index, DIM=1))

        Ver_1 = Msh%LstObj(Common_Edg,1); Ver_2 = Msh%LstObj(Common_Edg,2)

        write(*,fmt='("Fusing cells", i3," with ",i3)')
        print *, Msh%LstObj(Common_Edg,:)

        ! CRef (?)

        ! Sommets
        ! ALLOCATE ???
        print *,"iCell", iCell%Sommets(:)
        print *,"jCell", jCell%Sommets(:)
        print *, Ver_1, Ver_2

        ! index_edg splits the Sommets list of iCell in two [1,..,index_edg] // [index_edg+1,...,NbSom]
        ! find the edge in the Receiving cell
        index_edg =     FINDLOC(iCell%Sommets(:),Ver_1,DIM=1)
        index_edg = min(FINDLOC(iCell%Sommets(:),Ver_2,DIM=1),index_edg)

        IF(index_edg ==1 .AND. (iCell%Sommets(iCell%NbSom) == Ver_1 .OR. iCell%Sommets(iCell%NbSom) == Ver_2) ) index_edg = iCell%NbSom


        print *,"separation edge :", index_edg

        Sommets_temp(1:index_edg) = iCell%Sommets(1:index_edg)
        size_list = index_edg
        index_temp = FINDLOC(jCell%Sommets(:), iCell%Sommets(index_edg), DIM=1)
        ! IF(jCell%Sommets(mod(index_temp, jCell%NbSom)+1) == iCell%Sommets(mod(index_edg, jCell%NbSom)+1)) index_temp = index_temp +1

        ! Add the deleted cells after the edge
        DO WHILE(( Ver_temp .NE. Ver_1 .AND. Ver_temp .NE. Ver_2) .AND. (size_list .LT.  iCell%NbSom + jCell%NbSom -2))
            index_temp = mod(index_temp, jCell%NbSom)+1
            Ver_temp = jCell%Sommets( index_temp)
            print *,"ver temp", Ver_temp
            ! print *, Ver_1,Ver_2, Ver_temp .NE. Ver_1 ,Ver_temp .NE. Ver_2
            
            ! IF( Ver_temp == Ver_1 .OR. Ver_temp == Ver_2) EXIT

            size_list = size_list + 1
            Sommets_temp(size_list) = Ver_temp

            print *,"iter", Sommets_temp
            print *,"added", Ver_temp, "at position", index_temp
        END DO

        print *,"deleted cell received in full"

        ! Add the Receiving cell's missing vertex to the list
        IF(size_list .LT.  iCell%NbSom + jCell%NbSom -2) THEN 
            index_temp = mod(index_edg+1, iCell%NbSom)+1
            Ver_temp = iCell%Sommets( index_temp) 
            size_list = size_list + 1
            Sommets_temp(size_list) = Ver_temp
        END IF
            
        DO WHILE(( Ver_temp .NE. Ver_1 .AND. Ver_temp .NE. Ver_2) .AND. (size_list .LT.  iCell%NbSom + jCell%NbSom -2))
            index_temp = mod(index_temp, iCell%NbSom)+1
            Ver_temp = iCell%Sommets( index_temp)
            print *,"ver temp", Ver_temp

            IF( Ver_temp == Ver_1 .OR. Ver_temp == Ver_2) EXIT
            size_list = size_list + 1
            Sommets_temp(size_list) = Ver_temp

            print *,"iter", Sommets_temp
            print *,"added", Ver_temp, "at position", index_temp
        END DO


        ! Sommets_temp(size_list:) = iCell%Sommets(index_edg:)

        DEALLOCATE(iCell%Sommets)
        ALLOCATE(iCell%Sommets(size_list))
        
        print *, size_list,',',sommets_temp
        iCell%Sommets(:) = Sommets_Temp(1:size_list)


        ! NbSom
        iCell%NbSom = iCell%NbSom + jCell%NbSom -2
        Msh%max_neigh_Cell = max(Msh%max_neigh_Cell, iCell%NbSom)

        print *,'-------------------------'

        Msh%Cells(iCell%index) = iCell
        Msh%Cells(jCell%index) = Void_cell


    END SUBROUTINE


END MODULE mod_divers
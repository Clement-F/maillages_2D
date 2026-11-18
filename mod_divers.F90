!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE CONTAINING DIFFERENT SECONDARY FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mod_divers
  USE mod_var
  IMPLICIT NONE

    ! REAL(prec), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230

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
        REAL(prec), DIMENSION(2) :: subs1,subs2

        subs1 = B-A
        subs2 = C-A
        Signed_Surface = 0.5_prec*(det(subs1,subs2))
    END FUNCTION Signed_Surface

    FUNCTION Coord_Bary(A,B,C,P)  ! Calcul les coordonnées barycentriques du point P dans le triangle ABC
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C,P
        REAL(prec), DIMENSION(3) :: Coord_Bary

        Coord_Bary(1) = Signed_Surface(B,C,P)/Signed_Surface(A,B,C)
        Coord_Bary(2) = Signed_Surface(C,A,P)/Signed_Surface(A,B,C)
        Coord_Bary(3) = Signed_Surface(A,B,P)/Signed_Surface(A,B,C)

    END FUNCTION Coord_Bary


    FUNCTION Coord_Bary_FromTri(iTri,Msh,P)  ! Calcul les coordonnées barycentriques du point P dans le triangle ABC
        IMPLICIT NONE
        TYPE(Cell) :: iTri
        TYPE(Mesh) :: Msh
        REAL(prec), DIMENSION(2), INTENT(IN) :: P
        REAL(prec), DIMENSION(3) :: Coord_Bary_FromTri
        REAL(prec), DIMENSION(2) :: A,B,C

        A = Msh%Coord(iTri%Sommets(1),:)
        B = Msh%Coord(iTri%Sommets(2),:)
        C = Msh%Coord(iTri%Sommets(3),:)

        Coord_Bary_FromTri(1) = Signed_Surface(B,C,P)/Signed_Surface(A,B,C)
        Coord_Bary_FromTri(2) = Signed_Surface(C,A,P)/Signed_Surface(A,B,C)
        Coord_Bary_FromTri(3) = Signed_Surface(A,B,P)/Signed_Surface(A,B,C)

    END FUNCTION Coord_Bary_FromTri


    FUNCTION Is_inside_triangle(A,B,C,P)        ! Check if the point P is inside the triangle ABC
        IMPLICIT NONE
        REAL(prec), DIMENSION(2), INTENT(IN) :: A,B,C,P
        LOGICAL :: Is_inside_triangle

        Is_inside_triangle = .FALSE.
        IF(ALL(Coord_Bary(A,B,C,P) >= 0.0_prec)) Is_inside_triangle = .TRUE.
        IF(ALL(Coord_Bary(A,B,C,P) <= 0.0_prec)) STOP "Point undefined"

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

    SUBROUTINE Calc_Surface(Msh)
        TYPE(Mesh), INTENT(INOUT) :: Msh
        TYPE(Cell) :: nCell
        INTEGER :: n,i
        REAL(prec) :: surf
        REAL(prec), DIMENSION(2) :: A,B,C

        DO n=1,Msh%Nb_Cell
            nCell = Msh%Cells(n)
            C = nCell%Center
            surf =0._prec

            DO i =1,nCell%NbSom
                A = Msh%Coord(nCell%Sommets(mod(i   ,nCell%NbSom)+1),:)
                B = Msh%Coord(nCell%Sommets(mod(i +1,nCell%NbSom)+1),:)
                surf = surf + 0.5_prec*abs(det(B-A,C-A))
            END DO
            nCell%Surface = surf
            Msh%Cells(n) = nCell
        END DO

    END SUBROUTINE Calc_Surface

    SUBROUTINE Calc_Centroid(Msh)
        TYPE(MESH), INTENT(INOUT) :: Msh 
        TYPE(CELL) :: nCell
        INTEGER :: n,i 
        REAL(prec), DIMENSION(2) :: centroid, Ver1,Ver2, normal

        DO n=1,Msh%Nb_Cell
            nCell = Msh%Cells(n)

            DO i =1,nCell%NbSom
                Ver1 = Msh%Coord(nCell%Sommets(mod(i   ,nCell%NbSom)+1),:)
                Ver2 = Msh%Coord(nCell%Sommets(mod(i +1,nCell%NbSom)+1),:)
                normal(1) = Ver2(2) - Ver1(2); normal(2) = Ver1(1) - Ver2(1)

                centroid = centroid + 0.5_prec*(DOT_Product(Ver1**2,normal) +  DOT_Product(Ver2**2,normal) + 4*DOT_Product( 0.5_prec*(Ver1+Ver2),normal))
            END DO
            nCell%centroid = centroid/nCell%Surface
            Msh%Cells(n) = nCell
        END DO

    END SUBROUTINE Calc_Centroid

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
        INTEGER :: Ver_1, Ver_2, Ver_temp, index_temp, size_fused_cell
        INTEGER, DIMENSION(iCell%NbSom + jCell%NbSom -2) :: Sommets_temp

        IF(iCell%index == jCell%index) STOP "FUSE CELLS ERROR : SAME CELLS"
        
        size_fused_cell = iCell%NbSom + jCell%NbSom -2

        ! CNeigh(jEdg) = Cell_Nghb ; Cell2Edg(jEdg) = j_hsh
        Common_Edg = iCell%Cell2Edg(FINDLOC(iCell%CNeigh(:), jCell%index, DIM=1))

        Ver_1 = Msh%LstObj(Common_Edg,1); Ver_2 = Msh%LstObj(Common_Edg,2)

        write(*,fmt='("Fusing cells", i3," with ",i3)') iCell%index, jCell%index
        ! CRef (?)

        ! Sommets
        ! ALLOCATE ???

        ! index_edg splits the Sommets list of iCell in two [1,..,index_edg] // [index_edg+1,...,NbSom]
        ! find the edge in the Receiving cell
        index_edg =     FINDLOC(iCell%Sommets(:),Ver_1,DIM=1)
        index_edg = min(FINDLOC(iCell%Sommets(:),Ver_2,DIM=1),index_edg)

        IF(index_edg ==1 .AND. (iCell%Sommets(iCell%NbSom) == Ver_1 .OR. iCell%Sommets(iCell%NbSom) == Ver_2) ) index_edg = iCell%NbSom

        Sommets_temp(1:index_edg) = iCell%Sommets(1:index_edg)
        size_list = index_edg
        index_temp = FINDLOC(jCell%Sommets(:), iCell%Sommets(index_edg), DIM=1)

        ! Add the deleted cells after the edge
        DO WHILE(( Ver_temp .NE. Ver_1 .AND. Ver_temp .NE. Ver_2) .AND. (size_list .LT.  size_fused_cell))
            index_temp = mod(index_temp, jCell%NbSom)+1
            Ver_temp = jCell%Sommets( index_temp)
            size_list = size_list + 1
            Sommets_temp(size_list) = Ver_temp
        END DO

        ! Add the Receiving cell's missing vertex to the list
        IF(size_list .LT.  iCell%NbSom + jCell%NbSom -2) THEN 
            index_temp = mod(index_edg+1, iCell%NbSom)+1
            Ver_temp = iCell%Sommets( index_temp) 
            size_list = size_list + 1
            Sommets_temp(size_list) = Ver_temp
        END IF
            
        DO WHILE(( Ver_temp .NE. Ver_1 .AND. Ver_temp .NE. Ver_2) .AND. (size_list .LT.  size_fused_cell))
            index_temp = mod(index_temp, iCell%NbSom)+1
            Ver_temp = iCell%Sommets( index_temp)
            print *,"ver temp", Ver_temp

            IF( Ver_temp == Ver_1 .OR. Ver_temp == Ver_2) EXIT
            size_list = size_list + 1
            Sommets_temp(size_list) = Ver_temp
        END DO


        DEALLOCATE(iCell%Sommets)
        ALLOCATE(iCell%Sommets(size_list))
        
        ! print *, size_list,',',sommets_temp
        iCell%Sommets(:) = Sommets_Temp(1:size_list)


        ! NbSom
        iCell%NbSom = size_list
        Msh%max_neigh_Cell = max(Msh%max_neigh_Cell, iCell%NbSom)

        ! 


        Msh%Cells(iCell%index) = iCell
        Msh%Cells(jCell%index) = Void_cell



    END SUBROUTINE

    SUBROUTINE Cell_Delete(iCell,Msh)
        INTEGER, INTENT(IN) :: iCell
        TYPE(MESH), INTENT(INOUT):: Msh
        TYPE(Cell) :: Cell_Neigh
        INTEGER :: iCell_Neigh
        INTEGER :: j,k


        Do j=1,3
            iCell_Neigh = Msh%Cells(iCell)%CNeigh(j); 
            IF(iCell_Neigh .NE. 0) THEN 
                Cell_Neigh  = Msh%Cells(iCell_Neigh)
                DO k=1,3
                IF(Cell_Neigh%CNeigh(k) == iCell ) Msh%Cells(iCell_Neigh)%CNeigh(k) = 0 ! Removing from the list of neighbours
                END DO
            END IF
        END DO

    END SUBROUTINE Cell_Delete


    SUBROUTINE REALLOC_Mesh(Msh, Nb_Cell, Nb_Ver, Nb_Edg)
        IMPLICIT NONE
        TYPE(Mesh), INTENT(INOUT) :: Msh
        INTEGER, INTENT(IN) :: Nb_Cell, Nb_Ver, Nb_Edg
        TYPE(Mesh) :: Msh_Cpy
        TYPE(Cell) :: nCell
        INTEGER :: mCell, mVer, mEdg
        INTEGER :: n

        print *,"before after       :",Msh%Nb_Cell,"Cells, ",Msh%Nb_Ver,"Vertex",Msh%Nb_Edg,"Edges"
        Print *,"Realloc mesh under :", Nb_Cell,"Cells, ",Nb_Ver,"Vertex", Nb_Edg,"Edges"

        CALL Copy_Mesh(Msh,Msh_Cpy)

        mCell = min(Nb_Cell, Msh%Nb_Cell)
        mVer  = min(Nb_Ver , Msh%Nb_Ver )
        mEdg  = min(Nb_Edg , Msh%Nb_Edg )
        ! Msh%Nb_Cell = Nb_Cell; 
        ! Msh%Nb_Ver = Nb_Ver; 
        ! Msh%Nb_Edg = Nb_Edg

        DEALLOCATE(Msh%Cells)
        DEALLOCATE(Msh%VRef,  Msh%NbNeigh, Msh%VNeigh, Msh%Ver2Cell, Msh%Coord)
        DEALLOCATE(Msh%Edg,   Msh%EdgNeigh)
        DEALLOCATE(Msh%Bdr,   Msh%VBdr,     Msh%VIsOnBdr, Msh%EIsOnBdr)
        ! DEALLOCATE(Msh%Head,  Msh%LstObj)


        ! Cell
        ALLOCATE(Msh%Cells(Nb_Cell))
        Msh%Cells = Void_Cell
        DO n =1,mCell
            nCell = Msh_Cpy%Cells(n)
            Msh%Cells(n) = nCell
        END DO

        ! Vertex
        ALLOCATE(Msh%VRef(Nb_Ver),  Msh%NbNeigh(Nb_Ver),  Msh%Coord(Nb_Ver,2), Msh%VNeigh(Nb_Ver,Msh%max_neigh_Ver))
        ALLOCATE(Msh%Ver2Cell(Nb_Ver,Msh%max_neigh_Ver))
        Msh%VRef(1:mVer) = Msh_Cpy%VRef(1:mVer);                    Msh%NbNeigh(1:mVer)=  Msh_Cpy%NbNeigh(1:mVer);      Msh%Coord(1:mVer,:) = Msh_Cpy%Coord(1:mVer,:)
        Msh%VNeigh(1:mVer,:) = Msh_Cpy%VNeigh(1:mVer,:);            Msh%Ver2Cell(1:mVer,:) = Msh_Cpy%Ver2Cell(1:mVer,:)

        ALLOCATE(Msh%Edg(Nb_Edg,2),     Msh%EdgNeigh(Nb_Edg,2))
        Msh%Edg(1:mEdg,:) = Msh_Cpy%Edg(1:mEdg,:);                  Msh%EdgNeigh(1:mEdg,:) = Msh_Cpy%EdgNeigh(1:mEdg,:)

        ALLOCATE(Msh%Bdr(Msh%Nb_EBdr),       Msh%VBdr(Msh%Nb_VBdr),   Msh%VIsOnBdr(Nb_Ver),   Msh%EIsOnBdr(Nb_Edg)) 
        Msh%Bdr = Msh_Cpy%Bdr;               Msh%VBdr = Msh_Cpy%VBdr;        
        Msh%VIsOnBdr(1:mVer) = Msh_Cpy%VIsOnBdr(1:mVer);     Msh%EIsOnBdr(1:mEdg) = Msh_Cpy%EIsOnBdr(1:mEdg)

        ! ALLOCATE(Msh%Head( Msh%SizHead), Msh%LstObj(Msh%Nb_Edg,5))
        ! Msh_Cpy%Head = Msh%Head;                    Msh_Cpy%LstObj = Msh%LstObj
    END SUBROUTINE

    SUBROUTINE ALLOCATE_Mesh(Msh,Nb_Cell, Nb_Ver, Nb_Edg)
        IMPLICIT NONE
        TYPE(Mesh), INTENT(INOUT) :: Msh
        INTEGER, INTENT(IN) :: Nb_Cell, Nb_Ver, Nb_Edg
        
        ALLOCATE(Msh%Cells(Nb_Cell))
        ALLOCATE(Msh%VRef(Nb_Ver),  Msh%NbNeigh(Nb_Ver),  Msh%Coord(Nb_Ver,2), Msh%VNeigh(Nb_Ver,Msh%max_neigh_Ver))
        ALLOCATE(Msh%Ver2Cell(Nb_Ver,Msh%max_neigh_Ver))
        ALLOCATE(Msh%Edg(Nb_Edg,2),     Msh%EdgNeigh(Nb_Edg,2))
        ALLOCATE(Msh%Bdr(Msh%Nb_EBdr),       Msh%VBdr(Msh%Nb_VBdr),   Msh%VIsOnBdr(Nb_Ver),   Msh%EIsOnBdr(Nb_Edg)) 

    END SUBROUTINE

    SUBROUTINE Copy_Mesh(Msh,Msh_Cpy)
        TYPE(Mesh), INTENT(IN) :: Msh
        TYPE(Mesh), INTENT(OUT):: Msh_Cpy
        TYPE(Cell) :: nCell
        INTEGER :: n

        ! Mesh
        Msh_Cpy%DIM     = Msh%DIM
        Msh_Cpy%Nb_Cell = Msh%Nb_Cell
        Msh_Cpy%Nb_Ver = Msh%Nb_Ver
        Msh_Cpy%Nb_Edg = Msh%Nb_Edg
        Msh_Cpy%max_neigh_Cell = Msh%max_neigh_Cell
        Msh_Cpy%max_neigh_Ver = Msh%max_neigh_Ver
        Msh_Cpy%x_deb = Msh%x_deb; Msh_Cpy%x_fin = Msh%x_fin; Msh_Cpy%y_deb = Msh%y_deb; Msh_Cpy%y_fin = Msh%y_fin; 

        Msh_Cpy%SizHead = 2*Msh_Cpy%Nb_Ver
        Msh_Cpy%NbrObj = 0;

        Msh_Cpy%Nb_CBdr = Msh%Nb_CBdr
        Msh_Cpy%Nb_VBdr = Msh%Nb_VBdr
        Msh_Cpy%Nb_EBdr = Msh%Nb_EBdr

        ! Cell
        ALLOCATE(Msh_Cpy%Cells(Msh_Cpy%Nb_Cell))
        ! Msh_Cpy%Cells = Msh%Cells
        DO n =1,Msh%Nb_Cell
            nCell = Msh%Cells(n)
            Msh_Cpy%Cells(n) = nCell
        END DO

        ! Vertex
        ALLOCATE(Msh_Cpy%VRef(Msh_Cpy%Nb_Ver),      Msh_Cpy%NbNeigh(Msh_Cpy%Nb_Ver),    Msh_Cpy%Coord(Msh_Cpy%Nb_Ver,2))
        ALLOCATE(Msh_Cpy%VNeigh(Msh_Cpy%Nb_Ver,Msh_Cpy%max_neigh_Ver),   Msh_Cpy%Ver2Cell(Msh_Cpy%Nb_Ver,Msh_Cpy%max_neigh_Ver))
        Msh_Cpy%VRef = Msh%VRef;                    Msh_Cpy%NbNeigh=  Msh%NbNeigh;      Msh_Cpy%Coord = Msh%Coord
        Msh_Cpy%VNeigh = Msh%VNeigh;                Msh_Cpy%Ver2Cell = Msh%Ver2Cell

        ALLOCATE(Msh_Cpy%Edg(Msh_Cpy%Nb_Edg,2),     Msh_Cpy%EdgNeigh(Msh_Cpy%Nb_Edg,2))
        Msh_Cpy%Edg = Msh%Edg;                      Msh_Cpy%EdgNeigh = Msh%EdgNeigh

        ! print *, Msh_Cpy%Nb_Edg
        ALLOCATE(Msh_Cpy%Bdr(Msh_Cpy%Nb_EBdr),       Msh_Cpy%VBdr(Msh_Cpy%Nb_VBdr),   Msh_Cpy%VIsOnBdr(Msh_Cpy%Nb_Ver),   Msh_Cpy%EIsOnBdr(Msh_Cpy%Nb_Edg)) 
        Msh_Cpy%Bdr= Msh%Bdr;                           Msh_Cpy%VBdr = Msh%VBdr;        Msh_Cpy%VIsOnBdr = Msh%VIsOnBdr;     Msh_Cpy%EIsOnBdr = Msh%EIsOnBdr

        ! ALLOCATE(Msh_Cpy%Head( Msh_Cpy%SizHead),    Msh_Cpy%LstObj(Msh%Nb_Edg,5))
        ! Msh_Cpy%Head = Msh%Head;                    Msh_Cpy%LstObj = Msh%LstObj

    END SUBROUTINE Copy_Mesh

    SUBROUTINE Out_Cell(inCell)
        IMPLICIT NONE
        TYPE(Cell), INTENT(IN) :: inCell
        CHARACTER (len=63):: fmt_
        INTEGER :: i

        write(*,fmt='("----------- cell ",i3," -----------")') inCell%index
        write(*,fmt='("Ref of Cell :", i3)') inCell%CRef
        write(*,fmt='("Nb of vertex:", i3)') inCell%NbSom

        fmt_ = '("Sommets :"'//Repeat(",i4",inCell%NbSom)//")"
        write(*,fmt=fmt_) inCell%Sommets(:)
        fmt_ = '("Voisins :"'//Repeat(",i4",inCell%NbSom)//")"
        write(*,fmt=fmt_) inCell%CNeigh(:)
        fmt_ = '("Edges :"'//Repeat(",i4",inCell%NbSom)//")"
        write(*,fmt=fmt_) inCell%Cell2Edg(:)

        write(*,fmt='("Centre :",f8.3,";",f8.3)') inCell%Center
        write(*,fmt='("Surface:",f8.3)') inCell%Surface




    END SUBROUTINE 

END MODULE mod_divers  
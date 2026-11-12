MODULE mod_mesh

    use mod_var
    use mod_divers

    IMPLICIT NONE
    
    INTEGER, DIMENSION(0:2,2), PARAMETER :: tri2edg = reshape([1,2, 2,3, 3,1], [3,2])

    Contains


    SUBROUTINE ReadMesh(namefile, Msh)
        ! ----------------------------------------------------------------!
        !  READ MESH FROM .MESH FILE
        !  Msh%max_neigh_Cell = 3 (triangles)
        !  besoin de fixé Msh%max_neigh_Cell avant l'appel de cette fonction
        ! ----------------------------------------------------------------!


        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        CHARACTER(LEN=15) :: blank
        REAL(prec) :: Rtemp
        TYPE(Cell) :: nCell
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER :: io
        INTEGER :: i, dim, itemp

        print *, "Reading mesh file:", TRIM(namefile)

        OPEN(unit=10, file=TRIM(namefile), status='old', action='read', iostat=io)
        IF (io /= 0) THEN
            print *, "Error opening file. IOSTAT =", io
            STOP
        END IF
        
        IF (ALLOCATED(Msh%Vref)) DEALLOCATE(Msh%Vref)
        IF (ALLOCATED(Msh%Coord)) DEALLOCATE(Msh%Coord)
        IF (ALLOCATED(Msh%Cells)) DEALLOCATE(Msh%Cells)

        CALL JUMP(2,10)
        READ(10,*) dim 
        CALL JUMP(1,10)

        READ(10,*) Msh%Nb_Ver
        print *, "Number of vertices:", Msh%Nb_Ver
        ALLOCATE(Msh%Vref(Msh%Nb_Ver), Msh%Coord(Msh%Nb_Ver, 2), Msh%NbNeigh(Msh%Nb_Ver))
        IF(dim==2)READ(10,*) (Msh%Coord(i,1:2),         Msh%Vref(i),i=1,Msh%Nb_Ver)  
        IF(dim==3)READ(10,*) (Msh%Coord(i,1:2),Rtemp,   Msh%Vref(i),i=1,Msh%Nb_Ver)  

        READ(10,*) blank
        IF(TRIM(blank)=="Edges") THEN
            READ(10,*) Msh%Nb_Ebdr
            ALLOCATE(Msh%Bdr(Msh%Nb_Ebdr), Msh%Edg(Msh%Nb_Ebdr,2))
            CALL JUMP(Msh%Nb_Ebdr,10)
            ! DO i = 1, Msh%Nb_Ebdr
            !     READ(10,*) itemp
            !     ! Msh%Bdr(i)=i
            ! END DO
            
            READ(10,*) blank
        END IF

        IF(TRIM(blank)=="Cells" .OR. TRIM(blank)=='Triangles') THEN 
            READ(10,*) Msh%Nb_Cell
            print *, "Number of cells:", Msh%Nb_Cell

            ALLOCATE(Msh%Cells(Msh%Nb_Cell))

            Msh%Cells(:) = Void_cell

            DO i = 1, Msh%Nb_Cell
                nCell =  Msh%Cells(i)

                ALLOCATE(nCell%Sommets(3), nCell%CNeigh(3), nCell%Cell2Edg(3))
                nCell%NbSom =3; nCell%index = i

                READ(10,*) nCell%Sommets(1:3),nCell%cref
                nCell%CNeigh(:) =0; nCell%Cell2Edg(:) =0
                Msh%Cells(i) = nCell

            END DO
        END IF

        print *, Msh%Cells(1)%Sommets(:)

        ! CALCUL DES VOISINS
        print *, "Calculating neighbors..."
        CALL msh_neighbors_Tri(Msh)

        ! CALCUL DES CENTRES DES MAILLES
        print *, "Calculating cell centers..."
        CALL Calc_center(Msh)

    END SUBROUTINE ReadMesh

    SUBROUTINE Readvtk(namefile,Msh)
        ! ----------------------------------------------------------------!
        !  READ MESH FROM .VTK FILE en format ASCII, POLYDATA
        !  Msh%max_neigh_Cell = 3 (triangles)
        !  besoin de fixé Msh%max_neigh_Cell avant l'appel de cette fonction
        ! ----------------------------------------------------------------!
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell) :: jCell
        CHARACTER(LEN=1) :: blank
        INTEGER :: io
        INTEGER :: i, j, nb_Ele, j_temp
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: POLY, Som_
        INTEGER, ALLOCATABLE, DIMENSION(:)   :: NbSom_
        LOGICAL, ALLOCATABLE, DIMENSION(:)   :: isbdr
        

        print *, "Reading VTK file:", TRIM(namefile)
        OPEN(unit=10, file=TRIM(namefile), status='old', action='read', iostat=io)
        IF (io /= 0) THEN
            print *, "Error opening file. IOSTAT =", io
            STOP
        END IF
        
        IF (ALLOCATED(Msh%Vref)) DEALLOCATE(Msh%Vref)
        IF (ALLOCATED(Msh%Coord)) DEALLOCATE(Msh%Coord)
        IF (ALLOCATED(Msh%NbNeigh)) DEALLOCATE(Msh%NbNeigh)
        IF (ALLOCATED(Msh%Ver2Cell)) DEALLOCATE(Msh%Ver2Cell)
        IF (ALLOCATED(Msh%Cells))   DEALLOCATE(Msh%Cells)

        CALL JUMP(4,10)  
        
        READ(10,*) blank, Msh%Nb_Ver, blank
        print *,Msh%Nb_Ver
        ALLOCATE(Msh%Vref(Msh%Nb_Ver), Msh%Coord(Msh%Nb_Ver, 2), Msh%NbNeigh(Msh%Nb_Ver))
        DO i = 1, Msh%Nb_Ver
            READ(10,*) Msh%Coord(i,1:2)
        END DO

        ! CALL JUMP(1,10)
        READ(10,*) blank, nb_Ele, j_temp
        print *, nb_Ele, j_temp
        ALLOCATE(POLY(nb_Ele, Msh%max_neigh_Cell), isbdr(nb_Ele))
        ALLOCATE(Som_(nb_Ele,Msh%max_neigh_Cell))
        ALLOCATE(NbSom_(nb_Ele))

        NbSom_ = 0
        Msh%Nb_EBdr =0

        DO i = 1, nb_Ele
            READ(10,*) NbSom_(i), Som_(i,1:NbSom_(i))
            IF(NbSom_(i) .GT. 1)  Som_(i,1:NbSom_(i)) = Som_(i,1:NbSom_(i)) +1

            IF(NbSom_(i) == 2) THEN
                isbdr(i) = .TRUE.
            ELSE
                isbdr(i) = .FALSE.
            END IF
        END DO

        Msh%Nb_EBdr = COUNT(isbdr)
        Msh%Nb_Cell = nb_Ele - Msh%Nb_EBdr

        ALLOCATE(Msh%Cells(Msh%nb_Cell))

        j = 1
        DO i = 1,nb_Ele; IF(.not. isbdr(i)) THEN
            jCell = Msh%Cells(j)
            jCell%NbSom = NbSom_(i); 
            ALLOCATE(jCell%Sommets(jCell%NbSom))
            jCell%Sommets(:) = Som_(i,:)
            j = j + 1
        END IF; END DO

        ! CALCUL DES VOISINS
        print *, "Calculating neighbors..."
        CALL msh_neighbors_Tri(Msh)

        ! CALCUL DES CENTRES DES MAILLES
        print *, "Calculating cell centers..."
        CALL Calc_center(Msh)       

    END SUBROUTINE Readvtk

    SUBROUTINE WriteMesh(namefile, Msh)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell) :: iCell

        INTEGER :: i
        CHARACTER(LEN=1) :: blank

        blank =""
        OPEN(unit=20, file=TRIM(namefile), status='unknown', action='write')

        WRITE(20,'(A)') "MeshVersionFormatted 2"
        WRITE(20,'(A)') blank
        WRITE(20,'(A)') "Dimension 2"
        WRITE(20,'(A)') blank

        WRITE(20,'(A)') "Vertices"
        WRITE(20,'(I10)') Msh%NB_Ver
        DO i = 1, Msh%NB_Ver
            WRITE(20,'(2F15.6,I10)') Msh%Coord(i,:), Msh%Vref(i)
        END DO
        WRITE(20,'(A)') blank

        WRITE(20,'(A)') "Triangles"
        WRITE(20,'(I10)') Msh%Nb_Cell
        DO i = 1, Msh%Nb_Cell
            iCell = Msh%Cells(i)
            WRITE(20,'(3I10,I10)') iCell%Sommets(:), iCell%Cref
        END DO
        WRITE(20,'(A)') blank
        

    END SUBROUTINE WriteMesh

    SUBROUTINE Write_vtk(fichier,Msh)
        !******** Creation du fichier paraview qui represente la grille ************
        !** Sens des variables : 
        !       fichier = nom du fichier dans lequel mettre le maillage
        !       Msh%NB_Ver = nombre de Msh%Sommets
        !       xs,ys = coordonnées des Msh%Sommets
        !       Msh%Nb_Cell = nombre de mailles
        !       Msh%Sommets = tableau des Msh%Sommets de chaque maille
        !       Msh%NbBdr = nombre d'aretes sur le bord du maillage
        !       Msh%Bdr = tableau des aretes sur le bord du maillage
        !** Sortie : fichier paraview en format .vtk (ASCII, POLYDATA)
        !**************************************

        CHARACTER(LEN=*) :: fichier
        Character(LEN=1) :: blank

        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell) :: iCell

        INTEGER :: i, j

        blank = ""
        OPEN(unit=21, File=fichier, status="unknown")
        WRITE(21,"(A)") "# vtk DataFile Version 2.0"
        WRITE(21,"(A)") "Grille"
        WRITE(21,"(A)") "ASCII"

        WRITE(21,"(A)") "DATASET POLYDATA"
        WRITE(21,"(A,I10,X,A)") "POINTS",Msh%NB_Ver,"float"
        DO i=1, Msh%NB_Ver
        WRITE(21,"(F16.8,F16.8,F16.8)") Msh%Coord(i,1:2),0D0
        ENDDO
        write(21,"(A)") blank

        ! WRITE(21,"(A,I10,X,I10)") "POLYGONS",Msh%Nb_Cell+Msh%Nb_EBdr,Msh%Nb_Cell*4+Msh%Nb_EBdr*3
        ! WRITE(21,"(A,I10,X,I10)") "POLYGONS",Msh%Nb_Cell,Msh%Nb_Cell*4
        WRITE(21,"(A,I10,X,I10)") "POLYGONS",Msh%Nb_Cell,Sum(Msh%Cells%NbSom)+ Msh%Nb_Cell
        
        
        ! bcl_edg :DO i=1, Msh%Nb_EBdr
        ! IF(Msh%Bdr(i) == 0) Exit bcl_edg
        ! IF(Msh%Edg(Msh%Bdr(i),1) == 0 .AND. Msh%Edg(Msh%Bdr(i),2) ==0) Exit bcl_edg
        ! WRITE(21,*) 2,(Msh%Edg(Msh%Bdr(i),j)-1,j=1,2)
        ! ENDDO bcl_edg

        DO i=1, Msh%Nb_Cell
            iCell = Msh%Cells(i)
        WRITE(21,*) iCell%NbSom,(iCell%Sommets(j)-1,j=1,iCell%NbSom)
        ENDDO

        CLOSE(21)

    END SUBROUTINE Write_vtk

    SUBROUTINE msh_neighbors_Tri(Msh)

        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh

        TYPE(Cell) :: iTri, jTri
        INTEGER :: iEdg, jEdg
        INTEGER :: iVer1, iVer2, jVer1, jVer2
        INTEGER :: i,j,k,j_hsh
        INTEGER, DIMENSION(0:2,2), PARAMETER :: tri2edg = reshape([1,2, 2,3, 3,1], [3,2])

        Msh%SizHead = 2*Msh%NB_Ver
        Msh%Nb_Edg = Msh%Nb_Ver + Msh%Nb_Cell; Msh%NbrObj = 0;

        print *, "Check allocation"

        IF(ALLOCATED(Msh%Edg))      DEALLOCATE(Msh%Edg)
        IF(ALLOCATED(Msh%EdgNeigh)) DEALLOCATE(Msh%EdgNeigh)
        IF(ALLOCATED(Msh%EIsOnBdr)) DEALLOCATE(Msh%EIsOnBdr)

        ALLOCATE(Msh%LstObj  (Msh%Nb_Edg+1,5),   Msh%Edg     (Msh%Nb_Edg,2), Msh%EdgNeigh(Msh%Nb_Edg,2), Msh%Head(Msh%SizHead+1))
        ALLOCATE(Msh%EIsOnBdr(Msh%Nb_Edg),       Msh%VIsOnBdr(Msh%Nb_Ver))


        Msh%Head = 0;   Msh%EdgNeigh = 0; Msh%LstObj = 0;

        print *, "Calculating neighbors"
        DO i = 1, Msh%Nb_Cell;   
            iTri = Msh%Cells(i)
            DO iEdg = 0, 2

                iVer1 = iTri%Sommets(tri2edg(iEdg,1))
                iVer2 = iTri%Sommets(tri2edg(iEdg,2))

                j_hsh = hash_find(iVer1,iVer2,Msh)

                IF(j_hsh /= 0) THEN !Edge found in hashtable (update)

                    CALL hash_add(iVer1,iVer2,iTri%index,Msh,j_hsh)
                    jTri = Msh%Cells(Msh%LstObj(j_hsh,3))
                    iTri%CNeigh(iEdg+1) = jTri%index
                    iTri%Cell2Edg(iEdg+1) = j_hsh

                    DO jEdg = 0, 2
                        jVer1 = jTri%Sommets(tri2edg(jEdg,1))
                        jVer2 = jTri%Sommets(tri2edg(jEdg,2))
                        IF((jVer1 == iVer1 .AND. jVer2 == iVer2) .OR. (jVer1 == iVer2 .AND. jVer2 == iVer1)) THEN
                            jTri%CNeigh(jEdg+1) = iTri%index
                            jTri%Cell2Edg(jEdg+1) = j_hsh
                        END IF
                    END DO
                    Msh%Cells(jTri%index) = jTri

                ELSE    !Edge not found in hashtable (add)

                    CALL hash_add(iVer1,iVer2,iTri%index,Msh,0)
                END IF
            END DO 
            Msh%Cells(i) = iTri
        END DO

        Msh%EIsOnBdr = Msh%LstObj(:,4) == 0
        Msh%Nb_EBdr = COUNT(Msh%EIsOnBdr)
        
        Msh%max_neigh_Ver = maxval(Msh%NbNeigh)
        print *, "Maximum number of neighbors for a vertex:", Msh%max_neigh_Ver

        ! print *,Msh%Bdr
        
        IF(ALLOCATED(Msh%VNeigh))   DEALLOCATE(Msh%VNeigh)
        IF(ALLOCATED(Msh%Bdr))      DEALLOCATE(Msh%Bdr)

        ALLOCATE(Msh%VNeigh(Msh%Nb_Ver, Msh%max_neigh_Ver), Msh%Bdr(Msh%Nb_EBdr))
        ALLOCATE(Msh%Ver2Cell(Msh%Nb_Ver, Msh%max_neigh_Ver))

        Msh%Ver2Cell = 0;        Msh%NbNeigh = 0
        DO i = 1, Msh%Nb_Cell
            iTri = Msh%Cells(i)
            Msh%NbNeigh(iTri%Sommets(:)) = Msh%NbNeigh(iTri%Sommets(:)) + 1
            DO k = 1,3; Msh%Ver2Cell(iTri%Sommets(k),Msh%NbNeigh(iTri%Sommets(k))) = iTri%index;            END DO
        END DO

        Msh%VIsOnBdr(:) = .FALSE.
        Msh%EIsOnBdr(:) = .FALSE.
        Msh%Cells%CIsOnBdr = .FALSE.
        Msh%NbNeigh = 0;        Msh%Nb_EBdr = 0

        DO iEdg = 1, Msh%Nb_Edg-1

            iVer1 = Msh%Edg(iEdg,1)
            iVer2 = Msh%Edg(iEdg,2)

            Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
            Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1

            Msh%VNeigh(iVer1,Msh%NbNeigh(iVer1)) = iVer2
            Msh%VNeigh(iVer2,Msh%NbNeigh(iVer2)) = iVer1
            IF(Msh%EdgNeigh(iEdg,2) == 0) THEN
                Msh%EIsOnBdr(iEdg) = .TRUE.

                iTri  = Msh%Cells(Msh%EdgNeigh(iEdg,1))
                iTri%CIsOnBdr  = .TRUE.
                Msh%VIsOnBdr(iVer1) = .TRUE.
                Msh%VIsOnBdr(iVer2) = .TRUE.

                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
                Msh%Bdr(Msh%Nb_EBdr) = iEdg
            END IF
        END DO

        Msh%Nb_EBdr = COUNT(Msh%EIsOnBdr)
        Msh%Nb_Vbdr = COUNT(Msh%VIsOnBdr)
        Msh%Nb_CBdr = COUNT(Msh%Cells%CIsOnBdr)

        print *,Msh%Nb_Edg, "edges and", Msh%Nb_EBdr, "boundary edges found."

        IF(ALLOCATED(Msh%VBdr)) DEALLOCATE(Msh%VBdr)
        ALLOCATE(Msh%VBdr(Msh%Nb_Vbdr))

        i = 0
        DO iVer1 = 1,Msh%Nb_Ver
            IF(Msh%VIsOnBdr(iVer1)) THEN 
                i = i + 1
                Msh%Vbdr(i) = iVer1                
                Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) - 1
            END IF
        END DO

        print *, "Neighbors calculated successfully."


    END SUBROUTINE msh_neighbors_Tri

    SUBROUTINE msh_neighbors(Msh)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell) :: iCell,jCell

        ! FUNCTION
        INTEGER :: iEdg, jEdg
        INTEGER :: iVer1, iVer2, jVer1, jVer2
        INTEGER :: i,j

        print *,"using msh_neighbors function"

        Msh%SizHead = 2*Msh%NB_Ver
        Msh%NbrObj = 0;

        print *, "Check allocation"

        IF(ALLOCATED(Msh%Edg))      DEALLOCATE(Msh%Edg)
        IF(ALLOCATED(Msh%EdgNeigh)) DEALLOCATE(Msh%EdgNeigh)
        IF(ALLOCATED(Msh%EIsOnBdr)) DEALLOCATE(Msh%EIsOnBdr)

        ALLOCATE(Msh%LstObj  (Msh%Nb_Edg+1,5),   Msh%Edg     (Msh%Nb_Edg,2), Msh%EdgNeigh(Msh%Nb_Edg,2), Msh%Head(Msh%SizHead+1))
        ALLOCATE(Msh%EIsOnBdr(Msh%Nb_Edg),       Msh%VIsOnBdr(Msh%Nb_Ver))
        Msh%Head = 0;   Msh%EdgNeigh = 0; Msh%LstObj = 0; Msh%EIsOnBdr = .FALSE.

        print *, "Calculating neighbors"
        DO i = 1, Msh%Nb_Cell
            iCell = Msh%Cells(i)
            ! print *, i, iCell%NbSom,'|', iCell%Sommets
            DO iEdg = 1, iCell%NbSom

                iVer1 = iCell%Sommets(mod(iEdg,   iCell%NbSom)+1)
                iVer2 = iCell%Sommets(mod(iEdg+1, iCell%NbSom)+1)

                j = hash_find(iVer1,iVer2,Msh)

                IF(j /= 0) THEN
                    CALL hash_add(iVer1,iVer2,iCell%index,Msh,j)
                    jCell = Msh%Cells(Msh%LstObj(j,3))
                    iCell%CNeigh(iEdg) = jCell%index
                    iCell%Cell2Edg(iEdg) = j

                    ! print *, jCell%index, jCell%NbSom,'|', jCell%Sommets

                    DO jEdg = 1, jCell%NbSom
                        jVer1 = jCell%Sommets( mod(jEdg,   jCell%NbSom)+1)
                        jVer2 = jCell%Sommets( mod(jEdg+1, jCell%NbSom)+1)
                        IF((jVer1 == iVer1 .AND. jVer2 == iVer2) .OR. (jVer1 == iVer2 .AND. jVer2 == iVer1)) THEN
                            jCell%CNeigh(jEdg) = iCell%index
                            jCell%Cell2Edg(jEdg) = j
                        END IF
                    END DO
                    Msh%Cells(jCell%index) = jCell
                ELSE
                    CALL hash_add(iVer1,iVer2,iCell%index,Msh,0)
                END IF
            END DO
            Msh%Cells(i) = iCell
        END DO

        ! OUTS the whole hashtable
        ! Do i =1,Msh%Nb_Edg
        !     print *,Msh%LstObj(i,:)
        ! END DO


        Msh%max_neigh_Ver = maxval(Msh%NbNeigh)  ! Each edge contributes to two vertices, so divide by 2
        Msh%EIsOnBdr = Msh%LstObj(:,4) == 0
        Msh%Nb_EBdr = COUNT(Msh%EIsOnBdr)
        
        
        IF(ALLOCATED(Msh%VNeigh))   DEALLOCATE(Msh%VNeigh)
        IF(ALLOCATED(Msh%Bdr))      DEALLOCATE(Msh%Bdr)
        IF(ALLOCATED(Msh%Ver2Cell)) DEALLOCATE(Msh%Ver2Cell)

        ALLOCATE(Msh%VNeigh(Msh%Nb_Ver, Msh%max_neigh_Ver), Msh%Bdr(Msh%Nb_EBdr))
        ALLOCATE(Msh%Ver2Cell(Msh%Nb_Ver, Msh%max_neigh_Ver))

        Msh%Ver2Cell = 0;        Msh%NbNeigh = 0
        DO i = 1, Msh%Nb_Cell
            iCell = Msh%Cells(i)
            Msh%NbNeigh(iCell%Sommets(1:iCell%NbSom)) = Msh%NbNeigh(iCell%Sommets(1:iCell%NbSom)) + 1
            DO j=1, iCell%NbSom
                Msh%Ver2Cell(iCell%Sommets(j),Msh%NbNeigh(iCell%Sommets(j))) = iCell%index
            END DO
        END DO

        Msh%VIsOnBdr(:) = .FALSE.
        Msh%Cells%CIsOnBdr = .FALSE.

        Msh%NbNeigh = 0;        Msh%Nb_EBdr = 0
        DO iEdg = 1, Msh%Nb_Edg

            iVer1 = Msh%Edg(iEdg,1)
            iVer2 = Msh%Edg(iEdg,2)
            IF(iVer1 == 0 .OR. iVer2 == 0) EXIT
            Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
            Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1

            Msh%VNeigh(iVer1,Msh%NbNeigh(iVer1)) = iVer2
            Msh%VNeigh(iVer2,Msh%NbNeigh(iVer2)) = iVer1
            IF(Msh%EdgNeigh(iEdg,2) == 0) THEN

                iCell  = Msh%Cells(Msh%EdgNeigh(iEdg,1))
                iCell%CIsOnBdr  = .TRUE.
                Msh%VIsOnBdr(iVer1) = .TRUE.
                Msh%VIsOnBdr(iVer2) = .TRUE.

                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
                Msh%Bdr(Msh%Nb_EBdr) = iEdg
            END IF
        END DO

        Msh%Nb_CBdr = COUNT(Msh%Cells%CIsOnBdr)
        Msh%Nb_Vbdr = COUNT(Msh%VIsOnBdr)

        print *,Msh%Nb_Edg, "edges and", Msh%Nb_EBdr, "boundary edges found."

        IF(ALLOCATED(Msh%VBdr)) DEALLOCATE(Msh%VBdr)
        ALLOCATE(Msh%VBdr(Msh%Nb_Vbdr))

        i = 0
        DO iVer1 = 1,Msh%Nb_Ver
            IF(Msh%VIsOnBdr(iVer1)) THEN 
                i = i + 1
                Msh%Vbdr(i) = iVer1                
                Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) - 1
            END IF
        END DO

        print *, "Neighbors calculated successfully."


    END SUBROUTINE msh_neighbors

    SUBROUTINE OUTMESH(Msh)
        IMPLICIT NONE

        TYPE(MESH) :: Msh
        TYPE(Cell) :: iCell
        INTEGER :: i
        INTEGER ::iVer,iEdg

        print *,"Vertex"
        DO iVer=1,Msh%Nb_Ver
            print *,'-----------------------------------------'
            print *,iVer,'on Bdr:',Msh%VIsOnBdr(iVer),'Neighbors :',Msh%VNeigh(iVer,:)
            print *,"between cells :",Msh%Ver2Cell(iVer,:)
            write(*,'("at coord", 2f8.4 )')Msh%Coord(iVer,:)
        END DO

        print *,"========================================================="
        print *,"========================================================="
        print *,"Cells"
        DO i=1,Msh%Nb_Cell
            iCell = Msh%Cells(i)
            print *,'--------------'
            print *,iCell%index,'on Bdr:',iCell%CIsOnBdr,'Neighbors :',iCell%CNeigh(:)
            print *,"made of :",iCell%Sommets(:)
            write(*,'("at coord", 2f8.4 )')iCell%Center(:)
        END DO

    END SUBROUTINE

    SUBROUTINE Voronoi(Msh_P,Msh)
        IMPLICIT NONE
        Type(MESH), INTENT(IN)  :: Msh_P
        TYPE(MESH), INTENT(OUT) :: Msh

        TYPE(Cell) :: iCell,jCell

        INTEGER :: iEdg, ind
        INTEGER ::  Ver_Pi, Ver_Pj, Ver_Di, Ver_Dj      !Ver_Pi = Cell_Di
        INTEGER :: Cell_Pi,Cell_Pj,Cell_Di,Cell_Dj      !Ver_Di = Cell_Pi   Sauf sur le bord

        INTEGER, DIMENSION(3) :: ID
        REAL(prec), DIMENSION(2) :: i_To_j,j_To_k, Coord_j, i_To_jd
        INTEGER :: i,j

        Msh%Nb_Cell         = Msh_P%Nb_Ver
        Msh%NB_Ver          = Msh_P%Nb_Cell + Msh_P%Nb_EBdr + Msh_P%Nb_VBdr
        Msh%Nb_Edg          = Msh%NB_Ver + Msh%Nb_Cell -2 +1 +50 !Caractéristique -2 + Vide + Survive
        Msh%max_neigh_Cell  = Msh_P%max_neigh_Ver
        Msh%max_neigh_Ver   = Msh_P%max_neigh_Cell   
        Msh%Nb_Ebdr         = 2*Msh_P%Nb_Ebdr   ! Chaque arête du bord est divisé en deux autres arête du bord

        ALLOCATE(Msh%Coord(Msh%NB_Ver,2))
        ALLOCATE(Msh%Vref(Msh%NB_Ver),      Msh%NbNeigh(Msh%NB_Ver), Msh%VNeigh(Msh%NB_Ver, Msh%max_neigh_Ver), Msh%Ver2Cell(Msh%NB_Ver, Msh%max_neigh_Ver))
        ALLOCATE(Msh%Edg(Msh%Nb_Edg,2),     Msh%EdgNeigh(Msh%Nb_Edg,2))
        ALLOCATE(Msh%Cells(Msh%Nb_Cell))

        ! assign dual vertex 
        DO i =1,Msh_P%Nb_Cell
            iCell = Msh_P%Cells(i)
            Msh%Vref(i)         = iCell%Cref
            Msh%NbNeigh(i)      = iCell%NbSom
            Msh%VNeigh(i,:)     = iCell%CNeigh(:)
            Msh%Coord(i,:)      = iCell%Center(:)
            Msh%Ver2Cell(i,:)   = iCell%Sommets(:)
        END DO

        ! assign dual Cells
        DO i=1,Msh_P%Nb_Ver
            iCell       = Msh%Cells(i)
            iCell%index = i
            iCell%NbSom = Msh_P%NbNeigh(i)
            iCell%CIsOnBdr = Msh_P%VIsOnBdr(i)
            iCell%CRef      = Msh_P%VRef(i)
            
            IF(iCell%CIsOnBdr) ALLOCATE(iCell%Sommets(iCell%NbSom+3),    iCell%CNeigh(iCell%NbSom+3),  iCell%Cell2Edg(iCell%NbSom+3))
            IF(.NOT. iCell%CIsOnBdr) ALLOCATE(iCell%Sommets(iCell%NbSom),    iCell%CNeigh(iCell%NbSom),  iCell%Cell2Edg(iCell%NbSom))

            iCell%Sommets(1:iCell%NbSom)= Msh_P%Ver2Cell(i,1:iCell%NbSom)
            iCell%CNeigh(1:iCell%NbSom) = Msh_P%VNeigh(i,1:iCell%NbSom)
            iCell%Center(:) = Msh_P%Coord(i,:)
            Msh%Cells(i) = iCell
        END DO

        ! Msh%Vref   (1:Msh_P%Nb_Cell  )  = Msh_P%cref;           Msh%Cref   (1:Msh%Nb_Cell)  = Msh_P%Vref
        ! Msh%VNeigh (1:Msh_P%Nb_Cell,:)  = Msh_P%CNeigh;         Msh%CNeigh (1:Msh%Nb_Cell,:)= Msh_P%VNeigh
        ! Msh%NbNeigh(1:Msh_P%Nb_Cell  )  = Msh_P%NbSom;          Msh%NbSom  (1:Msh%Nb_Cell)  = Msh_P%NbNeigh
        ! Msh%Coord  (1:Msh_P%Nb_Cell,:)  = Msh_P%Center;         Msh%Center (1:Msh%Nb_Cell,:)= Msh_P%Coord       
        ! Msh%Ver2Cell                    = 0   ;                 Msh%Sommets(1:Msh%Nb_Cell,:)= 0
        ! Msh%Ver2Cell(1:Msh_P%Nb_Cell,:) = Msh_P%Sommets;        Msh%Sommets                 = Msh_P%Ver2Cell

        DO i=1,Msh%Nb_Cell
            iCell = Msh%Cells(i)
            CALL Gift_wrapping(iCell,Msh)
        END DO

        DO i=1,Msh_P%Nb_Vbdr
            Ver_Pi = Msh_P%VBdr(i)
            ind = Msh_P%Nb_Cell + i
            Msh%Vref    (ind)   = 0
            Msh%Coord   (ind,:) = Msh_P%Coord(Ver_Pi,:)

            Cell_Di = Ver_Pi;   
            Msh%Cells(Cell_Di)%NbSom =  Msh%Cells(Cell_Di)%NbSom +1 
            Msh%Cells(Cell_Di)%Sommets(Msh%Cells(Cell_Di)%NbSom)= ind
        END DO

        DO i=1,Msh_P%Nb_Ebdr
            iEdg = Msh_P%Bdr(i)
            ind = Msh_P%Nb_Cell + Msh_P%Nb_Vbdr + i

            Msh%Vref    (ind)   = 0
            Msh%Coord   (ind,1) = Sum(Msh_P%Coord(Msh_P%Edg(iEdg,1:2),1))/2._prec
            Msh%Coord   (ind,2) = Sum(Msh_P%Coord(Msh_P%Edg(iEdg,1:2),2))/2._prec

            Ver_Pi = Msh_P%Edg(iEdg,1); 
            Cell_Di = Ver_Pi
            Msh%Cells(Cell_Di)%NbSom =  Msh%Cells(Cell_Di)%NbSom +1 
            Msh%Cells(Cell_Di)%Sommets(Msh%Cells(Cell_Di)%NbSom)= ind

            Ver_Pi = Msh_P%Edg(iEdg,2); 
            Cell_Di = Ver_Pi
            Msh%Cells(Cell_Di)%NbSom =  Msh%Cells(Cell_Di)%NbSom +1 
            Msh%Cells(Cell_Di)%Sommets(Msh%Cells(Cell_Di)%NbSom)= ind
        END DO

        ! -----------------------------
        ! En essence, on essaye d'arranger les sommets dans un sens trigo
        ! Tout en rajoutant à chaque vertex du bord du maillage primaire, ce sommets et les milieu de ces arêtes.
        ! -----------------------------
        ! 
        !   k             kd        jd------i------j    | où i est un sommet du maillage primaire sur le bord
        !   |              |        |              |    |    j/jd est le milieu de l'arête entre i et un autre vertex du bord
        !   j----- i -----jd        kd             k    |    k/kd est un sommet du maillage dual
        ! 
        ! -----------------------------

        DO i=1,Msh_P%Nb_Vbdr

            Ver_Pi = Msh_P%VBdr(i); Cell_Di= Ver_Pi
            iCell = Msh%Cells(Cell_Di)
            ! print *,i,iCell%index,'|', iCell%Sommets
            iCell%NbSom =  iCell%NbSom -3 
            ID(1:3) = iCell%Sommets(iCell%NbSom+1:iCell%NbSom+3)

            j=1
            Ver_Pj = Msh_P%VNeigh(Ver_Pi,j)
            DO WHILE( .NOT. Msh_P%VIsOnBdr(Ver_Pj) )
                Ver_Pj = Msh_P%VNeigh(Ver_Pi,j)
                j = j + 1
            END DO
            iEdg = hash_find(Ver_Pi,Ver_Pj,Msh_P) 
            Cell_Pj = Msh_P%LstObj(iEdg,3); Ver_Dj = Cell_Pj

            Coord_j = (Msh_P%Coord(Ver_Pj,:) + Msh_P%Coord(Ver_Pi,:))/2._prec
            i_To_j  = Coord_j            - Msh_P%Coord(Ver_Pi,:)
            i_To_jd = Msh%Coord(ID(3),:) - Msh_P%Coord(Ver_Pi,:)
            j_To_k = Msh_P%Cells(Cell_Pj)%Center(:) - Coord_j

            iCell%Sommets(iCell%NbSom+1:iCell%NbSom+3) = 0

            IF (det(i_To_j,j_To_k) > eps0 ) THEN
                ! print *, "det(i_To_j,j_To_k) > eps0"
                IF(norm2(i_To_j-i_To_jd) .LT. eps0 )THEN
                    ! print *,"j is jd", ID(3)
                    CALL Insert_Vertex_onCell(ID(2),Ver_Dj,iCell,'aft') ! on place jd 
                    CALL Insert_Vertex_onCell(ID(1),Ver_Dj,iCell,'aft') ! on place i
                    CALL Insert_Vertex_onCell(ID(3),Ver_Dj,iCell,'aft') ! on place j
                ELSE 
                    ! print *,"j is not jd", ID(3)
                    CALL Insert_Vertex_onCell(ID(3),Ver_Dj,iCell,'aft') ! on place j
                    CALL Insert_Vertex_onCell(ID(1),Ver_Dj,iCell,'aft') ! on place i
                    CALL Insert_Vertex_onCell(ID(2),Ver_Dj,iCell,'aft') ! on place jd 
                END IF
            ELSE
                ! print *, "det(i_To_j,j_To_k) <= eps0"
                IF(norm2(i_To_j-i_To_jd) .LT. eps0 )THEN
                    print *,"j is jd", ID(3)
                    CALL Insert_Vertex_onCell(ID(3),Ver_Dj,iCell,'aft') ! on place j
                    CALL Insert_Vertex_onCell(ID(1),Ver_Dj,iCell,'aft') ! on place i
                    CALL Insert_Vertex_onCell(ID(2),Ver_Dj,iCell,'aft') ! on place jd 
                ELSE
                    ! print *,"j is not jd", ID(3)
                    CALL Insert_Vertex_onCell(ID(3),Ver_Dj,iCell,'bef') ! on place j
                    CALL Insert_Vertex_onCell(ID(1),Ver_Dj,iCell,'bef') ! on place i
                    CALL Insert_Vertex_onCell(ID(2),Ver_Dj,iCell,'bef') ! on place jd 
                END IF
            END IF
            msh%Cells(iCell%index) = iCell
        END DO

        print *, "Calculating neighbors for the dual mesh"
        CALL msh_neighbors(Msh)

    END SUBROUTINE

    SUBROUTINE Gift_Wrapping(nCell,Msh)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell), INTENT(INOUT) :: nCell
        REAL(prec),DIMENSION(nCell%NbSom,2) :: Coord_Cell
        REAL(prec),DIMENSION(2) :: Coord_i,Coord_j, Coord_start, vect, vect_next,vect_prev
        INTEGER, DIMENSION(nCell%NbSom) :: Convex_hull
        INTEGER :: Ver_i,Ver_j, Start_Point,Next_Point, Start_ind, Next_ind
        INTEGER :: i,j, NbPoints

        ! print *,"==============="

        IF(nCell%NbSom == 0) THEN
            print *, "Error: Cell", nCell%index, "has no vertices."
            STOP " ERROR Gift_Wrapping"
        END IF

        IF(nCell%NbSom == 1) THEN
            RETURN
        END IF

        IF(nCell%NbSom == 2) THEN  
            vect      = nCell%Center(:) - Msh%Coord(nCell%Sommets(1),:)
            vect_next = Msh%Coord(nCell%Sommets(2),:) - nCell%Center(:)
            IF(det(vect,vect_next) < 0._prec) THEN
                Convex_hull(1:2) = nCell%Sommets(1:2)
            ELSE
                Convex_hull(1:2) = nCell%Sommets(2:1:-1)
            END IF
        END IF

        IF(nCell%NbSom > 2) THEN
        Convex_hull = 0
        Coord_Cell = Msh%Coord(nCell%Sommets(1:nCell%NbSom),1:2)
        
        Coord_start(1) = Minval(Coord_Cell(:,1))
        Start_ind   = FINDLOC(Coord_Cell(:,1),Coord_start(1),DIM=1)
        Start_Point = nCell%Sommets(Start_ind)

        Ver_i           = Start_Point
        Convex_hull(1)  = Start_Point
        Coord_Start     = Msh%Coord(Start_Point,:) 

        vect_prev   = 0._prec
        Next_ind    = mod(FINDLOC(nCell%Sommets(:),Ver_i,DIM=1),nCell%NbSom) + 1
        Next_Point  = nCell%Sommets(Next_ind)

        i = 2;      NbPoints = 1
        ! print *,Start_Point, nCell%NbSom
        ! print *, nCell%Sommets(:)
        ! print *,'a',Msh%Coord(nCell%Sommets(1:nCell%NbSom),:)
        DO WHILE(Next_Point .NE. Start_Point)
            ! print *,'----------------', nCell%index, '--------------'
            Coord_i = Msh%Coord(Ver_i,:)
            
            Next_ind = mod(FINDLOC(nCell%Sommets(:),Ver_i,DIM=1),nCell%NbSom) + 1
            Next_Point = nCell%Sommets(Next_ind)

            Coord_j = Msh%Coord(Next_Point,:)
            vect_Next = Coord_j- Coord_i
            DO j=1,nCell%NbSom
                Ver_j   = nCell%Sommets(j) 
                Coord_j = Msh%Coord(Ver_j,:)
                vect = Coord_j- Coord_i


                ! print *,Ver_i,ver_j, Next_Point,  norm2(vect-vect_prev), norm2(vect-vect_prev) .GT. eps0
                ! IF it is to the left of the previous guess and it's not the same vector as the last
                IF((det(vect,vect_next)<eps0 .AND. Ver_j .NE. Ver_i) .AND.  norm2(vect-vect_prev) .GT. eps0) THEN
                    Next_Point = Ver_j
                    vect_next  = vect
                END IF
            END DO

            vect_prev = vect_next
            Ver_i = Next_Point

            IF(Next_Point == Start_Point) Exit

            Convex_hull(i)= Ver_i
            i = i+1;    NbPoints = NbPoints + 1
        END DO
        
        nCell%Sommets(1:NbPoints) = Convex_hull(1:NbPoints)
        nCell%NbSom = NbPoints
        END IF
        Msh%Cells(nCell%index) = nCell

    END SUBROUTINE Gift_Wrapping

    SUBROUTINE Correct_Mesh(Msh)
        IMPLICIT NONE

        TYPE(MESH), INTENT(INOUT) :: Msh
        TYPE(Cell) :: iCell, jCell
        INTEGER :: i,j,k
        REAL(prec), DIMENSION(2) :: iCenter, jCenter 

        print *,"CLEANING MESH"

        ! Check If two centers are overlapping
        
        Do i =1,Msh%Nb_Cell
            iCell = Msh%Cells(i)
            iCenter = iCell%Center(:)

            DO j=1, iCell%NbSom
                ! print *,i,j,iCell%CNeigh(j)
                IF(iCell%CNeigh(j) .NE. 0) THEN
                jCell = Msh%Cells(iCell%CNeigh(j))
                    jCenter = jCell%Center(:)

                    ! print *, iCell%index,jCell%index

                    ! IF(norm2(jCenter- iCenter) .LT. eps0) THEN
                    !     CALL Fuse_cells(iCell,jCell,Msh)
                    ! END IF
                END IF
            END DO
        END DO

    END SUBROUTINE Correct_Mesh


    
END MODULE mod_mesh
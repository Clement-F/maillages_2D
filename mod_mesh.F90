MODULE mod_mesh

    use mod_var
    use mod_divers

    IMPLICIT NONE
    
    INTEGER, DIMENSION(0:2,2), PARAMETER :: tri2edg = reshape([1,2, 2,3, 3,1], [3,2])

    Contains


    SUBROUTINE ReadMesh(namefile, Msh)
        ! ----------------------------------------------------------------!
        !  READ MESH FROM .MESH FILE
        !  Msh%maxneigh = 3 (triangles)
        !  besoin de fixé Msh%maxneigh avant l'appel de cette fonction
        ! ----------------------------------------------------------------!


        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        CHARACTER(LEN=15) :: blank
        REAL(prec) :: Rtemp
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
        IF (ALLOCATED(Msh%cref)) DEALLOCATE(Msh%cref)
        IF (ALLOCATED(Msh%Sommets)) DEALLOCATE(Msh%Sommets)
        IF (ALLOCATED(Msh%CNeigh)) DEALLOCATE(Msh%CNeigh)

        CALL JUMP(2,10)
        READ(10,*) dim 
        CALL JUMP(1,10)

        READ(10,*) Msh%Nb_Ver
        print *, "Number of vertices:", Msh%Nb_Ver
        ALLOCATE(Msh%Vref(Msh%Nb_Ver), Msh%Coord(Msh%Nb_Ver, 2), Msh%NbNeigh(Msh%Nb_Ver))
        IF(dim==2)READ(10,*) (Msh%Coord(i,1:2),Msh%Vref(i),i=1,Msh%Nb_Ver)  
        IF(dim==3)READ(10,*) (Msh%Coord(i,1:2),Rtemp,Msh%Vref(i),i=1,Msh%Nb_Ver)  

        READ(10,*) blank
        IF(TRIM(blank)=="Edges") THEN
            READ(10,*) Msh%Nb_Ebdr
            ALLOCATE(Msh%Bdr(Msh%Nb_Ebdr), Msh%Edg(Msh%Nb_Ebdr,2))
            DO i = 1, Msh%Nb_Ebdr
                READ(10,*) itemp
                ! Msh%Bdr(i)=i
            END DO
            READ(10,*) blank
        END IF
        IF(TRIM(blank)=="Cells" .OR. TRIM(blank)=='Triangles') THEN 
            READ(10,*) Msh%Nb_Cell
            print *, "Number of cells:", Msh%Nb_Cell

            ALLOCATE(Msh%cref(Msh%Nb_Cell), Msh%Sommets(Msh%Nb_Cell, Msh%maxneigh), Msh%CNeigh(Msh%Nb_Cell, Msh%maxneigh))
            ALLOCATE(Msh%NbSom(Msh%Nb_Cell))
            Msh%CNeigh = 0; Msh%Sommets =0

            DO i = 1, Msh%Nb_Cell
                READ(10,*) Msh%Sommets(i,1:3),Msh%cref(i)
            END DO
        END IF

        Msh%NbSom = 3
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
        !  Msh%maxneigh = 3 (triangles)
        !  besoin de fixé Msh%maxneigh avant l'appel de cette fonction
        ! ----------------------------------------------------------------!
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        TYPE(MESH), INTENT(INOUT) :: Msh
        CHARACTER(LEN=1) :: blank
        INTEGER :: io
        INTEGER :: i, j, i_temp, j_temp
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
        IF (ALLOCATED(Msh%cref)) DEALLOCATE(Msh%cref)
        IF (ALLOCATED(Msh%Sommets)) DEALLOCATE(Msh%Sommets)
        IF (ALLOCATED(Msh%CNeigh)) DEALLOCATE(Msh%CNeigh)
        IF (ALLOCATED(Msh%NbNeigh)) DEALLOCATE(Msh%NbNeigh)
        IF (ALLOCATED(Msh%Ver2Cell)) DEALLOCATE(Msh%Ver2Cell)
        IF (ALLOCATED(Msh%NbSom)) DEALLOCATE(Msh%NbSom)

        CALL JUMP(4,10)  
        
        READ(10,*) blank, Msh%Nb_Ver, blank
        print *,Msh%Nb_Ver
        ALLOCATE(Msh%Vref(Msh%Nb_Ver), Msh%Coord(Msh%Nb_Ver, 2), Msh%NbNeigh(Msh%Nb_Ver))
        DO i = 1, Msh%Nb_Ver
            READ(10,*) Msh%Coord(i,1:2)
        END DO

        ! CALL JUMP(1,10)
        READ(10,*) blank, i_temp, j_temp
        print *, i_temp, j_temp
        ALLOCATE(POLY(i_temp, Msh%maxneigh), isbdr(i_temp))
        ALLOCATE(Som_(i_temp,Msh%maxneigh))
        ALLOCATE(NbSom_(i_temp))

        NbSom_ = 0
        Msh%Nb_EBdr =0
        DO i = 1, i_temp
            READ(10,*) NbSom_(i), Som_(i,1:NbSom_(i))
            IF(NbSom_(i) .GT. 1)  Som_(i,1:NbSom_(i)) = Som_(i,1:NbSom_(i)) +1

            IF(NbSom_(i) == 2) THEN
                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
                isbdr(i) = .TRUE.
            ELSE
                isbdr(i) = .FALSE.
            END IF
        END DO

        Msh%Nb_Cell = i_temp - Msh%Nb_EBdr
        ALLOCATE(Msh%cref(Msh%Nb_Cell), Msh%NbSom(Msh%Nb_Cell))
        ALLOCATE(Msh%Sommets(Msh%Nb_Cell,Msh%maxneigh), Msh%CNeigh(Msh%Nb_Cell,Msh%maxneigh))
        j = 1
        Msh%NbSom = 0
        DO i = 1,i_temp; IF(.not. isbdr(i)) THEN
            Msh%NbSom(j) = NbSom_(i)
            Msh%Sommets(j,:) = Som_(i,:)
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
            WRITE(20,'(3I10,I10)') Msh%Sommets(i,:), Msh%Cref(i)
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

        INTEGER :: i, j

        blank = ""
        OPEN(unit=21, File=fichier, status="unknown")
        WRITE(21,"(A)") "# vtk DataFile Version 2.0"
        WRITE(21,"(A)") "Grille"
        WRITE(21,"(A)") "ASCII"

        WRITE(21,"(A)") "DATASET POLYDATA"
        WRITE(21,"(A,I10,X,A)") "POINTS",Msh%NB_Ver,"float"
        DO i=1, Msh%NB_Ver
        WRITE(21,"(F12.5,F12.5,F12.5)") Msh%Coord(i,1:2),0D0
        ENDDO
        write(21,"(A)") blank
        WRITE(21,"(A,I10,X,I10)") "POLYGONS",Msh%Nb_Cell+Msh%Nb_EBdr,Msh%Nb_Cell*3+Msh%Nb_EBdr*2+Msh%Nb_Cell+Msh%Nb_EBdr
        
        bcl_edg :DO i=1, Msh%Nb_EBdr
        IF(Msh%Bdr(i) == 0) Exit bcl_edg
        IF(Msh%Edg(Msh%Bdr(i),1) == 0 .AND. Msh%Edg(Msh%Bdr(i),2) ==0) Exit bcl_edg
        WRITE(21,*) 2,(Msh%Edg(Msh%Bdr(i),j)-1,j=1,2)
        ENDDO bcl_edg

        DO i=1, Msh%Nb_Cell
        WRITE(21,*) Msh%NbSom(i),(Msh%Sommets(i,j)-1,j=1,Msh%NbSom(i))
        ENDDO

        CLOSE(21)

    END SUBROUTINE Write_vtk

    SUBROUTINE msh_neighbors_Tri(Msh)

        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh

        INTEGER :: iTri, jTri, iEdg, jEdg
        INTEGER :: iVer1, iVer2, jVer1, jVer2
        INTEGER :: i,j_hsh

        INTEGER, DIMENSION(0:2,2), PARAMETER :: tri2edg = reshape([1,2, 2,3, 3,1], [3,2])

        ! HASH TABLE
        INTEGER :: SizHead                                          ! Maximum key value, the key is in [0,SizHead-1]
        INTEGER :: NbrObj                                           ! Number of objects in the hash table
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  Head                 ! Head[key%(SizHead)] = link to the first object having this key in the LstObj list
        INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  LstObj       ! List of objects in the hash table


        SizHead = 2*Msh%NB_Ver
        Msh%Nb_Edg = Msh%Nb_Ver + Msh%Nb_Cell; NbrObj = 0;

        print *, "Check allocation"

        IF(ALLOCATED(Msh%Edg)) DEALLOCATE(Msh%Edg)
        IF(ALLOCATED(Msh%EdgNeigh)) DEALLOCATE(Msh%EdgNeigh)
        IF(ALLOCATED(Msh%EIsOnBdr)) DEALLOCATE(Msh%EIsOnBdr)

        ALLOCATE(LstObj(Msh%Nb_Edg+1,5), Msh%Edg(Msh%Nb_Edg,2), Msh%EdgNeigh(Msh%Nb_Edg,2))
        ALLOCATE(Msh%EIsOnBdr(Msh%Nb_Edg), Msh%VIsOnBdr(Msh%Nb_Ver), Msh%CIsOnBdr(Msh%Nb_Cell))
        ALLOCATE(Head(SizHead+1))
        Head = 0;   Msh%EdgNeigh = 0; LstObj = 0; Msh%EIsOnBdr = .FALSE.

        print *, "Calculating neighbors"

        DO iTri = 1, Msh%Nb_Cell
            ! print *, "Processing triangle:", iTri, "of", Msh%Nb_Cell
            DO iEdg = 0, 2
                iVer1 = Msh%Sommets(iTri, tri2edg(iEdg,1))
                iVer2 = Msh%Sommets(iTri, tri2edg(iEdg,2))

                j_hsh = hash_find(iVer1,iVer2)
                ! print *, "Hash find result for edge (", iVer1, ",", iVer2, "):", j_hsh
                IF(j_hsh /= 0) THEN
                    CALL hash_add(iVer1,iVer2,iTri,j_hsh)

                    jTri = LstObj(j_hsh,3)
                    Msh%CNeigh(iTri,iEdg+1) = jTri
                    DO jEdg = 0, 2
                        jVer1 = Msh%Sommets(jTri, tri2edg(jEdg,1))
                        jVer2 = Msh%Sommets(jTri, tri2edg(jEdg,2))
                        IF((jVer1 == iVer1 .AND. jVer2 == iVer2) .OR. (jVer1 == iVer2 .AND. jVer2 == iVer1)) THEN
                            Msh%CNeigh(jTri,jEdg+1) = iTri
                        END IF
                    END DO
                ELSE
                    CALL hash_add(iVer1,iVer2,iTri,0)
                END IF
            END DO
        END DO

        
        Msh%max_neigh_Ver = maxval(Msh%NbNeigh)  ! Each edge contributes to two vertices, so divide by 2

        print *, "Maximum number of neighbors for a vertex:", Msh%max_neigh_Ver
        Msh%Nb_EBdr = 0
        DO iEdg = 1, Msh%Nb_Edg-1
            ! print *, "Edge", iEdg, "vertex", Msh%Edg(iEdg,1), "to", Msh%Edg(iEdg,2), "neighbors:", Msh%EdgNeigh(iEdg,1), Msh%EdgNeigh(iEdg,2)
            IF(Msh%EdgNeigh(iEdg,2) == 0) THEN
                Msh%EIsOnBdr(iEdg) = .TRUE.
                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
            ELSE
                Msh%EIsOnBdr(iEdg) = .FALSE.
            END IF
        END DO

        print *,Msh%Nb_Edg, "edges and", Msh%Nb_EBdr, "boundary edges found."

        IF(ALLOCATED(Msh%VNeigh)) DEALLOCATE(Msh%VNeigh)
        IF(ALLOCATED(Msh%Bdr)) DEALLOCATE(Msh%Bdr)

        ALLOCATE(Msh%VNeigh(Msh%Nb_Ver, Msh%max_neigh_Ver), Msh%Bdr(Msh%Nb_EBdr))
        ALLOCATE(Msh%Ver2Cell(Msh%Nb_Ver, Msh%max_neigh_Ver))

        Msh%Ver2Cell = 0;        Msh%NbNeigh = 0
        DO iTri = 1, Msh%Nb_Cell
            Msh%NbNeigh(Msh%Sommets(iTri,1:3)) = Msh%NbNeigh(Msh%Sommets(iTri,1:3)) + 1
            ! print *,"tri :",iTri,"composed of ",Msh%Sommets(iTri,:), Msh%NbNeigh(Msh%Sommets(iTri,1:3))
            ! write(*,'(13i3)')Msh%NbNeigh(:)
            print *,"Ver",Msh%Sommets(iTri,1:3),"added ",Msh%NbNeigh(Msh%Sommets(iTri,1:3))
            DO i = 1,3
                Msh%Ver2Cell(Msh%Sommets(iTri,i),Msh%NbNeigh(Msh%Sommets(iTri,i))) = iTri
            END DO

        END DO

        Msh%NbNeigh = 0; Msh%Nb_EBdr = 0
        Msh%VIsOnBdr(:) = .FALSE.
        Msh%CIsOnBdr(:) = .FALSE.

        DO iEdg = 1, Msh%Nb_Edg-1
            iVer1 = Msh%Edg(iEdg,1)
            iVer2 = Msh%Edg(iEdg,2)

            Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
            Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1

            Msh%VNeigh(iVer1,Msh%NbNeigh(iVer1)) = iVer2
            Msh%VNeigh(iVer2,Msh%NbNeigh(iVer2)) = iVer1

            IF(Msh%EIsOnBdr(iEdg)) THEN
                iTri  = Msh%EdgNeigh(iEdg,1)
                Msh%CIsOnBdr(iTri)  = .TRUE.
                Msh%VIsOnBdr(iVer1) = .TRUE.
                Msh%VIsOnBdr(iVer2) = .TRUE.
                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
                write(*,fmt='("(",i3,",",i3,") is on the edge" )') iVer1, iVer2

                Msh%Bdr(Msh%Nb_EBdr) = iEdg
            END IF

        END DO

        Msh%Nb_EBdr = COUNT(Msh%EIsOnBdr)
        Msh%Nb_CBdr = COUNT(Msh%CIsOnBdr)
        Msh%Nb_Vbdr = COUNT(Msh%VIsOnBdr)

        IF(ALLOCATED(Msh%VBdr)) DEALLOCATE(Msh%VBdr)
        ALLOCATE(Msh%VBdr(Msh%Nb_Vbdr))

        iVer2 = 0
        DO iVer1 = 1,Msh%Nb_Ver
            IF(Msh%VIsOnBdr(iVer1)) THEN 
                print *, iVer1," Is On Boundary"
                iVer2 = iVer2 + 1
                Msh%Vbdr(iVer2) = iVer1
            END IF
        END DO

        print *, "Neighbors calculated successfully."

        Contains
        FUNCTION hash_find(iVer1, iVer2) RESULT(j_hsh)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: iVer1, iVer2
            INTEGER :: j_hsh, key

            key = iVer1 + iVer2
            ! print *,key
            j_hsh = Head(key)

            DO WHILE (j_hsh /= 0)
                ! print *,"looping through hash list, j_hsh:", j_hsh
                IF ((LstObj(j_hsh,1) == iVer1 .AND. LstObj(j_hsh,2) == iVer2) .OR. &
                    (LstObj(j_hsh,1) == iVer2 .AND. LstObj(j_hsh,2) == iVer1)) THEN
                    RETURN
                END IF
                j_hsh = LstObj(j_hsh,5)
            END DO

            j_hsh = 0

        END FUNCTION hash_find

        SUBROUTINE hash_add(iVer1, iVer2, iTri, i_hsh)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: iVer1, iVer2, iTri, i_hsh
            INTEGER :: key


            ! i_hsh is a initial guess. it should be 0 by default, 
            ! but if you have already run hash_find, you can bypass the find to imput the element at i_hsh

            ! if(i_hsh==0)  i_hsh = hash_find(hsh,iVer1,iVer2); ! check if the element isn't in the hash_list already

            IF (i_hsh == 0) THEN

                IF(NbrObj >= Msh%Nb_Edg) THEN
                    print *, "Error: Exceeded maximum number of edges."
                    STOP
                END IF
                key = iVer1 + iVer2
                NbrObj = NbrObj + 1
                LstObj(NbrObj,1) = iVer1
                LstObj(NbrObj,2) = iVer2
                LstObj(NbrObj,3) = iTri
                LstObj(NbrObj,5) = Head(key)
                Head(key) = NbrObj

                Msh%Edg(NbrObj,1) = iVer1    
                Msh%Edg(NbrObj,2) = iVer2
                Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
                Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1
                Msh%EdgNeigh(NbrObj,1) = iTri

            ELSE
                LstObj(i_hsh,4) = iTri
                Msh%EdgNeigh(i_hsh,2) = iTri
            END IF

        END SUBROUTINE hash_add


    END SUBROUTINE msh_neighbors_Tri

    SUBROUTINE msh_neighbors(Msh)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh

        ! FUNCTION
        INTEGER :: iPol, jPol, iEdg, jEdg
        INTEGER :: iVer1, iVer2, jVer1, jVer2
        INTEGER :: i,j

        INTEGER :: MaxSom
        INTEGER,DIMENSION(0:maxval(Msh%NbSom(:))-1,2) :: Cell2Edg

        ! HASH TABLE
        INTEGER :: SizHead                                          ! Maximum key value, the key is in [0,SizHead-1]
        INTEGER :: NbrObj                                           ! Number of objects in the hash table
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  Head                 ! Head[key%(SizHead)] = link to the first object having this key in the LstObj list
        INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  LstObj       ! List of objects in the hash table


        print *,"using msh_neighbors function"

        SizHead = 3*Msh%NB_Ver
        Msh%Nb_Edg = Msh%NB_Ver + Msh%Nb_Cell +22 ; NbrObj = 0;

        print *, "Check allocation"

        IF(ALLOCATED(Msh%Edg)) DEALLOCATE(Msh%Edg)
        IF(ALLOCATED(Msh%EdgNeigh)) DEALLOCATE(Msh%EdgNeigh)
        IF(ALLOCATED(Msh%EIsOnBdr)) DEALLOCATE(Msh%EIsOnBdr)

        ALLOCATE(LstObj(Msh%Nb_Edg+1,5), Msh%Edg(Msh%Nb_Edg,2), Msh%EdgNeigh(Msh%Nb_Edg,2))
        ALLOCATE(Msh%EIsOnBdr(Msh%Nb_Edg), Msh%VIsOnBdr(Msh%NB_Ver), Msh%CIsOnBdr(Msh%Nb_Cell))
        ALLOCATE(Head(SizHead+1))
        Head = 0;   Msh%EdgNeigh = 0; LstObj = 0; Msh%EIsOnBdr = .FALSE.

        print *, "Calculating neighbors"
        MaxSom = maxval(Msh%NbSom(:))
        print *, "max number of Vertex on a cell :",MaxSom

        ! IF(ALLOCATED(Cell2Edg)) DEALLOCATE(Cell2Edg)
        ! ALLOCATE(Cell2Edg(maxval(Msh%NbSom(:)),2))
        Cell2Edg = 0
        DO iPol = 1, Msh%Nb_Cell
            print *, "Processing Cell:", iPol, "of", Msh%Nb_Cell, Msh%NbSom(iPol)
            ! print *,[(i,i+1,i=1,Msh%NbSom(iPol)-1),Msh%NbSom(iPol),1,(i*0,i*0, i=Msh%NbSom(iPol),MaxSom)]
            ! Cell2Edg =  reshape([(i,i+1,i=1,Msh%NbSom(iPol)-1),Msh%NbSom(iPol),1,(i*0,i*0, i=Msh%NbSom(iPol),MaxSom)], [MaxSom,2])
            Cell2Edg(:,1) = [(i,i=1,Msh%NbSom(iPol)-1),  Msh%NbSom(iPol),(0, i=Msh%NbSom(iPol),MaxSom) ]
            Cell2Edg(:,2) = [(i+1,i=1,Msh%NbSom(iPol)-1),      1        ,(0, i=Msh%NbSom(iPol),MaxSom) ]
            ! print *, CELL2EDG
            ! print *,Cell2Edg(:,1)
            ! print *,Cell2Edg(:,2)
            DO iEdg = 0, Msh%NbSom(iPol)-1
                iVer1 = Msh%Sommets(iPol, Cell2Edg(iEdg,1))
                iVer2 = Msh%Sommets(iPol, Cell2Edg(iEdg,2))
                ! print *, Cell2Edg(iEdg,:)
                j = hash_find(iVer1,iVer2)
                ! print *, "Hash find result for edge (", iVer1, ",", iVer2, "):", j
                IF(j /= 0) THEN
                    CALL hash_add(iVer1,iVer2,iPol,j)

                    jPol = LstObj(j,3)
                    Msh%CNeigh(iPol,iEdg+1) = jPol
                    DO jEdg = 0, Msh%NbSom(iPol)-1
                        jVer1 = Msh%Sommets(jPol, Cell2Edg(jEdg,1))
                        jVer2 = Msh%Sommets(jPol, Cell2Edg(jEdg,2))
                        IF((jVer1 == iVer1 .AND. jVer2 == iVer2) .OR. (jVer1 == iVer2 .AND. jVer2 == iVer1)) THEN
                            Msh%CNeigh(jPol,jEdg+1) = iPol
                        END IF
                    END DO
                ELSE
                    CALL hash_add(iVer1,iVer2,iPol,0)
                END IF
            END DO
        END DO

        Msh%max_neigh_Ver = maxval(Msh%NbNeigh)  ! Each edge contributes to two vertices, so divide by 2

        print *, "Maximum number of neighbors for a vertex:", Msh%max_neigh_Ver
        Msh%Nb_EBdr = 0
        DO iEdg = 1, Msh%Nb_Edg-1
            ! print *, "Edge", iEdg, "vertex", Msh%Edg(iEdg,1), "to", Msh%Edg(iEdg,2), "neighbors:", Msh%EdgNeigh(iEdg,1), Msh%EdgNeigh(iEdg,2)
            IF(Msh%EdgNeigh(iEdg,2) == 0) THEN
                Msh%EIsOnBdr(iEdg) = .TRUE.
                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
            ELSE
                Msh%EIsOnBdr(iEdg) = .FALSE.
            END IF
        END DO

        print *,Msh%Nb_Edg, "edges and", Msh%Nb_EBdr, "boundary edges found."

        IF(ALLOCATED(Msh%VNeigh)) DEALLOCATE(Msh%VNeigh)
        IF(ALLOCATED(Msh%Bdr)) DEALLOCATE(Msh%Bdr)
        IF(ALLOCATED(Msh%Ver2Cell)) DEALLOCATE(Msh%Ver2Cell)


        ALLOCATE(Msh%VNeigh(Msh%Nb_Ver, Msh%max_neigh_Ver), Msh%Bdr(Msh%Nb_EBdr))
        ALLOCATE(Msh%Ver2Cell(Msh%Nb_Ver, Msh%max_neigh_Ver))

        Msh%NbNeigh = 0

        DO iPol = 1, Msh%Nb_Cell
            ! print *,Msh%Sommets(iPol,:)
            Msh%NbNeigh(Msh%Sommets(iPol,1:Msh%NbSom(iPol))) = Msh%NbNeigh(Msh%Sommets(iPol,1:Msh%NbSom(iPol))) + 1
            Msh%Ver2Cell(Msh%Sommets(iPol,1:Msh%NbSom(iPol)),Msh%NbNeigh(Msh%Sommets(iPol,1:Msh%NbSom(iPol)))) = iPol
        END DO

        Msh%NbNeigh = 0; Msh%Nb_EBdr = 0

        bcl_edg:DO iEdg = 1, Msh%Nb_Edg-1
            iVer1 = Msh%Edg(iEdg,1)
            iVer2 = Msh%Edg(iEdg,2)
            IF(iVer1 == 0 .AND. iVer2 ==0) Exit bcl_edg
            ! print *,iEdg,iVer1,iVer2
            Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
            Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1
            Msh%VNeigh(iVer1,Msh%NbNeigh(iVer1)) = iVer2
            Msh%VNeigh(iVer2,Msh%NbNeigh(iVer2)) = iVer1
            IF(Msh%EIsOnBdr(iEdg)) THEN
                iPol  = Msh%EdgNeigh(iEdg,1)
                Msh%CIsOnBdr(iPol)  = .TRUE.
                Msh%VIsOnBdr(iVer1) = .TRUE.
                Msh%VIsOnBdr(iVer2) = .TRUE.
                Msh%Nb_EBdr = Msh%Nb_EBdr + 1
                Msh%Bdr(Msh%Nb_EBdr) = iEdg
            END IF
        END DO bcl_edg

        Msh%Nb_EBdr = COUNT(Msh%EIsOnBdr)
        Msh%Nb_CBdr = COUNT(Msh%CIsOnBdr)
        Msh%Nb_Vbdr = COUNT(Msh%VIsOnBdr)


        IF(ALLOCATED(Msh%VBdr)) DEALLOCATE(Msh%VBdr)
        ALLOCATE(Msh%VBdr(Msh%Nb_Vbdr))

        iVer2 = 0
        DO iVer1 = 1,Msh%Nb_Ver
            IF(Msh%VIsOnBdr(iVer1)) THEN 
                iVer2 = iVer2 + 1
                Msh%Vbdr(iVer2) = iVer1
            END IF
        END DO

        print *, "Neighbors calculated successfully."

        Contains
        FUNCTION hash_find(iVer1, iVer2) RESULT(j_hsh)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: iVer1, iVer2
            INTEGER :: j_hsh, key

            key = iVer1 + iVer2
            ! print *,key
            j_hsh = Head(key)

            DO WHILE (j_hsh /= 0)
                ! print *,"looping through hash list, j_hsh:", j_hsh
                IF ((LstObj(j_hsh,1) == iVer1 .AND. LstObj(j_hsh,2) == iVer2) .OR. &
                    (LstObj(j_hsh,1) == iVer2 .AND. LstObj(j_hsh,2) == iVer1)) THEN
                    RETURN
                END IF
                j_hsh = LstObj(j_hsh,5)
            END DO

            j_hsh = 0

        END FUNCTION hash_find

        SUBROUTINE hash_add(iVer1, iVer2, iPol, i_hsh)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: iVer1, iVer2, iPol, i_hsh
            INTEGER :: key


            ! i_hsh is a initial guess. it should be 0 by default, 
            ! but if you have already run hash_find, you can bypass the find to imput the element at i_hsh

            ! if(i_hsh==0)  i_hsh = hash_find(hsh,iVer1,iVer2); ! check if the element isn't in the hash_list already

            IF (i_hsh == 0) THEN

                IF(NbrObj >= Msh%Nb_Edg) THEN
                    print *, "Error: Exceeded maximum number of edges."
                    STOP
                END IF
                key = iVer1 + iVer2
                NbrObj = NbrObj + 1
                LstObj(NbrObj,1) = iVer1
                LstObj(NbrObj,2) = iVer2
                LstObj(NbrObj,3) = iPol
                LstObj(NbrObj,5) = Head(key)
                Head(key) = NbrObj

                Msh%Edg(NbrObj,1) = iVer1    
                Msh%Edg(NbrObj,2) = iVer2
                Msh%NbNeigh(iVer1) = Msh%NbNeigh(iVer1) + 1
                Msh%NbNeigh(iVer2) = Msh%NbNeigh(iVer2) + 1
                Msh%EdgNeigh(NbrObj,1) = iPol

            ELSE
                LstObj(i_hsh,4) = iPol
                Msh%EdgNeigh(i_hsh,2) = iPol
            END IF

        END SUBROUTINE hash_add


    END SUBROUTINE msh_neighbors

    SUBROUTINE Calc_center(Msh)
        IMPLICIT NONE
        TYPE(MESH), INTENT(INOUT) :: Msh
        INTEGER :: iTri
        REAL(prec), DIMENSION(2) :: A, B, C
        
        IF(ALLOCATED(Msh%Center)) DEALLOCATE(Msh%Center)
        ALLOCATE(Msh%Center(Msh%Nb_Cell,2))

        ! Barycenter
        DO iTri = 1, Msh%Nb_Cell
            Msh%Center(iTri,1) = (Msh%Coord(Msh%Sommets(iTri,1),1) + Msh%Coord(Msh%Sommets(iTri,2),1) + Msh%Coord(Msh%Sommets(iTri,3),1)) / 3.0_prec
            Msh%Center(iTri,2) = (Msh%Coord(Msh%Sommets(iTri,1),2) + Msh%Coord(Msh%Sommets(iTri,2),2) + Msh%Coord(Msh%Sommets(iTri,3),2)) / 3.0_prec
        END DO

        ! Circumscribed circle Msh%Center
        DO iTri = 1, Msh%Nb_Cell
            A = Msh%Coord(Msh%Sommets(iTri,1),:)
            B = Msh%Coord(Msh%Sommets(iTri,2),:)
            C = Msh%Coord(Msh%Sommets(iTri,3),:)
            Msh%Center(iTri,:) = Circumcenter(A, B, C)
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

    SUBROUTINE Voronoi(Msh_P,Msh)
        IMPLICIT NONE
        Type(MESH), INTENT(IN)  :: Msh_P
        TYPE(MESH), INTENT(OUT) :: Msh

        INTEGER :: i
        INTEGER :: iPol,iVer,iEdg
        INTEGER :: dPol,dVer,dEdg
        INTEGER :: size

        Msh%Nb_Cell = Msh_P%Nb_Ver
        Msh%NB_Ver  = Msh_P%Nb_Cell + Msh_P%Nb_EBdr + Msh_P%Nb_VBdr
        Msh%Nb_Edg = Msh%NB_Ver + Msh%Nb_Cell

        Msh%maxneigh        = Msh_P%max_neigh_Ver
        Msh%max_neigh_Ver   = Msh_P%maxneigh   

        ALLOCATE(Msh%Cref(Msh%Nb_Cell), Msh%NbSom(Msh%Nb_Cell))
        ALLOCATE(Msh%Sommets(Msh%Nb_Cell,Msh%max_neigh_Ver + 2 ), Msh%CNeigh(Msh%Nb_Cell,Msh%maxneigh))
        ALLOCATE(Msh%Center(Msh%Nb_Cell,2), Msh%Coord(Msh%NB_Ver,2))
        ALLOCATE(Msh%Vref(Msh%NB_Ver), Msh%NbNeigh(Msh%NB_Ver), Msh%VNeigh(Msh%NB_Ver, Msh%max_neigh_Ver), Msh%Ver2Cell(Msh%NB_Ver, Msh%max_neigh_Ver))
        ALLOCATE(Msh%Edg(Msh%Nb_Edg,2),Msh%EdgNeigh(Msh%Nb_Edg,2))

        ! print *,Msh_P%max_neigh_Ver , Msh_P%maxneigh
        ! print *,"nb Edgbdr",Msh_P%Nb_EBdr ,"Nb Vbdr", Msh_P%Nb_VBdr
        ! print *,Msh_P%VNeigh (1,:)
        ! print *,Msh_P%CNeigh(1,:)
        
        Msh%Vref   (1:Msh_P%Nb_Cell  )  = Msh_P%cref;           Msh%Cref   = Msh_P%Vref
        Msh%VNeigh (1:Msh_P%Nb_Cell,:)  = Msh_P%CNeigh;         Msh%CNeigh = Msh_P%VNeigh
        Msh%NbNeigh(1:Msh_P%Nb_Cell  )  = Msh_P%NbSom;          
        Msh%NbSom  = Msh_P%NbNeigh
        Msh%Coord  (1:Msh_P%Nb_Cell,:)  = Msh_P%Center;         Msh%Center = Msh_P%Coord       
        Msh%Ver2Cell = 0   ;    Msh%Sommets =0
        Msh%Ver2Cell(1:Msh_P%Nb_Cell,:) = Msh_P%Sommets;        
        Msh%Sommets= Msh_P%Ver2Cell

        DO i=1,Msh_P%Nb_Ebdr
            iEdg = Msh_P%Bdr(i)

            Msh%Vref    (Msh_P%Nb_Cell+i)   = 0
            ! Msh%VNeigh  (Msh_P%Nb_Cell+i,:) = 0
            ! Msh%NbNeigh (Msh_P%Nb_Cell+i)   = 2
            Msh%Coord   (Msh_P%Nb_Cell+i,1) = Sum(Msh_P%Coord(Msh_P%Edg(iEdg,1:2),1))/2._prec
            Msh%Coord   (Msh_P%Nb_Cell+i,2) = Sum(Msh_P%Coord(Msh_P%Edg(iEdg,1:2),2))/2._prec

            ! write(*,'(2f8.4)') Msh_P%Coord(Msh_P%Edg(iEdg,1:2),:)
            iPol = Msh_P%Edg(iEdg,1)
            Msh%Sommets(iPol,Msh%NbSom(iPol))   = Msh_P%Nb_Cell+i
            Msh%NbSom(iPol) = Msh%NbSom(iPol) +1 

            iPol = Msh_P%Edg(iEdg,2)
            Msh%Sommets(iPol,Msh%NbSom(iPol))   = Msh_P%Nb_Cell+i
            Msh%NbSom(iPol) = Msh%NbSom(iPol) +1 
            
            ! print *,iEdg, Msh_P%Edg(iEdg,:), Msh_P%Nb_Cell+i
            ! print *,"Edge", iPol, "Som", Msh%Sommets(iPol,:)
            ! print *,"Edge", iPol, "Som", Msh%Sommets(iPol,:)
        END DO

        DO i=1,Msh_P%Nb_Vbdr
            iVer = Msh_P%VBdr(i)
            size = Msh_P%Nb_Cell +Msh_P%Nb_Ebdr 
            Msh%Vref    (size+i)   = 0
            ! Msh%VNeigh  (size+i,:) = 0
            ! Msh%NbNeigh (size+i)   = 2
            Msh%Coord   (size+i,:) = Msh_P%Coord(iVer,:)

            iPol = iVer;
            ! print *,"Ver", iPol, "Som", Msh%Sommets(iPol,:)
            Msh%Sommets(iPol,Msh%NbSom(iPol))   = size+i
            ! Msh%NbSom(iPol) = Msh%NbSom(iPol) +1 
            print *,"origine : ",iVer,"Ver", iPol, "Som", Msh%Sommets(iPol,:)
        END DO

        ! ADD COMMENT 
        ! Msh%NbSom(:) = Msh%NbSom(:) -1



        print *,"Cleaning the Mesh"
        CALL Correct_Mesh(Msh)

        print *, "adding the neighbors"
        CALL msh_neighbors(Msh)
        

        ! rajouter les points au bords du domaine.



    END SUBROUTINE

    SUBROUTINE OUTMESH(Msh)
        IMPLICIT NONE

        TYPE(MESH) :: Msh
        INTEGER ::iVer,iEdg,iCell

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
        DO iCell=1,Msh%Nb_Cell
            print *,'--------------'
            print *,iCell,'on Bdr:',Msh%CIsOnBdr(iCell),'Neighbors :',Msh%CNeigh(iCell,:)
            print *,"made of :",Msh%Sommets(iCell,:)
            write(*,'("at coord", 2f8.4 )')Msh%Center(iCell,:)
        END DO

    END SUBROUTINE

    SUBROUTINE Correct_Mesh(Msh)
        IMPLICIT NONE

        TYPE(MESH) :: Msh
        INTEGER :: iVer, jVer,iCell
        INTEGER :: to_Del, to_Rec
        INTEGER, DIMENSION(:), ALLOCATABLE :: temp_Neigh, temp_Cell
        INTEGER :: i,j,k,p, size
        REAL(prec), DIMENSION(2) :: c1,c2

        print *,"CLEANING MESH"

        ! Check if two vertexes are at the same coord (epsilon close)
        ! IF it is relegate the data from the higher indexed one to the lower indexed one
        DO iCell=1,Msh%Nb_Cell
            DO i = 1,Msh%NbSom(iCell)
                iVer = Msh%Sommets(iCell,i)
                c1 = Msh%Coord(iVer,:)

                DO j= i+1, Msh%NbSom(iCell)

                    jVer = Msh%Sommets(iCell,j)
                    ! print *, iVer,jVer, iCell
                    c2 = Msh%Coord(jVer,:)

                    ! IF(maxval(abs(c1-c2)) .LT. eps0) CALL Replace_Vertex(iVer,jVer)
                END DO
            END DO
        END DO

        DO iCell=1,Msh%Nb_Cell
            CALL Convex_Check(iCell)
        END DO

        Contains

        SUBROUTINE Replace_Vertex(iVer,jVer)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: iVer,jVer
            INTEGER :: to_Del,to_Rec,nCell
            INTEGER :: i,j, loc_Del

            to_Del = max(iVer,jVer); to_Rec = min(iVer,jVer)
            print *,'-------- DEL --------(',to_Del,'-',to_Rec,')'
            print *,c1
            print *,c2

            ! what to do with Ref ?
            ! Msh%VIsOnBdr(to_Rec) = Msh%VIsOnBdr(to_Rec) .OR. Msh%VIsOnBdr(to_Del)  !updates the bdr logic

            ! k = FINDLOC(Msh%VBdr,to_Del,DIM=1)
            ! p = FINDLOC(Msh%VBdr,to_Rec,DIM=1)
            ! Deals with Vbdr array
            ! IF(k .NE. 0 .AND. p)  THEN
            !     MSH%VBdr
            ! END IF

            ! update the neighbors vertex without duplicates
            ! CALL unique_sort(size,2*Msh%max_neigh_Ver,temp_Neigh,[Msh%VNeigh(to_Rec,:), Msh%VNeigh(to_Del,:)])
            ! Msh%VNeigh(to_Rec,:) = 0
            ! Msh%VNeigh(to_Rec,1:size) = temp_Neigh

            ! update the neighbors cells without duplicates
            ! CALL unique_sort(size,2*Msh%max_neigh_Ver,temp_Neigh,[Msh%Ver2Cell(to_Rec,:), Msh%Ver2Cell(to_Del,:)])
            ! Msh%Ver2Cell(to_Rec,:) = 0
            ! Msh%Ver2Cell(to_Rec,1:size) = temp_Neigh

            ! update the cell's vertex
            DO i=1,Msh%NbNeigh(to_Rec)
                nCell = Msh%Ver2Cell(to_Rec,i)

                ! print *,"updating cell", nCell
                loc_Del = FINDLOC(Msh%Sommets(nCell,:),to_Rec,DIM=1)
                Msh%Sommets(nCell,loc_Del) = Msh%Sommets(nCell,Msh%NbSom(nCell))
                ! Msh%Sommets(nCell,Msh%NbSom(nCell)) = 0
                Msh%NbSom(nCell) = Msh%NbSom(nCell) -1
            END DO

            ! Msh%NbSom(iCell) = Msh%NbSom(iCell) -1

            print *, temp_cell

            ! Msh%Nb_Ver = Msh%Nb_Ver -1
        END SUBROUTINE Replace_Vertex

        SUBROUTINE Convex_Check(nCell)
            ! Regarde si la liste de Sommet le passage d'un vertex à un autre 
            ! se déroule dans le sens horaire.
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nCell
            REAL(prec),DIMENSION(2) :: Coord_i, Coord_j, Coord_k
            REAL(prec),DIMENSION(2) :: i_To_j, j_To_k
            REAL(prec) :: signed
            LOGICAL :: IsConvex
            INTEGER, DIMENSION(Msh%NbSom(nCell),3) :: Cell2Edg
            INTEGER :: iVer,jVer,kVer
            INTEGER :: i,j,k

            print *,Msh%NbSom(nCell)
            print *,Msh%Sommets(nCell,:)
            Cell2Edg(:,1) = [(i,i=1,Msh%NbSom(nCell)-1),  Msh%NbSom(nCell)]
            Cell2Edg(:,2) = [(i+1,i=1,Msh%NbSom(nCell)-1),      1         ]
            Cell2Edg(:,3) = [(i+2,i=1,Msh%NbSom(nCell)-2),      1,2         ]

            iVer = Msh%Sommets(nCell,Cell2Edg(1,1))
            jVer = Msh%Sommets(nCell,Cell2Edg(1,2))
            kVer = Msh%Sommets(nCell,Cell2Edg(1,3))
            Coord_i = Msh%Coord(iVer,:)
            Coord_j = Msh%Coord(jVer,:)
            Coord_k = Msh%Coord(kVer,:)
            i_To_j = Coord_j-Coord_i
            j_To_k = Coord_k-Coord_j
            signed = Sign(1._prec,det(i_To_j,j_To_k))

            IsConvex = .FALSE.

            i =1; j=1
            DO WHILE(.NOT. ISConvex)
            
            iVer = Msh%Sommets(nCell,Cell2Edg(i,1))
            jVer = Msh%Sommets(nCell,Cell2Edg(i,2))
            kVer = Msh%Sommets(nCell,Cell2Edg(i,3))
            Coord_i = Msh%Coord(iVer,:)
            Coord_j = Msh%Coord(jVer,:)
            Coord_k = Msh%Coord(kVer,:)
            i_To_j = Coord_j-Coord_i
            j_To_k = Coord_k-Coord_j

            ! write(*,fmt='(i3,1x,i3,1x,"in",i3,1x," de det ",f6.3,1x,"avec vect :",2f6.3,1x,2f6.3)')iVer,jVer,nCell,det(i_To_j,j_To_k), i_To_j,j_To_k
            write(*,fmt='("Vertex :",i3,1x,2f6.3)') iVer, Coord_i
            write(*,fmt='("Vertex :",i3,1x,2f6.3)') jVer, Coord_j
            write(*,fmt='("Vertex :",i3,1x,2f6.3)') kVer, Coord_k
            write(*,fmt='(i3,1x,i3,1x,i3,1x,"in",i3,1x," de det ",f6.3,1x,"avec vect :",2f6.3,1x,2f6.3," sign Cell :",f6.3)')iVer,jVer,kVer,nCell,det(i_To_j,j_To_k), i_To_j,j_To_k,signed
                
            IF(det(i_To_j,j_To_k) < -eps0) THEN
                CALL Permute_Vertex_onCell(jVer,kVer,nCell)
                iVer = Msh%Sommets(nCell,Cell2Edg(1,1))
                jVer = Msh%Sommets(nCell,Cell2Edg(1,2))
                kVer = Msh%Sommets(nCell,Cell2Edg(1,3))
                Coord_i = Msh%Coord(iVer,:)
                Coord_j = Msh%Coord(jVer,:)
                Coord_k = Msh%Coord(kVer,:)
                i_To_j = Coord_j-Coord_i
                j_To_k = Coord_k-Coord_j
                signed = Sign(1._prec,det(i_To_j,j_To_k))
                j = 0; i=0                    
            END IF

            IF(j>Msh%NbSom(nCell)+1) ISConvex = .TRUE.

            j = j +1;
            i = i +1
            IF(i>Msh%NbSom(nCell)) i=1

            END DO


        END SUBROUTINE Convex_Check

        SUBROUTINE Gift_Wrapping(nCell)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nCell
            REAL(prec),DIMENSION(:,:), ALLOCATABLE :: Coord_Cell
            REAL(prec),DIMENSION(2) :: Coord_i, Coord_j, Coord_start, vect
            REAL(prec) :: Min_angle
            INTEGER, DIMENSION(Msh%NbSom(nCell)) :: Convex_hull
            INTEGER :: Ver_i,Ver_j, Start_Point,Next_Point
            INTEGER :: i,j

            print *, "Gift wrapping ",nCell
            print *, Msh%Sommets(nCell,:)

            ALLOCATE(Coord_Cell(Msh%NbSom(nCell),2))
            Coord_Cell = Msh%Coord(Msh%Sommets(nCell,1:Msh%NbSom(nCell)),:)
            

            print *,Coord_Cell(:,1)

            ! Find point leftmost of nCell
            Coord_start(1) = Minval(Coord_Cell(:,1))
            Start_Point = FINDLOC(Coord_Cell(:,1),Coord_start(1),DIM=1)
            Start_Point = Msh%Sommets(nCell,Start_Point)

            Ver_i= Start_Point;
            Convex_hull(1) = Start_Point
            i = 2
            DO WHILE(Next_Point .NE. Start_Point)
                print *,Ver_i 
                Coord_i = Msh%Coord(Ver_i,:)
                Min_angle = 2*pi

                DO j=1,Msh%NbSom(nCell)
                    Ver_j = Msh%Sommets(nCell,j)
                    Coord_j = Msh%Coord(Ver_j,:)
                    vect = Coord_j-Coord_i
                    ! print *,'vect',vect
                    print *,'angle', pi/2._prec -atan2(vect(2),vect(1)), Min_Angle
                    IF(pi/2._prec - atan2(vect(2),vect(1)) < Min_angle .AND. Ver_j .NE. Ver_i) THEN 
                        Next_Point = Ver_j
                        Min_angle = pi/2._prec - atan2(vect(2),vect(1))
                    END IF
                END DO
                print *,"To", Next_Point
                Ver_i = Next_Point
                Convex_hull(i)= Ver_i
                i = i+1
            END DO
            
            Msh%Sommets(nCell,1:Msh%NbSom(nCell)) = Convex_hull

        END SUBROUTINE Gift_Wrapping

        SUBROUTINE Permute_Vertex_onCell(iVer,jVer,nCell)
            IMPLICIT NONE
            INTEGER,INTENT(IN) :: iVer,jVer,nCell
            INTEGER :: loc_i,loc_j

            print *,"permute ",iVer,jVer,'in cell',nCell

            loc_i = FINDLOC(Msh%Sommets(nCell,:),iVer,DIM=1)
            loc_j = FINDLOC(Msh%Sommets(nCell,:),jVer,DIM=1)
            Msh%Sommets(nCell,loc_i) = jVer
            MSh%Sommets(nCell,loc_j) = iVer
            print *,Msh%Sommets(nCell,:)

        END SUBROUTINE

    END SUBROUTINE Correct_Mesh




    
END MODULE mod_mesh
MODULE mod_mesh

    use mod_var
    use mod_divers

    IMPLICIT NONE
    
    INTEGER, DIMENSION(0:2,2), PARAMETER :: tri2edg = reshape([1,2, 2,3, 3,1], [3,2])

    Contains


    SUBROUTINE ReadMesh(namefile)
        ! ----------------------------------------------------------------!
        !  READ MESH FROM .MESH FILE
        !  maxneigh = 3 (triangles)
        !  besoin de fixé maxneigh avant l'appel de cette fonction
        ! ----------------------------------------------------------------!


        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        INTEGER :: io
        INTEGER :: i

        print *, "Reading mesh file:", TRIM(namefile)

        OPEN(unit=10, file=TRIM(namefile), status='old', action='read', iostat=io)
        IF (io /= 0) THEN
            print *, "Error opening file. IOSTAT =", io
            STOP
        END IF
        
        IF (ALLOCATED(Vref0)) DEALLOCATE(Vref0)
        IF (ALLOCATED(Coord0)) DEALLOCATE(Coord0)
        IF (ALLOCATED(cref0)) DEALLOCATE(cref0)
        IF (ALLOCATED(Sommets0)) DEALLOCATE(Sommets0)
        IF (ALLOCATED(CNeigh0)) DEALLOCATE(CNeigh0)

        CALL JUMP(3,10)
        READ(10,*) nb_Ver0
        print *, "Number of vertices:", nb_Ver0
        ALLOCATE(Vref0(nb_Ver0), Coord0(nb_Ver0, 2), NbNeigh0(nb_Ver0))
        READ(10,*) (Coord0(i,1:2),Vref0(i),i=1,nb_Ver0)  

        CALL JUMP(1,10)
        READ(10,*) nb_cell0
        print *, "Number of cells:", nb_cell0

        ALLOCATE(cref0(nb_cell0), Sommets0(nb_cell0, maxneigh), CNeigh0(nb_cell0, maxneigh))
        ALLOCATE(Ver2Cell0(nb_Ver0, maxneigh), NbSom0(nb_cell0))

        Ver2Cell0 = 0
        DO i = 1, nb_cell0
            READ(10,*) Sommets0(i,1:3),cref0(i)
            Ver2Cell0(Sommets0(i,1:3),:) = i
        END DO

        ! CALCUL DES VOISINS
        print *, "Calculating neighbors..."
        CALL msh_neighbors

        ! CALCUL DES CENTRES DES MAILLES
        print *, "Calculating cell centers..."
        CALL Calc_center

        Nbsom0 = 3

    END SUBROUTINE ReadMesh

    SUBROUTINE Readvtk(namefile)
        ! ----------------------------------------------------------------!
        !  READ MESH FROM .VTK FILE en format ASCII, POLYDATA
        !  maxneigh = 3 (triangles)
        !  besoin de fixé maxneigh avant l'appel de cette fonction
        ! ----------------------------------------------------------------!
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
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
        
        IF (ALLOCATED(Vref0)) DEALLOCATE(Vref0)
        IF (ALLOCATED(Coord0)) DEALLOCATE(Coord0)
        IF (ALLOCATED(cref0)) DEALLOCATE(cref0)
        IF (ALLOCATED(Sommets0)) DEALLOCATE(Sommets0)
        IF (ALLOCATED(CNeigh0)) DEALLOCATE(CNeigh0)
        IF (ALLOCATED(NbNeigh0)) DEALLOCATE(NbNeigh0)
        IF (ALLOCATED(Ver2Cell0)) DEALLOCATE(Ver2Cell0)
        IF (ALLOCATED(NbSom0)) DEALLOCATE(NbSom0)

        CALL JUMP(4,10)  
        
        READ(10,*) blank, Nb_Ver0, blank
        print *,Nb_Ver0
        ALLOCATE(Vref0(Nb_Ver0), Coord0(Nb_Ver0, 2), NbNeigh0(Nb_Ver0))
        DO i = 1, Nb_Ver0
            READ(10,*) Coord0(i,1:2)
        END DO

        ! CALL JUMP(1,10)
        READ(10,*) blank, i_temp, j_temp
        print *, i_temp, j_temp
        ALLOCATE(POLY(i_temp, maxneigh), isbdr(i_temp))
        ALLOCATE(Som_(i_temp,maxneigh))
        ALLOCATE(Ver2Cell0(Nb_Ver0, maxneigh), NbSom_(i_temp))

        Ver2Cell0 = 0
        Nb_bdr0 =0
        DO i = 1, i_temp
            READ(10,*) NbSom_(i), Som_(i,1:NbSom_(i))
            Som_(i,1:NbSom_(i)) = Som_(i,1:NbSom_(i)) +1

            IF(NbSom_(i) == 2) THEN
                Nb_bdr0 = Nb_bdr0 + 1
                isbdr(i) = .TRUE.
            ELSE
                isbdr(i) = .FALSE.
            END IF
        END DO

        Nb_Cell0 = i_temp - Nb_bdr0
        ALLOCATE(cref0(Nb_Cell0), NbSom0(Nb_Cell0))
        ALLOCATE(Sommets0(Nb_Cell0,maxneigh), Cneigh0(Nb_cell0,maxneigh))
        j = 1
        DO i = 1,i_temp; IF(.not. isbdr(i)) THEN
            NbSom0(j) = NbSom_(i)
            Sommets0(j,:) = Som_(i,:)
            Ver2Cell0(Som_(i,:),:)=j
            j = j + 1
        END IF; END DO

        ! CALCUL DES VOISINS
        print *, "Calculating neighbors..."
        CALL msh_neighbors

        ! CALCUL DES CENTRES DES MAILLES
        print *, "Calculating cell centers..."
        CALL Calc_center        

    END SUBROUTINE Readvtk

    SUBROUTINE msh_neighbors
        IMPLICIT NONE
        INTEGER :: iTri, jTri, iEdg, jEdg
        INTEGER :: iVer1, iVer2, jVer1, jVer2
        INTEGER :: j_hsh


        ! HASH TABLE
        INTEGER :: SizHead                                          ! Maximum key value, the key is in [0,SizHead-1]
        INTEGER :: NbrObj                                           ! Number of objects in the hash table
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  Head                 ! Head[key%(SizHead)] = link to the first object having this key in the LstObj list
        INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  LstObj       ! List of objects in the hash table


        SizHead = 2*NB_Ver0
        Nb_Edg0 = NB_Ver0 + NB_Cell0; NbrObj = 1;

        print *, "Check allocation"

        IF(ALLOCATED(Edg0)) DEALLOCATE(Edg0)
        IF(ALLOCATED(EdgNeigh0)) DEALLOCATE(EdgNeigh0)
        IF(ALLOCATED(is_bdr0)) DEALLOCATE(is_bdr0)

        ALLOCATE(LstObj(Nb_Edg0+1,5), Edg0(Nb_Edg0,2), EdgNeigh0(Nb_Edg0,2),is_bdr0(Nb_Edg0))
        ALLOCATE(Head(SizHead+1))
        Head = 0;   EdgNeigh0 = 0; LstObj = 0; is_bdr0 = .FALSE.

        print *, "Calculating neighbors"

        DO iTri = 1, NB_Cell0
            ! print *, "Processing triangle:", iTri, "of", NB_Cell0
            DO iEdg = 0, 2
                iVer1 = Sommets0(iTri, tri2edg(iEdg,1))
                iVer2 = Sommets0(iTri, tri2edg(iEdg,2))

                Edg0(NbrObj,1) = iVer1    
                Edg0(NbrObj,2) = iVer2
                NbNeigh0(iVer1) = NbNeigh0(iVer1) + 1
                NbNeigh0(iVer2) = NbNeigh0(iVer2) + 1
                EdgNeigh0(NbrObj,1) = iTri

                j_hsh = hash_find(iVer1,iVer2)
                ! print *, "Hash find result for edge (", iVer1, ",", iVer2, "):", j_hsh
                IF(j_hsh /= 0) THEN
                    CALL hash_add(iVer1,iVer2,iTri,j_hsh)

                    jTri = LstObj(j_hsh,3)
                    CNeigh0(iTri,iEdg+1) = jTri
                    DO jEdg = 0, 2
                        jVer1 = Sommets0(jTri, tri2edg(jEdg,1))
                        jVer2 = Sommets0(jTri, tri2edg(jEdg,2))
                        IF((jVer1 == iVer1 .AND. jVer2 == iVer2) .OR. (jVer1 == iVer2 .AND. jVer2 == iVer1)) THEN
                            CNeigh0(jTri,jEdg+1) = iTri
                        END IF
                    END DO
                ELSE
                    CALL hash_add(iVer1,iVer2,iTri,0)
                END IF
            END DO
        END DO

        max_neigh_Ver0 = maxval(NbNeigh0/2)  ! Each edge contributes to two vertices, so divide by 2

        print *, "Maximum number of neighbors for a vertex:", max_neigh_Ver0
        Nb_bdr0 = 0
        DO iEdg = 1, Nb_Edg0-1
            ! print *, "Edge", iEdg, "vertex", Edg0(iEdg,1), "to", Edg0(iEdg,2), "neighbors:", EdgNeigh0(iEdg,1), EdgNeigh0(iEdg,2)
            IF(EdgNeigh0(iEdg,2) == 0) THEN
                is_bdr0(iEdg) = .TRUE.
                Nb_bdr0 = Nb_bdr0 + 1
            ELSE
                is_bdr0(iEdg) = .FALSE.
            END IF
        END DO

        print *,Nb_Edg0, "edges and", Nb_bdr0, "boundary edges found."

        IF(ALLOCATED(VNeigh0)) DEALLOCATE(VNeigh0)
        IF(ALLOCATED(Bdr0)) DEALLOCATE(Bdr0)

        ALLOCATE(VNeigh0(Nb_Ver0, max_neigh_Ver0), Bdr0(Nb_bdr0, 2))

        NbNeigh0 = 0; Nb_bdr0 = 0

        DO iEdg = 1, Nb_Edg0-1
            iVer1 = Edg0(iEdg,1)
            iVer2 = Edg0(iEdg,2)
            NbNeigh0(iVer1) = NbNeigh0(iVer1) + 1
            NbNeigh0(iVer2) = NbNeigh0(iVer2) + 1
            VNeigh0(iVer1,NbNeigh0(iVer1)) = iVer2
            VNeigh0(iVer2,NbNeigh0(iVer2)) = iVer1
            IF(is_bdr0(iEdg)) THEN
                is_bdr0(iEdg) = .TRUE.
                Nb_bdr0 = Nb_bdr0 + 1
                Bdr0(Nb_bdr0,1) = iVer1
                Bdr0(Nb_bdr0,2) = iVer2
            END IF
        END DO

        Nb_bdr0 = Nb_bdr0 -1
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

                IF(NbrObj >= Nb_Edg0) THEN
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
            ELSE
                LstObj(i_hsh,4) = iTri
                EdgNeigh0(i_hsh,2) = iTri
            END IF

        END SUBROUTINE hash_add


    END SUBROUTINE msh_neighbors

    SUBROUTINE Calc_center
        IMPLICIT NONE
        INTEGER :: iTri
        REAL(prec), DIMENSION(2) :: A, B, C
        
        IF(ALLOCATED(Center0)) DEALLOCATE(Center0)
        ALLOCATE(Center0(nb_cell0,2))

        ! Barycenter
        DO iTri = 1, nb_cell0
            Center0(iTri,1) = (Coord0(Sommets0(iTri,1),1) + Coord0(Sommets0(iTri,2),1) + Coord0(Sommets0(iTri,3),1)) / 3.0_prec
            Center0(iTri,2) = (Coord0(Sommets0(iTri,1),2) + Coord0(Sommets0(iTri,2),2) + Coord0(Sommets0(iTri,3),2)) / 3.0_prec
        END DO

        ! Circumscribed circle center
        DO iTri = 1, nb_cell0
            A = Coord0(Sommets0(iTri,1),:)
            B = Coord0(Sommets0(iTri,2),:)
            C = Coord0(Sommets0(iTri,3),:)
            Center0(iTri,:) = Circumcenter(A, B, C)
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

    SUBROUTINE WriteMesh(namefile, nb_cell, cref, Sommets, nb_Ver, Vref, Coord)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: namefile
        INTEGER, INTENT(IN) :: nb_cell, nb_Ver
        INTEGER, DIMENSION(:), INTENT(IN) :: cref, Vref
        INTEGER, DIMENSION(:,:), INTENT(IN) :: Sommets
        REAL(prec), DIMENSION(:,:) , INTENT(IN) ::  Coord

        INTEGER :: i
        CHARACTER(LEN=1) :: blank

        blank =""
        OPEN(unit=20, file=TRIM(namefile), status='unknown', action='write')

        WRITE(20,'(A)') "MeshVersionFormatted 2"
        WRITE(20,'(A)') blank
        WRITE(20,'(A)') "Dimension 2"
        WRITE(20,'(A)') blank

        WRITE(20,'(A)') "Vertices"
        WRITE(20,'(I10)') nb_Ver
        DO i = 1, nb_Ver
            WRITE(20,'(2F15.6,I10)') Coord(i,:), Vref(i)
        END DO
        WRITE(20,'(A)') blank

        WRITE(20,'(A)') "Triangles"
        WRITE(20,'(I10)') nb_cell
        DO i = 1, nb_cell
            WRITE(20,'(3I10,I10)') Sommets(i,:), cref(i)
        END DO
        WRITE(20,'(A)') blank
        

    END SUBROUTINE WriteMesh

    SUBROUTINE Write_vtk(fichier, Nb_Ver, xs, ys, Nbcell, Sommets, NbBdr,Bdr)
        !******** Creation du fichier paraview qui represente la grille ************
        !** Sens des variables : 
        !       fichier = nom du fichier dans lequel mettre le maillage
        !       Nb_Ver = nombre de sommets
        !       xs,ys = coordonnées des sommets
        !       Nbcell = nombre de mailles
        !       Sommets = tableau des sommets de chaque maille
        !       NbBdr = nombre d'aretes sur le bord du maillage
        !       Bdr = tableau des aretes sur le bord du maillage
        !** Sortie : fichier paraview en format .vtk (ASCII, POLYDATA)
        !**************************************

        CHARACTER(LEN=*) :: fichier
        Character(LEN=1) :: blank
        INTEGER :: Nb_Ver, Nbcell, NbBdr
        REAL(prec), DIMENSION(:) :: xs, ys
        INTEGER, DIMENSION(:,:) ::  Sommets
        INTEGER, DIMENSION(:,:) :: Bdr

        INTEGER :: i, j

        blank = ""
        OPEN(unit=21, File=fichier, status="unknown")
        WRITE(21,"(A)") "# vtk DataFile Version 2.0"
        WRITE(21,"(A)") "Grille"
        WRITE(21,"(A)") "ASCII"

        WRITE(21,"(A)") "DATASET POLYDATA"
        WRITE(21,"(A,I10,X,A)") "POINTS",Nb_Ver,"float"
        DO i=1, Nb_Ver
        WRITE(21,"(F12.5,F12.5,F12.5)") xs(i),ys(i),0D0
        ENDDO
        write(21,"(A)") blank
        WRITE(21,"(A,I10,X,I10)") "POLYGONS",Nbcell+NbBdr,Nbcell*3+NbBdr*2+Nbcell+NbBdr
        DO i=1, NbBdr
        WRITE(21,*) 2,(Bdr(i,j)-1,j=1,2)
        ENDDO
        DO i=1, Nbcell
        WRITE(21,*) NbSom0(i),(Sommets(i,j)-1,j=1,NbSom0(i))
        ENDDO

        CLOSE(21)

    END SUBROUTINE Write_vtk


END MODULE mod_mesh


MODULE Mod_Delaunay

    use mod_mesh

    IMPLICIT NONE

    INTEGER ::  sizeof_cavity = 0;  ! size of the cavity/index of the last element added to the cavity 
    INTEGER ::  sizeof_boundary =0; ! var global


    INTEGER, DIMENSION(3,2), PARAMETER :: tri2edg = reshape([2,3 , 1,3 , 1,2], [3,2])

CONTAINS


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
                      IsDelaun = .NOT. IsInCircum(nCell,Msh%Coord(tVer,:),Msh)
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

  FUNCTION IsInCircum(nCell,Point,Msh)         ! Check if the vertex tVer is inside the circumcircle of the triangle nCell
      IMPLICIT NONE
      REAL(Prec), DIMENSION(2) :: Point
      TYPE(Mesh), INTENT(IN) :: Msh
      TYPE(Cell), INTENT(IN) :: nCell
      REAL(prec), DIMENSION(2) :: c1,c2,c3
      REAL(prec), DIMENSION(3) :: v1,v2,v3
      LOGICAL :: IsInCircum

      c1 = Msh%Coord(nCell%Sommets(1),:)
      c2 = Msh%Coord(nCell%Sommets(2),:)
      c3 = Msh%Coord(nCell%Sommets(3),:)

      v1(1:2) = c1-Point; v1(3) = (norm2(c1-Point))**2
      v2(1:2) = c2-Point; v2(3) = (norm2(c2-Point))**2
      v3(1:2) = c3-Point; v3(3) = (norm2(c3-Point))**2

      IF(det(v1,v2,v3) >=0) THEN; IsInCircum = .TRUE.
      ELSE ; IsInCircum = .FALSE.
      END IF

  END FUNCTION IsInCircum

  ! ----------------------------------------------------------------

  FUNCTION Localise_Tri(Msh, Point)

    ! ----------------------------------
    ! Rends l'indice du triangle dans lequel se trouve le Point
    ! La fonction n'affecte pas le maillage
    !
    ! Compass_check nous donne si la prochaine direction a été choisi
    !
    ! On assume que le maillage est triangulaire (Msh%max_neigh_Cell = 3)
    ! On assume que le point n'a pas déjà été ajouté (onPoint == False) (à rajouter)
    ! ----------------------------------


    IMPLICIT NONE
    INTEGER, DIMENSION(4)  :: Localise_Tri
    TYPE(Mesh),               INTENT(IN):: Msh
    REAL(PREC), DIMENSION(2), INTENT(IN) :: Point

    TYPE(Cell) :: Tri
    REAL(PREC), DIMENSION(3) :: Coord_bary
    LOGICAL :: inTri,onEdg,onPoint, Compass_check
    INTEGER :: i, iTri, eTri, edge_case
    INTEGER, DIMENSION(0:Msh%Nb_Cell) :: Mrk

    iTri = 1
    Mrk =0; Mrk(0) = 1

    inTri = .False. ; onEdg = .False.

      edge_case = 1
    DO WHILE(.NOT. inTri .AND. .NOT. onEdg)
      Compass_check = .FALSE.
      ! traitement des boucles ou des maillages non connexes

      DO WHILE(Mrk(iTri) == Mrk(0) +1 .OR. iTri == 0)
        ! print *,"skip",iTri
        iTri = edge_case
        edge_case = edge_case +1
      END DO 

      Tri = Msh%Cells(iTri)
      IF(iTri .NE. 0) THEN 
        ! CALL Out_Cell(Tri)
        ! print *,Mrk(Tri%CNeigh)

        Mrk(iTri) = Mrk(0)+1
        Coord_bary = Coord_Bary_FromTri(Tri, Msh, Point)

        IF(all(Coord_bary .LT. -eps0) ) STOP "ERROR POINT OR TRIANGLE WRONGLY DEFINED"


        IF(abs(Coord_bary(1)) .LT. eps0 .AND. (Coord_bary(2) .GT. 0._prec .AND. Coord_bary(3) .GT. 0._prec)) THEN; Compass_check= .true.; onEdg = .true.; eTri = Tri%CNeigh(1); ENDIF
        IF(abs(Coord_bary(2)) .LT. eps0 .AND. (Coord_bary(3) .GT. 0._prec .AND. Coord_bary(1) .GT. 0._prec)) THEN; Compass_check= .true.; onEdg = .true.; eTri = Tri%CNeigh(2); ENDIF
        IF(abs(Coord_bary(3)) .LT. eps0 .AND. (Coord_bary(1) .GT. 0._prec .AND. Coord_bary(2) .GT. 0._prec)) THEN; Compass_check= .true.; onEdg = .true.; eTri = Tri%CNeigh(3); ENDIF

        ! print *, Coord_bary
        ! print *, Coord_bary .LT. -eps0

        IF(Coord_bary(1) .LT. 0._prec .AND. .NOT. Compass_check) THEN; Compass_check= .true.;  iTri= Tri%CNeigh(1); ENDIF; 
        IF(Coord_bary(2) .LT. 0._prec .AND. .NOT. Compass_check) THEN; Compass_check= .true.;  iTri= Tri%CNeigh(2); ENDIF; 
        IF(Coord_bary(3) .LT. 0._prec .AND. .NOT. Compass_check) THEN; Compass_check= .true.;  iTri= Tri%CNeigh(3); ENDIF; 

        IF(all(Coord_bary .GT. eps0) )  THEN; Compass_check = .TRUE.; inTri = .True.;  EXIT; END IF;

        ! if(iTri == 0 .AND. Mrk(iTri) == Mrk(0)+1) THEN; edge_case = edge_case + 1; iTri = edge_case; END IF
        
        ! print *,"Itri",iTri
      END IF
      ! print *,"a ",iTri
    END DO

    ! print *, "found in :", iTri

    Localise_Tri = 0
    IF(OnEdg) inTri = .True.

    IF(inTri) THEN; Localise_Tri(1) = 1; Localise_Tri(2) = iTri; ENDIF
    IF(OnEdg) THEN; Localise_Tri(3) = 1; Localise_Tri(4) = eTri; ENDIF
      
  END FUNCTION Localise_Tri

  SUBROUTINE Digging_a_hole(Msh,Point, cavity) ! diggy diggy hole
    IMPLICIT NONE
    TYPE(Mesh),                       INTENT(IN) :: Msh
    REAL(prec), DIMENSION(2),         INTENT(IN) :: Point
    INTEGER, DIMENSION(Msh%Nb_Cell),  INTENT(OUT) :: cavity  !besoin de limiter la taille de ses listes 
    INTEGER, DIMENSION(Msh%Nb_Cell) :: pile                 !besoin de limiter la taille de ses listes 

    TYPE(Cell) :: iCell
    REAL(prec), DIMENSION(3,2) :: Sommets_Crd
    INTEGER, DIMENSION(3) :: Sommets
    INTEGER, DIMENSION(4) :: Loc_Tri
    INTEGER :: sizeof_pile
    INTEGER :: iTri, jTri, jEdg
    LOGICAL :: IsInCavity, IsInPile
    INTEGER :: k


    Loc_Tri = Localise_Tri(Msh,Point)
    IF(Loc_Tri(1) .NE. 1) STOP "Point Not Found"

    ! initializing the pile
    pile = 0; sizeof_pile = 0 ;sizeof_cavity = 0
    if(Loc_Tri(1) ==1) THEN;  pile(1) = Loc_Tri(2); sizeof_pile = sizeof_pile +1; END IF
    if(Loc_Tri(3) ==1) THEN;  pile(2) = Loc_Tri(4); sizeof_pile = sizeof_pile +1; END IF

    ! print *,Point
    ! print *, Loc_Tri
    DO WHILE(sizeof_pile .GT. 0)
      ! add the element to the cavity and remove it from the pile
      iTri = pile(sizeof_pile);   sizeof_pile   = sizeof_pile -1
      sizeof_cavity = sizeof_cavity +1
      cavity(sizeof_cavity) = iTri; 
      iCell = Msh%Cells(iTri)

      ! check if the neighbours are in the cavity
      DO jEdg = 1,3
        IsInCavity =.False. ; IsInPile = .False.
        jTri = iCell%CNeigh(jEdg)
        Sommets = iCell%Sommets(1:3)
        Sommets_Crd= Msh%Coord(Sommets,:)

        ! Check if the element hasn't been added
        Do k=1,sizeof_cavity; IF(cavity(k) == jTri) THEN; IsInCavity = .True.;  Exit; END IF; END DO;
        Do k=1,sizeof_pile;   IF(pile(k) == jTri)   THEN; IsInPile = .True.;    Exit; END IF; END DO;

        IF(.NOT. IsInCavity .AND. (.NOT. IsInPile .AND. jTri .NE. 0 ))THEN
          ! If it is in the circumcircle
          IF(IsInCircum(Msh%Cells(jTri),Point,Msh) ) THEN 
            Sizeof_pile = Sizeof_pile +1
            pile(Sizeof_pile) = jTri
          END IF
        END IF

      END DO

    ! print *,'size of cavity', Sizeof_Cavity
    ! print *,cavity(1:Sizeof_Cavity)
    END DO


  END SUBROUTINE Digging_a_hole

  SUBROUTINE Deleting_cavity(Msh,Cavity)
    IMPLICIT NONE
    Type(Mesh), INTENT(INOUT) :: Msh
    INTEGER, DIMENSION(Msh%Nb_Cell), INTENT(IN) :: cavity

    TYPE(Cell) :: Tri_Neigh, Tri_Cavity
    INTEGER :: iTri_Neigh, iTri_Cavity
    INTEGER :: i,j,k

    DO i=1,Sizeof_Cavity
      ! Deleting it from the list of neighbours
      iTri_Cavity = Cavity(i);
      ! print *,'Delete', iTri_Cavity
      CALL Cell_Delete(iTri_Cavity, Msh)
    END DO

  END SUBROUTINE Deleting_cavity

  SUBROUTINE boundary_hole(Msh, Cavity, boundary )
    IMPLICIT NONE
    Type(Mesh), INTENT(INOUT) :: Msh
    INTEGER, DIMENSION(Msh%Nb_Cell), INTENT(IN) :: cavity
    INTEGER, DIMENSION(3*Sizeof_Cavity,3) :: boundary

    TYPE(Cell) :: Cell_
    LOGICAL :: IsOnBdr
    INTEGER :: iVer1,iVer2, NbZero
    INTEGER :: i,j,k

    boundary = 1

    DO i=1,Sizeof_Cavity
      ! print *, 'boundary size', sizeof_boundary
      DO j=1,3
        Cell_ = Msh%Cells(cavity(i))
        iVer1 = Cell_%Sommets(Tri2Edg(j,1))
        iVer2 = Cell_%Sommets(Tri2Edg(j,2))
        IsOnBdr = .False.

        ! print *, iVer1,iVer2, Cell_%index

        ! on regarde si elle a déjà été ajouté à la liste frontière
        Do k=1,sizeof_boundary
          ! print *, 'test against ', boundary(k,:)
          IF((iVer1 == boundary(k,1) .AND. iVer2 == boundary(k,2)) .OR. &
          &  (iVer2 == boundary(k,1) .AND. iVer1 == boundary(k,2))) THEN  !  if it is in the boundary -> remove it 

            IsOnBdr = .True. 
            boundary(k,1) = boundary(sizeof_boundary,1)
            boundary(k,2) = boundary(sizeof_boundary,2)
            boundary(k,3) = boundary(sizeof_boundary,3)
            ! print *,"removed", boundary(sizeof_boundary,:)
            sizeof_boundary = sizeof_boundary - 1;
            ! EXIT
          END IF
        END DO

      ! if it is not in the boundary -> adding it
        IF(.NOT. IsOnBdr) THEN 
        sizeof_boundary = sizeof_boundary + 1;
        boundary(sizeof_boundary,1)=iVer1;
        boundary(sizeof_boundary,2)=iVer2;
        boundary(sizeof_boundary,3)=Cell_%CNeigh(j);
        ! print *,'Cell :', Cell_%index,' : ', Cell_%CNeigh(:)
        ! print *,"added", boundary(sizeof_boundary,:)
        END IF
      END DO
    END DO

    ! Do i=1,sizeof_Boundary
      ! print *, i,"e bdry", boundary(i,:)
    ! END DO
    ! print *,boundary(:,3)

    NbZero = 0
    Do i=1,Sizeof_Boundary
      IF(boundary(i,1)==0) THEN 
        boundary(i,:) = boundary(sizeof_boundary-NbZero-1,:);  
        boundary(sizeof_boundary-NbZero-1,:)=0;
        NbZero = NbZero +1;
      END IF
    END DO

    sizeof_boundary = sizeof_boundary - NbZero

  END SUBROUTINE boundary_hole

  SUBROUTINE Starring_Point(Msh,cavity, boundary, Point)
    IMPLICIT NONE

    TYPE(Mesh),                         INTENT(INOUT) :: Msh
    INTEGER, DIMENSION(sizeof_boundary),    INTENT(INOUT) :: cavity
    INTEGER, DIMENSION(3*Sizeof_Cavity,3),INTENT(IN) :: boundary
    REAL(prec), DIMENSION(2),           INTENT(IN) :: Point

    TYPE(Cell) :: Tri,Tri_bdry, Tri2
    INTEGER :: iTri, iTri_add, iTri_bdry, jTri
    INTEGER :: iVer1,iVer2,jVer1,jVer2
    REAL(prec) :: surface
    REAL(prec), DIMENSION(2) :: Ver1,Ver2
    INTEGER :: i,j,k,n


    ! print *,"---------------------------------"
    ! print *,'         starring                '
    ! print *,"---------------------------------"

    ! starring the point to the cavity's boundary
    i = 0
    DO While (i .LT. sizeof_boundary)
      i = i +1
      ! print *, boundary(i,1:2)
      iVer1 = boundary(i,1); iVer2 = boundary(i,2)
      Ver1 = Msh%Coord(iVer1,:); Ver2= Msh%Coord(iVer2,:)


      ! toutes les arrêtes (qui existent) dont le triangle généré n'est pas dégénéré
      ! print *,Ver1,Ver2
      DO while( ((boundary(i,1)==0 .OR. boundary(i,2)==0) .AND.  i<sizeof_boundary) .OR. abs(Signed_Surface(Ver1,Ver2,Point)) <eps0)
        ! print *,"skip",i
        i =i + 1;
        iVer1 = boundary(i,1); 
        iVer2 = boundary(i,2); 
        Ver1 = Msh%Coord(iVer1,:); Ver2= Msh%Coord(iVer2,:)
        
        if(i > sizeof_boundary) Exit;
      END DO
      if(i > sizeof_boundary) Exit;


      if(sizeof_cavity>0) THEN; iTri_add = cavity(sizeof_cavity);
      else ; iTri_add =Msh%Nb_Cell +1; cavity(i)= iTri_add; END IF

      ! print *, "index to add",iTri_add, 'with vertexes :', iVer1,iVer2

      ! add to Msh
      iVer1 = boundary(i,1); Ver1 = Msh%Coord(iVer1,:)
      iVer2 = boundary(i,2); Ver2 = Msh%Coord(iVer2,:)
      surface = Signed_Surface(Ver1,Ver2,Point);

      Tri = Void_cell
      ALLOCATE(Tri%Sommets(3), Tri%CNeigh(3), Tri%Cell2Edg(3))
      if(surface< -eps0)  THEN ! triangle mal défini (mauvais sens)
        Tri%index = iTri_add
        Tri%NbSom = 3
        Tri%Sommets(1) = iVer2
        Tri%Sommets(2) = iVer1
        Tri%Sommets(3) = Msh%Nb_Ver
        Tri%Surface = -surface
        Tri%CRef = 1
        Tri%CNeigh(:) = 0
      else
        Tri%index = iTri_add
        Tri%NbSom = 3
        Tri%Sommets(1) = iVer1
        Tri%Sommets(2) = iVer2
        Tri%Sommets(3) = Msh%Nb_Ver
        Tri%Surface = surface
        Tri%CRef = 1
        Tri%CNeigh(:) = 0
      END IF

      ! linking the edge of the cavity to the outside
      iTri_bdry = boundary(i,3);  ! triangle de l'autre côté de la barrière voisin à Tri
      ! print *,"bdry", iTri_Bdry
      Tri%CNeigh(3) =iTri_bdry;
      IF(iTri_Bdry .NE. 0) THEN;      
      DO j=1,3
        Tri_bdry = Msh%Cells(iTri_bdry)
        jVer1 = Tri_bdry%Sommets(Tri2Edg(j,1))
        jVer2 = Tri_bdry%Sommets(Tri2Edg(j,2))
        if((iVer1==jVer1 .AND. iVer2 == jVer2) .OR. (iVer2 == jVer1 .AND. iVer1 == jVer2)) THEN
          ! print *,'bbbb',boundary(i,:)
          Msh%Cells(iTri_bdry)%CNeigh(j)=iTri_add;
        END IF
      END DO; 
      END IF

      sizeof_cavity = sizeof_cavity - 1;
      Msh%Nb_Cell = Msh%Nb_Cell +1 

      Msh%Cells(iTri_add) = Tri

      ! print *," addded new triangle : ", iTri_add
      ! print *," de sommets : ",Tri%Sommets(:)
    END DO

    ! print *,"size of boundary :",sizeof_boundary
    ! updating the neighbours of the star using a quadratic algorithm.
    DO i=1,sizeof_boundary
      iTri = Cavity(i)
      IF(iTri .NE. 0) THEN 
      Tri = Msh%Cells(iTri)
      iVer1 = Tri%Sommets(1); 
      iVer2 = Tri%Sommets(2)

      ! or just stock it inside an array
      DO j=1,sizeof_boundary

        jTri = cavity(j); 
        IF(jTri .NE. 0) THEN 
        Tri2 = Msh%Cells(jTri)

        ! print *, Msh%Nb_Ver
        jVer1 = Tri2%Sommets(1)
        jVer2 = Tri2%Sommets(2)

        IF(iTri .NE. jTri) THEN 
          ! print *,'------------',  iTri, jTri,'--------------'
          ! print *, iVer1,iVer2
          ! print *, jVer1,jVer2
          if(jVer2 == iVer1) THEN; Tri2%CNeigh(1) = iTri; Tri%CNeigh(2) = jTri; END IF
          if(jVer1 == iVer2) THEN; Tri2%CNeigh(2) = iTri; Tri%CNeigh(1) = jTri; END IF
          if(jVer1 == iVer1 .OR. jVer2 == iVer2) STOP "ERROR TRIANGLE NOT DEF CORRECTLY \n"
          Msh%Cells(jTri) = Tri2
        END IF
        END IF

      END DO

      ! CALL Out_Cell(Tri)
      
      Msh%Cells(iTri) = Tri
      END IF
      ! print *,"updated", iTri
      ! print *,Tri%index," sommest :", Tri%Sommets
    END DO

  END SUBROUTINE Starring_Point

  SUBROUTINE add_Point(Msh, Point)
    !---------------------------------------------
    ! Ajoute un point au maillage msh aux coordonnées données
    !
    ! Certaines données (VNeigh,Edg,Bdr,..) ne sont pas encore mis à jour, 
    ! il est conseiller de recalculer les voisins après une ou plusieurs itération de cette routine
    ! il est cependant possible de réajouter un point sans avoir à recalculer les voisins.
    !
    ! Présuppose que le point n'a pas déjà été ajouté
    ! Présuppose que le point à ajouté se trouve dans un triangle du maillage
    ! Présuppose que le maillage est triangulaire
    !---------------------------------------------

    IMPLICIT NONE
    Type(Mesh), INTENT(INOUT) :: Msh
    REAL(prec), DIMENSION(2) :: Point
    INTEGER, DIMENSION(Msh%Nb_Cell) :: cavity
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: boundary

    sizeof_cavity=0; sizeof_boundary =0;
    CALL digging_a_hole(Msh,Point, cavity);

    if(sizeof_cavity == 0) STOP "ERROR : point already there \n"

    ALLOCATE(boundary(3*sizeof_cavity,3))
    boundary = 0
    CALL Boundary_hole(Msh, Cavity, boundary ) 
    ! print *,"sizes :", sizeof_Cavity, sizeof_Boundary
    ! print *,"cavity :", Cavity(1:sizeof_Cavity)
    ! print *,"boundary:",boundary(1:sizeof_Boundary,3)
    CALL Realloc_Point(Msh)
    CALL deleting_cavity(Msh,cavity)
    Msh%Coord(Msh%Nb_Ver,:) = Point

    CALL Starring_Point(Msh,cavity,boundary,Point)

    Deallocate(boundary)


  END SUBROUTINE add_Point

  SUBROUTINE Realloc_Point(Msh)
    IMPLICIT NONE
    TYPE(Mesh) :: Msh, Msh_cpy
    TYPE(Cell) :: nCell
    INTEGER :: security
    INTEGER :: NbCellMax
    INTEGER :: Nb_Cell, Nb_Ver, Nb_Edg

    ! Msh_cpy = Msh
    security = 0
    NbCellMax= Msh%Nb_Cell + sizeof_boundary + security
    Nb_Ver  = Msh%Nb_Ver +1
    Nb_Edg  = Msh%Nb_Edg + sizeof_boundary  +1

    CALL REALLOC_Mesh(Msh, NbCellMax, Nb_Ver, Nb_Edg)
    Msh%Nb_Cell = NbCellMax
    Msh%Nb_Ver  = Nb_Ver
    Msh%Nb_Edg  = Nb_Edg

    Msh%Nb_Cell = Msh%Nb_Cell -( sizeof_boundary  + security) 
    Msh%max_neigh_Ver = max(Msh%max_neigh_Ver, sizeof_boundary)

    ! print *, 'realloc cells to', NbCellMax
    ! print *, Msh%Nb_Cell, Sizeof_cavity

  END SUBROUTINE Realloc_Point

  SUBROUTINE Mesh_Creation(Msh,Msh_D)
    !--------------------------------------------
    !
    ! prends un maillage avec uniquement des noeuds et un bords
    ! Calcul le maillage de Delaunay défini par cela
    ! Si le maillage possède d'autres propriété, elles sont ignorée
    !
    !-------------------------------------------
    IMPLICIT NONE
    TYPE(Mesh), INTENT(INOUT) :: Msh
    TYPE(Mesh) :: Msh_D
    TYPE(Cell) :: Cell_, Cello
    REAL(prec) :: epsi
    INTEGER, DIMENSION(2) :: Edg
    INTEGER :: iVer1,iVer2, i_Edg_d
    INTEGER :: i, iBdr, node
    INTEGER :: Nb_Cellf,Nb_Edgf
    INTEGER :: Color_suppr

    ! INTEGER, DIMENSION(Msh%Nb_Cell) :: cavity
    ! INTEGER, DIMENSION(:,:), ALLOCATABLE :: boundary

    ! Alloue une estimée du nombre de point/Cellules/Edges du maillage final
    Nb_Cellf = 4*Msh%Nb_Ver; Nb_Edgf = 4*Msh%Nb_Ver

    Msh_D%max_neigh_Cell=3
    Msh_D%Nb_EBdr = 0; Msh_D%Nb_VBdr =0; Msh_D%max_neigh_Ver = 0
    CALL ALLOCATE_Mesh(Msh_D,Nb_Ver = Msh%Nb_Ver +4 ,Nb_Cell = Nb_Cellf, Nb_Edg= Nb_Edgf)

    Msh_D%Cells(:) = Void_Cell
    ! crée 4 points en plus  qui seront la boite du maillage (+ epsilon )(!= eps0)
    ! ils seront enlevée du maillage final
    Msh_D%Nb_Ver = 4;

    epsi = 1._prec

    Msh_D%Coord(1,:) = (/Msh%x_deb -epsi, Msh%y_deb -epsi /)
    Msh_D%Coord(2,:) = (/Msh%x_deb -epsi, Msh%y_fin +epsi /)
    Msh_D%Coord(3,:) = (/Msh%x_fin +epsi, Msh%y_deb -epsi /)
    Msh_D%Coord(4,:) = (/Msh%x_fin +epsi, Msh%y_fin +epsi /)

    print *, Msh%x_deb, Msh%y_deb
    print *, Msh%x_fin, Msh%y_fin
    

    ALLOCATE(Cell_%Sommets(3), Cell_%CNeigh(3), Cell_%Cell2Edg(3))
    ALLOCATE(Cello%Sommets(3), Cello%CNeigh(3), Cello%Cell2Edg(3))

    Cell_%index = 1                   ; Cello%index = 2
    Cell_%NbSom = 3                   ; Cello%NbSom = 3
    Cell_%Sommets(1) = 3              ; Cello%Sommets(1) = 2
    Cell_%Sommets(2) = 2              ; Cello%Sommets(2) = 3
    Cell_%Sommets(3) = 1              ; Cello%Sommets(3) = 4
    Cell_%Surface = 0                 ; Cello%Surface = 0
    Cell_%CRef = 1                    ; Cello%CRef = 1  
    Cell_%CNeigh(:) = (/0,0,2/)       ; Cello%CNeigh(:) = (/0,0,1/)
    
    Msh_D%Cells(1) = Cell_             
    Msh_D%Cells(2) = Cello
    Msh_D%Nb_Cell =2
    Msh_D%NbNeigh = 0

    ! Ajoute tout les points du bord
    print *,"Nb Ver :", Msh_D%Nb_Ver
    CALL msh_neighbors_Tri(Msh_D)

    DO i=1,Msh%Nb_Ver
      write(*,fmt='("======================= adding ",i3,"e Point ", f8.3,f8.3," =================")') i, Msh%Coord(i,:)
      CALL add_Point(Msh_D, Msh%Coord(i,:))
      CALL msh_neighbors_Tri(Msh_D)

      Call Write_vtk("carre_last_.vtk", Msh_d)
      ! IF(i== 27) STOP
    END DO

    ! Remove the first 4 vertex and cells connected to them
    print *, Msh_D%Nb_Cell
    
    ! TODO : ADD enforce boundary

    ! inside edges are now also boundary 

    DO i=1,Msh%Nb_Edg
      Edg = Msh%Edg(i,:)
      i_Edg_d = Hash_find(Edg(1)+4,Edg(2)+4, Msh_D)
      IF(i_Edg_d == 0) STOP  'edge not found, Mesh Creation'
      Msh_d%EIsOnBdr(i_Edg_d) = Msh%EIsOnBdr(i)
      ! print *,Edg+4, Msh_d%EIsOnBdr(i_Edg_d)
    END DO

    CALL Connex_Comp(Msh_D)

    i = hash_find(1,2,Msh_D)
    i = Msh_D%LstObj(i,3)
    Color_suppr = Msh_D%Cells(i)%Color

    print *, Color_suppr
    Write(*, '("There is ",i2," Connex Component")') maxval(Msh_D%Cells%Color)

    Msh_d%Cells(:)%cref = Msh_d%Cells(:)%color
    Call WriteMesh("carre_last_.mesh", Msh_d)

    DO i=1,Msh_d%Nb_Cell
      Cell_ = Msh_d%Cells(i)
      ! print *,Cell_%color
      IF( Cell_%color == Color_Suppr) THEN 
        CALL Cell_Delete(i,Msh_D)
        ! Cell_ = Void_Cell
        ! Msh_d%Cells(i) = Void_Cell
      END IF
    END DO

    print *, Msh_D%Nb_Cell
    




  END SUBROUTINE Mesh_Creation

  SUBROUTINE Connex_Comp(Msh)
    IMPLICIT NONE
    TYPE(Mesh), INTENT(INOUT) :: Msh
    TYPE(Cell) :: Cell_
    INTEGER,    DIMENSION(Msh%Nb_Cell) :: influence_ring
    INTEGER :: nbr_influence, nbr_colored, id_last_seed, color
    INTEGER :: iCell, i_edg
    LOGICAL :: isvalid

    print *,"================== Coloring ====================="

    Msh%Cells(:)%Color = 0    ! stored color of the triangles
    influence_ring = 0        ! neighbohood of the seed, the pile
    nbr_influence = 1         ! top of the pile 
    nbr_colored   = 0         ! number of colored triangles
    id_last_seed  = 1         ! last seed
    color = 1                 ! color of the CC

    DO WHILE(id_last_seed .LT. Msh%Nb_Cell)
      ! INIT THE SEED
      influence_ring(1)= id_last_seed;  nbr_influence =1;
      iCell = id_last_seed

      ! loop on the rings of neighbors  
      DO WHILE(iCell .NE. 0 .AND. nbr_influence .GT. 0)
        
        iCell = influence_ring(nbr_influence)
        ! print *,iCell, Color, Msh%Cells(iCell)%CNeigh(:)
        ! print *,influence_ring(1:Nbr_influence)
        ! print *, "---------------------------"
        Cell_ = Msh%Cells(iCell)
        ! if the element isn't colored, color it and add it's valid edges

        if(Cell_%Color==0) THEN 
      
        Cell_%Color = color; 
        DO i_Edg= 1,3
        
          ! an edge is valid if the edge isn't a border and the neighboring triangle isn't colored
          ! and it is not in the pile
          isvalid = valid_edge(Msh,Cell_,i_edg)
          ! print *, isvalid
          if(isvalid .AND. .NOT. any(influence_ring(1:Nbr_influence-1) == Cell_%CNeigh(i_edg))) THEN
            ! if the edge is valid, add on top of the pile
            influence_ring(nbr_influence) = Cell_%CNeigh(i_edg)
            nbr_influence = nbr_influence + 1;
            ! print *,"added ", Cell_%CNeigh(i_edg)
          END IF
        END DO
        nbr_colored = nbr_colored + 1;
        nbr_influence = nbr_influence - 1;
        Msh%Cells(iCell) = Cell_
        ! if the element is colored, skip to the next element
        ! the pile is one less high.

        else;  nbr_influence = nbr_influence - 1;

        END IF

      END DO
      Color = Color + 1
      id_last_seed = 0
      ! print *, Msh%Cells(1:Msh%Nb_Cell)%Color
      DO iCell = 1,Msh%Nb_Cell
        IF(Msh%Cells(iCell)%Color ==0 ) THEN; id_last_seed = iCell; EXIT; END IF
      END DO
      IF(id_last_seed == 0) EXIT
      ! print *,id_last_seed

      ! IF(Color == 4) STOP
    END DO

  END SUBROUTINE Connex_Comp

  FUNCTION valid_edge(Msh,Cell_,i_edg)
    IMPLICIT NONE
    TYPE(Mesh), INTENT(IN) :: Msh
    TYPE(Cell), INTENT(IN) :: Cell_
    INTEGER,    INTENT(IN) :: i_edg

    TYPE(Cell)  :: Cell_Neigh
    LOGICAL :: valid_edge
    INTEGER :: iVer1, iVer2, edg

    valid_edge = .TRUE.

    ! is not a adjacent to void
    IF(Cell_%CNeigh(i_edg) == 0) THEN 
      valid_edge = .FALSE. 
      Return
    END IF

    ! is not a boundary edge
    iVer1 = Cell_%Sommets(tri2edg(i_edg,1))
    iVer2 = Cell_%Sommets(tri2edg(i_edg,2))
    edg = hash_find(iVer1,iVer2,Msh)

    IF(edg == 0) STOP "edg not found, Valid_edge"

    IF(Msh%EIsOnbdr(edg)) THEN 
      valid_edge = .FALSE. 
      Return
    END IF

    ! is not colored
    Cell_Neigh = Msh%Cells(Cell_%CNeigh(i_edg))
    IF(Cell_Neigh%color .NE. 0) THEN 
      valid_edge = .FALSE. 
      Return
    END IF

  END FUNCTION

END MODULE Mod_Delaunay
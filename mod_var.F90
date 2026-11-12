module mod_Var

  use precis
  !************************
  ! Variables globales du programme, accessibles dans tous les modules
  !************************
  IMPLICIT NONE

  TYPE Cell
    INTEGER :: index
    INTEGER :: cref                 ! Ref des mailles
    INTEGER :: NbSom               ! Nombre de voisins de chaque maille
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Sommets, CNeigh   ! Sommets assignée aux mailles et voisins des mailles
                                                                ! Sommets0 et Cneigh0 sont de dimension (nb_cell0, maxneigh) où la 2e coordonée est le stockage
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Cell2Edg            ! Pointer to EdgeTable
    REAL(prec), DIMENSION(2) :: Center          ! Coordonnées des centres des mailles du maillage , de dimension (nb_cell0, 2)
    LOGICAL :: CIsOnBdr

  END TYPE Cell

  TYPE MESH

    ! Mesh
    INTEGER :: DIM  ! restera 2
    REAL(prec) :: x_deb,x_fin,y_deb,y_fin ! bounding box

    INTEGER :: nb_cell                                         ! Nombre de mailles dans le maillage
    INTEGER :: Nb_Ver                                         ! nombre de sommets dans le maillage 
    INTEGER :: Nb_Edg                                           ! nombre d'aretes dans le maillage  = Nb_Ver0 + nb_cell0 - 2
    INTEGER :: max_neigh_Cell                                         ! nombre maximum de voisins d'une maille
    INTEGER :: max_neigh_Ver                                         ! nombre maximum de voisins d'un sommet

    ! Cell
    TYPE(CELL), ALLOCATABLE, DIMENSION(:) :: Cells

    ! Vertex
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Vref                ! Ref des sommets
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NbNeigh               ! Nombre de vertex voisins de chaque sommet
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VNeigh             ! Sommets voisins des sommets
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Ver2Cell            ! Mailles auquel appartient chaque sommet 
    REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: Coord           ! Coordonnées des sommets du maillage , de dimension (Nb_Ver0, 2)

    ! Edges 
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Edg                 ! Indices des deux sommets composant l'arete
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: EdgNeigh            ! Indices des deux cellules voisines de l'arete 

    ! Boundary
    INTEGER :: Nb_CBdr                                            ! nombre de cellules sur le bord du maillage 
    INTEGER :: Nb_VBdr                                            ! nombre de sommets sur le bord du maillage 
    INTEGER :: Nb_EBdr                                            ! nombre d'aretes sur le bord du maillage 
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Bdr                     ! Indices de l'arête sur le bord du maillage
    INTEGER, ALLOCATABLE, DIMENSION(:) :: VBdr                    ! Lst of Bdr vertex's index ex : (5,14,20,36,...)
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: VIsOnBdr
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: EIsOnbdr              

    ! HASH TABLE
    INTEGER :: SizHead                                          ! Maximum key value, the key is in [0,SizHead-1]
    INTEGER :: NbrObj                                           ! Number of objects in the hash table
    INTEGER, ALLOCATABLE, DIMENSION(:) ::  Head                 ! Head[key%(SizHead)] = link to the first object having this key in the LstObj list
    INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  LstObj             ! List of edges in the hash table (ver1,ver2,cell1,cell2) (cell2 !=0 if the edge is not on the border)


  END TYPE MESH

  REAL(prec), parameter :: eps0=0.1_prec**(2*prec-3)


  TYPE(Cell) :: Void_cell
  

  END MODULE mod_Var


module mod_Var

  use precis
  !************************
  ! Variables globales du programme, accessibles dans tous les modules
  !************************
  IMPLICIT NONE

  TYPE MESH

    ! Mesh
    INTEGER :: DIM  ! restera 2
    REAL(prec) :: x_deb,x_fin,y_deb,y_fin ! bounding box

    ! Cell
    INTEGER :: nb_cell                                         ! Nombre de mailles dans le maillage
    INTEGER :: maxneigh                                         ! nombre maximum de voisins d'une maille
    INTEGER :: Nb_CBdr                                           ! nombre d'aretes sur le bord du maillage 
    INTEGER, ALLOCATABLE, DIMENSION(:) :: cref                 ! Ref des mailles
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NbSom               ! Nombre de voisins de chaque maille
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Sommets, CNeigh   ! Sommets assignée aux mailles et voisins des mailles
                                                                ! Sommets0 et Cneigh0 sont de dimension (nb_cell0, maxneigh) où la 2e coordonée est le stockage
    REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: Center          ! Coordonnées des centres des mailles du maillage , de dimension (nb_cell0, 2)
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: CIsOnBdr

    ! Vertex
    INTEGER :: Nb_Ver                                         ! nombre de sommets dans le maillage 
    INTEGER :: Nb_VBdr                                           ! nombre d'aretes sur le bord du maillage 
    INTEGER :: max_neigh_Ver                                         ! nombre maximum de voisins d'un sommet
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: VIsOnBdr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: VBdr                ! Lst of Bdr vertex's index ex : (5,14,20,36,...)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Vref                ! Ref des sommets
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NbNeigh               ! Nombre de vertex voisins de chaque sommet
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VNeigh             ! Sommets voisins des sommets
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Ver2Cell            ! Mailles auquel appartient chaque sommet 
    REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: Coord           ! Coordonnées des sommets du maillage , de dimension (Nb_Ver0, 2)

    ! Edges 
    INTEGER :: Nb_Edg                                           ! nombre d'aretes dans le maillage  = Nb_Ver0 + nb_cell0 - 1
    INTEGER :: Nb_EBdr                                           ! nombre d'aretes sur le bord du maillage 
    INTEGER, ALLOCATABLE, DIMENSION(:) :: Bdr                 ! Indices des deux sommets composant l'arete
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: EIsOnbdr                 ! tableau de booléens pour savoir si une arete est sur le bord ou non
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Edg                 ! Indices des deux sommets composant l'arete
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: EdgNeigh            ! Indices des deux cellules voisines de l'arete 

  END TYPE MESH

  REAL(prec), parameter :: eps0=0.1_prec**(2*prec-3)


  END MODULE mod_Var


module mod_Var

  use precis
  !************************
  ! Variables globales du programme, accessibles dans tous les modules
  !************************
  IMPLICIT NONE

  ! Mesh
  INTEGER :: DIM  ! restera 2
  REAL(prec) :: x_deb,x_fin,y_deb,y_fin ! bounding box

  ! Cell
  INTEGER :: nb_cell0, nb_cell                                         ! Nombre de mailles dans le maillage
  INTEGER :: maxneigh0, maxneigh                                         ! nombre maximum de voisins d'une maille
  INTEGER, ALLOCATABLE, DIMENSION(:) :: cref0, cref                 ! Ref des mailles
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NbSom0, NbSom               ! Nombre de voisins de chaque maille
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Sommets0, Sommets, CNeigh0, CNeigh   ! Sommets assignée aux mailles et voisins des mailles
                                                              ! Sommets0 et Cneigh0 sont de dimension (nb_cell0, maxneigh) où la 2e coordonée est le stockage
  REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: Center0, Center          ! Coordonnées des centres des mailles du maillage , de dimension (nb_cell0, 2)

  ! Vertex
  INTEGER :: Nb_Ver0, Nb_Ver                                         ! nombre de sommets dans le maillage 
  INTEGER :: max_neigh_Ver0, max_neigh_Ver                                         ! nombre maximum de voisins d'un sommet
  INTEGER, ALLOCATABLE, DIMENSION(:) :: Vref0, Vref                ! Ref des sommets
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NbNeigh0, NbNeigh               ! Nombre de voisins de chaque sommet
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VNeigh0, VNeigh             ! Sommets voisins des sommets
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Ver2Cell0, Ver2Cell            ! Mailles auquel appartient chaque sommet 
  REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: Coord0, Coord           ! Coordonnées des sommets du maillage , de dimension (Nb_Ver0, 2)

  ! Edges 
  INTEGER :: Nb_Edg0                                           ! nombre d'aretes dans le maillage  = Nb_Ver0 + nb_cell0 - 1
  INTEGER :: Nb_bdr0                                           ! nombre d'aretes sur le bord du maillage 
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: is_bdr0, is_bdr                 ! tableau de booléens pour savoir si une arete est sur le bord ou non
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Edg0, Edg                 ! Indices des deux sommets composant l'arete
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Bdr0, Bdr                 ! Indices des deux sommets composant l'arete
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: EdgNeigh0, EdgNeigh            ! Indices des deux cellules voisines de l'arete 


  
END MODULE mod_Var


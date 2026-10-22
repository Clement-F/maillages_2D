MODULE mod_variables
  !
  !** Module contenant les declarations de variables generales, utilisees dans le programme
	!	principal et d'autres modules. Surtout pour la grille.
  !
  !** Sens des variables principales de structure de la grille:
  !     nbsom = nombre de sommets
  !     ncv = nombre de volumes de controle
  !     nptvois = plus grande valeur du pointeur ptvois, taille de nusom
  !     ptvois = pointeur des "voisins" d'une maille, cf ci-dessous
  !     nusom = sommets de chaque maille, a lire avec le pointeur ptvois
  !               nusom(ptvois(i)+1) -> nusom(ptvois(i+1)) : numero des sommets de la maille i,
  !                                                         parcourus dans le sens trigonometrique.
  !     nuvois = lu avec ptvois, indique le numero de la maille de l'autre cote.
  !               nuvois(ptvois(i)+1) -> nuvois(ptvois(i+1)) : par rapport a la maille i, donne
  !                  le numero de la maille de l'autre cote de l'arete [nusom(k-1),nusom(k)]
  !                  avec k=ptvois(i)+1,ptvois(i+1) (et k-1=ptvois(i+1) si k=ptvois(i)+1: on boucle
  !                  sur les aretes de la maille i uniquement)
  !     numar = lu avec ptvois, apres numerotation des aretes, indique les numeros des aretes de chaque maille.
  !               numar(ptvois(i)+1) -> numar(ptvois(i+1)) : numero des aretes de la maille i.
  !                             si k=ptvois(i)+1,ptvois(i+1), numar(k) est l'arete entre les
  !                             sommets nusom(k-1) et nusom(k) (avec k-1=ptvois(i+1) si k=ptvois(i)+1: on ne
  !                             boucle que sur les sommets de la maille i)
  !
  !     nar = nombre d'aretes
  !     xs et ys = coordonnees x et y des sommets
  !
  !     Exemple: *si ptvois = (0,4,7,11,...) alors :
  !                   a) les donnees correspondant au sommet 1 dans nusom et numar sont a
  !                       lire entre les cases 1 et 4 de ces tableaux
  !                   b) les donnees correspondant au sommet 2 dans nusom et numar sont a
  !                       lire entre les cases 5 et 7 de ces tableaux
  !                   c) etc...
  !          *si nusom = (14,54,21,35,85,96,5,4,7,...) alors:
  !                   a) les numeros des sommets de la maille 1, dans le sens trigo, sont
  !                      14,54,21,35
  !                   b) les numeros des sommets de la maille 2, dans le sens trigo, sont
  !                       85,96,5
  !                   c) etc...
  !          *si numar = (14,6,4,85,123,2,74,85,96,...) alors:
  !                   a) les numeros des aretes de la maille 1 sont 14,6,4,85; 
  !                               l'arete 14 correspond aux sommets 35->14; 
  !                               l'arete 6 correspond aux sommets 14->54;
  !                               l'arete 4 correspond aux sommets 54->21;
  !                               l'arete 85 correspond aux sommets 21->35;
  !                   b) etc...
  !
  !
  !************************

  IMPLICIT NONE

  ! Variables generales de grille, affectees lors de la lecture du fichier d.dat
  ! ou des calculs suivants.

  ! Variables lues dans le fichier
  INTEGER :: paradigme, solveur, type_sorties, type_decentre, nutest, nucbord
  DOUBLE PRECISION :: nu, alpha, theta, epsil

  ! Variables affectees par les procedures de construction/division du maillage
  INTEGER :: ncv, nptvois, nbsom
  INTEGER, ALLOCATABLE, DIMENSION (:) :: ptvois, nusom
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xs, ys
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: xcv, ycv

  ! Variables affectees par calcul_pointeurs
  INTEGER :: nar, nar_internes
  INTEGER, ALLOCATABLE, DIMENSION (:) :: nuvois, numar

  ! Variables affectees par cal_mcv_reg
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: mcv, diam
  DOUBLE PRECISION :: sizem, reg_maillage, airetot

	! Pour les valeurs une maille
	!		%p: valeur de la pression
	!		%px, %py: gradient de la pression
	!		%s: source
	!		%Dxx, Dxy, Dyy: diffusion
  TYPE datamaille
    DOUBLE PRECISION :: p, px, py, s, Dxx, Dxy, Dyy
  END TYPE datamaille

	! Pour la convection dans une maille
	!		%x, %y: composantes Vx, Vy au centre de maille
	!		%div: divergence au centre de maille
	!		%flux: alloue de taille ptvois(i+1)-ptvois(i),
	!				%flux(ja) = flux msVKs=\int_s V.n_Ks pour l'arete
	!										j=ptvois(i)+1,ptvois(i+1), avec ja=j-ptvois(i)
	!		%linf: max de |Vx|, |Vy| au centre de maille et aux sommets
	!		%lip: estimation de la constant de lipschiz, max des derivees
	!					de Vx et Vy au centre et sommets
  TYPE convection
    DOUBLE PRECISION :: x, y, div, linf, lip
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: flux
  END TYPE convection

	! variable generale pour calculer les donnees de convection dans une
	! maille. Doit etre declare ici car ne peut pas etre passe en transfert
	! de donnee vu qu'il a une composante allouable.
	TYPE(convection) :: conv

	! varsigma = 0 correspond aux conditions Dirichlet imposees en fort,
	!	varsigma = 1 a celle imposees faiblement (article Ern-Di Pietro)
	DOUBLE PRECISION :: varsigma


END MODULE mod_variables

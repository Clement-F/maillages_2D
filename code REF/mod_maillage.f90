MODULE mod_maillage
  !
  !** Module contenant les routines pour generer des grilles par subdivision
  !   a partir d'une grille de base
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
  ! ** Variables annexes:
  !       type_sorties: indique le type de fichiers de sortie qui doit etre genere
  !             pour la grille et la solution (-1=aucun, 0=tous, 1=xfig, 2=vtk)
  !
  !************************

	USE mod_variables
  USE mod_dessin

  IMPLICIT NONE

  ! Variables generales de grille, affectees lors de la lecture du fichier d.dat
  ! ou des calculs suivants. Doivent etre connues ici pour pouvoir etre affectees
  ! ou allouees, et seront connues dans le programme principal grace au USE du module
  !
  ! Variables lues dans le fichier
  INTEGER :: typepoly, ndiv, ncv0, nbsom0, ndiv1, ndiv2, typediv, ndivx, ndivy
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xs0, ys0, xcv0, ycv0
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: sommets0
  INTEGER, ALLOCATABLE, DIMENSION(:) :: na0
  DOUBLE PRECISION :: anisotropy, layer

  ! Pointeurs pour le maillage initial lu dans le fichier
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ptvois0, nusom0, nuvois0

  ! Variables affectees par les procedures de construction/division du maillage
  INTEGER :: maxna0

  CONTAINS


!---------------------------------------------------------------!
!-----------------PROCEDURES DE LECTURE ------------------------!
!---------------------------------------------------------------!

!----------------------------------------------------------------!
!*****************lec_cons_grille****************************************!

  SUBROUTINE lec_cons_grille
  !**********************************
  ! Lit les donnees du test et du maillage dans le fichier d.dat
  ! et construit la grille.
  ! Les numeros des sommets de chaque maille doivent etre indiques
  ! en ordre trigonometrique dans d.dat
  !************************
    INTEGER :: i, j
    INTEGER :: i0, i1, i2
    CHARACTER (LEN=15) :: texte_paradigme, texte_solveur, texte_type_sorties,&
    					texte_type_decentre
    CHARACTER (LEN=50) :: fichmaill
    DOUBLE PRECISION :: pi, xs_ori, ys_ori, xcv_ori, ycv_ori, r

    !** On ouvre le fichier de donnees
    OPEN(unit=10, File="d.dat",status="old")

    !** Lecture du test
		READ(10,*) type_sorties
		READ(10,*) nutest, epsil

		! paradigme=1 est HFV, paradigme=0 est MFV, paradigme=-1 est VF2
		paradigme=1

    !** Nombre de sommets du maillage de base, nombre de divisions,
    ! nombre de volumes de controle initial
    READ(10,*) typepoly

    !** On ne suppose ici que triangles ou quadrangles
    IF (typepoly == -3) THEN
      ! Triangles
      READ(10,*) ndiv
      READ(10,*) ncv0

      ALLOCATE(sommets0(ncv0,3))

      !** Calcul du nombre exact de triangles apres subdivision
      !** et variables associees aux centres de ces triangles
      ncv = ncv0 * ndiv * ndiv
      nptvois = 3*ncv

      !** Nombre de sommets du maillage de base, et lecture des coordonnees
      !** de ces sommets
      READ(10,*) nbsom0
      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      READ(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      !** Pour chaque triangle les numeros de ses sommets
      READ(10,*) ((sommets0(i,j),j=1,3),i=1,ncv0)

      CLOSE (10)

      !** On va ensuite diviser ce maillage de base ndiv fois.
      !** On prepare cette division en allouant les variables pour les
      !** coordonnees des sommets, et on appelle la procedure de division
      !** depuis le module 'mod_grille'
      nbsom = ncv0 * (ndiv+1) * (ndiv+2) / 2

      ! Division du maillage initial
      CALL divtri

      DEALLOCATE(xs0,ys0,sommets0)

      !** Calcul des centres des triangles dans le maillage complet
      CALL calxcvtri

    ELSEIF (typepoly == -9) THEN
      ! Quadrangles
      READ(10,*) ndiv
      READ(10,*) ncv0

      ALLOCATE(sommets0(ncv0,4))

      ! calcul du nombre exact de rectangles apres subdivision
      ncv = ncv0 * ndiv * ndiv
      nptvois = 4*ncv

      ! nombre de sommets du maillage de base
      read(10,*) nbsom0

      ALLOCATE(xs0(nbsom0),ys0(nbsom0))

      ! coordonnees des sommets
      read(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      ! pour chaque rectangle, les numeros de sommet
      read(10,*) ((sommets0(i,j),j=1,4),i=1,ncv0)

      CLOSE (10)

      ! nbsom est un majorant du nombre de sommets (certains comptes deux
      ! fois), qui sera remis a jour dans divquad
      nbsom = ncv0 * (ndiv+1) * (ndiv+1)

      ! Division du maillage initial
      CALL divquad

      DEALLOCATE(xs0,ys0,sommets0)

      ! calcul centres de maille
      CALL calxcvquad
      

    ELSEIF (typepoly == -20) THEN
      ! Mailles generales divises inegalement (les quadrangles et triangles sont divises,
      ! pas les mailles plus compliquees pour l'instant).
      READ(10,*) ndiv1, ndiv2, typediv
      READ(10,*) ncv0

      ! Sommets de base
      READ(10,*) nbsom0
      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      READ(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      ! Le fichier d.dat doit indiquer le nombre d'aretes
      ! de chaque maille initiale
      ALLOCATE(na0(ncv0))
      READ(10,*) (na0(i),i=1,ncv0)
      ! On cherche le max des na afin d'allouer sommets0
      maxna0 = 0
      DO i=1, ncv0
        IF (maxna0<na0(i)) maxna0 = na0(i)
      ENDDO

      ALLOCATE(sommets0(ncv0,maxna0))
      ! Lecture des sommets de chaque maille, dans l'ordre trigo
      DO i=1, ncv0
        READ(10,*) (sommets0(i,j),j=1,na0(i))
      ENDDO

      CLOSE (10)

      ! On fait la division, qui affectera ncv, nbsom, nptvois, ptvois, nusom,
      ! xs, ys. La routine divheterogene tourne a partir de pointeurs ptvois0/nusom0/nuvois0
      ! du maillage initial, de sorte a pouvoir tourner potentiellement sur de gros
      ! maillages (pour lesquels les donnees de type arete0(i,j) ne sont pas manipulables).
      ! Il faut donc commencer par generer ces pointeurs pour le maillage que l'on vient
      ! de lire dans le fichier d.dat
      CALL pointeurs_initial
      CALL divheterogene

      DEALLOCATE(sommets0,na0)

      ! calcul centres de maille
      ALLOCATE(xcv(ncv),ycv(ncv))
   		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,2)

    ELSEIF (typepoly == -30 .OR. typepoly == -31 .OR. typepoly == -40 .OR. typepoly == -41) THEN
      ! Lecture d'un maillage initial puis reproduction dans des subdivisions de l'ouvert,
      ! et enfin deformation eventuelle (si -40 ou -41)
      !
      ! On lit le maillage initial, puis on decoupe l'ouvert en ndivx selon x
      ! et ndivy selon y et on reproduit le maillage initial dans les subdivisions
      ! ainsi crees.
      ! Si typepoly == -30 ou -40, la reproduction se fait par simple translation (il faut donc que
      ! chaque sommet de x=1 corresponde a un sommet de x=0, et idem pour y)
      ! Si typeoly == -31 ou -41, la reproduction se fait par symetrie selon les axes et il
      ! n'y a donc pas de condition sur les sommets au bord 
      READ(10,*) ndivx, ndivy
      READ(10,*) ncv0

      ! Sommets de base
      READ(10,*) nbsom0
      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      READ(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      ! Le fichier d.dat doit indiquer le nombre d'aretes
      ! de chaque maille initiale
      ALLOCATE(na0(ncv0))
      READ(10,*) (na0(i),i=1,ncv0)
      ! On cherche le max des na afin d'allouer sommets0
      maxna0 = 0
      DO i=1, ncv0
        IF (maxna0<na0(i)) maxna0 = na0(i)
      ENDDO

      ALLOCATE(sommets0(ncv0,maxna0))
      ! Lecture des sommets de chaque maille, dans l'ordre trigo
      DO i=1, ncv0
        READ(10,*) (sommets0(i,j),j=1,na0(i))
      ENDDO

      CLOSE (10)

      ! On genere ptvois0 et nusom0 a partir du maillage initial lu
      CALL pointeurs_initial

			!
      ! CALCUL DES CENTRES DE MAILLES INITIAUX
      !	permet de modifier ces centres avant la reproduction, soit en changeant la formule ci-dessous soit en les
      ! faisant lire dans le fichier de donnees.
      ALLOCATE(xcv0(ncv0),ycv0(ncv0))
   		CALL calcul_centres(ncv0,ptvois0,nusom0,xs0,ys0,xcv0,ycv0,1)

			!
			! REPRODUCTION DU MAILLAGE INITIAL
      ! On reproduit le maillage dans les subdivision du maillage initial
      CALL reproduction_maillage

      DEALLOCATE(sommets0,na0,xcv0,ycv0)
      
      !
      ! DEFORMATION EVENTUELLE du maillage
      !
      IF (typepoly == -40 .OR. typepoly == -41) THEN
	      pi = 4D0*atan(1D0)
	      DO i=1, nbsom
	      	xs_ori = xs(i)
	      	ys_ori = ys(i)
  !	      	xs(i) = xs_ori + 0.05D0*sin(2D0*pi*xs_ori)*sin(2*pi*ys_ori)
  !	      	ys(i) = ys_ori + 0.05D0*sin(2D0*pi*xs_ori)*sin(2*pi*ys_ori)

          IF (xs_ori > 0. .AND. xs_ori<1. .AND. ys_ori>0. .AND. ys_ori<1.) THEN
            CALL RANDOM_NUMBER(r)   
	        	xs(i) = xs_ori + 0.05D0*r
            CALL RANDOM_NUMBER(r)   
	        	ys(i) = ys_ori + 0.05D0*r
          ENDIF
	      ENDDO
	      DO i=1, ncv
	      	xcv_ori = xcv(i)
	      	ycv_ori = ycv(i)
	      	xcv(i) = xcv_ori + 0.05D0*sin(2D0*pi*xcv_ori)*sin(2*pi*ycv_ori)
	      	ycv(i) = ycv_ori + 0.05D0*sin(2D0*pi*xcv_ori)*sin(2*pi*ycv_ori)
	      ENDDO
	      
	      !
	      ! RE-CALCUL EVENTUEL DES CENTRES
	      !
	  		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,5)
			ENDIF      

    ELSEIF (typepoly == -50) THEN
      ! Lecture d'un maillage initial puis reproduction (par translation) dans des subdivisions de l'ouvert,
      ! de maniere anisotrope
      READ(10,*) ndiv, anisotropy
      READ(10,*) ncv0

      ! Sommets de base
      READ(10,*) nbsom0
      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      READ(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      ! Le fichier d.dat doit indiquer le nombre d'aretes
      ! de chaque maille initiale
      ALLOCATE(na0(ncv0))
      READ(10,*) (na0(i),i=1,ncv0)
      ! On cherche le max des na afin d'allouer sommets0
      maxna0 = 0
      DO i=1, ncv0
        IF (maxna0<na0(i)) maxna0 = na0(i)
      ENDDO

      ALLOCATE(sommets0(ncv0,maxna0))
      ! Lecture des sommets de chaque maille, dans l'ordre trigo
      DO i=1, ncv0
        READ(10,*) (sommets0(i,j),j=1,na0(i))
      ENDDO

      CLOSE (10)

      ! On genere ptvois0 et nusom0 a partir du maillage initial lu
      CALL pointeurs_initial

			!
      ! CALCUL DES CENTRES DE MAILLES INITIAUX
      !	permet de modifier ces centres avant la reproduction, soit en changeant la formule ci-dessous soit en les
      ! faisant lire dans le fichier de donnees.
      ALLOCATE(xcv0(ncv0),ycv0(ncv0))
   		CALL calcul_centres(ncv0,ptvois0,nusom0,xs0,ys0,xcv0,ycv0,1)

			!
			! REPRODUCTION DU MAILLAGE INITIAL
      ! On reproduit le maillage dans les subdivision du maillage initial
      CALL reproduction_maillage_anisotrope

      DEALLOCATE(sommets0,na0,xcv0,ycv0)
      

    ELSEIF (typepoly == -100) THEN
      ! Lecture dans un fichier non formatte et deux divisions
      READ(10,*) fichmaill
      READ(10,*) ndiv1, ndiv2, typediv
      READ(10,*) layer
      CLOSE (10)

      ! Lecture des donnees du maillage initial puis divisions
      CALL acces_fichier_maillage("lire",fichmaill)
      CALL divheterogene

      ! calcul centres de maille
      ALLOCATE(xcv(ncv),ycv(ncv))
   		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,2)

		ELSEIF (typepoly == -150) THEN
			! Lecture d'un fichier de maillage "typ2"
			READ(10,"(A)") fichmaill			
			CALL acces_maillage_typ2("lire",fichmaill)

		ELSEIF (typepoly == -230 .OR. typepoly == -231) THEN
			! Creation d'un maillage dual
      ! 	1. Lecture dans le fichier du maillage initial
      !		2. Creation du maillage primal par reproduction de ce maillage initial
			!		3. Deformation eventuelle du maillage primal
      !		4. Construction d'un maillage dual en reliant les centres des mailles primales
      !			autour de chaque sommet primal.

			!
			! LECTURE MAILLAGE INITIAL
			!
      READ(10,*) ndivx, ndivy
      READ(10,*) ncv0

      ! Sommets de base
      READ(10,*) nbsom0
      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      READ(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      ! Le fichier d.dat doit indiquer le nombre d'aretes
      ! de chaque maille initiale
      ALLOCATE(na0(ncv0))
      READ(10,*) (na0(i),i=1,ncv0)
      ! On cherche le max des na afin d'allouer sommets0
      maxna0 = 0
      DO i=1, ncv0
        IF (maxna0<na0(i)) maxna0 = na0(i)
      ENDDO

      ALLOCATE(sommets0(ncv0,maxna0))
      ! Lecture des sommets de chaque maille, dans l'ordre trigo
      DO i=1, ncv0
        READ(10,*) (sommets0(i,j),j=1,na0(i))
      ENDDO

      CLOSE (10)

			!
			!	CREATION MAILLAGE PRIMAL PAR REPRODUCTION
			!
      ! On genere ptvois0 et nusom0 a partir du maillage initial lu
      CALL pointeurs_initial

      ! Calcul des centres de maille initiaux (on peut choisir plusieurs
      !	methodes de calcul, voir calcul_centres)
      ALLOCATE(xcv0(ncv0),ycv0(ncv0))
   		CALL calcul_centres(ncv0,ptvois0,nusom0,xs0,ys0,xcv0,ycv0,0)

      ! Reproduction du maillage initial
			typepoly = typepoly + 200
      CALL reproduction_maillage
			typepoly = typepoly - 200

      DEALLOCATE(sommets0,na0,xcv0,ycv0)

      ! Calcul des pointeurs du maillage primal (les donnees
      ! ncv, nbsom, ptvois, nusom, nuvois, numar, xs, ys, xcv, ycv correspondent
      ! ici au maillage primal)
      CALL calcul_pointeurs(0)

  
      !
      ! DEFORMATION EVENTUELLE du maillage primal
      !
      pi = 4D0*atan(1D0)
      DO i=1, nbsom
      	xs_ori = xs(i)
      	ys_ori = ys(i)
      	xs(i) = xs_ori + 0.1D0*sin(2D0*pi*xs_ori)*sin(2*pi*ys_ori)
      	ys(i) = ys_ori + 0.1D0*sin(2D0*pi*xs_ori)*sin(2*pi*ys_ori)
      ENDDO
  		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,1)
  !			CALL dessinergrille(type_sorties,"gridp","pasnumeros")


			!
			! CREATION MAILLAGE DUAL
			!
			!	nbsom, ncv, nptvois, ptvois, nusom, xcv, ycv, xs, ys deviennent ceux du
			! maillage dual
			CALL calcul_maillage_dual
			ALLOCATE(xcv(ncv),ycv(ncv))
  		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,0)

		ELSEIF (typepoly == -430) THEN
			! Lecture maillage triangulaire comme typepoly=-3, avec divisions, puis maillage dual
      ! Triangles
      READ(10,*) ndiv
      READ(10,*) ncv0

      ALLOCATE(sommets0(ncv0,3))

      !** Calcul du nombre exact de triangles apres subdivision
      !** et variables associees aux centres de ces triangles
      ncv = ncv0 * ndiv * ndiv
      nptvois = 3*ncv

      !** Nombre de sommets du maillage de base, et lecture des coordonnees
      !** de ces sommets
      READ(10,*) nbsom0
      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      READ(10,*) (xs0(i),ys0(i),i=1,nbsom0)

      !** Pour chaque triangle les numeros de ses sommets
      READ(10,*) ((sommets0(i,j),j=1,3),i=1,ncv0)

      CLOSE (10)

      !** On va ensuite diviser ce maillage de base ndiv fois.
      !** On prepare cette division en allouant les variables pour les
      !** coordonnees des sommets, et on appelle la procedure de division
      !** depuis le module 'mod_grille'
      nbsom = ncv0 * (ndiv+1) * (ndiv+2) / 2

      ! Division du maillage initial
      CALL divtri

      ! Calcul des centres de maille initiaux (on peut choisir plusieurs
      !	methodes de calcul, voir calcul_centres)
      ALLOCATE(xcv(ncv),ycv(ncv))
   		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,1)

      ! Calcul des pointeurs du maillage primal (les donnees
      ! ncv, nbsom, ptvois, nusom, nuvois, numar, xs, ys, xcv, ycv correspondent
      ! ici au maillage primal)
      CALL calcul_pointeurs(0)

      !
      ! DEFORMATION EVENTUELLE du maillage primal
      !
      pi = 4D0*atan(1D0)
      DO i=1, nbsom
      	xs_ori = xs(i)
      	ys_ori = ys(i)
      	xs(i) = xs_ori + 0.05D0*sin(2D0*pi*xs_ori)*sin(2*pi*ys_ori)
      	ys(i) = ys_ori + 0.05D0*sin(2D0*pi*xs_ori)*sin(2*pi*ys_ori)
      ENDDO
  		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,1)
  !			CALL dessinergrille(type_sorties,"gridp","pasnumeros")


			!
			! CREATION MAILLAGE DUAL
			!
			!	nbsom, ncv, nptvois, ptvois, nusom, xcv, ycv, xs, ys deviennent ceux du
			! maillage dual
			CALL calcul_maillage_dual
			ALLOCATE(xcv(ncv),ycv(ncv))
  		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,0)
			
      DEALLOCATE(xs0,ys0,sommets0)
  !stop

    ELSE
      print *,"typepoly inconnu",typepoly
      STOP

    ENDIF

		! Affichages
		IF (paradigme == 1) THEN
			write(texte_paradigme,"(A)") "HFV"
		ELSEIF (paradigme == 0) THEN
			write(texte_paradigme,"(A)") "MFV"
    ELSEIF (paradigme == -1) THEN
      write(texte_paradigme,"(A)") "VF2"
    ELSE
      print *,"Paradigme inconnu:", paradigme
      STOP
		ENDIF
		IF (solveur == 0) THEN
			write(texte_solveur,"(A)") "perso"
		ELSEIF (solveur == 1) THEN
			write(texte_solveur,"(A)") "ma41ad"
		ELSEIF (solveur == 2) THEN
			write(texte_solveur,"(A)") "hsl_ma48"
		ELSEIF (solveur == 3) THEN
			write(texte_solveur,"(A)") "BiCGStab"
		ELSE
			write(texte_solveur,"(A)") "GMRes"
		ENDIF
		IF (type_sorties == 0) THEN
			write(texte_type_sorties,"(A)") "tous"
		ELSEIF (type_sorties == 1) THEN
			write(texte_type_sorties,"(A)") "xfig"
		ELSEIF (type_sorties == 2) THEN
			write(texte_type_sorties,"(A)") "vtk"
		ELSE
			write(texte_type_sorties,"(A)") "aucun"
		ENDIF
		IF (type_decentre == 0) THEN
			write(texte_type_decentre,"(A)") "aretes"
		ELSE
			write(texte_type_decentre,"(A)") "mailles"
		ENDIF
		print "(2X,A,A,A,A,2X,A,A,A,A,A)","Paradigme ",trim(texte_paradigme),", decentrement ",trim(texte_type_decentre),&
							" (solveur: ",trim(texte_solveur),", figures: ",trim(texte_type_sorties),")"

    print "(A,I4,6X,A,I4,6X,A,I5)",'nutest=',nutest,'nucbord=',nucbord,'typepoly=',typepoly
    IF (typepoly == -3) THEN
      print "(A,I4,6X,A,I9,6X,A,I9)","ndiv=",ndiv,"nb sommets=",nbsom,"nb triangles=", ncv
    ELSEIF (typepoly == -9) THEN
      print "(A,I4,6X,A,I9,6X,A,I9)","ndiv=",ndiv,"nb sommets=",nbsom,"nb rectangles=", ncv
    ELSEIF (typepoly == -20) THEN
      print "(A,I4,3X,A,I4,6X,A,I9,6X,A,I9)","ndiv1=",ndiv1,"ndiv2=",ndiv2,"nb sommets=",nbsom,"nb mailles=", ncv
    ELSEIF (typepoly == -30 .OR. typepoly == -31) THEN
      print "(A,I4,3X,A,I4,6X,A,I9,6X,A,I9)","ndivx=",ndivx,"ndivy=",ndivy,"nb sommets=",nbsom,"nb mailles=", ncv
    ELSEIF (typepoly == -150) THEN
      print "(A,X,A)","Lecture du fichier maillage",fichmaill
    ELSEIF (typepoly == -230 .OR. typepoly == -231) THEN
      print "(A,I4,3X,A,I4,6X,A,I9,6X,A,I9)","ndivx primal=",ndivx,"ndivy primal=",ndivy,&
      		"nb sommets=",nbsom,"nb mailles=", ncv
    ENDIF

  !			CALL acces_maillage_typ2("ecrire","mesh.typ2")


  END SUBROUTINE lec_cons_grille


!----------------------------------------------------------------!
!***************** acces_fichier_maillage ***********************************!

  SUBROUTINE acces_fichier_maillage(action,nomfichier,nbsom,ncv,xs,ys,ptvois,nusom,nuvois)
    !*************************************************************!
    ! Ecrit ou lit un maillage (donnees nbsom, ncv, xs, ys,
    ! ptvois, nusom, nuvois) dans un fichier, sous forme non-formattee
    !
    !   action = chaine de caractere "lire" ou "ecrire"
    !   nomfichier = nom du fichier a lire/ecrire
    !
    !   les donnees "nbsom, ncv, xs, ys, ptvois, nusom, nuvois" sont celles qui
    !   doivent etre ecrites dans le fichier, et ne sont passees que si action="ecrire"
    !
    !   Dans le cas action="lire", les donnees nbsom0, ncv0, xs0, ys0,
    !         ptvois0, nusom0, nuvois0 seront affectes.
    !
    !**************************************************************!
    CHARACTER(LEN=*) :: action, nomfichier
    ! Optionnels
    INTEGER, OPTIONAL :: nbsom, ncv
    DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xs, ys
    INTEGER, DIMENSION(:), OPTIONAL :: ptvois, nusom, nuvois

    INTEGER :: i

    IF (action == "ecrire") THEN
      OPEN(unit=44, file=nomfichier, form="unformatted", status="unknown")
      write(44) nbsom, ncv
      write(44) (xs(i),ys(i),i=1,nbsom)
      write(44) (ptvois(i),i=1,ncv+1)
      write(44) (nusom(i),i=1,ptvois(ncv+1))
      write(44) (nuvois(i),i=1,ptvois(ncv+1))
      CLOSE(44)

    ELSEIF (action == "lire") THEN

      OPEN(unit=44, file=nomfichier, form="unformatted", status="unknown")

      read(44) nbsom0, ncv0

      ALLOCATE(xs0(nbsom0),ys0(nbsom0))
      read(44) (xs0(i),ys0(i),i=1,nbsom0)

      ALLOCATE(ptvois0(ncv0+1))
      read(44) (ptvois0(i),i=1,ncv0+1)

      ALLOCATE(nusom0(ptvois0(ncv0+1)),nuvois0(ptvois0(ncv0+1)))
      read(44) (nusom0(i),i=1,ptvois0(ncv0+1))
      read(44) (nuvois0(i),i=1,ptvois0(ncv0+1))

      CLOSE(44)

    ELSE
      print *,"Action inconnue dans acces_fichier_maillage"
      STOP
    ENDIF

    END SUBROUTINE acces_fichier_maillage


  ! SUBROUTINE read_msh(nomfichier,nbsom,ncv,xs,ys,ptvois,nusom,nuvois)
  !   IMPLICIT NONE

  !   CHARACTER(LEN=*) :: nomfichier
  !   ! Optionnels
  !   INTEGER, OPTIONAL :: nbsom, ncv
  !   DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xs, ys
  !   INTEGER, DIMENSION(:), OPTIONAL :: ptvois, nusom, nuvois

  !   INTEGER :: nelem,pt,npt,ncell,nface,ios_back,ntype

  !   INTEGER,DIMENSION(9) :: temp
  !   DOUBLE PRECISION :: tempr
  !   DOUBLE PRECISION,DIMENSION(2) :: XX1

  !   OPEN(11,file=TRIM(nomfichier),status="old")
  !   CALL JUMP(4,11)

  !   READ(11,*) nbsom0
  !   ALLOCATE(xs0(nbsom0),ys0(nbsom0))
    
  !   !! Grid points allocation !!
  !   DO pt=1,nbsom0
  !       READ(11,*) temp(1),XX1,tempr

  !       point(pt)%Coord=XX1
  !   END DO
  !   CALL JUMP(2,11)


  !   !! Determine the number of tri !!
  !   READ(11,*) nelem

  !   ALLOCATE(bdr_edge_nodes(nelem,2),bdr_edge_ref(nelem),bdr_edge_used(nelem))
  !   nb_bdr_edges=0
  !   bdr_edge_used=0

  !   bcl_face:DO nface=1,nelem
  !      READ(11,*) temp(1:7)

  !      IF (temp(2)/=1) EXIT bcl_face

  !      nb_bdr_edges=nb_bdr_edges+1
  !      bdr_edge_nodes(nb_bdr_edges,:)=(/temp(6),temp(7)/)
  !      bdr_edge_ref(nb_bdr_edges)=temp(4)
       
  !      point(temp(6))%ref=temp(4)
  !   END DO bcl_face

  !   nb_face_bdr_gmsh=nb_bdr_edges     ! Number of bdry faces !
  !   nb_cell=nelem-nb_face_bdr_gmsh    ! Number of cells      !
  !   ALLOCATE(mesh(nb_cell))

  !   CLOSE(11)
  !   OPEN(11,file=TRIM(name_mesh),status="old")
  !   CALL JUMP(8+nb_pts+nb_face_bdr_gmsh,11)

  !   !! Mesh nodes and number of neighboring cells of the node !!
  !   DO ncell=1,nb_cell
  !      READ(11,*) temp(1:8)
  !      BACKSPACE(11, IOSTAT=ios_back)
  !      IF (ios_back /= 0) THEN
  !         WRITE(*,*) "ERROR IN BACKSPACE IN  SUBROUTINE read_msh"
  !         STOP
  !      END IF
       
  !      ntype=temp(2)+1
  !      mesh(ncell)%nb_face=ntype
  !      ALLOCATE(mesh(ncell)%node(ntype),mesh(ncell)%neigh(ntype))
  !      ALLOCATE(mesh(ncell)%length(ntype),mesh(ncell)%Normal(ntype,2))
       
  !      READ(11,*) temp(1:5+ntype)
  !      mesh(ncell)%ref=temp(4)
  !      mesh(ncell)%node(1:ntype)=temp(6:5+ntype)

  !      DO npt=1,ntype
  !         pt=mesh(ncell)%node(npt)
  !         point(pt)%nb_neigh=point(pt)%nb_neigh+1
  !      END DO
  !   END DO

  !   CLOSE(11)

  !   END SUBROUTINE read_msh
  ! SUBROUTINE acces_fichier_maillage_paraview(action,nomfichier,nbsom,ncv,xs,ys,ptvois,nusom,nuvois)
    !*************************************************************!
    ! Ecrit ou lit un maillage (donnees nbsom, ncv, xs, ys,
    ! ptvois, nusom, nuvois) dans un fichier, sous format .vtk pour Paraview
    !
    !   action = chaine de caractere "lire" ou "ecrire"
    !   nomfichier = nom du fichier a lire/ecrire
    !
    !   les donnees "nbsom, ncv, xs, ys, ptvois, nusom, nuvois" sont celles qui
    !   doivent etre ecrites dans le fichier, et ne sont passees que si action="ecrire"
    !
    !   Dans le cas action="lire", les donnees nbsom0, ncv0, xs0, ys0,
    !         ptvois0, nusom0, nuvois0 seront affectes.
    !
    !**************************************************************!
    ! CHARACTER(LEN=*) :: action, nomfichier
    ! ! Optionnels
    ! INTEGER, OPTIONAL :: nbsom, ncv
    ! DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: xs, ys
    ! INTEGER, DIMENSION(:), OPTIONAL :: ptvois, nusom, nuvois

    ! INTEGER :: i,j, temp_i
    ! DOUBLE PRECISION :: temp_d
    ! CHARACTER(1) :: temp_c

    ! IF (action == "ecrire") THEN
      
    !   OPEN(unit=44, File=nomfichier, status="unknown")
    !   WRITE(44,"(A)") "# vtk DataFile Version 2.0"
    !   WRITE(44,"(A)") "Grille"
    !   WRITE(44,"(A)") "ASCII"

    !   WRITE(44,"(A)") "DATASET POLYDATA"
    !   WRITE(44,"(A,I10,X,A)") "POINTS",nbsom,"float"
    !   DO i=1, nbsom
    !     WRITE(44,"(E12.5,E12.5,E12.5)") xs(i),ys(i),0D0
    !   ENDDO
    !   WRITE(44,"(A,I10,X,I10)") "POLYGONS",ncv,ptvois(ncv+1)+ncv
    !   DO i=1, ncv
    !     WRITE(44,*) ptvois(i+1)-ptvois(i),(nusom(j)-1,j=ptvois(i)+1,ptvois(i+1))
    !   ENDDO

    !   CLOSE(44)

    ! ELSEIF (action == "lire") THEN

    !   OPEN(unit=44, File=nomfichier, status="unknown")
    !   READ(44,"(A)") temp_c
    !   READ(44,"(A)") temp_c
    !   READ(44,"(A)") temp_c
    !   READ(44,"(A)") temp_c
    !   READ(44,"(A,I10,X,A)") temp_c,nbsom,temp_c

    !   ALLOCATE(xs0(nbsom0),ys0(nbsom0))
    !   DO i=1, nbsom
    !     READ(44,"(E12.5,E12.5,E12.5)") xs(i),ys(i),temp_d
    !   ENDDO

    !   READ(44,"(A,I10,X,I10)") temp_c,ncv,temp_i

    !   ALLOCATE(ptvois0(ncv0+1))
    !   ptvois(ncv+1) = temp_i-ncv

    !   DO i=1, ncv
    !     READ(44,*) temp_i
    !     ptvois(i) = ptvois(i-1)+temp_i
    !     READ(44,*) (nusom(j),j= ptvois(i)+1,temp_i)
    !   ENDDO

    !   CLOSE(44)

    ! ELSE
    !   print *,"Action inconnue dans acces_fichier_maillage_paraview"
    !   STOP
    ! ENDIF

    ! END SUBROUTINE acces_fichier_maillage_paraview


!----------------------------------------------------------------!
!***************** acces_maillage_typ2 ***********************************!

  SUBROUTINE acces_maillage_typ2(afaire,nomfichier)
  !*************************************************************!
  ! Ecrit ou lit un maillage (utilisant les donnees globales "nbsom, xs, ys, ncv,
  !	xcv, ycv, nptvois, ptvois, nusom")	dans un fichier, sous forme "Bench FVCA5 - typ2".
  !	Lors de la lecture, les variables "xs, ys, xcv, ycv, nusom" sont aussi
  !	allouees.
  !
  !   afaire = chaine de caractere "lire" ou "ecrire"
  !   nomfichier = nom du fichier a lire/ecrire
  !
	!	Format typ2 (avec centres)
	!			"Vertices"
	!     -> the  total number of vertices
	!     -> the coordinates of the vertices
	!			"cells"
	!     -> the number of control volumes
	!     -> the control volumes with
	!               + the total number of vertices on the boundary of volume
	!               + the number of each of their vertices ordered clockwise
	!			"centers"
	!			-> coordinates of the "centers" of the cells, given in the
	!				same order as the cells
	!
	!
  !**************************************************************!
    CHARACTER(LEN=*) :: afaire, nomfichier

    INTEGER :: i, j, na
		CHARACTER(1) :: ligneinutile

    IF (afaire == "ecrire") THEN

      OPEN(unit=55, file=nomfichier, form="formatted", status="unknown")
			write(55,"(A)") "Vertices"
      write(55,"(I10)") nbsom
      DO i=1,nbsom
      	write(55,*) xs(i),ys(i)
      ENDDO
			write(55,"(A)") "cells"
			write(55,"(I10)") ncv
			DO i=1,ncv
				na = ptvois(i+1)-ptvois(i)
				write(55,*) na,(nusom(j),j=ptvois(i)+1,ptvois(i+1))
			ENDDO
			write(55,"(A)") "centers"
			DO i=1,ncv
				write(55,*) xcv(i),ycv(i)
			ENDDO

      CLOSE(55)

    ELSEIF (afaire == "lire") THEN

      OPEN(unit=55, file=nomfichier, form="formatted", status="unknown")

			! Affectation des sommets
			read(55,"(A)") ligneinutile
      read(55,"(I1000)") nbsom
      ALLOCATE(xs(nbsom),ys(nbsom))
      DO i=1,nbsom
      	read(55,*) xs(i),ys(i)
      ENDDO
			read(55,"(A)") ligneinutile
			read(55,"(I100)") ncv

			! Premier passage pour allouer ptvois en lisant le nombre de
			! sommets de chaque maille
			ALLOCATE(ptvois(ncv+1))
			ptvois(1)=0
			DO i=1,ncv
				read(55,*) na
				ptvois(i+1) = ptvois(i)+na
			ENDDO
			nptvois = ptvois(ncv+1)

			ALLOCATE(xcv(ncv),ycv(ncv))
			! Affectation des centres de mailles / certain typ2 uniquement
			! read(55,"(A)") ligneinutile
			! DO i=1,ncv
			! 	read(55,*) xcv(i),ycv(i)
			! ENDDO

      CLOSE(55)

			! Deuxieme passage pour lire les sommets de chaque maille
      OPEN(unit=55, file=nomfichier, form="formatted", status="unknown")

			DO i=1,4+nbsom
				read(55,"(A)") ligneinutile
			ENDDO
			ALLOCATE(nusom(nptvois))
			DO i=1,ncv
				read(55,*) na,(nusom(j),j=ptvois(i)+1,ptvois(i+1))
			ENDDO

			! Creation centre de mailles, si pas lus en haut
			CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,0)

      CLOSE(55)

    ELSE
      print *,"Action inconnue dans acces_maillage_typ2"
      STOP
    ENDIF

    END SUBROUTINE acces_maillage_typ2


!----------------------------------------------------------------!
!***************** acces_maillage_typ1 ***********************************!

  SUBROUTINE acces_maillage_typ1(afaire,nomfichier)
  !*************************************************************!
  ! Ecrit ou lit un maillage (utilisant les donnees globales "nbsom, xs, ys, ncv,
  !	xcv, ycv, nptvois, ptvois, nusom")	dans un fichier, sous forme "Bench FVCA5 - typ1".
  !	Lors de la lecture, les variables "xs, ys, xcv, ycv, nusom" sont aussi
  !	allouees.
  !
  !   afaire = chaine de caractere "lire" ou "ecrire"
  !   nomfichier = nom du fichier a lire/ecrire
  !
	!	Format typ1 (avec centres)
	!			"Vertices"
	!     -> the  total number of vertices
	!     -> the coordinates of the vertices
	!			"Triangles"
	!			-> nb of triangles
	!			-> vertices of the triangles (clockwise)
	!			"quadrangles"
	!			-> nb of quads
	!			-> vertices of the quads (clockwise)
	!			"edges of the boundary"
	!			-> nb of edges on the boundary
	!			-> list of edges by vertices (nb of vertices of the edges)
	!			"all edges"
	!			-> nb of edges
	!			-> vertices of edges + cells (with 0 if boundary cell)
	!
	!
  !**************************************************************!
    CHARACTER(LEN=*) :: afaire, nomfichier

    INTEGER :: i, j, na, j1, j2, narb
		CHARACTER(1) :: ligneinutile

    IF (afaire == "ecrire") THEN

      OPEN(unit=55, file=nomfichier, form="formatted", status="unknown")
			write(55,"(A)") "Vertices"
      write(55,"(I10)") nbsom
      DO i=1,nbsom
      	write(55,*) xs(i),ys(i)
      ENDDO
			write(55,"(A)") "triangles"
			write(55,"(I10)") 0
			write(55,"(A)") "quadrangles"
			write(55,"(I10)") ncv
			DO i=1,ncv
				write(55,*) (nusom(j),j=ptvois(i)+1,ptvois(i+1))
			ENDDO
			!
			narb = 0
			DO i=1,ncv
				DO j=ptvois(i)+1,ptvois(i+1)
					IF (nuvois(j)==0) narb=narb+1
				ENDDO
			ENDDO
			write(55,"(A)") "edges of the boundary"
			write(55,"(I10)") narb
			DO i=1,ncv
				j1=nusom(ptvois(i+1))
				DO j=ptvois(i)+1,ptvois(i+1)
					j2=nusom(j)
					IF (nuvois(j)==0) THEN
						write(55,"(I10,I10)") j1,j2
					ENDIF
					j1=j2
				ENDDO
			ENDDO
			!
			write(55,"(A)") "all edges"
			write(55,"(I10)") nar
			DO i=1,ncv
				j1=nusom(ptvois(i+1))
				DO j=ptvois(i)+1,ptvois(i+1)
					j2=nusom(j)
					IF (i>nuvois(j)) THEN
						write(55,"(I10,I10,I10,I10)") j1,j2,i,nuvois(j)
					ENDIF
					j1=j2
				ENDDO
			ENDDO
	

      CLOSE(55)

    ELSEIF (afaire == "lire") THEN


    ELSE
      print *,"Action inconnue dans acces_maillage_typ1"
      STOP
    ENDIF

    END SUBROUTINE acces_maillage_typ1



!---------------------------------------------------------------!
!---------------- PROCEDURES DE DIVISION -----------------------!
!---------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------!
!******************** TRIANGLES ***************************************************************************!

  SUBROUTINE divtri
  !************** Division du maillage triangulaire initial *********************!
  !
  !     sommets0(i,:) = numeros des sommets de la maille i dans le maillage initial (avant subdivision)
  !     arete0(i,j,k) = numero du k-eme sommet ajoute par subdivision de l'arete entre les
  !                     sommets initiaux i et j
  !
  ! Un "0" indique que cela se rapporte au maillage initial, et l'absence de "0" au maillage final
  !
  !****************
  ! Utilise: nbsom0,ncv0,sommets0,xs0,ys0,ndiv,nbsom,ncv,nusom,xs,ys,ptvois
  ! Affecte: xs, ys, nusom
  !***********************************!
    INTEGER :: i, j, k, is, icv, i0, i1, i2
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:)::  arete0
    INTEGER, ALLOCATABLE, DIMENSION(:) :: lim, lip
    DOUBLE PRECISION ::  a0, a1, a2

    ALLOCATE(arete0(nbsom0,nbsom0,ndiv+1),lim(ndiv+1),lip(ndiv+1))

    !- initialisation par les anciens sommets -!
    is = nbsom0
    icv = ncv0

    ! initialisation du tableau portant les numeros des sommets crees sur les aretes !
    do i=1, nbsom0
       do j=1, nbsom0
          do k=1, ndiv+1
             arete0(i,j,k) =0
          end do
       end do
    end do

    ! changement de donnees dans lecture
    ALLOCATE(xs(nbsom),ys(nbsom))
    do i=1, nbsom0
       xs(i)=xs0(i)
       ys(i)=ys0(i)
    end do

    ! creation des nouveaux sommets sur les aretes principales par subdivision

    do i=1, ncv0

       do j = 1, 3

          k = j+1
          if (k>3) k=1

          i0 = sommets0(i,j)
          i1 = sommets0(i,k)

          ! ici on recherche si l'arete i0 i1 a deja ete rencontree

          if (arete0(i0,i1,1)==0) then
             ! l'arete est nouvelle

             arete0(i0,i1,1)      = i0
             arete0(i0,i1,ndiv+1) = i1
             arete0(i1,i0,1)      = i1
             arete0(i1,i0,ndiv+1) = i0

             do k =1, ndiv-1

                a1= dble(k)/ dble(ndiv)
                a2=1. - a1

                is = is+1
                xs(is) = a1*xs0(i0)+a2*xs0(i1)
                ys(is) = a1*ys0(i0)+a2*ys0(i1)
                arete0(i1,i0,k+1) = is
                arete0(i0,i1,ndiv-k+1) = is

             end do

          end if

       end do

    END DO

    ! creation des nouveaux triangles par subdivision
    ALLOCATE(nusom(nptvois),ptvois(ncv+1))
    do k=1,ncv+1
       ptvois(k) = 3*k-3
    END DO


    do i =1, ncv0

       i0 = sommets0(i,1)
       i1 = sommets0(i,2)
       i2 = sommets0(i,3)


       ! parcours des bandes de triangles  entre deux lignes paralleles

       ! le tableau lim contient les numeros des sommets de la ligne du bas
       ! le tableau lip contient les numeros des sommets de la ligne du haut

       ! initialisation de lim par les sommets de l'arete initiale i2 - i0

       do j = 1,ndiv+1
          lim(j) = arete0(i2,i1,j)
       end do

       do j =1, ndiv-1

          ! initialisation des deux extremites de lip

          lip(1) = arete0(i2,i0,j+1)
          lip(ndiv-j+1) = arete0(i1,i0,j+1)


          a0=dble(j)/dble(ndiv)

          do k = 1,  ndiv -j -1
             a1=dble(k)/dble(ndiv)
             a2=1.- a0 - a1
             ! nouveau sommet interne
             is=is+1
             xs(is) = a0*xs0(i0) + a1*xs0(i1) + a2*xs0(i2)
             ys(is) = a0*ys0(i0) + a1*ys0(i1) + a2*ys0(i2)
             lip(k+1) = is
          end do


          ! creation des triangles pointe en haut
          do k=1, ndiv-j +1
             icv = icv+1
             nusom(3*icv-2) = lim(k)
             nusom(3*icv-1) = lip(k)
             nusom(3*icv) = lim(k+1)
          END DO

          ! creation des triangles pointe en bas

          do k = 1, ndiv-j
             icv=icv+1
             nusom(3*icv-2) = lip(k)
             nusom(3*icv-1) = lip(k+1)
             nusom(3*icv) = lim(k+1)
          END DO

          ! mise de lip dans lim (car la ligne du haut devient la ligne du bas de la bande suivante)

          do k =1, ndiv - j +1
             lim(k)=lip(k)
          END DO

       end do

       ! je redonne au dernier triangle le numero initial du gros triangle
       nusom(3*i-2) = lim(1)
       nusom(3*i-1) = i0
       nusom(3*i) = lim(2)
    END DO

    if (icv/=ncv) then
       print *," la division a donne ",icv," triangles au lieu de ",ncv
       stop
    ENDIF

    ! voila le nouveau nombre de sommets et de triangles

    nbsom = is
    ncv = icv


  !!$    do i=1, nbsom
  !!$     print*, " i",i," xs", xs(i)," ys", ys(i)
  !!$    end do
  !!$    do i=1, ncv
  !!$        print *,"triangle:", i
  !!$       do j=1, 3
  !!$              print *, j,"   nusom: ",nusom(i,j)
  !!$       end do
  !!$    end do
  !!$   print*, "fin division"

    DEALLOCATE(arete0,lim,lip)


  END SUBROUTINE divtri

!---------------------------------------------------------------!
!******************** calxcvtri ********************************!

  SUBROUTINE calxcvtri
  !*********** Calcul des centres des triangles ***********************
  ! Calcule les coordonnees xcv et ycv des centres des triangles
  !   d'un maillage.
  !   Voir le sens des variables dans la routine "divtri".
  !
  !********
  ! Utilise: ncv,nusom,xs,ys,xcv,ycv
  ! Affecte: xcv, ycv
  !*******************************************************
    DOUBLE PRECISION :: a11, a12, b1, a21, a22, b2, deter, aire
    INTEGER :: i,i0,i1,i2

    ALLOCATE(xcv(ncv),ycv(ncv))

    DO i=1, ncv

      ! Version "centre du cercle circonscrit" (i.e. intersection
      ! des mediatrices). Si tous les angles sont aigus, ce point est
      ! dans le triangle et rend le maillage admissible 
      ! pour les VF standards.
      ! Pourquoi les calculs suivants marchent: avec a,b,c les sommets
      ! du triangle, on cherche m tel que |a-m|^2=|b-m|^2=|c-m|^2
      ! ce qui mene a |b|^2-|a|^2=2(b-a).m et |c|^2-|a|^2=2(c-a).m, soit B=2Am
      ! avec A=(a_ij) et B=(b1,b2) ci-dessous.
      i0 = nusom (3*i-2)
      i1 = nusom (3*i-1)
      i2 = nusom (3*i)
      a11 = xs(i1) - xs(i0)
      a12 = ys(i1) - ys(i0)
      b1 = a11*(xs(i1)+xs(i0))+a12*(ys(i1)+ys(i0))
      a21 = xs(i2) - xs(i0)
      a22 = ys(i2) - ys(i0)
      b2 = a21*(xs(i2)+xs(i0))+a22*(ys(i2)+ys(i0))
      ! Resolution de 2Am=B
      aire = (a11*a22-a12*a21)
      deter = 0.5D0/aire
      xcv(i) = deter*(b1*a22-b2*a12)
      ycv(i) = deter*(a11*b2-a21*b1)

      ! Version "centre de gravite"
      ! xcv(i) = (xs(i0) +xs(i1) +xs(i2) )/3.
      ! ycv(i) = (ys(i0) +ys(i1) +ys(i2) )/3.

      IF (aire<0) THEN
        print *, "aire < 0 ",i
      ENDIF

    ENDDO

  END SUBROUTINE calxcvtri

!--------------------------------------------------------------------------------------------------------------!
!**************************** QUADRANGLES *********************************************************************!

  SUBROUTINE divquad
  !************ Division du maillage quandrangulaire initial
  !
  !     sommets0(i,:) = numeros des sommets de la maille i dans le maillage initial (avant subdivision)
  !     arete0(i,j,k) = numero du k-eme sommet ajoute par subdivision de l'arete entre les 
  !                     sommets initiaux i et j
  !
  ! Un "0" indique que cela se rapporte au maillage initial, et l'absence de "0" au maillage final
  !
  !****************
  ! Utilise: nbsom0,ncv0,sommets0,xs0,ys0,ndiv,nbsom,ncv,nusom,xs,ys,ptvois
  ! Affecte: xs, ys, nusom
  !***********************************************************************!
    INTEGER :: i, j, k, i0, i1, i2, i3, is, icv
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) ::  arete0
    INTEGER, ALLOCATABLE, DIMENSION(:) :: lim, lip
    DOUBLE PRECISION :: a1, a2

    ALLOCATE(arete0(nbsom0,nbsom0,ndiv+1),lim(ndiv+1),lip(ndiv+1))

    ! initialisation du tableau portant les numeros des sommets crees sur les aretes
    do i=1, nbsom0
       do j=1, nbsom0
          do k=1, ndiv+1
             arete0(i,j,k) =0
          end do
       end do
    end do

    ! on enregistre les sommets originels
    ALLOCATE(xs(nbsom),ys(nbsom))
    do i=1, nbsom0
       xs(i)=xs0(i)
       ys(i)=ys0(i)
    end do

    ! Creation des nouveaux sommets sur les aretes principales par subdivision
    ! On demarre a l'indice correspondant au dernier sommet initial
    is = nbsom0
    do i=1, ncv0
       do j = 1, 4
          k = j+1
          if (k>4) k=1

          i0 = sommets0(i,j)
          i1 = sommets0(i,k)

          ! ici on recherche si l'arete [i0 i1] a deja ete rencontree
          if (arete0(i0,i1,1)==0) then
             ! l'arete est nouvelle

             arete0(i0,i1,1)      = i0
             arete0(i0,i1,ndiv+1) = i1
             arete0(i1,i0,1)      = i1
             arete0(i1,i0,ndiv+1) = i0

             do k =1, ndiv-1
                ! on ajoute un nouveau sommet
                a1= dble(k)/ dble(ndiv)
                a2= 1D0 - a1

                is = is+1
                xs(is) = a1*xs0(i0)+a2*xs0(i1)
                ys(is) = a1*ys0(i0)+a2*ys0(i1)
                arete0(i1,i0,k+1) = is
                arete0(i0,i1,ndiv-k+1) = is
             end do
          end if
       end do
    end do

    ! creation des nouveaux rectangles par subdivision
    ALLOCATE(nusom(nptvois),ptvois(ncv+1))
    do k=1,ncv+1
       ptvois(k) = 4*k-4
    enddo

    ! On demarre avec l'indice du dernier rectangle du maillage initial
    icv = ncv0

    do i =1, ncv0

       i0 = sommets0(i,1)
       i1 = sommets0(i,2)
       i2 = sommets0(i,3)
       i3 = sommets0(i,4)

       ! parcours des bandes de rectangles  entre deux lignes paralleles
       ! le tableau lim contient les numeros des sommets de la ligne du bas
       ! le tableau lip contient les numeros des sommets de la ligne du haut

       ! initialisation de lim par les sommets de l'arete initiale i1 - i0

       do j = 1,ndiv+1
          lim(j) = arete0(i1,i0,j)
       end do

       do j =1, ndiv-1

          ! initialisation des deux extremites de lip

          lip(1) = arete0(i1,i2,j+1)
          lip(ndiv+1) = arete0(i0,i3,j+1)

          do k = 1,  ndiv -1
             a1=dble(k)/dble(ndiv)
             is=is+1
             lip(k+1) = is
             ! nouveau sommet interne: en notant A=lip(1) et B=lip(ndiv+1)
             ! les extremites de lip, ce sommet est en A+(k/ndiv).AB
             xs(is) = xs(lip(1))+a1*(xs(lip(ndiv+1))-xs(lip(1)))
             ys(is) = ys(lip(1))+a1*(ys(lip(ndiv+1))-ys(lip(1)))
          end do

          ! creation des rectangles
          do k=1, ndiv
             icv = icv+1
             nusom(4*icv-3) = lim(k+1)
             nusom(4*icv-2) = lim(k)
             nusom(4*icv-1) = lip(k)
             nusom(4*icv) = lip(k+1)
          end do


          ! mise de lip dans lim (car la ligne du haut devient la ligne du bas de la bande suivante)

          do k =1, ndiv+1
             lim(k)=lip(k)
          end do

       end do

       ! la derniere ligne de rectangles

       do j = 1,ndiv+1
          lip(j) = arete0(i2,i3,j)
       end do

       do k =1, ndiv-1
          ! creation des rectangles
          icv = icv+1
          nusom(4*icv-3) = lim(k+1)
          nusom(4*icv-2) = lim(k)
          nusom(4*icv-1) = lip(k)
          nusom(4*icv) = lip(k+1)
       end do

       ! je redonne au dernier rectangle le numero initial du gros triangle
       nusom(4*i-3) = lim(ndiv+1)
       nusom(4*i-2) = lim(ndiv)
       nusom(4*i-1) = lip(ndiv)
       nusom(4*i) = lip(ndiv+1)
    end do

    if (icv/=ncv) then
       print *," la division a donne ",icv," rectangles au lieu de ",ncv
       stop
    endif

    ! voila le nouveau nombre de sommets et de rectangles
    nbsom = is
    ncv = icv

    !   print*, "fin division"
    deallocate(arete0,lim,lip)

  END SUBROUTINE divquad

!------------------------------------------------------------------------!
!****************************calcxcvquad*************************************!

  SUBROUTINE calxcvquad
  !******************************************!
  ! Calcul des centres des quadrilateres
  !
  !*****
  ! Utilise: ncv,nusom,xs,ys,xcv,ycv
  ! Affecte: xcv, ycv
  !******************************************!
    INTEGER :: i, i0, i1, i2, i3

    ALLOCATE(xcv(ncv),ycv(ncv))

    DO i=1, ncv
       i0 = nusom(4*i-3)
       i1 = nusom(4*i-2)
       i2 = nusom(4*i-1)
       i3 = nusom(4*i)
       xcv(i) = 0.25D0*(xs(i0)+xs(i1)+xs(i2)+xs(i3))
       ycv(i) = 0.25D0*(ys(i0)+ys(i1)+ys(i2)+ys(i3))

  !       ! Barycentre decale 1
  !        xcv(i) = (5D0*xs(i0)+xs(i1)+xs(i2)+5D0*xs(i3))/12D0
  !        ycv(i) = (5D0*ys(i0)+ys(i1)+ys(i2)+5D0*ys(i3))/12D0
  !       ! Barycentre decale 2
  !        xcv(i) = (5D0*xs(i0)+xs(i1)+xs(i2)+5D0*xs(i3))/12D0
  !        ycv(i) = (5D0*ys(i0)+ys(i1)+ys(i2)+5D0*ys(i3))/12D0
  !				IF (xcv(i)>.5D0) THEN
  !	        xcv(i) = (xs(i0)+5D0*xs(i1)+5D0*xs(i2)+xs(i3))/12D0
  !	        ycv(i) = (ys(i0)+5D0*ys(i1)+5D0*ys(i2)+ys(i3))/12D0
  !				ENDIF
          ! Barycentre decale 3
  !        xcv(i) = (1/8D0)*xs(i0)+(1/8D0)*xs(i1)+(1/8D0+1/2D0)*xs(i2)+(1/8D0)*xs(i3)
  !        ycv(i) = (1/8D0)*ys(i0)+(1/8D0)*ys(i1)+(1/8D0+1/2D0)*ys(i2)+(1/8D0)*ys(i3)
          ! Barycentre decale 4
  !				IF ( MOD(i,2) == 0) THEN
  !					xcv(i) = (4D0*xs(i0)+xs(i1)+xs(i2)+xs(i3))/7D0
  !	        ycv(i) = (4D0*ys(i0)+ys(i1)+ys(i2)+ys(i3))/7D0
  !				ELSE
  !					xcv(i) = (xs(i0)+xs(i1)+4D0*xs(i2)+xs(i3))/7D0
  !	        ycv(i) = (ys(i0)+ys(i1)+4D0*ys(i2)+ys(i3))/7D0		
  !				ENDIF
      ENDDO

  END SUBROUTINE calxcvquad


!------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------!
!******************* DEUX DIVISIONS DIFFERENTES *************************************************************!

!-------------------------------------------------------------------!
!************** pointeurs_initial **********************************!

  SUBROUTINE pointeurs_initial
  !******************************************************************!
  ! Transforme les donnees na0, sommets0 du maillage lu dans le fichier en forme
  ! de pointeurs standards ptvois0, nusom0, nuvois0 (pour pouvoir etre utilisees dans
  ! divheterogene).
  ! On fait ici des manips un peu bourrines pour calculer ces pointeurs
  ! car, de toutes facons, sommets0 vient d'un fichier lu et ne peut pas
  ! etre monstrueusement gros. Ca evite de faire appel a la procedure
  ! calcul_pointeurs ci-dessous qui est plus lourde que ce dont on a besoin,
  ! et qui demanderait a davantage localiser ses arguments.
  !
  !******************************************!
  ! Utilise: ncv0, nbsom0, na0, sommets0
  ! Affecte: ptvois0, nusom0, nuvois0
  !*******************************************************************!
    INTEGER :: i, j, k, l, j0, j1, l0, trouve

    ! Creation de ptvois0
    ALLOCATE(ptvois0(ncv0+1))
    ptvois0(1) = 0
    DO i=2, ncv0+1
      ptvois0(i) = ptvois0(i-1) + na0(i-1)
    ENDDO

    ! Creation de nusom0
    ALLOCATE(nusom0(ptvois0(ncv0+1)))
    DO i=1, ncv0
      DO j=ptvois0(i)+1, ptvois0(i+1)
        nusom0(j) = sommets0(i,j-ptvois0(i))
      ENDDO
    ENDDO

    ! Creation de nuvois0
    ALLOCATE(nuvois0(ptvois0(ncv0+1)))
    DO i=1, ncv0
      ! Pour chaque arete de la maille i, on reparcours les autres mailles
      ! pour voir laquelle trouve les sommets dans l'ordre inverse. Si on n'en
      ! trouve pas, nuvois0 sera nul (arete de bord)
      j0 = nusom0(ptvois0(i+1))
      DO j=ptvois0(i)+1, ptvois0(i+1)
        j1 = nusom0(j)

        k = 0
        trouve = 0
        DO WHILE (trouve == 0 .AND. k<ncv0)
          k = k+1
          l0 = nusom0(ptvois0(k+1))
          l = ptvois0(k)+1
          DO WHILE (l<=ptvois0(k+1) .AND. (l0 /= j1 .OR. nusom0(l) /= j0))
            l0 = nusom0(l)
            l = l+1
          ENDDO
          IF (l<= ptvois0(k+1) .AND. l0 == j1 .AND. nusom0(l) == j0) trouve = 1
        ENDDO

        IF (trouve == 1) THEN
          nuvois0(j) = k
        ELSE
          nuvois0(j) = 0
        ENDIF

        j0 = j1
      ENDDO

    ENDDO


  END SUBROUTINE pointeurs_initial


!-------------------------------------------------------------------!
!************** divheterogene **************************************!

  SUBROUTINE divheterogene
  !******************************************************************!
  ! Divise des mailles de maniere heterogene: ndiv1 dans une zone du domaine,
  ! ndiv2 dans l'autre (selon typediv). Les mailles non quadrangulaires ou triangulaires
  ! ne sont pas divisees (mais voir le commentaire un peu plus bas)
  !
  !***************
  ! Utilise: ncv0, nbsom0, xs0, ys0, ndiv1, ndiv2, ptvois0, nusom0, nuvois0
  ! Affecte: nptvois, ptvois, ncv, nbsom, nusom, xs, ys
  !******************************************************************!
    INTEGER :: i, ii, j, k, l, j0, j1, i0, i1, i2, i3, k0, pos, pos1, pos2, nd1, nd2, nd, is, icv
    INTEGER :: ncv_temp, nbsom_temp, nptvois_temp
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ptvois_temp, nusom_temp, lim, lip
    DOUBLE PRECISION :: a1
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xs_temp, ys_temp, xcv0, ycv0
    INTEGER, ALLOCATABLE, DIMENSION(:) :: arete0_nbs
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: arete0_sub, arete0_som
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: arete0_posdiv


    !!!!!!
    ! Majorants du nombre de mailles, de sommets, et de la taille de nusom (nptvois)
    !
    ! Au pire, chaque maille est divisee en max des deux divisions
    ncv_temp = ncv0 * (max(ndiv1,ndiv2)) * (max(ndiv1,ndiv2))

    ! Il faut ensuite majorer la somme sur les mailles du nombre de sommets de chaque
    ! maille (pour la longueur de nusom).
    ! La plupart des mailles finales seront a l'interieur d'une maille initiale
    ! (divisee au pire par le max des deux divisions) et, dans ce cas, auront
    ! exactement le meme nombre de sommets que la maille initiale: on somme donc sur
    ! i les (ptvois0(i+1)-ptvois0(i))*max(divisions)*max(divisions), ce qui donne
    ! ptvois0(ncv0+1)*max(divisions)*max(divisions)
    ! Se rajoutent a ca les aretes que l'on cree par subdivision des aretes
    ! initiales, chacune de ces aretes subdivisees rajoutant aux mailles qui la bordent
    ! autant de sommets que le nombre de sommets de l'autre cote: le rajout global
    ! pour ces aretes est donc ndiv1+ndiv2 (ndiv1 d'un cote, ndiv2 de l'autre). Il faut
    ! donc a nouveau sommer sur i la quantite (ptvois0(i+1)-ptvois0(i))*(ndiv1+ndiv2)
    ! (nombre d'aretes de la maille i fois le rajout precedent).
    nptvois_temp = ptvois0(ncv0+1)*(max(ndiv1,ndiv2))*(max(ndiv1,ndiv2)) + ptvois0(ncv0+1)*(ndiv1+ndiv2)

    ! Nombre max de sommets: comme si on avait divise chaque maille par le max des deux
    nbsom_temp = ncv0 * (max(ndiv1,ndiv2)+1) * (max(ndiv1,ndiv2)+1)

    ! On aura besoin des centres des mailles initiales
    ALLOCATE(xcv0(ncv0),ycv0(ncv0))
    DO i=1, ncv0
      xcv0(i) = 0D0
      ycv0(i) = 0D0
      DO j=ptvois0(i)+1, ptvois0(i+1)
        xcv0(i) = xcv0(i) + xs0(nusom0(j))
        ycv0(i) = ycv0(i) + ys0(nusom0(j))
      ENDDO
      xcv0(i) = xcv0(i) / dble(ptvois0(i+1)-ptvois0(i))
      ycv0(i) = ycv0(i) / dble(ptvois0(i+1)-ptvois0(i))
    ENDDO

    !!!!!! RECHERCHE DES SUBDIVISIONS DE CHAQUE ARETE DE BASE
    ! Les tableaux arete0 doivent etre lus avec le pointeur ptvois0: pour une maille i
    ! et j=ptvois0(i)+1, ptvois0(i+1), arete0(j,...) donne des infos sur l'arete correspondante
    ! de i.
    !
    ! Pour chaque arete j, on cherche d'abord quelles subdivisions
    ! elle subira; on note ces subdivisions dans arete0_sub(j,1) et arete0_sub(j,2)
    ALLOCATE(arete0_sub(ptvois0(ncv0+1),2))
    DO j=1, ptvois0(ncv0+1)
      arete0_sub(j,1)=0
      arete0_sub(j,2)=0
    ENDDO

    DO i=1, ncv0
      ! On ne divise que les quadrangles et triangles, pas les autres mailles, peu importe
      ! leur position dans le domaine (donc division a 1 pour elles).
      ! Si on souhaite diviser quand meme les aretes de ces mailles, commenter
      ! ce "IF" et le "ELSEIF" de la boucle sur j (reperant le nombre de divisions
      ! de la maille de l'autre cote de j) pour toujours mettre nd = niv_div,
      ! sauf si nuvois0(j)=0
      IF ( (ptvois0(i+1)-ptvois0(i) /= 4) .AND. (ptvois0(i+1)-ptvois0(i) /= 3) ) THEN
        nd1 = 1
      ELSE
        nd1 = niv_div(typediv,xcv0(i),ycv0(i),ndiv1,ndiv2)
      ENDIF

      DO j=ptvois0(i)+1, ptvois0(i+1)
        ! Subdivision de la maille de l'autre cote de j?
        IF (nuvois0(j) == 0) THEN
          nd2 = 1
        ELSEIF ( (ptvois0(nuvois0(j)+1)-ptvois0(nuvois0(j)) /= 4) &
                    .AND. (ptvois0(nuvois0(j)+1)-ptvois0(nuvois0(j)) /= 3) ) THEN
          nd2 = 1
        ELSE
          nd2 = niv_div(typediv,xcv0(nuvois0(j)),ycv0(nuvois0(j)),ndiv1,ndiv2)
        ENDIF

        ! Subdivisions subies par j
        arete0_sub(j,1) = nd1
        arete0_sub(j,2) = nd2
      ENDDO

    ENDDO


    !!!!!!! SUBDIVISION DES ARETES DE BASE
    ! On va ensuite subdiviser chaque arete selon ce qu'il doit etre.
    ! arete0_som(j,:) indique les numeros des sommets (que l'on cree simultanement)
    !   crees par subdivison de l'arete j, dans l'ordre correspondant a la maille i correspondante
    !   (on a j=ptvois0(i)+1, ptvois0(i+1), on va donc de nusom0(j-1) a nusom0(j))
    ! arete0_nbs(j) est le nombre de sommets finaux sur l'arete j
    ! arete0_posdiv(j,m,:) indique les m+1 indices dans arete0_som(j,:) correspondant
    !   aux sommets crees par division m (non-nul uniquement pour m=ndiv1, m=ndiv2 et m=1)
    ALLOCATE(arete0_nbs(ptvois0(ncv0+1)),arete0_som(ptvois0(ncv0+1),ndiv1+ndiv2))
    ALLOCATE(arete0_posdiv(ptvois0(ncv0+1),max(ndiv1,ndiv2),max(ndiv1,ndiv2)+1))
    DO j=1, ptvois0(ncv0+1)
      DO k=1,ndiv1+ndiv2
        arete0_som(j,k)=0
      ENDDO
      arete0_nbs(j)=0
      DO k=1,max(ndiv1,ndiv2)
        DO l=1,max(ndiv1,ndiv2)+1
          arete0_posdiv(j,k,l)=0
        ENDDO
      ENDDO
    ENDDO

    ! Initialisation des sommets par les sommets du maillage initial
    ALLOCATE(xs_temp(nbsom_temp),ys_temp(nbsom_temp))
    xs_temp(:) = 0D0
    ys_temp(:) = 0D0
    DO i=1, nbsom0
      xs_temp(i) = xs0(i)
      ys_temp(i) = ys0(i)
    ENDDO
    ! Indiquera le sommet en cours de creation
    is = nbsom0

    ! On divise les aretes en creant les sommets dessus. On ne traite que les
    ! aretes qui n'ont pas deja ete vues dans une autre maille (arete0_som affecte)
    DO i=1, ncv0
      j0 = nusom0(ptvois0(i+1))
      DO j=ptvois0(i)+1, ptvois0(i+1)
        j1 = nusom0(j)

        IF (arete0_som(j,1) == 0) THEN
          ! Les niveaux de subdivision de l'arete
          nd1 = arete0_sub(j,1)
          nd2 = arete0_sub(j,2)
          IF (nd1 == 0 .OR. nd2 == 0) THEN
            print *,"Erreur divisions arete"
            STOP
          ENDIF
          ! Les premiers termes
          arete0_som(j,1) = j0
          arete0_posdiv(j,nd1,1) = 1
          arete0_posdiv(j,nd2,1) = 1

          ! Positions en cours dans arete0_som, arete0_posdiv(m=nd1), arete0_posdiv(m=nd2)
          pos = 1
          pos1 = 1
          pos2 = 1
          DO WHILE (pos1<nd1 .OR. pos2<nd2)
            ! On se demande quel type de sommet (division nd1 ou nd2) apparait ensuite.
            ! C'est selon que pos1/nd1 est plus petit ou non que pos2/nd2 (la position
            ! du point est donne par ces quantites, cf a1)
            IF ( pos1*nd2 < pos2*nd1 ) THEN
              ! Sommet de type nd1
              ! Creation du sommet et le mettre dans arete0_som et arete0_posdiv
              pos = pos+1
              is = is+1
              arete0_som(j,pos) = is

              a1= dble(pos1)/ dble(nd1)
              xs_temp(is) = a1*xs0(j1)+(1D0-a1)*xs0(j0)
              ys_temp(is) = a1*ys0(j1)+(1D0-a1)*ys0(j0)

              pos1 = pos1 + 1
              arete0_posdiv(j,nd1,pos1) = pos

            ELSEIF ( pos1*nd2 == pos2*nd1 ) THEN
              ! Sommet commun nd1 et nd2
              pos = pos+1
              is = is+1
              arete0_som(j,pos) = is

              a1= dble(pos1)/ dble(nd1)
              xs_temp(is) = a1*xs0(j1)+(1D0-a1)*xs0(j0)
              ys_temp(is) = a1*ys0(j1)+(1D0-a1)*ys0(j0)

              pos1 = pos1 + 1
              arete0_posdiv(j,nd1,pos1) = pos
              pos2 = pos2 + 1
              arete0_posdiv(j,nd2,pos2) = pos

            ELSE
              ! Sommet de type nd2
              pos = pos+1
              is = is+1
              arete0_som(j,pos) = is

              a1= dble(pos2)/ dble(nd2)
              xs_temp(is) = a1*xs0(j1)+(1D0-a1)*xs0(j0)
              ys_temp(is) = a1*ys0(j1)+(1D0-a1)*ys0(j0)

              pos2 = pos2 + 1
              arete0_posdiv(j,nd2,pos2) = pos

            ENDIF
          ENDDO
          IF (pos>ndiv1+ndiv2 .OR. pos1/=nd1 .OR. pos2/=nd2 .OR. is>nbsom_temp) THEN
            print *,"Erreur position"
            STOP
          ENDIF

          ! Derniers sommets
          pos = pos+1
          arete0_nbs(j) = pos
          arete0_som(j,pos) = j1
          arete0_posdiv(j,nd1,nd1+1) = pos
          arete0_posdiv(j,nd2,nd2+1) = pos


          ! S'il y a une maille de l'autre cote, on indique les sommets qui viennent d'etre
          ! crees sur son arete.
          IF (nuvois0(j) /= 0) THEN
            ii = nuvois0(j)
            ! On cherche a quel endroit dans la maille ii l'arete j=j0j1 apparait (sous forme
            ! j1j0)
            k0 = nusom0(ptvois0(ii+1))
            k = ptvois0(ii)+1
            DO WHILE ( k<=ptvois0(ii+1) .AND. (k0/=j1 .OR. nusom0(k)/=j0) )
              k0 = nusom0(k)
              k=k+1
            ENDDO
            IF (k>ptvois0(ii+1) .OR. k0/=j1 .OR. nusom0(k)/=j0) THEN
              print *,"Erreur recherche position dans voisin"
              STOP
            ENDIF

            arete0_nbs(k) = arete0_nbs(j)
            DO l=1, arete0_nbs(j)
              arete0_som(k,l) = arete0_som(j,arete0_nbs(j)+1-l)
            ENDDO
            DO l=1, nd1+1
              arete0_posdiv(k,nd1,l) = arete0_nbs(j) + 1 - arete0_posdiv(j,nd1,nd1+2-l)
            ENDDO
            DO l=1, nd2+1
              arete0_posdiv(k,nd2,l) = arete0_nbs(j) + 1 -  arete0_posdiv(j,nd2,nd2+2-l)
            ENDDO
          ENDIF

        ENDIF

        j0 = j1
      ENDDO
    ENDDO



  ! DO i=1, ncv0
  !   print *,"*****************",i
  !   print *,(nusom0(j),j=ptvois0(i)+1,ptvois0(i+1))
  !   print *,(nuvois0(j),j=ptvois0(i)+1,ptvois0(i+1))
  !   print *,(arete0_sub(j,1),j=ptvois0(i)+1,ptvois0(i+1))
  !   print *,(arete0_sub(j,2),j=ptvois0(i)+1,ptvois0(i+1))
  !   DO j=ptvois0(i)+1, ptvois0(i+1)
  !     print *,"-----------------------------",nusom0(j)
  !     print *,arete0_nbs(j)
  !     print *,arete0_som(j,:)
  !     print *,arete0_posdiv(j,arete0_sub(j,1),:)
  !     print *,arete0_posdiv(j,arete0_sub(j,2),:)
  !   ENDDO
  ! ENDDO
  ! 
  ! DO i=1,is
  !   print *,i,xs_temp(i),ys_temp(i)
  ! ENDDO


    !!!! CREATION DES MAILLES, ET DES SOMMETS INTERIEURS
    !
    ! On cree les mailles en divisant les mailles de base
    ALLOCATE(nusom_temp(nptvois_temp),ptvois_temp(ncv_temp+1))
    nusom_temp(:) = 0
    ptvois_temp(:) = 0

    ! Indice de la nouvelle maille en cours (is est le sommet en cours)
    icv = 0
    ! Position dans nusom
    pos = 0

    DO i=1, ncv0
      ! On ne divise que les mailles quadrangulaires ou triangulaires. Les autres, on
      ! ne fait que rajouter les sommets necessaires au bord
      IF ( (ptvois0(i+1)-ptvois0(i) /= 4)  .AND. (ptvois0(i+1)-ptvois0(i) /= 3) ) THEN
        icv = icv+1
        ptvois_temp(icv) = pos

        ! On rajoute les sommets crees sur les aretes initiales bordant la maille
        ! en question

        DO j=ptvois0(i)+1, ptvois0(i+1)
          ! On ne met pas le dernier de chaque arete, pour eviter de le compter
          ! deux fois
          DO k=1, arete0_nbs(j)-1
            pos = pos + 1
            nusom_temp(pos) = arete0_som(j,k)
          ENDDO

        ENDDO

      ELSEIF ( ptvois0(i+1)-ptvois0(i) == 4 ) THEN
        ! Maille quandrangulaire
        nd = niv_div(typediv,xcv0(i),ycv0(i),ndiv1,ndiv2)
        ALLOCATE(lim(nd+1),lip(nd+1))

        ! Sommets de la maille, partant d'en bas a gauche et tournant dans
        ! le sens trigo. On imagine un carre avec i0 en bas a gauche, puis les
        ! autres sommets en sens trigo
        i0 = ptvois0(i)+1
        i1 = ptvois0(i)+2
        i2 = ptvois0(i)+3
        i3 = ptvois0(i)+4

        ! Initialisation de lim par la bande du bas (sommets venant de la division
        ! nd de i0i1). On parcours cette bande de gauche a droite (sens naturel)
        DO k=1, nd+1
          l = arete0_posdiv(i1,nd,k)
          lim(k) = arete0_som(i1,l)
        ENDDO

        ! On monte ensuite etage par etage le long des aretes i1i2 et i0i3
        ! lim est le plancher de l'etage (gauche a droite), lip le plafond (droite a gauche)
        DO j=1, nd
          ! Les termes extremes de lip. Facile pour lip(1) car i1i2 est dans le bon sens,
          ! mais i0i3 est dans le mauvais sens et lip(nd+1) doit etre le j+1-eme sommet
          ! correspondant a la division nd dans le sens i0i3.
          l = arete0_posdiv(i2,nd,j+1)
          lip(1) = arete0_som(i2,l)
          l = arete0_posdiv(i0,nd,nd+2-(j+1))
          lip(nd+1) = arete0_som(i0,l)

          ! Creation des sommets internes sur lip, uniquement si j<nd (i.e. lip n'est
          ! pas sur i2i3). Sinon on charge lip avec les sommets de i2i3 venant de la division n
          IF (j<nd) THEN
            DO k=1, nd-1
              a1=dble(k)/dble(nd)
              is=is+1
              lip(k+1) = is
              xs_temp(is) = xs_temp(lip(1))+a1*(xs_temp(lip(nd+1))-xs_temp(lip(1)))
              ys_temp(is) = ys_temp(lip(1))+a1*(ys_temp(lip(nd+1))-ys_temp(lip(1)))
            ENDDO
          ELSE
            DO k=2, nd
              l = arete0_posdiv(i3,nd,k)
              lip(k) = arete0_som(i3,l)
            ENDDO
          ENDIF
          ! Creation des mailles, de sommets de base lim(k),lim(k+1),lip(nd+1-k),lip(nd+2-k)
          ! Il faut faire attention aux cas d'une maille
          ! qui borde une arete de base: j=1 (bas), j=nd (haut) et dans
          ! chaque bande les 1er (droite) et nd-eme (gauche) terme
          DO k=1, nd
            icv = icv + 1
            ptvois_temp(icv) = pos

            ! Premier sommet en bas a gauche
            pos = pos + 1
            nusom_temp(pos) = lim(k)

            ! Sommets en bas, de lim(k) a lim(k+1) (droite a gauche). Si j=1, on est
            ! sur l'arete i0i1 et il faut donc compter les eventuels sommets
            ! supplementaires pas vus dans lim (venant de l'autre division de l'arete i0i1 entre
            ! les sommets de position k et k+1 pour la division nd)
            IF (j == 1) THEN
              DO l=arete0_posdiv(i1,nd,k)+1,arete0_posdiv(i1,nd,k+1)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i1,l)
              ENDDO
            ENDIF
            ! Dernier sommet du bas
            pos = pos + 1
            nusom_temp(pos) = lim(k+1)

            ! Sommets a droite, de lim(k+1) a lip(nd+1-k) (de bas en haut). Si k=nd, on est
            ! sur l'arete i1i2 (compter les sommets en plus)
            IF (k == nd) THEN
              DO l=arete0_posdiv(i2,nd,j)+1,arete0_posdiv(i2,nd,j+1)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i2,l)
              ENDDO
            ENDIF
            ! Dernier sommet a droite
            pos = pos + 1
            nusom_temp(pos) = lip(nd+1-k)

            ! Sommets en haut, de lip(nd+1-k) a lip(nd+2-k) (droite a gauche). Si j=nd, on est
            ! sur l'arete i2i3 (compter les sommets en plus)
            IF (j == nd) THEN
              DO l=arete0_posdiv(i3,nd,nd+1-k)+1,arete0_posdiv(i3,nd,nd+2-k)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i3,l)
              ENDDO
            ENDIF
            ! Dernier sommet en haut
            pos = pos + 1
            nusom_temp(pos) = lip(nd+2-k)

            ! Les eventuels sommets ente lip(nd+2-k) et lim(k), si on est sur une arete
            ! de base i3i0 a gauche (i.e. k=1).
            ! Rappel: lip(nd+2-k)=lip(nd+1) est sur l'arete i0 (=i3i0) le sommet en
            ! position arete0_posdiv(i0,nd,nd+2-(j+1)), et de lip(nd+1) a lim(1)
            ! on parcours l'arete i3i0, autrement dit dans le bon sens.
            IF (k == 1) THEN
              DO l=arete0_posdiv(i0,nd,nd+2-(j+1))+1,arete0_posdiv(i0,nd,nd+2-j)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i0,l)
              ENDDO
            ENDIF

          ENDDO

          ! On met lip dans lim (la ligne du haut devient celle du bas a l'etage
          ! suivant), mais attention a l'ordre: lip etait de droite a gauche et on veut lim
          ! de gauche a droite.
          DO k=1, nd+1
            lim(k) = lip(nd+2-k)
          ENDDO
        ENDDO

        DEALLOCATE(lim,lip)

      ELSEIF ( ptvois0(i+1)-ptvois0(i) == 3 ) THEN
        ! Maille triangulaire
        nd = niv_div(typediv,xcv0(i),ycv0(i),ndiv1,ndiv2)
        ALLOCATE(lim(nd+1),lip(nd))

        ! Sommets de la maille, partant d'en bas a gauche et tournant dans
        ! le sens trigo. On imagine un triangle pointe en haut avec i0 en bas a gauche, puis les
        ! autres sommets en sens trigo
        i0 = ptvois0(i)+1
        i1 = ptvois0(i)+2
        i2 = ptvois0(i)+3

        ! Initialisation de lim par la bande du bas (sommets venant de la division
        ! nd de i0i1). On parcours cette bande de gauche a droite (sens naturel)
        DO k=1, nd+1
          l = arete0_posdiv(i1,nd,k)
          lim(k) = arete0_som(i1,l)
        ENDDO

        ! On monte ensuite etage par etage le long des aretes i1i2 et i0i3
        ! lim est le plancher de l'etage (gauche a droite), lip le plafond (droite a gauche)
        ! Pour un etage j donne, lim aura donc une longueur nd+2-j et lip un de moins (nd+1-j)
        DO j=1, nd
          ! Les termes extremes de lip. Facile pour lip(1) car i1i2 est dans le bon sens,
          ! mais i0i2 est dans le mauvais sens et lip(nd+1-j) doit etre le j+1-eme sommet
          ! correspondant a la division nd dans le sens i0i2.
          l = arete0_posdiv(i2,nd,j+1)
          lip(1) = arete0_som(i2,l)
          l = arete0_posdiv(i0,nd,nd+2-(j+1))
          lip(nd+1-j) = arete0_som(i0,l)

          ! Creation des sommets internes sur lip. A l'etape j, l'etage plafond lip
          ! est divise en nd-j morceaux (et les premier 1 et dernier nd+1-j indices de lip
          ! ont deja ete rentres)
          DO k=1, nd-j-1
            a1=dble(k)/dble(nd-j)
            is=is+1
            lip(k+1) = is
            xs_temp(is) = xs_temp(lip(1))+a1*(xs_temp(lip(nd+1-j))-xs_temp(lip(1)))
            ys_temp(is) = ys_temp(lip(1))+a1*(ys_temp(lip(nd+1-j))-ys_temp(lip(1)))
          ENDDO

          ! Creation des mailles.
          ! La longueur de lim etant nd+2-j, k va de 1 a nd+1-j.
          ! On cree d'abord un triangle pointe en haut lim(k),lim(k+1),lip(nd+2-j-k), puis
          ! un triangle pointe en bas lim(k+1),lip(nd+2-j-k-1),lip(nd+2-j-k) [uniquement si
          ! k n'est pas au bout de la boucle, i.e. k<nd+1-j].
          ! Pour le triangle pointe en haut, il faut faire attention aux situations ou
          ! l'une de ses aretes est sur une arete initiale i0i1 (j=1), i1i2 (k=nd+1-j)
          ! ou i2i0 (k=1).
          DO k=1, nd+1-j

            ! TRIANGLE POINTE EN HAUT
            icv = icv + 1
            ptvois_temp(icv) = pos

            ! Premier sommet en bas a gauche
            pos = pos + 1
            nusom_temp(pos) = lim(k)

            ! Sommets en bas, de lim(k) a lim(k+1) (droite a gauche). Si j=1, on est
            ! sur l'arete i0i1 et il faut donc compter les eventuels sommets
            ! supplementaires pas vus dans lim (venant de l'autre division de l'arete i0i1 entre
            ! les sommets de position k et k+1 pour la division nd)
            IF (j == 1) THEN
              DO l=arete0_posdiv(i1,nd,k)+1,arete0_posdiv(i1,nd,k+1)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i1,l)
              ENDDO
            ENDIF
            ! Dernier sommet du bas
            pos = pos + 1
            nusom_temp(pos) = lim(k+1)

            ! Sommets a droite, de lim(k+1) a lip(nd+2-j-k) (de bas en haut). Si k=nd+1-j, on est
            ! sur l'arete i1i2 (compter les sommets en plus)
            IF (k == nd+1-j) THEN
              DO l=arete0_posdiv(i2,nd,j)+1,arete0_posdiv(i2,nd,j+1)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i2,l)
              ENDDO
            ENDIF
            ! Dernier sommet en haut
            pos = pos + 1
            nusom_temp(pos) = lip(nd+2-j-k)

            ! Eventuels sommets a gauche entre lip(nd+2-j-k) et lim(k), si on est sur l'arete
            ! i2i0 (i.e. k=1).
            ! Rappel: lip(nd+2-j-k)=lip(nd+1-j) est sur l'arete i0 (=i2i0) le sommet en
            ! position arete0_posdiv(i0,nd,nd+2-(j+1)), et de lip(nd+1-j) a lim(1)
            ! on parcours l'arete i2i0, autrement dit dans le bon sens.
            IF (k == 1) THEN
              DO l=arete0_posdiv(i0,nd,nd+2-(j+1))+1,arete0_posdiv(i0,nd,nd+2-j)-1
                pos = pos+1
                nusom_temp(pos) = arete0_som(i0,l)
              ENDDO
            ENDIF

            ! TRIANGLE POINTE EN HAUT, si k n'est pas au bout.
            IF (k<nd+1-j) THEN
              icv = icv + 1
              ptvois_temp(icv) = pos

              ! Triangle lim(k+1),lip(nd+2-j-k-1),lip(nd+2-j-k), aucun risque d'etre sur
              ! une arete initiale
              pos = pos + 1
              nusom_temp(pos) = lim(k+1)
              pos = pos + 1
              nusom_temp(pos) = lip(nd+2-j-k-1)
              pos = pos + 1
              nusom_temp(pos) = lip(nd+2-j-k)

            ENDIF

          ENDDO

          ! On met lip dans lim (la ligne du haut devient celle du bas a l'etage
          ! suivant), mais attention a l'ordre: lip etait de droite a gauche et on veut lim
          ! de gauche a droite.
          DO k=1, nd+1-j
            lim(k) = lip(nd+2-j-k)
          ENDDO
        ENDDO

        DEALLOCATE(lim,lip)

      ENDIF
    ENDDO


  ! DO i=1,is
  !   print *,i,xs_temp(i),ys_temp(i)
  ! ENDDO
  ! DO i=1, icv
  !   print *,"************",i
  !   print *,(nusom_temp(j),j=ptvois_temp(i)+1,ptvois_temp(i+1))
  ! ENDDO

    IF (icv>ncv_temp .OR. is>nbsom_temp .OR. pos>nptvois_temp) THEN
      print *,"Erreur evaluation initiale",icv,ncv_temp,is,nbsom_temp,pos,nptvois_temp
      STOP
    ENDIF

    ! Nombre de mailles, de sommets, et nptvois finaux
    ncv = icv
    nbsom = is
    nptvois = pos

    ! On charge dans les vecteurs definitifs
    ALLOCATE(ptvois(ncv+1),nusom(nptvois),xs(nbsom),ys(nbsom))
    DO i=1,ncv
      ptvois(i) = ptvois_temp(i)
    ENDDO
    ptvois(ncv+1) = nptvois
    DO i=1,nptvois
      nusom(i) = nusom_temp(i)
    ENDDO
    DO i=1,nbsom
      xs(i) = xs_temp(i)
      ys(i) = ys_temp(i)
    ENDDO


    DEALLOCATE(nusom_temp,ptvois_temp,xs_temp,ys_temp,arete0_som,arete0_nbs,arete0_posdiv)

    DEALLOCATE(xs0,ys0,ptvois0,nusom0,nuvois0)

  END SUBROUTINE divheterogene

!-------------------------------------------------------------------!
!************** niv_div ********************************************!

  FUNCTION niv_div(typediv,x,y,nd1,nd2)
    !********************************************************!
    ! Determine, en fonction de typediv, dans quelle zone (1 ou 2)
    ! le point x,y se trouve et renvoie donc le niveau de division
    ! nd1 ou nd2 que subira une maille ayant x,y comme centre
    !
    !********************************************************!
    INTEGER :: typediv, nd1, nd2
    DOUBLE PRECISION :: x, y

    INTEGER :: niv_div

    IF (typediv == 1) THEN
      ! Division selon x=0.5, nd1 a droite et nd2 a gauche

      IF (x <= 0.5D0) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSEIF (typediv == 2) THEN
      ! nd1 dans le quart inferieur gauche, nd2 ailleurs

      IF (x<=0.5D0 .AND. y<=0.5D0) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSEIF (typediv == 3) THEN
      ! nd1 dans le tiers superieur droit, nd2 ailleurs

      IF (x>=0.66D0 .AND. y>=0.66D0) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSEIF (typediv == 4) THEN
      ! nd1 dans le tiers inferieur, nd2 ailleurs

      IF (y<0.35D0) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSEIF (typediv == 5) THEN
      ! nd1 le long des bords, nd2 ailleurs

      IF (y<layer .OR. x<layer .OR. y>1.-layer .OR. x>1.-layer) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSEIF (typediv == 6) THEN
      ! nd1 le long des bords et du quart de domaine en haut a gauche, nd2 ailleurs

      IF (y<layer .OR. x<layer .OR. y>1.-layer .OR. x>1.-layer &
            .OR. (x<layer .AND. y > .5-layer) .OR. (y>1.-layer .AND. x < .5+layer) &
              .OR. (abs(x-.5)<layer .AND. y>.5-layer) .OR. (abs(y-.5)<layer .AND. x<.5+layer)) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSEIF (typediv == 7) THEN
      ! nd1 autour du point (0.25,1), nd2 ailleurs

      IF ( ABS(x-.25)<layer .AND. y>1-layer) THEN
        niv_div = nd1
      ELSE
        niv_div = nd2
      ENDIF

    ELSE
      print *,"Type de division inconnue",typediv
      STOP

    ENDIF


  END FUNCTION niv_div


!------------------------------------------------------------------------!
!****************************calcul_centres******************************!

  SUBROUTINE calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,typecentre)
  !******************************************!
  ! Calcul des centres de maille, a partir
  ! de ncv, nusom, xs, ys et ptvois.
  ! Remplit xcv et ycv, qui doivent etre alloues avant
  !
  !	typecentre indique le type de centres qu'on veut:
  !			1 : isobarycentre des sommets de la maille
  !			2 : centre de gravite du bord de la maille
  !			3 : centre des cercles circonscrits (uniquement triangles)
  !			4+ : manuel
  !
  !******************************************!
		INTEGER :: ncv, typecentre
		INTEGER, DIMENSION(:) :: ptvois, nusom
		DOUBLE PRECISION, DIMENSION(:) :: xs, ys, xcv, ycv
		
    INTEGER :: i, j, j0, j1, i0, i1, i2
    DOUBLE PRECISION :: coef, long
    DOUBLE PRECISION :: a11, a12, a21, a22, b1, b2, deter
    DOUBLE PRECISION :: coef1, coef2, coef3, coeft
		DOUBLE PRECISION :: xcv_temp, ycv_temp, airetri, airetot
		INTEGER, PARAMETER :: verbose=0
	
		IF (typecentre == 0) THEN
			! Vrais centres de gravite, calcules en coupant chaque
			! maille en triangles (centre des sommets-sommet1-sommet2)
			!	et en faisant la moyenne des centres de chaque triangle ponderee par leur aire
			IF (verbose==1) print *,"x_K = vrais centres de gravite\n"
			DO i=1,ncv
				xcv(i) = 0D0
				ycv(i) = 0D0
				
				! Calcul d'un point a l'interieur de la maille (centre des sommets ici)
				xcv_temp = 0D0
				ycv_temp = 0D0
				DO j=ptvois(i)+1, ptvois(i+1)
				  xcv_temp = xcv_temp + xs(nusom(j))
				  ycv_temp = ycv_temp + ys(nusom(j))
				ENDDO
				coef = 1D0/(dble(ptvois(i+1)-ptvois(i)))
				xcv_temp = xcv_temp * coef
				ycv_temp = ycv_temp * coef
				! Ici, juste un des sommets de la maille
  !				xcv_temp = xs(nusom(ptvois(i+1)))
  !				ycv_temp = ys(nusom(ptvois(i+1)))

				airetot = 0D0

				j0 = nusom(ptvois(i+1))
				DO j=ptvois(i)+1,ptvois(i+1)
					j1 = nusom(j)
					! Calcul de l'aire du triangle cv_temp-j0-j1 (moitie de l'aire du parallelogramme)
					airetri = (xs(j0)-xcv_temp)*(ys(j1)-ycv_temp) - (ys(j0)-ycv_temp)*(xs(j1)-xcv_temp)
					airetri = 0.5D0 * airetri
					IF (airetri<0D0) THEN
						print *,"Erreur aire triangle dans calcul centres",airetri,' ',i
						STOP
					ENDIF
					airetot = airetot+airetri

					! Ajout du centre de gravite du triangle pondere par son aire
					xcv(i) = xcv(i) + airetri * (xcv_temp+xs(j0)+xs(j1))/3D0
					ycv(i) = ycv(i) + airetri * (ycv_temp+ys(j0)+ys(j1))/3D0

					j0 = j1
				ENDDO
				! On pondere par l'aire totale
				xcv(i) = xcv(i) / airetot
				ycv(i) = ycv(i) / airetot
												
			ENDDO

		ELSEIF (typecentre == 1) THEN
			! Version "centre de gravite des sommets de la maille"
			IF (verbose==1) print *,"x_K = centres de gravite des sommets\n"
			DO i=1, ncv
				xcv(i) = 0D0
				ycv(i) = 0D0

				DO j=ptvois(i)+1, ptvois(i+1)
				  xcv(i) = xcv(i) + xs(nusom(j))
				  ycv(i) = ycv(i) + ys(nusom(j))
				ENDDO
				coef = 1D0/(dble(ptvois(i+1)-ptvois(i)))
				xcv(i) = xcv(i) * coef
				ycv(i) = ycv(i) * coef

			ENDDO

		ELSEIF (typecentre == 2) THEN
	    ! Version "centre de gravite des milieux des aretes,
	    ! ponderes par les longueurs des aretes". Il s'agit donc en fait
	    ! du centre de gravite du bord de la maille.
			print *,"x_K = moyenne des milieux des aretes, ponderes par longeur aretes\n"
	    DO i=1, ncv
	      xcv(i) = 0D0
	      ycv(i) = 0D0
	      coef = 0D0
	      j0 = nusom(ptvois(i+1))
	      DO j=ptvois(i)+1, ptvois(i+1)
	        j1 = nusom(j)

	        long = sqrt( (xs(j1)-xs(j0))**2 + (ys(j1)-ys(j0))**2 )
	        coef = coef + long
	        xcv(i) = xcv(i) + .5D0*(xs(j0)+xs(j1)) * long
	        ycv(i) = ycv(i) + .5D0*(ys(j0)+ys(j1)) * long

	        j0 = j1
	      ENDDO
	      xcv(i) = xcv(i) / coef
	      ycv(i) = ycv(i) / coef
	    ENDDO

		ELSEIF (typecentre == 3) THEN
			! Version centre des cercles circonscrits, uniquement pour triangles
			IF (verbose==1) print *,"x_K = centres cercles circonscrits\n"
	    DO i = 1, ncv
	    	IF (ptvois(i+1)-ptvois(i)/=3) THEN
	    		print *,"Calcul centre de cercle circonscrit pour une maille non triangulaire"
	    		STOP
	    	ENDIF
        i0 = nusom(ptvois(i)+1)
        i1 = nusom(ptvois(i)+2)
        i2 = nusom(ptvois(i)+3)
        a11 = xs(i1) - xs(i0)
        a12 = ys(i1) - ys(i0)
        b1 = a11*(xs(i1)+xs(i0))+a12*(ys(i1)+ys(i0))
        a21 = xs(i2) - xs(i0)
        a22 = ys(i2) - ys(i0)
        b2 = a21*(xs(i2)+xs(i0))+a22*(ys(i2)+ys(i0))
        deter = (a11*a22-a12*a21)
        deter = 0.5D0/deter
        xcv(i) = deter*(b1*a22-b2*a12)
        ycv(i) = deter*(a11*b2-a21*b1)
      ENDDO
    
    ELSEIF (typecentre == 4) THEN
      !
      ! Centres suivants uniquement pour triangles...
      !
			IF (verbose==1) print *,"x_K = bizarroides\n"
      DO i = 1, ncv
	      IF (ptvois(i+1)-ptvois(i)/=3) THEN
	    		print *,"Calcul centre de type triangles pour une maille non triangulaire"
	    		STOP
	    	ENDIF
        i0 = nusom(ptvois(i)+1)
        i1 = nusom(ptvois(i)+2)
        i2 = nusom(ptvois(i)+3)
        coef1 = 1.D0 + dble(modulo(i,5)+1)
        coef2 = 4.D0 + dble(modulo(i,5)+1)
        coef3 = 7.D0 + dble(modulo(i,5)+1)
        coeft = coef1+coef2+coef3
        coef1 = coef1/coeft
        coef2 = coef2/coeft
        coef3 = coef3/coeft
        xcv(i) = coef1*xs(i0)+coef2*xs(i1)+coef3*xs(i2)
        ycv(i) = coef1*ys(i0)+coef2*ys(i1)+coef3*ys(i2)
      ENDDO
		
		ELSEIF (typecentre == 5) THEN
			! Barycentre des trois premiers sommets uniquement
			IF (verbose==1) print *,"x_K = barycentre des trois premiers sommets\n"
			DO i=1, ncv
				xcv(i) = 0D0
				ycv(i) = 0D0

				DO j=ptvois(i)+1, ptvois(i)+3
				  xcv(i) = xcv(i) + xs(nusom(j))
				  ycv(i) = ycv(i) + ys(nusom(j))
				ENDDO
				coef = 1D0/3D0
				xcv(i) = xcv(i) * coef
				ycv(i) = ycv(i) * coef

			ENDDO

		ELSE
			print *,"typecentre",typecentre,"inconnu"
			STOP
		ENDIF

  END SUBROUTINE calcul_centres


!------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------!
!******************* REPRODUCTION D'UN MAILLAGE INITIAL DANS DES SUBDIVISIONS DE L'OUVERT *************************!


  SUBROUTINE reproduction_maillage
  !******************************************!
  ! Reproduit un maillage initial dans des subdivisions de l'ouvert, obtenues
  ! en coupant en ndivx selon x et en ndivy selon y
  !
  ! Ne marche que sur un carre (0,1)^2
  !
  ! Si typepoly == -30 ou -40, la reproduction se fait par simple translation (il faut donc que
  ! chaque sommet de x=1 corresponde a un sommet de x=0, et idem pour y)
  ! Si typeoly == -31 ou -41, la reproduction se fait par symetrie selon les axes et il
  ! n'y a donc pas de condition sur les sommets au bord
  !
  !*****
  ! Utilise: ncv0, nbsom0, ptvois0, nusom0, xcv0, ycv0, xs0, ys0
  ! Affecte: ncv, nbsom, nptvois, ptvois, nusom, xcv, ycv, xs, ys
  !******************************************!
    INTEGER :: i, j, k, is, ii
    INTEGER :: pos_ptvois, pos_nusom, pos_xsys
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nv_som0
    DOUBLE PRECISION :: xsnew, ysnew, bordi, bordj
    INTEGER :: debut, fin, increment

    REAL :: t1, t2, ttot
    ttot = 0.

  !   CALL CPU_TIME(t1)

      ALLOCATE(nv_som0(nbsom0))

      ! Calcul du nombre de mailles, de nptvois, et d'un majorant du nombre de sommets
      ncv = ncv0 * ndivx * ndivy
      nptvois = ptvois0(ncv0+1) * ndivx * ndivy
      nbsom = nbsom0 * ndivx * ndivy

      ALLOCATE(xcv(ncv),ycv(ncv),xs(nbsom),ys(nbsom),nusom(nptvois),ptvois(ncv+1))
      xcv(:) = 0D0
      ycv(:) = 0D0
      xs(:) = 0D0
      ys(:) = 0D0
      nusom(:) = 0
      ptvois(:) = 0

      ! On coupe l'ouvert en ndivx selon x et ndivy selon y puis on se ballade dans chaque
      ! ligne (reperee par j, la position dans la ligne etant reperee par i) ainsi cree
      ! pour creer les nouvelles mailles et nouveaux sommets
      !
      ! pos_ptvois repere la position actuelle dans ptvois, i.e. le numero de la maille en cours
      !   de creation
      ! pos_nusom est la position actuelle dans nusom
      ! pos_xsys est la position actuelle dans xsys, i.e. le numero du sommet en cours de creation
      pos_ptvois = 0
      pos_nusom = 0
      pos_xsys = 0


      DO j = 1, ndivy
        DO i = 1, ndivx
          ! Dans le petit ouvert correspondant a (i,j), on reproduit le maillage initial
          ! On commence par creer les nouveaux sommets: nv_som0 est un tableau qui contient
          ! les nouveaux numeros (lies aux xs,ys en cours de creation et au futur nusom)
          ! des sommets du maillage initial dans cette reproduction a l'interieur de l'ouvert (i,j)
          nv_som0(:) = 0

          DO k = 1, nbsom0
            ! Coordonnees du sommet que l'on va creer a partir du sommet k du maillage initial,
            ! soit par simple translation, soit par symetrie
            IF (typepoly == -30 .OR. typepoly == -40) THEN
              xsnew = (xs0(k) + i-1D0)/ndivx
              ysnew = (ys0(k) + j-1D0)/ndivy
            ELSE
              xsnew = modulo(i,2) * (xs0(k) + i-1D0)/ndivx + modulo(i+1,2) * (-xs0(k) + i)/ndivx
              ysnew = modulo(j,2) * (ys0(k) + j-1D0)/ndivy + modulo(j+1,2) * (-ys0(k) + j)/ndivy
            ENDIF

            ! Dans le cas d'une reproduction par translation, si le sommet k en question est sur le
            ! bord x=0 ou y=0, celui que l'on veut creer a deja ete cree (si j>1 ou i>); sinon on le cree
            ! Dans le cas d'une reproduction par symetrie, le bord en question (x=0 ou x=1, y=0 ou y=1)
            ! depend des parites de i, j
            IF (typepoly == -30 .OR. typepoly == -40) THEN
              bordi = 0D0
              bordj = 0D0
            ELSE
              bordi = dble(modulo(i+1,2))
              bordj = dble(modulo(j+1,2))
            ENDIF

            IF ( ( (abs(xs0(k)-bordi)<1D-13) .AND. (i>1) ) .OR. ( (abs(ys0(k)-bordj)<1D-13) .AND. (j>1) ) ) THEN
              ! Recherche du sommet deja cree ayant xsnew et ysnew comme coordonnees
              ! On part du dernier sommet cree et on remonte car le sommet que l'on recherche
              ! a ete cree au pire au rang i-1 ou j-1, donc fait partie globalement des derniers
              ! sommets crees, pas des premiers

              ii = pos_xsys
              DO WHILE ( (ii>=1) .AND. ( (abs(xs(ii)-xsnew)>1D-13) .OR. (abs(ys(ii)-ysnew)>1D-13) ) )
                ii = ii - 1
              ENDDO
              IF ( (abs(xs(ii)-xsnew)>1D-13) .OR. (abs(ys(ii)-ysnew)>1D-13) ) THEN
                print *,"Erreur sommet non-trouve",i,j,k,xsnew,ysnew
                STOP
              ENDIF
              nv_som0(k) = ii

            ELSE
              ! On cree le sommet
              pos_xsys = pos_xsys + 1
              xs(pos_xsys) = xsnew
              ys(pos_xsys) = ysnew
              nv_som0(k) = pos_xsys

            ENDIF

          ENDDO

          ! On cree ensuite les mailles a partir des mailles du maillage initial
          DO k = 1, ncv0
            ! Creation d'une maille
            pos_ptvois = pos_ptvois + 1
            IF (typepoly == -30 .OR. typepoly == -40) THEN
              xcv(pos_ptvois) = (xcv0(k) + i-1D0)/ndivx
              ycv(pos_ptvois) = (ycv0(k) + j-1D0)/ndivy
            ELSE
              xcv(pos_ptvois) = modulo(i,2) * (xcv0(k) + i-1D0)/ndivx + modulo(i+1,2) * (-xcv0(k) + i)/ndivx
              ycv(pos_ptvois) = modulo(j,2) * (ycv0(k) + j-1D0)/ndivy + modulo(j+1,2) * (-ycv0(k) + j)/ndivy
            ENDIF

            ! Reperage des sommets de la maille, a partir de ceux de la maille k

            ! Dans le cas d'une reproduction par translation, il suffit de remettre les sommets sur la
            ! nouvelle maille dans l'ordre dans lequel ils apparaissent dans la maille initiale
            ! Dans le cas d'une reproduction par symetrie, il faut parfois prendre l'ordre inverse
            ! (selon la parite du nombre de symetries effectuees) pour etre sur d'ecrire les sommets
            ! de la nouvelle maille dans l'ordre trigonometrique

            ptvois(pos_ptvois+1) = ptvois(pos_ptvois) + ptvois0(k+1)-ptvois0(k)
            IF (typepoly == -30  .OR. typepoly == -40 .OR. modulo(i+j,2)==0) THEN
              debut = ptvois0(k)+1
              fin = ptvois0(k+1)
              increment = 1
            ELSE
              debut = ptvois0(k+1)
              fin = ptvois0(k)+1
              increment = -1
            ENDIF

            DO is = debut,fin,increment
              IF (nv_som0(nusom0(is)) == 0) THEN
                print *,"Erreur sommet non repere"
                STOP
              ENDIF

              pos_nusom = pos_nusom + 1
              nusom(pos_nusom) = nv_som0(nusom0(is))
            ENDDO

          ENDDO
        ENDDO
      ENDDO


      nbsom = pos_xsys
      ! pos_ptvois et pos_nusom doivent ici etre a ncv et nptvois
      IF ( (pos_ptvois /= ncv) .OR. (pos_nusom /= nptvois) ) THEN
        print *,"Erreur pos_ptvois/=ncv ou pos_nusom/=nptvois",pos_ptvois,ncv,pos_nusom,nptvois
        STOP
      ENDIF


      DEALLOCATE(nv_som0)

  !   CALL CPU_TIME(t2)
  !   ttot=ttot+t2-t1
  !   print *,"Temps reproduction maillage",ttot


  END SUBROUTINE reproduction_maillage


!------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------!
!******************* REPRODUCTION ANISOTROPIQUE D'UN MAILLAGE INITIAL DANS DES SUBDIVISIONS DE L'OUVERT *************************!


  SUBROUTINE reproduction_maillage_anisotrope
  !******************************************!
  ! Reproduit par translation un maillage initial dans des subdivisions de l'ouvert, obtenues
  ! en coupant en ndiv, et de maniere anisotrope
  !
  ! Ne marche que sur un carre (0,1)^2
  !
  ! Il faut que chaque sommet de x=1 corresponde a un sommet de x=0, et idem pour y
  !
  ! L'anisotropie est faite ainsi: avec a=anisotropie, on veut que si les tranches en y pour y>1/2
  ! ont une epaisseur de c, celles pour y<1/2 ont une epaisseur a*c. Le nombre total de tranches sera
  ! ndiv. Donc avec n1 le nb de tranches superieures et n2 celui des tranches inferieures, on a
  !   c*n1=1/2, a*c*n2=1/2 donc n1=a*n2 (a doit donc etre rationel!)
  !   n1+n2=ndiv
  ! Donc n2= 1/(a+1) * ndivy et n1=a/(a+1) * ndiv.
  ! La division en x est fixee par 2*n1 (les mailles en haut seront isotropes)
  !
  !
  !*****
  ! Utilise: ncv0, nbsom0, ptvois0, nusom0, xcv0, ycv0, xs0, ys0, anisotropy
  ! Affecte: ncv, nbsom, nptvois, ptvois, nusom, xcv, ycv, xs, ys
  !******************************************!
    INTEGER :: i, j, k, is, ii
    INTEGER :: pos_ptvois, pos_nusom, pos_xsys, n1, n2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nv_som0
    DOUBLE PRECISION :: xsnew, ysnew, bordi, bordj
    INTEGER :: debut, fin, increment

    REAL :: t1, t2, ttot
    ttot = 0.

    ! Compute n1 and n2
    n2 = floor(ndiv/(anisotropy+1))
    n1 = ndiv - n2
    IF (n1 == 0 .OR. n2 == 0) THEN 
      print *,"Nombre de subdivisions nulles, on arete", n1, n2
      STOP
    ENDIF

    ndivx = 2*n1
    ndivy = ndiv

  !   CALL CPU_TIME(t1)

    ALLOCATE(nv_som0(nbsom0))

    ! Calcul du nombre de mailles, de nptvois, et d'un majorant du nombre de sommets
    ncv = ncv0 * ndivx * ndivy
    nptvois = ptvois0(ncv0+1) * ndivx * ndivy
    nbsom = nbsom0 * ndivx * ndivy

    ALLOCATE(xcv(ncv),ycv(ncv),xs(nbsom),ys(nbsom),nusom(nptvois),ptvois(ncv+1))
    xcv(:) = 0D0
    ycv(:) = 0D0
    xs(:) = 0D0
    ys(:) = 0D0
    nusom(:) = 0
    ptvois(:) = 0

    ! On coupe l'ouvert en ndivx selon x et ndivy selon y puis on se ballade dans chaque
    ! ligne (reperee par j, la position dans la ligne etant reperee par i) ainsi cree
    ! pour creer les nouvelles mailles et nouveaux sommets
    !
    ! pos_ptvois repere la position actuelle dans ptvois, i.e. le numero de la maille en cours
    !   de creation
    ! pos_nusom est la position actuelle dans nusom
    ! pos_xsys est la position actuelle dans xsys, i.e. le numero du sommet en cours de creation
    pos_ptvois = 0
    pos_nusom = 0
    pos_xsys = 0


    DO j = 1, ndivy
      DO i = 1, ndivx
        ! Dans le petit ouvert correspondant a (i,j), on reproduit le maillage initial
        ! On commence par creer les nouveaux sommets: nv_som0 est un tableau qui contient
        ! les nouveaux numeros (lies aux xs,ys en cours de creation et au futur nusom)
        ! des sommets du maillage initial dans cette reproduction a l'interieur de l'ouvert (i,j)
        nv_som0(:) = 0

        DO k = 1, nbsom0
          ! Coordonnees du sommet que l'on va creer a partir du sommet k du maillage initial
          xsnew = (xs0(k) + i-1D0)/ndivx
          IF ( j<=n2 ) THEN
            ysnew = (ys0(k) + j-1D0)/(2*n2)
          ELSE
            ysnew = 0.5D0 + (ys0(k) + (j-n2-1D0))/(2*n1)
          ENDIF

          ! Dans le cas d'une reproduction par translation, si le sommet k en question est sur le
          ! bord x=0 ou y=0, celui que l'on veut creer a deja ete cree (si j>1 ou i>); sinon on le cree
          bordi = 0D0
          bordj = 0D0

          IF ( ( (abs(xs0(k)-bordi)<1D-13) .AND. (i>1) ) .OR. ( (abs(ys0(k)-bordj)<1D-13) .AND. (j>1) ) ) THEN
            ! Recherche du sommet deja cree ayant xsnew et ysnew comme coordonnees
            ! On part du dernier sommet cree et on remonte car le sommet que l'on recherche
            ! a ete cree au pire au rang i-1 ou j-1, donc fait partie globalement des derniers
            ! sommets crees, pas des premiers

            ii = pos_xsys
            DO WHILE ( (ii>=1) .AND. ( (abs(xs(ii)-xsnew)>1D-13) .OR. (abs(ys(ii)-ysnew)>1D-13) ) )
              ii = ii - 1
            ENDDO
            IF ( (abs(xs(ii)-xsnew)>1D-13) .OR. (abs(ys(ii)-ysnew)>1D-13) ) THEN
              print *,"Erreur sommet non-trouve",i,j,k,xsnew,ysnew
              STOP
            ENDIF
            nv_som0(k) = ii

          ELSE
            ! On cree le sommet
            pos_xsys = pos_xsys + 1
            xs(pos_xsys) = xsnew
            ys(pos_xsys) = ysnew
            nv_som0(k) = pos_xsys

          ENDIF

        ENDDO

        ! On cree ensuite les mailles a partir des mailles du maillage initial
        DO k = 1, ncv0
          ! Creation d'une maille
          pos_ptvois = pos_ptvois + 1
          xcv(pos_ptvois) = (xcv0(k) + i-1D0)/ndivx
          IF ( j<=n2 ) THEN
            ycv(pos_ptvois) = (ycv0(k) + j-1D0)/(2*n2)
          ELSE
            ycv(pos_ptvois) = 1/2 + (ycv0(k) + (j-n2-1D0))/(2*n1)
          ENDIF

          ! Reperage des sommets de la maille, a partir de ceux de la maille k

          ! Dans le cas d'une reproduction par translation, il suffit de remettre les sommets sur la
          ! nouvelle maille dans l'ordre dans lequel ils apparaissent dans la maille initiale
          ! Dans le cas d'une reproduction par symetrie, il faut parfois prendre l'ordre inverse
          ! (selon la parite du nombre de symetries effectuees) pour etre sur d'ecrire les sommets
          ! de la nouvelle maille dans l'ordre trigonometrique

          ptvois(pos_ptvois+1) = ptvois(pos_ptvois) + ptvois0(k+1)-ptvois0(k)
          debut = ptvois0(k)+1
          fin = ptvois0(k+1)
          increment = 1

          DO is = debut,fin,increment
            IF (nv_som0(nusom0(is)) == 0) THEN
              print *,"Erreur sommet non repere"
              STOP
            ENDIF

            pos_nusom = pos_nusom + 1
            nusom(pos_nusom) = nv_som0(nusom0(is))
          ENDDO

        ENDDO
      ENDDO
    ENDDO

    nbsom = pos_xsys
    ! pos_ptvois et pos_nusom doivent ici etre a ncv et nptvois
    IF ( (pos_ptvois /= ncv) .OR. (pos_nusom /= nptvois) ) THEN
      print *,"Erreur pos_ptvois/=ncv ou pos_nusom/=nptvois",pos_ptvois,ncv,pos_nusom,nptvois
      STOP
    ENDIF

    DEALLOCATE(nv_som0)


  END SUBROUTINE reproduction_maillage_anisotrope


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!******************* CALCUL MAILLAGE DUAL *********************************************************!


  SUBROUTINE calcul_maillage_dual
  !******************************************!
  ! A partir d'un maillage primal complet (avec tous ses pointeurs: nuvois, etc.), 
  !	calcule le maillage dual obtenu en reliant les centres des mailles autour de chaque
  !	sommet (et les milieux des aretes de bords et sommets de bord quand au bord)
  !
  !*****
  ! Utilise: pointeurs du maillage primal
  ! Affecte: nbsom, ncv, nptvois, ptvois, nusom, xcv, ycv, xs, ys
  !******************************************!
		INTEGER :: ncv_p, nbsom_p, nar_p, nptvois_p
		INTEGER, ALLOCATABLE, DIMENSION(:) :: ptvois_p, nusom_p, numar_p, nuvois_p
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xs_p, ys_p, xcv_p, ycv_p
		INTEGER :: i, j, j0, j1, jj, k, jp1
		INTEGER, ALLOCATABLE, DIMENSION(:) :: nusomd_ar, nusomd_so, somtrait
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xmil_ar, ymil_ar
		INTEGER :: nb_somd_ar, nb_somd_so, traiter_sommet, voisin_suivant, arete_cours
		INTEGER, ALLOCATABLE, DIMENSION(:) :: indices

		!
		! On sauve les donnees du maillage primal
		!
		ncv_p = ncv
		nbsom_p = nbsom
		nar_p = nar
		nptvois_p = ptvois(ncv+1)
		ALLOCATE(ptvois_p(ncv_p+1),nusom_p(nptvois_p),numar_p(nptvois_p),nuvois_p(nptvois_p))
		DO i=1, ncv_p+1
			ptvois_p(i) = ptvois(i)
		ENDDO
		DO i=1, nptvois_p
			nusom_p(i) = nusom(i)
			numar_p(i) = numar(i)
			nuvois_p(i) = nuvois(i)
		ENDDO
		ALLOCATE(xs_p(nbsom_p),ys_p(nbsom_p),xcv_p(ncv_p),ycv_p(ncv_p))
		DO i=1, nbsom_p
			xs_p(i) = xs(i)
			ys_p(i) = ys(i)
		ENDDO
		DO i=1, ncv_p
			xcv_p(i) = xcv(i)
			ycv_p(i) = ycv(i)
		ENDDO
		DEALLOCATE(ptvois,nusom,numar,nuvois,xs,ys,xcv,ycv)
		
		!
		! Les aretes et sommets primaux au bord seront des sommets duaux: on les
		!	repere en notant leurs numeros (nb_XXX) de nouveaux sommets duaux.
		!	nusomd_ar(i) ou nusomd_so(i) reperent ce numero (relatif) pour l'arete/sommet i (avec 0
		! si l'arete/sommet n'est pas au bord)
		!	On repere aussi les milieux des aretes de bord
		!
		nb_somd_ar = 0
		nb_somd_so = 0
		ALLOCATE(nusomd_ar(nar_p),nusomd_so(nbsom_p),xmil_ar(nar_p),ymil_ar(nar_p))
		nusomd_ar(:) = 0
		nusomd_so(:) = 0
		DO i=1, ncv_p
			j0 = nusom_p(ptvois_p(i+1))
			DO j=ptvois_p(i)+1,ptvois_p(i+1)
				j1 = nusom_p(j)

				IF (nuvois_p(j)==0) THEN
					! Le sommet nusom_p(j) et l'arete numar_p(j) sont au bord, on les repere
					! si ca n'etait pas deja fait
					IF (nusomd_so(nusom_p(j)) == 0) THEN
						nb_somd_so = nb_somd_so+1
						nusomd_so(nusom_p(j)) = nb_somd_so
					ENDIF
					IF (nusomd_ar(numar_p(j)) == 0) THEN
						nb_somd_ar = nb_somd_ar+1
						nusomd_ar(numar_p(j)) = nb_somd_ar
						xmil_ar(numar_p(j)) = 0.5D0*(xs_p(j0)+xs_p(j1))
						ymil_ar(numar_p(j)) = 0.5D0*(ys_p(j0)+ys_p(j1))
					ENDIF
				ENDIF

				j0=j1
			ENDDO	! Fin boucle j=ptvois_p(i)+1,ptvois_p(i+1)
		ENDDO	! Fin boucle i=1, ncv_p

		!
		! SOMMETS DU MAILLAGE DUAL
		!		1 -> ncv_p : centres de mailles
		!		ncv_p+1 -> ncv_p+nb_somd_ar : milieux d'aretes de bord
		!		ncv_p+nb_somd_ar -> ncv_p+nb_somd_ar+nb_somd_so : sommets primaux de bord
		!
		! On en profite pour mettre les numeros absolus dans nusomd_ar et nusomd_so
		nbsom = ncv_p+nb_somd_ar+nb_somd_so
		ALLOCATE(xs(nbsom),ys(nbsom))
		DO i=1, ncv_p
			xs(i) = xcv_p(i)
			ys(i) = ycv_p(i)
		ENDDO
		DO i=1, nar_p
			IF (nusomd_ar(i) /= 0) THEN
				nusomd_ar(i) = nusomd_ar(i) + ncv_p
				xs(nusomd_ar(i)) = xmil_ar(i)
				ys(nusomd_ar(i)) = ymil_ar(i)
			ENDIF
		ENDDO
		DO i=1, nbsom_p
			IF (nusomd_so(i) /= 0) THEN
				nusomd_so(i) = nusomd_so(i) + ncv_p + nb_somd_ar
				xs(nusomd_so(i)) = xs_p(i)
				ys(nusomd_so(i)) = ys_p(i)
			ENDIF
		ENDDO
		DEALLOCATE(xmil_ar,ymil_ar)
		
		!
		!	MAILLES DUALES
		!
		!	Autant que sommets primaux: la numerotation des sommets primaux sera
		! celle des mailles duales
		
		ncv = nbsom_p
		ALLOCATE(ptvois(ncv+1))
		ptvois(:) = 0
		! On repere le nombre de sommets pour chaque maille, i.e.
		! le nombre de mailles primale auquel chaque sommet primal est lie
		DO i=1, ncv_p
			DO j=ptvois_p(i)+1,ptvois_p(i+1)
				ptvois(nusom_p(j)) = ptvois(nusom_p(j))+1
			ENDDO
		ENDDO	! Fin boucle i=1,ncv_p
		
		! Les sommets primaux de bord correspondent a des mailles duales ayant
		! 3 sommets duaux de plus (les deux centres d'aretes de bord, et le sommet
		! de bord lui-meme)
		DO i=1, nbsom_p
			IF (nusomd_so(i) /= 0) THEN
				ptvois(i) = ptvois(i) + 3
			ENDIF
		ENDDO
		
		! On transforme ptvois en un pointeur croissant
		DO i=2, ncv+1
			ptvois(i) = ptvois(i)+ptvois(i-1)
		ENDDO
		ptvois(ncv+1) = ptvois(ncv)
		nptvois = ptvois(ncv+1)

		! On cree les mailles duales en remplissant nusom (sommets de chaque maille)
		ALLOCATE(nusom(nptvois),somtrait(nbsom_p))
		nusom(:) = 0
		somtrait(:) = 0
		! On se ballade dans les mailles et sommets primaux et, pour chaque sommet
		! n'ayant pas encore ete traite (somtrait=0) on parcours les mailles de voisin
		! en voisin autour du sommet en question pour creer une maille duale
		DO i=1, ncv_p
			DO j=ptvois_p(i)+1, ptvois_p(i+1)
				jj = nusom_p(j)

				! On regarde si on doit creer la maille duale autour du sommet primal jj				
				traiter_sommet = 0
				IF (somtrait(jj) == 0) THEN
					IF (nusomd_so(jj) /= 0) THEN
						! Un sommet de bord n'est traite que si le voisin de l'arete suivante est 0
						! Cela assure que, en prenant les voisins successifs lorsqu'on tourne autour du
						! sommet, on attrape toutes les mailles
						jp1 = j+1
						IF (jp1>ptvois_p(i+1)) jp1 = ptvois_p(i)+1
						IF (nuvois_p(jp1) == 0) traiter_sommet = 1

					ELSE
						traiter_sommet = 1
					ENDIF
				ENDIF
				
				IF (traiter_sommet == 1) THEN
					somtrait(jj) = 1
					! On traite le sommet primal jj (maille duale jj)
					nusom(ptvois(jj)) = i
					ptvois(jj) = ptvois(jj)-1
	
					voisin_suivant = nuvois_p(j)
					arete_cours = numar_p(j)
					! On parcours les voisins tant qu'on ne tombe pas sur le bord ou
					! sur la maille initiale
					DO WHILE (voisin_suivant /= 0 .AND. voisin_suivant /= i)
						! On met le voisin dans nusom
						nusom(ptvois(jj)) = voisin_suivant
						ptvois(jj) = ptvois(jj)-1
						
						! On cherche ou apparait le sommet jj dans la maille qu'on vient d'ajouter
						k = ptvois_p(voisin_suivant)+1
						DO WHILE (k<=ptvois_p(voisin_suivant+1) .AND. nusom_p(k) /= jj)
							k = k+1
						ENDDO
						IF (k>ptvois_p(voisin_suivant+1) .OR. nusom_p(k) /= jj) THEN
							print *,"Erreur construction maillage dual",i,voisin_suivant,jj
							STOP
						ENDIF

						voisin_suivant = nuvois_p(k)
						arete_cours = numar_p(k)
					ENDDO

					! Si on a un sommet de bord, on a du tomber sur 0 comme voisin_suivant: verification
					IF (nusomd_so(jj) /= 0 .AND. voisin_suivant /= 0) THEN
						print *,"Erreur sommet de bord, pas voisin_suivant nul",i,voisin_suivant,jj
						STOP
					ENDIF
					
					! Si le voisin suivant etait 0, on est tombe au bord: on ajoute
					! le milieu de la maille de bord comme nouveau sommet de la maille duale,
					! ainsi que le sommet de bord lui-meme et la deuxieme arete de bord qui le contient
					IF (voisin_suivant == 0) THEN
						nusom(ptvois(jj)) = nusomd_ar(arete_cours)
						ptvois(jj) = ptvois(jj)-1
						
						IF (nusomd_so(jj)==0) THEN
							print *,"Erreur: sommet au bord mal repere",i,jj
							STOP
						ENDIF
						nusom(ptvois(jj)) = nusomd_so(jj)
						ptvois(jj) = ptvois(jj)-1
						
						! La deuxieme arete de bord contenant jj est la suivante dans la maille i
						jp1 = j+1
						IF (jp1>ptvois_p(i+1)) jp1=ptvois_p(i)+1
						nusom(ptvois(jj)) = nusomd_ar(numar_p(jp1))
						ptvois(jj) = ptvois(jj)-1
						
					ENDIF
					
				ENDIF
			ENDDO	!	Fin boucle j=ptvois_p(i)+1,ptvois_p(i+1)
		ENDDO	!	Fin boucle i=1,ncv_p

		IF (ptvois(1)/=0) THEN
			print *,"Erreur pointeur maillage dual",ptvois(1)
			STOP
		ENDIF
		
		! Les sommets duaux ont ete entres en ordre trigo inverse, il faut les
		! remettre dans l'ordre trigonometrique
		ALLOCATE(indices(nbsom))
		indices(:) = 0
		DO i=1, ncv
			DO j=ptvois(i)+1,ptvois(i+1)
				indices(j-ptvois(i)) = nusom(j)				
			ENDDO
			DO j=ptvois(i)+1,ptvois(i+1)
				nusom(j) = indices(ptvois(i+1)+1-j)
			ENDDO
		ENDDO
		DEALLOCATE(indices)
		
		! Affichage pour verification
  !		DO i=1, ncv
  !			print *,"Maille duale",i
  !			DO j=ptvois(i)+1,ptvois(i+1)
  !				IF (nusom(j)<=ncv_p) THEN
  !					print *,"     maille",nusom(j)
  !				ELSEIF (nusom(j)<=ncv_p+nb_somd_ar) THEN
  !					! On cherche le numero d'arete
  !					k=1
  !					DO WHILE (k<=nar_p .AND. nusomd_ar(k) /= nusom(j))
  !						k=k+1
  !					ENDDO
  !					print *,"     arete",k
  !				ELSE
  !					! On cherche le numero de sommet
  !					k=1
  !					DO WHILE (k<=nbsom_p .AND. nusomd_so(k) /= nusom(j))
  !						k=k+1
  !					ENDDO
  !					print *,"     sommet",k
  !				ENDIF
  !			ENDDO
  !		ENDDO
		

		DEALLOCATE(nusomd_ar,nusomd_so)
		DEALLOCATE(ptvois_p,nusom_p,numar_p,nuvois_p,xs_p,ys_p,xcv_p,ycv_p)


	END SUBROUTINE calcul_maillage_dual



!-----------------------------------------------------------------------!
!---------------- CALCULS DE POINTEURS, NUMEROTATIONS ------------------!
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!*********************calcul_pointeurs**********************************!

  SUBROUTINE calcul_pointeurs(verbose)
    !**** Calcul des tableaux de donnees du maillage, avec leurs pointeurs ****
    !
    ! Sens des variables:
    !     ptsc = pointeur servant a lire les donnees d'un sommet i dans nusc et voisind/g
    !     nusc = lu avec le pointeur ptsc, indique les numeros de sommets connectes
    !           par une arete a chaque sommet.
    !               nusc(ptsc(i)+1) -> nusc(ptsc(i+1)) : numeros des sommets connectes
    !                                         au sommet i, dans l'ordre croissant
    !     voisind/g = lu avec le pointeur ptsc, indique pour chaque sommet les numeros
    !                   des mailles a droite et a gauche des aretes partant de ce sommet
    !               voising(ptsc(i)+1) -> voising(ptsc(i+1)) : numeros des mailles "a gauche"
    !                                         des aretes [i,k] avec k=nusc(ptsc(i)+1),nusc(ptsc(i+1))
    !                    "a gauche" signifie que l'arete [i,k] est parcouru dans le sens trigo dans
    !                               la maille concernee;
    !                    "a droite" l'arete est parcouru en sens trigo inverse dans la maille concernee.
    !
    ! les variables en "new" (ptvoisnew, nusomnew, nuvoisnew...) servent quand on renumerote les mailles,
    !     et une fois cette renumerotation faite on met ces tableaux dans ptvois, nusom, nuvois et numar.
    !
    !	verbose = 1 si on veut des sorties dans la console, 0 sinon
    !
    !***************************
    ! Utilise les variables globales: nbsom,ncv,nptvois,ptvois,nusom,xcv,ycv,xs,ys
    ! Affecte: numar,nar,nuvois
    !**************************************************************************
    INTEGER :: verbose
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ptsc, nusc, voisind, voising, renum, denum
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ptvoisnew, nusomnew, nuvoisnew, numarnew
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xcvnew, ycvnew
    INTEGER :: i, j, k, is, js, iar, icv, jj, jjj, j0, j1, jj0, i0, nbc
    CHARACTER(LEN=10) :: formataffichage


    ! Allocation de variables globales
    ALLOCATE(nuvois(nptvois),numar(nptvois))

    IF (verbose==1) print "(10X,A,I9,5X,A,I9)","nb sommets=",nbsom,"nb elements=",ncv

    ! Donnees extraites des sommets
    ! Initialisation du pointeur des sommets a 0
    ALLOCATE (ptsc(nbsom+1))

    DO i=1,nbsom
       ptsc(i) = 0
    ENDDO

    ! On met dans ptsc(is) un majorant du nombre de sommets relies
    ! au sommet is (on se ballade dans les mailles: chaque fois que is
    ! apparait dans une maille, il est relie a 2 autres sommets de cette maille.
    DO i=1,ncv
      DO j=ptvois(i)+1,ptvois(i+1)
        is = nusom(j)
        ptsc(is) = ptsc(is) + 2
      ENDDO
    ENDDO

    ! On translate ptsc pour en faire un tableau d'indices strictement
    ! croissant plutot qu'un tableau de nombre de sommets lies (on fait la
    ! somme de tous les precedents)
    DO i=2,nbsom
      ptsc(i) = ptsc(i-1)+ptsc(i)
    ENDDO
    ! Voir tout de suite apres pourquoi cette derniere affectation
    ptsc(nbsom+1) = ptsc(nbsom)
    nbc = ptsc(nbsom)

    ! Creation du tableau des sommets lies a chaque sommet. A chaque
    ! fois qu'on rencontre un sommet connecte a is, on le met a la position
    ! ptsc(is) et on reduit cette position de 1 (les sommets sont donc
    ! rentres de droite a gauche dans chaque sequence de nusc).
    ! ptsc(is) finit donc un cran avant le 1er sommet connecte a is.
    !
    ! On a affecte ptsc(nbsom+1) a ptsc(nbsom) car c'est bien le dernier
    ! indice du dernier sommet, avant que ptsc(nbsom) n'ait ete diminue par les
    ! affectation de nusc; ou encore: ptsc(nbsom+1) doit etre egal a ptsc(nbsom)
    ! plus le nombre de sommets relies au sommet nbsom... or ptsc(nbsom) a justement
    ! ete diminue de ce nombre de sommets, donc sa valeur initiale etait ce qu'on
    ! veut pour ptsc(nbsom+1)
    ALLOCATE (nusc(nbc))

    DO i=1,ncv
      is = nusom(ptvois(i+1))
      DO j=ptvois(i)+1,ptvois(i+1)
        js = nusom(j)

        nusc(ptsc(is)) = js
        ptsc(is) = ptsc(is) - 1
        nusc(ptsc(js)) = is
        ptsc(js) = ptsc(js) - 1

        is = js
      ENDDO
    ENDDO

  !!$    print *,(ptsc(i),i=1,nbsom+1)
  !!$
  !!$    print *,(nusc(i),i=1,nbc)


    ! Remise dans l'ordre: les sommets connectes a un sommet donne sont places
    ! dans nusc en ordre croissant
    DO i = 1,nbsom
      CALL tri(nusc(ptsc(i)+1),ptsc(i+1) - ptsc(i))
    ENDDO

  !!$    print *,' apres remise en ordre '
  !!$
  !!$    print *,(ptsc(i),i=1,nbsom+1)
  !!$
  !!$    print *,(nusc(i),i=1,nbc)

    ! Elimination des doublons
    nbc = 0
    DO i=1,nbsom
      i0 = ptsc(i) + 1
      ptsc(i) = nbc

      IF (ptsc(i+1)>i0-1) THEN
        k = nusc(i0)
        nbc = nbc + 1
        nusc(nbc) = k
        DO j = i0, ptsc(i+1)
          IF (nusc(j) /= k) THEN
            k = nusc(j)
            nbc = nbc + 1
            nusc(nbc) = k
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    ptsc(nbsom + 1)  = nbc
    nbc = ptsc(nbsom + 1)


    ! print *,' apres elimination des doublons '

    ! print *,(ptsc(i),i=1,nbsom+1)

    ! print *,(nusc(i),i=1,nbc)

    !****
    !**** Calcul des voisins
    !****
    ALLOCATE (voisind(nbc),voising(nbc))
    DO i = 1, nbc
      voisind(i) = 0
      voising(i) = 0
    ENDDO

    ! Parcourir les mailles
    DO i =1, ncv
      is = nusom(ptvois(i+1))
      ! On parcours les sommets de la maille i
      DO j=ptvois(i)+1,ptvois(i+1)
        js = nusom(j)

        ! On cherche ou, dans nusc, le sommet js est relie
        ! au sommet is (en utilisant le fait que ces sommets
        ! sont classes en ordre croissant)
        k = ptsc(is)+1
        DO WHILE (k<ptsc(is+1) .AND. nusc(k)<js)
          k = k + 1
        ENDDO
        IF (nusc(k)/=js) THEN
          print *," probleme dans les voisins "
          STOP
        ENDIF
        ! L'arete [is js] etant orientee dans le sens trigo,
        ! la maille i est celle qui est a gauche de cette arete.
        voising(k) = i

        ! Idem pour le sommet js
        k = ptsc(js)+1
        DO WHILE (k<ptsc(js+1) .AND. nusc(k)<is)
          k = k + 1
        ENDDO
        IF (nusc(k)/=is) THEN
          print *," probleme dans les voisins "
          STOP
        ENDIF
        voisind(k) = i

        is = js
      ENDDO
    ENDDO

  ! do i=1,nbsom
  ! print *,"---------------",i
  ! print *,(nusc(j),j=ptsc(i)+1,ptsc(i+1))
  ! print *,(voisind(j),j=ptsc(i)+1,ptsc(i+1))
  ! print *,(voising(j),j=ptsc(i)+1,ptsc(i+1))
  ! enddo


    !***
    !*** Fabrication du tableau des voisins de chaque maille
    !***
    DO i = 1,nptvois
      nuvois(i) = 0
    ENDDO

    ! On se ballade dans les mailles et leurs aretes
    DO i = 1,ncv
      is = nusom(ptvois(i+1))

      DO j=ptvois(i)+1,ptvois(i+1)
        js = nusom(j)

        ! On cherche a quel endroit de nusc le sommet js est relie
        ! au sommet is (ca indique aussi l'endroit dans voisind ou
        ! on lira le voisin droit de [is js], c'est a dire la maille
        ! de l'autre cote de i)
        k = ptsc(is)+1
        DO WHILE (k<ptsc(is+1) .AND. nusc(k) /= js)
          k = k + 1
        ENDDO
        nuvois(j) = voisind(k)

        is = js
      ENDDO
    ENDDO

  ! do i=1,ncv
  ! print *,"---------------",i
  ! print *,(nusom(j),j=ptvois(i)+1,ptvois(i+1))
  ! print *,(nuvois(j),j=ptvois(i)+1,ptvois(i+1))
  ! enddo

    !***
    !*** Numerotation des aretes
    !***

    DO i=1,nptvois
       numar(i) = 0
    ENDDO

    ! calcul du nombre d'aretes et d'aretes internes
    nar = 0
    nar_internes = 0

    DO i=1, ncv
      DO j=ptvois(i)+1,ptvois(i+1)
        jj = nuvois(j)
        IF (jj<i) THEN
          nar = nar + 1
          IF (jj /= 0) nar_internes = nar_internes + 1
        ENDIF
      ENDDO
    ENDDO

    IF (verbose==1) print "(10X,A,I9,6X,A,I9)","nb aretes=",nar,"nb aretes internes=",nar_internes


    !!  On numerote de maniere intelligente:
    ! On parcours les mailles de proche en proche en les renumerotant:
    ! a chaque fois que l'on arrive sur une nouvelle maille, on numerote ses
    ! aretes (si elles ne l'ont pas deja ete)
    ! "renum" indique la nouvelle numerotation:
    !        renum(jj) = nouveau numero de la maille d'ancien numero jj
    ! "denum" est la pile reciproque:
    !        denum(k) est l'ancien numero de la maille de nouveau numero k

    ALLOCATE(renum(ncv),denum(ncv))

    DO i = 1,ncv
       renum(i) = 0
       denum(i) = 0
    ENDDO

    renum(1) = 1
    denum(1) = 1
    ! Numero d'arete precedente
    iar = 0
    ! Numero de maille precedente
    k = 1

    DO i = 1,ncv
      icv = denum(i)

      IF (icv==0) THEN
          print *,' du rate dans denum ',i,k
          STOP
      ENDIF

      ! On parcours les aretes de la maille icv et on les numerote,
      ! si elles ne l'ont pas deja ete
      j0 = nusom(ptvois(icv+1))
      DO j = ptvois(icv)+1,ptvois(icv+1)
        j1 = nusom(j)

        IF (numar(j)==0) THEN
          iar = iar + 1
          numar(j) = iar

          jj = nuvois(j)
          ! Si la maille de l'autre cote de l'arete n'a pas ete
          ! mise dans la pile renum, on la rajoute
          IF (jj>0) THEN
            IF (renum(jj)==0) THEN
              k = k+1
              renum(jj) = k
              denum(k) = jj
            ENDIF

            ! On cherche ensuite, pour la maille jj, a quelle position
            ! de pointeur doit apparaitre l'arete iar (entre icv et jj)
            ! que l'on vient de numeroter, pour entrer cette info. Comme il peut
            ! y avoir plusieurs aretes entre les mailles icv et jj, on doit reperer
            ! cette arete a l'aide de ses sommets (qui doivent etre ceux vus de l'arete
            ! icv en cours mais en sens inverse)
            jj0 = nusom(ptvois(jj+1))
            jjj = ptvois(jj)+1
            DO WHILE ( jjj<=ptvois(jj+1) .AND. (jj0 /= j1 .OR. nusom(jjj) /= j0 .OR. nuvois(jjj)/=icv) )
              jj0 = nusom(jjj)
              jjj = jjj + 1
            ENDDO
            IF (jjj>ptvois(jj+1)) THEN
              print *,"Erreur num arete"
              STOP
            ENDIF
            numar(jjj) = iar

          ENDIF
        ENDIF

        j0 = j1
      ENDDO

    ENDDO

  ! do i=1,ncv
  ! print *,"---------------",i
  ! print *,(nusom(jjj),jjj=ptvois(i)+1,ptvois(i+1))
  ! print *,(numar(jjj),jjj=ptvois(i)+1,ptvois(i+1))
  ! enddo

    IF (iar/=nar) THEN
       print *,' du rate dans les aretes ',iar,nar
       STOP
    ENDIF
    IF (k/=ncv) THEN
       print *,' du rate dans les cv ',k,ncv
       STOP
    ENDIF

  ! Afficher ptvois a l'ecran, avec un format variable
  ! dependant de la taille de ptvois
  !     print *,"\n ptvois: "
  !     write(formataffichage, '(a,i2,a)') '(',ncv,'i3)'
  !     print formataffichage, (ptvois(k),k=1,ncv)
  !     print *,"\n"
  !
  ! Enregistrer ptvois et numar dans un fichier 
  !     open(unit=112, File="methode2.txt",status="unknown")
  !     write(112,"(a)") "ptvois"
  !     write(112,"(15i5)") (ptvois(k),k=1,ncv)
  !     write(112,"(a)") "numar"
  !     write(112,"(20i5)") (numar(k),k=1,nptvois)
  !

    !***
    !*** Renumerotation des mailles
    !***
    ! Il s'agit en fait de la renumerotation donnee par le renum precedent (i.e.
    ! on numerote les mailles au fur et a mesure qu'elles ont ete empilees)
    ALLOCATE (ptvoisnew(ncv+1),nusomnew(nptvois),nuvoisnew(nptvois),numarnew(nptvois))
    ALLOCATE (xcvnew(ncv),ycvnew(ncv))
    DO i = 1,ncv+1
      ptvoisnew(i) = 0
    ENDDO
    DO i = 1,nptvois
      nusomnew(i) = 0
      nuvoisnew(i) = 0
      numarnew(i) = 0
    ENDDO

    ! Position dans les nouveaux tableau nusomnew, nuvoisnew et numarnew
    k = 0
    ! On se ballade dans les mailles selon la nouvelle numerotation.
    ! Pour chaque maille, on va chercher les infos sur ses sommets/voisins/aretes
    ! dans les anciens tableaux et on les colle dans le nouveau tableau, en
    ! avancant la position a chaque fois
    DO i = 1,ncv
      ptvoisnew(i) = k
      DO j = ptvois(denum(i))+1,ptvois(denum(i)+1)
        k = k+1
        nusomnew(k) = nusom(j)
        IF (nuvois(j)/=0) THEN
          nuvoisnew(k) = renum(nuvois(j))
        ENDIF
        numarnew(k) = numar(j)
      ENDDO
    ENDDO

    ptvoisnew(ncv+1) = nptvois

    ! verifications
    IF (k/=nptvois) THEN
      print *, '\n\noups, nouvelle numerotation ratee\n\n'
      STOP
    ENDIF

    ! nouveaux tableaux pour les coordonnees de points
    DO i=1,ncv
      xcvnew(i) = xcv(denum(i))
      ycvnew(i) = ycv(denum(i))
    ENDDO

  !     print *,"indice, ptvois et ptvoisnew"
  !     print *,(k,ptvois(k),ptvoisnew(k),"\n",k=1,ncv+1)
  !
  !     print *,"indice, nuvois et nuvoisnew"
  !     print *,(k,nuvois(k),nuvoisnew(k),"\n",k=1,nptvois)

  !         print *,"anc. num cv, nouv. num cv"
  !         print *,(denum(i),i,"\t sommets:",(nusomnew(j),j=ptvoisnew(i)+1,ptvoisnew(i+1)),"\n",i=1,ncv)
  !
  !       print *,"anc. num cv, nouv. num cv"
  !       print *,(denum(i),i,"\t aretes:",(numarnew(j),j=ptvoisnew(i)+1,ptvoisnew(i+1)),"\n",i=1,ncv)
  !
  !       print *,"anc. num cv, nouv. num cv"
  !       print *,(i,renum(i),"\t anc. voisins:",(nuvois(j),j=ptvois(i)+1,ptvois(i+1)),"\n",i=1,ncv)
  !
  !       print *,"nouv. num cv"
  !       print *,(i,"\t nouv. voisins:",(nuvoisnew(j),j=ptvoisnew(i)+1,ptvoisnew(i+1)),"\n",i=1,ncv)
  !
  !       print *,"nouv. num cv"
  !       print *,(k,"\t aretes:",(numarnew(j),j=ptvoisnew(k)+1,ptvoisnew(k+1)),"\n",k=1,ncv)

    ! On installe la nouvelle numerotation, que l'on utilisera a partir de maintenant
    DO i=1,ncv
      ptvois(i) = ptvoisnew(i)
      xcv(i) = xcvnew(i)
      ycv(i) = ycvnew(i)
    ENDDO
    ptvois(ncv+1) = ptvoisnew(ncv+1)
    DO i=1,nptvois
      nusom(i) = nusomnew(i)
      nuvois(i) = nuvoisnew(i)
      numar(i) = numarnew(i)
    ENDDO

    !***
    !*** Renumerotation des aretes
    !***

    ! Renumerotation efficace pour les systemes hybrides en version sparse
		CALL renumeroter_aretes_dichotomie(xs,ys,ncv,xcv,ycv,nptvois,ptvois,nusom,nuvois,numar,nar,verbose)
		
    DEALLOCATE (ptsc,nusc,voisind,voising,renum,denum,ptvoisnew,nusomnew,nuvoisnew,numarnew,xcvnew,ycvnew)

  !      print *,"\n"
  !      print *,("num cv",i,"\t aretes:",(numar(j),j=ptvois(i)+1,ptvois(i+1)),"\n",i=1,ncv)

  END SUBROUTINE calcul_pointeurs



!------------------------------------------------------------------------!
!****************renumeroter_aretes_dichotomie***************************!

  SUBROUTINE renumeroter_aretes_dichotomie(xs,ys,ncv,xcv,ycv,nptvois,ptvois,nusom,nuvois,numar,nar,verbose)
    !*******************************************************************!
    ! Renumerote les aretes par un processus de dichotomie: on coupe le domaine en deux
    !   dans le sens des x, on numerote (en commencant par la fin) les aretes qui rencontrent
    !   de maniere "presque tangente" cette interface (i.e. les mailles autour de l'arete
    !   sont de part et d'autre de l'interface), puis la meme chose en y, et on
    !   recommence en coupant en 4, 8, etc...
    ! Les aretes de bord sont numerotees en dernier (i.e. dans les premiers nouveaux numeros)
    !
    ! nvnum(i) donne le nouveau numero de l'arete i
    !
    !*******************************************************************!
    INTEGER :: ncv, nptvois, nar, verbose
    DOUBLE PRECISION, DIMENSION(:) :: xs, ys, xcv, ycv
    INTEGER, DIMENSION(:) :: ptvois, nuvois, nusom, numar

    INTEGER :: i, j, iar, pass, nar_bord
    DOUBLE PRECISION :: xmax, xmin, ymax, ymin, coef
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nvnum
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xcvabs, ycvabs

    ! Variables pour dessiner les grilles
    CHARACTER(LEN=15) :: g_fichier
    INTEGER :: g_i, g_j, g_i0, g_i1
    DOUBLE PRECISION :: g_x, g_y
    CHARACTER(LEN=10) :: g_numero

    DOUBLE PRECISION :: T1, T2

    ALLOCATE(xcvabs(ncv),ycvabs(ncv))
    xcvabs(:) = 0D0
    ycvabs(:) = 0D0

    ! On commence par faire un changement d'echelle qui ramene au besoin
    ! tout le domaine dans le carre (0,1)^2. Les xcvabs et ycvabs seront, a la fin,
    ! ces coordonnees absolues des centres de gravite des sommets de chaque maille
    ! multiplies par 2 a chaque passage de renumerotation.
    ! On calcule aussi au passage le nombre d'arete de bord
    xmax = xs(1)
    xmin = xs(1)
    ymax = ys(1)
    ymin = ys(1)
    nar_bord = 0

    DO i=1, ncv
      ! xcvabs et ycvabs seront ici les centres de gravite des sommets de chaque
      ! maille
      DO j=ptvois(i)+1, ptvois(i+1)
        xcvabs(i) = xcvabs(i) + xs(nusom(j))
        ycvabs(i) = ycvabs(i) + ys(nusom(j))

        IF (xs(nusom(j)) > xmax) xmax=xs(nusom(j))
        IF (xs(nusom(j)) < xmin) xmin=xs(nusom(j))
        IF (ys(nusom(j)) > ymax) ymax=ys(nusom(j))
        IF (ys(nusom(j)) < ymin) ymin=ys(nusom(j))

        IF (nuvois(j) == 0) nar_bord = nar_bord+1

      ENDDO
      coef = 1D0 / (ptvois(i+1) - ptvois(i))
      xcvabs(i) = xcvabs(i) * coef
      ycvabs(i) = ycvabs(i) * coef
    ENDDO
    IF (xmax == xmin .OR. ymax == ymin) THEN
      print *,"La bonne blague..."
      STOP
    ENDIF

    ! On passe en coordonnees absolues, multipliees par 2 pour preparer le 1er passage
    ! de renumerotation
    xmax = 1D0/(xmax-xmin)
    ymax = 1D0/(ymax-ymin)
    DO i=1, ncv
      xcvabs(i) = 2D0*(xcvabs(i)-xmin)*xmax
      ycvabs(i) = 2D0*(ycvabs(i)-ymin)*ymax
    ENDDO

    ! On renumerote en partant de la fin: on coupe le domaine en 2 (en x puis y)
    ! et on numerote les aretes internes qui "rencontrent" cette decoupe (i.e. telles que
    ! les mailles que l'arete borde sont de part et d'autre de la decoupe en 2). On recommence
    ! ensuite en coupant en 4, 8, etc. tant qu'on trouve des aretes internes.
    ! Pour reperer si les mailles d'un cote et de l'autre d'une arete sont de part et
    ! d'autre d'une interface obtenue a la decoupe 2^N, on verifie si les parties
    ! entieres de 2^N fois leurs coordonnees de centre de gravite sont egales ou non
    ! (en fait, cela consiste a regarder les N premiers termes du developpement en base 2:
    ! la partie entiere de 2^N * a consiste a garder ces N premiers termes de a)

    ALLOCATE(nvnum(nar))
    nvnum(:) = 0

    ! iar est le numero d'arete que l'on est en train d'attribuer (on commence par
    ! la fin). pass le numero du passage (niveau di-adique ou l'on coupe)
    iar = nar+1
    pass = 0

    CALL CPU_TIME(T1)

    ! On itere tant qu'il reste d'autres aretes que les aretes de bord a numeroter.
    ! On met une securite: ne pas faire plus de 53 passages. La raison de ce nombre:
    ! a chaque passage N on doit differencier les xcvabs et ycvabs en fonction de leur
    ! N-eme coefficient di-adique; or 2^54>1E16, donc le 54e coef di-adique n'est pas
    ! distinguable a l'erreur machine... inutile donc de continuer.
    DO WHILE (iar > nar_bord + 1 .AND. pass < 53)
      pass = pass + 1

      ! Passage en x (on coupe en x=...)
      DO i=1, ncv
        DO j=ptvois(i)+1, ptvois(i+1)

          IF ( nuvois(j) /= 0 .AND. nvnum(numar(j)) == 0 .AND. int(xcvabs(i)) /= int(xcvabs(nuvois(j))) ) THEN
            iar = iar - 1
            nvnum(numar(j)) = iar
          ENDIF

        ENDDO
      ENDDO

      ! Passage en y (on coupe en y=...)
      DO i=1, ncv
        DO j=ptvois(i)+1, ptvois(i+1)
          IF ( nuvois(j) /= 0 .AND. nvnum(numar(j)) == 0 .AND. int(ycvabs(i)) /= int(ycvabs(nuvois(j))) ) THEN
            iar = iar - 1
            nvnum(numar(j)) = iar
          ENDIF
        ENDDO
      ENDDO

      ! Preparation des tableaux de centres pour l'iteration suivante (on multiplie par 2
      ! a chaque passage)
      DO i=1, ncv
        xcvabs(i) = 2D0*xcvabs(i)
        ycvabs(i) = 2D0*ycvabs(i)
      ENDDO

    !       print *,"Renumerotation: passage, niv aretes, nb aretes",pass,iar,nar

          !!!!!!!!!!!!
          ! Dessiner les grilles du passage, uniquement si le nombre de mailles est raisonnable
                IF (ncv <= 0) THEN
                  WRITE(g_fichier,'(A,I4.4,A)') "grid",pass,".fig"
                  OPEN(unit=21, File=g_fichier, status="unknown")
                  WRITE(21,"(A)")"#FIG 3.2"
                  WRITE(21,"(A)")"Landscape"
                  WRITE(21,"(A)")"Center"
                  WRITE(21,"(A)")"Metric"
                  WRITE(21,"(A)")"A4"
                  WRITE(21,"(A)")"100.00"
                  WRITE(21,"(A)")"Single"
                  WRITE(21,"(A)")"-2"
                  WRITE(21,"(A)")"1200 2"

                  DO g_i=1, ncv
                    ! On ecrit le numero de la maille
                    WRITE(g_numero,"(I5)") g_i
    !                     CALL positionner_texte(21,g_numero,xcv(g_i),ycv(g_i),"noir")

                    g_i0 = nusom(ptvois(g_i+1))
                    DO g_j=ptvois(g_i)+1,ptvois(g_i+1)
                      g_i1=nusom(g_j)
                      ! On ne trace que les aretes de bord ou celles venant d'etre numerotees,
                      ! puisqu'on numerote en sens inverse. En lisant les fichiers grid000X.fig
                      ! en sens inverse, on verra donc les aretes disparaitre dans l'ordre de leurs
                      ! nouveaux numeros.
                      IF (nvnum(numar(g_j))/=0 .OR. nuvois(g_j)==0) THEN
                        ! On trace l'arete
                        CALL line(21,xs(g_i0),ys(g_i0),xs(g_i1),ys(g_i1))
                        ! On ecrit son numero au centre de l'arete
                        g_x = (xs(g_i0)+xs(g_i1))/2
                        g_y = (ys(g_i0)+ys(g_i1))/2
                        WRITE(g_numero,"(I5)") numar(g_j)
    !                         CALL positionner_texte(21,g_numero,g_x,g_y,"rouge")
                      ENDIF
                      g_i0 = g_i1
                    ENDDO
                  ENDDO
                  CLOSE(21)
                ENDIF
          !!!!!!!!!!!!!!!!

    ENDDO

    ! Numerotation des aretes restantes (au bord)
    DO i=1, ncv
      DO j=ptvois(i)+1, ptvois(i+1)
        IF (nuvois(j) == 0) THEN
          iar = iar - 1
          nvnum(numar(j)) = iar
        ENDIF
      ENDDO
    ENDDO

    IF (iar /= 1) THEN
      print *,"Erreur renumerotation des aretes",iar
      STOP
    ENDIF

    CALL CPU_TIME(T2)
    !    IF (verbose==1) print "(A,I3,X,A,2X,ES9.2,X,A)","\nRenumerotation des aretes:",pass,"passages,",T2-T1,"secondes"

    ! Installation de la nouvelle numerotation
    DO i=1, nptvois
      ! print *,numar(i),nvnum(numar(i))
      numar(i) = nvnum(numar(i))
    ENDDO

    DEALLOCATE(xcvabs,ycvabs,nvnum)

  END SUBROUTINE renumeroter_aretes_dichotomie



!-----------------------------------------------------------------------------!
!------------------- TAILLES DU MAILLAGE -------------------------------------!
!-----------------------------------------------------------------------------!

!---------------------------------------------------------------------------!
!********************cal_mcv_reg*******************************************!

  SUBROUTINE cal_mcv_reg
    !***************************************************
    ! Calcule l'aire des mailles, leur diametre et la regularite du maillage
    !
    !*********
    ! Utilise les variables: ncv,ptvois,nusom,xs,ys,xcv,ycv,mcv,diam,sizem
    ! Affecte: sizem,mcv,diam
    !***************************************************
    INTEGER :: i, j, k, i1, i2
    DOUBLE PRECISION :: a11, a12, a21, a22,  ax, ay, aire, long, xx, diamtemp, maxlongar, &
    			xs1, ys1, dKs, dLs
		DOUBLE PRECISION :: tt
		tt=0D0
    ! Allocation de variables generales
    ALLOCATE(mcv(ncv),diam(ncv))

    reg_maillage = 1D0
    airetot = 0D0
    sizem = 0D0
    maxlongar = 0D0

    ! On balaye les mailles
    DO i=1, ncv
      mcv(i) = 0D0
      diam(i) = 0D0
      i1 = nusom(ptvois(i+1))
      ! vecteur centre-sommet
      a11 = xs(i1) - xcv(i)
      a12 = ys(i1) - ycv(i)

      ! Calcul de l'aire: on balaye toutes les aretes de la
      ! maille et on rajoute l'aire du triangle centre-sommet-sommet
      ! pour chaque arete (=moitie de l'aire du parallelogramme
      ! engendre par centre et arete)
      DO j=ptvois(i)+1,ptvois(i+1)

        i2 = nusom(j)
        ! vecteur centre-sommet
        a21 = xs(i2) - xcv(i)
        a22 = ys(i2) - ycv(i)

        ! aire du parallelogramme centre + arete
        aire = a11*a22 - a12*a21

        ! Si l'aire est negative (centre hors de la maille) on fait exploser
        ! la regularite
				IF (aire<0D0) reg_maillage=1D10

        ! longueur de l'arete
        ax = xs(i2) - xs(i1)
        ay = ys(i2) - ys(i1)
        long =ax*ax + ay*ay
        IF (maxlongar < sqrt(long)) maxlongar = sqrt(long)

				! Partie de la regularite correspondant a l'ecrasement du triangle
				! centre-arete
        xx = sqrt(long/abs(aire))
        IF (xx>reg_maillage) reg_maillage = xx
				xx = sqrt( ( (xs(i2)-xcv(i))**2+(ys(i2)-ycv(i))**2 ) / abs(aire))
        IF (xx>reg_maillage) reg_maillage = xx
				xx = sqrt( ( (xs(i1)-xcv(i))**2+(ys(i1)-ycv(i))**2 ) / abs(aire))
        IF (xx>reg_maillage) reg_maillage = xx

        ! on rajoute l'aire du triangle
        mcv(i) = mcv(i) + aire*.5D0
        i1 = i2
        a11 = a21
        a12 = a22

        ! Calcul du diametre = sup des distances entre 2 sommets de i,
        ! et du pas du maillage
        DO k=ptvois(i)+1,ptvois(i+1)
          diamtemp = sqrt( (xs(nusom(k))-xs(i2))**2 + (ys(nusom(k))-ys(i2))**2 )
          IF (diamtemp > diam(i) ) diam(i) = diamtemp
          IF (diam(i) > sizem ) sizem = diam(i)
        ENDDO
      ENDDO
      airetot = airetot + mcv(i)

			! Ajustement de la regularite pour prendre en compte le facteur
			!	diam(K)/dist(x_K,arete), le ratio long arete/diam maille
			! et le facteur dist(x_K,arete)/dist(x_L,arete)
			xs1 = xs(nusom(ptvois(i+1)))
			ys1 = ys(nusom(ptvois(i+1)))
			DO j=ptvois(i)+1,ptvois(i+1)
				dKs = distance_Ks(xcv(i),ycv(i),xs1,ys1,xs(nusom(j)),ys(nusom(j)))
				IF (diam(i)/dKs > reg_maillage) reg_maillage = diam(i)/dKs
				! Si le centre est de l'autre cote de l'arete, on fait exploser la
				! regularite
				IF ( -(ys(nusom(j))-ys1)*(xcv(i)-xs1) + (xs(nusom(j))-xs1)*(ycv(i)-ys1) <0D0) reg_maillage=1D10

				long = sqrt( (xs(nusom(j))-xs1)**2 + (ys(nusom(j))-ys1)**2 )
				IF (diam(i)/long > reg_maillage) reg_maillage = diam(i)/long
			
				IF (nuvois(j)>0) THEN
					dLs = distance_Ks(xcv(nuvois(j)),ycv(nuvois(j)),xs1,ys1,xs(nusom(j)),ys(nusom(j)))
					IF (dKs/dLs > reg_maillage) reg_maillage = dKs/dLs
				ENDIF
			
				xs1 = xs(nusom(j))
				ys1 = ys(nusom(j))
			ENDDO

    ENDDO

		! siezm = max des longueurs d'arete
    !		sizem = 0D0
    !		DO i=1, ncv
    !			i2 = nusom(ptvois(i+1))
    !			DO k=ptvois(i)+1,ptvois(i+1)
    !				i1 = nusom(k)
    !				sizem = max(sizem, sqrt( (xs(i1)-xs(i2))**2 + (ys(i1)-ys(i2))**2 ) )
    !				i2 = i1
    !			ENDDO
    !		ENDDO


    print '(A,X,ES22.15,5X,A,X,ES9.2,5X,A,X,ES9.2)',"\nreg maillage=",reg_maillage,&
    				"aire totale=",airetot,"pas=",sizem
    !		print *,"maxlongar",maxlongar

  END SUBROUTINE cal_mcv_reg



!**********************************************************!
!**********************************************************!
!----------------- Routines annexes -----------------------!
!----------------------------------------------------------!

  SUBROUTINE tri(numer,long)
    !*****************************************!
    ! Tri un tableau
    !
    !*****************************************!
    integer :: long, i, ip, j, k
    integer,dimension(long)::numer

    do j=long,2,-1
       do i=1,j-1
          ip = i+1
          if (numer(i)>numer(ip)) then
             ! permutation
             k = numer(ip)
             numer(ip) = numer(i)
             numer(i) = k

          endif
       end do
    end do
    END SUBROUTINE tri

  !--------------------------------------------------------!

  !-----------------------------------------!
  ! SKIP kn LINES IN THE READING OF FILE id !
  !-----------------------------------------!
  SUBROUTINE jump(kn,id)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: kn,id
    INTEGER :: ij
    CHARACTER(LEN=1) :: ab

    DO ij=1,kn
       READ(id,*) ab
    END DO

    END SUBROUTINE jump

!----------------------------------------------------!
!*************** distance_Ks ************************!

  FUNCTION distance_Ks(xm,ym,xs1,ys1,xs0,ys0)
    !******************************************************!
    ! Calcule la distance entre (xm,ym) et la droite generee
    ! par (xs1,ys1) et (xs0,ys0)
    !
    !******************************************************!
    DOUBLE PRECISION :: xm, ym, xs1, ys1, xs0, ys0

    DOUBLE PRECISION :: distance_Ks
    DOUBLE PRECISION :: l_arete, nx, ny

    ! Vecteur normal a la droite
    l_arete = sqrt( (xs1-xs0)**2 + (ys1-ys0)**2 )
    nx = - (ys1-ys0)/l_arete
    ny = (xs1-xs0)/l_arete

    ! Distance: valeur absolue du produit scalaire "pt de la droite->(xm,ym)" par le vecteur normal
    distance_Ks = abs( (xm - xs0)*nx + (ym - ys0)*ny )

  END FUNCTION distance_Ks




END MODULE mod_maillage

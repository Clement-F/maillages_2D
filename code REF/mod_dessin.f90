MODULE mod_dessin
  !
  !** Module de dessin de la grille et de la solution
  !

	USE mod_variables

  IMPLICIT NONE
  
  DOUBLE PRECISION, PARAMETER :: echelle=8D3, decalage_x=0D3, decalage_y=-1.6D4

  CONTAINS 


!---------------------------------------------------------------------!

  SUBROUTINE dessinergrille(type_sorties,fichier_base,detail_grille)
  !******** Dessiner la grille (pour xfig ou paraview) ************
  !** Sens des variables : 
  !       type_sorties = -1 pour aucun, 0 pour tous, 1 pour xfig, 2 pour vtk
  !       fichier_base = nom du fichier dans lequel mettre la grille, sans son extension
  !       detail_grille = "voirnumeros" si on veut voir les numeros
  !**************************************

    INTEGER :: type_sorties
    CHARACTER(LEN=*) :: fichier_base, detail_grille

    CHARACTER(LEN=20) :: fichier

    IF (type_sorties==0 .OR. type_sorties==1) THEN
      ! xfig
      WRITE(fichier,"(A,A)") fichier_base,".fig"
      CALL dessinergrille_xfig(fichier,detail_grille)
    ENDIF
    IF (type_sorties==0 .OR. type_sorties==2) THEN
      ! vtk
      WRITE(fichier,"(A,A)") fichier_base,".vtk"
      CALL dessinergrille_paraview(fichier,detail_grille)
    ENDIF

  END SUBROUTINE dessinergrille

!---------------------------------------------------------------------!

  SUBROUTINE dessinersolution(type_sorties,fichier_base,sol,solsom)
  !******** Creation du fichier xfig qui represente la solution ********
  !** Sens des variables : 
  !       type_sorties = -1 pour aucun, 0 pour tous, 1 pour xfig, 2 pour vtk
  !       fichier_base = nom du fichier dans lequel mettre la solution, sans son extension
  !       sol = solution aux centres de mailles
  !				solsom = solution aux sommets
  !**************************************

    INTEGER :: type_sorties
    CHARACTER(LEN=*) :: fichier_base
    DOUBLE PRECISION, DIMENSION(:) :: sol, solsom

    CHARACTER(LEN=20) fichier

    IF (type_sorties==0 .OR. type_sorties==1) THEN
      WRITE(fichier,"(A,A)") fichier_base,".fig"
      CALL dessinersolution_xfig(fichier,sol)
    ENDIF
    IF (type_sorties==0 .OR. type_sorties==2) THEN
      WRITE(fichier,"(A,A)") fichier_base,".vtk"
      CALL dessinersolution_paraview(fichier,sol,solsom)
    ENDIF


  END SUBROUTINE dessinersolution


!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!------------------------ XFIG ---------------------------------------!
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!

  SUBROUTINE dessinergrille_xfig(fichier,detail_grille)
  !******** Creation du fichier xfig qui represente la grille ************
  !** Sens des variables : 
  !         fichier = nom du fichier dans lequel mettre la grille
  !         detail_grille = "voirnumeros" si on veut voir les numeros
  !  Autres: voir mod_grille.f90
	!
	! On ramene a une echelle 1 dans la plus grande dimension
  !**************************************

    CHARACTER(LEN=*) :: fichier, detail_grille

    INTEGER :: i, j, i0, i1
    DOUBLE PRECISION :: x, y
    CHARACTER(LEN=10) :: numero
		DOUBLE PRECISION :: scaling
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xcv_s, ycv_s, xs_s, ys_s

		! calcul de la plus grande dimension, par laquelle on va mettre a l'echelle
		scaling = 0D0
		DO i=1,ncv
			DO j=i+1,ncv
				scaling = max(scaling,abs(xcv(i)-xcv(j))+abs(ycv(i)-ycv(j)))
			ENDDO
		ENDDO
		ALLOCATE(xcv_s(ncv),ycv_s(ncv),xs_s(nbsom),ys_s(nbsom))
		DO i=1,ncv
			xcv_s(i) = xcv(i) / scaling
			ycv_s(i) = ycv(i) / scaling
		ENDDO

		DO i=1,nbsom
			xs_s(i) = xs(i) / scaling
			ys_s(i) = ys(i) / scaling
		ENDDO

    OPEN(unit=21, File=fichier, status="unknown")
    WRITE(21,"(A)")"#FIG 3.2"
    WRITE(21,"(A)")"Landscape"
    WRITE(21,"(A)")"Center"
    WRITE(21,"(A)")"Metric"
    WRITE(21,"(A)")"A4"
    WRITE(21,"(A)")"100.00"
    WRITE(21,"(A)")"Single"
    WRITE(21,"(A)")"-2"
    WRITE(21,"(A)")"1200 2"

    DO i=1, ncv
      IF (detail_grille == "voirnumeros") THEN
        ! On ecrit le numero de la maille
        WRITE(numero,"(I5)") i
        CALL positionner_texte(21,numero,xcv_s(i),ycv_s(i),"noir")
			ELSEIF (detail_grille == "voircentres") THEN
				! On indique la position du centre
				CALL positionner_point(21,xcv_s(i),ycv_s(i),"noir")
      ENDIF
      i0 = nusom(ptvois(i+1))
      DO j=ptvois(i)+1,ptvois(i+1)
        i1=nusom(j)
        ! On trace l'arete
        CALL line(21,xs_s(i0),ys_s(i0),xs_s(i1),ys_s(i1))
        IF (detail_grille == "voirnumeros") THEN
          ! On ecrit son numero au centre de l'arete
          x = (xs_s(i0)+xs_s(i1))/2
          y = (ys_s(i0)+ys_s(i1))/2
          WRITE(numero,"(I5)") numar(j)
          CALL positionner_texte(21,numero,x,y,"rouge")
        ENDIF
        i0 = i1
      ENDDO
    ENDDO

    ! Ecrire les numeros des sommets
    IF (detail_grille == "voirnumeros") THEN
      DO i0=1, nbsom
        WRITE(numero,"(I5)") i0
        CALL positionner_texte(21,numero,xs_s(i0),ys_s(i0),"bleu")
      ENDDO
    ENDIF

    CLOSE(21)

  END SUBROUTINE dessinergrille_xfig

!---------------------------------------------------------------------!

  SUBROUTINE dessinersolution_xfig(fichier,sol)
  !******** Creation du fichier xfig qui represente la solution ********
  !** Sens des variables : fichier = nom du fichier dans lequel
  !                                      mettre la solution
  !                       sol = solution
	!
	! On ramene a une echelle 1 dans la plus grande dimension
  !**************************************

    CHARACTER(LEN=*) :: fichier
    DOUBLE PRECISION, DIMENSION(:) :: sol

    INTEGER :: icou, i, j, k, maxna
    DOUBLE PRECISION :: solmin, solmax
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)  :: xss,yss
		DOUBLE PRECISION :: scaling
		DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xcv_s, ycv_s, xs_s, ys_s

		! calcul de la plus grande dimension, par laquelle on va mettre a l'echelle
		scaling = 0D0
		DO i=1,ncv
			DO j=i+1,ncv
				scaling = max(scaling,abs(xcv(i)-xcv(j))+abs(ycv(i)-ycv(j)))
			ENDDO
		ENDDO
		ALLOCATE(xcv_s(ncv),ycv_s(ncv),xs_s(nbsom),ys_s(nbsom))
		DO i=1,ncv
			xcv_s(i) = xcv(i) / scaling
			ycv_s(i) = ycv(i) / scaling
		ENDDO
		DO i=1,nbsom
			xs_s(i) = xs(i) / scaling
			ys_s(i) = ys(i) / scaling
		ENDDO

    ! Calcul du nombre max d'aretes pour une maille donnee afin d'allouer xss et yss
    maxna = 0
    DO i=1, ncv
     IF (maxna < ptvois(i+1)-ptvois(i)) maxna = ptvois(i+1)-ptvois(i)
    ENDDO
    ALLOCATE(xss(maxna),yss(maxna))

    !** Calcul des max et min de la solution
    solmin = sol(1)
    solmax = sol(1)
    DO i = 2,ncv
       IF (sol(i)<solmin) solmin = sol(i)
       IF (sol(i)>solmax) solmax = sol(i)
    ENDDO

    OPEN(unit=22, File=fichier, status="unknown")
    WRITE(22,"(A)")"#FIG 3.2"
    WRITE(22,"(A)")"Landscape"
    WRITE(22,"(A)")"Center"
    WRITE(22,"(A)")"Metric"
    WRITE(22,"(A)")"A4"
    WRITE(22,"(A)")"100.00"
    WRITE(22,"(A)")"Single"
    WRITE(22,"(A)")"-2"
    WRITE(22,"(A)")"1200 2"

    icou = 32

!     DO i=0,127
!       write(22,"(A,I3,A,Z2.2,Z2.2,Z2.2)")"0 ",icou, " #",2*i,2*i,2*i
!       icou = icou + 1
!     ENDDO
        !!!!!!!!!!!!!
        ! VERSION COULEUR: commenter les 4 lignes precedentes et decommenter ce qui suit
        DO i=0,31
          write(22,"(A,I3,A,Z2.2,Z2.2,Z2.2)")"0 ",icou, " #",0,0,8*i
          icou = icou + 1
        ENDDO
        DO i=0,31
          write(22,"(A,I3,A,Z2.2,Z2.2,Z2.2)")"0 ",icou, " #",0,8*i,255
          icou = icou + 1
        ENDDO
        DO i=0,31
          write(22,"(A,I3,A,Z2.2,Z2.2,Z2.2)")"0 ",icou, " #",8*i,255,255-8*i
          icou = icou + 1
        ENDDO
        DO i=0,31
          write(22,"(A,I3,A,Z2.2,Z2.2,Z2.2)")"0 ",icou, " #",255,255-8*i,0
          icou = icou + 1
        ENDDO
        ! FIN VERSION COULEUR
        !!!!!!!!!!!!!

    DO i=1, ncv
      k = 0
      DO j=ptvois(i)+1,ptvois(i+1)
          k = k+1
          xss(k) = xs_s(nusom(j))
          yss(k) = ys_s(nusom(j))
      END DO
      CALL polyl(22,k,xss,yss,(sol(i)-solmin)/(solmax-solmin))
       !          call line(xcv(i),ycv(i),xcv(i)+vitech*vx(i),ycv(i)+vitech*vy(i))
    END DO
    CLOSE(22)

  END SUBROUTINE dessinersolution_xfig

!---------------------------------------------------------------------!
  !****************** line *******************************************!

  SUBROUTINE line(unit,x1,y1,x2,y2)
  !******** Tracage d'une ligne en commandes xfig *****************
  !** Sens des variables : unit = numero du fichier dans lequel on veut ecrire
  !                       x1,y1,x2,y2 = on trace une ligne du point (x1,y1)
  !                                    au point (x2,y2).
  !******************************************************

    INTEGER :: unit
    DOUBLE PRECISION :: x1, y1, x2, y2

    INTEGER :: i1, j1, i2, j2

    i1 = x1 * echelle + decalage_x
    j1 = (3 - y1) * echelle + decalage_y
    i2 = x2 * echelle + decalage_x
    j2 = (3 - y2) * echelle + decalage_y

    WRITE(unit,*)"2 1 0 1 0 7 100 0 -1 0.000 0 0 -1 0 0 2"
    WRITE(unit,*)i1, j1, i2, j2

  END SUBROUTINE line

!----------------------------------------------------------------------!
  !****************** positionner_texte *******************************!

  SUBROUTINE positionner_texte(unit,texte,x,y,couleur)
  !** Positionne dans le fichier .fig le texte aux positions x et y ***!
  ! Exemple de ligne:
  !  4 0 0 50 -1 0 8 0.0000 4 90 75 720 6930 5\001
  ! Explications:
  !     7e position (le 8): taille de caractere
  !     12e et 13e positions (720 et 6930): positions X et Y a partir du
  !                 coin haut gauche
  !     derniere position (5\001): le texte ('5' ici)
  !
  ! couleur = le code xfig est 0 pour noir, 1 pour bleu, 4 pour rouge
  !*******************************************************************!

    INTEGER :: unit
    CHARACTER(LEN=*) :: texte, couleur
    DOUBLE PRECISION :: x, y

    CHARACTER(LEN=60) :: texte_complet
    INTEGER :: i, j, c
		DOUBLE PRECISION :: decalage_local

    IF (couleur == "bleu") THEN
      c = 1
    ELSEIF (couleur == "rouge") THEN
      c = 4
    ELSE
      c = 0
    ENDIF

    decalage_local = - 20D0

    i = x * echelle + decalage_x + decalage_local
    j = (3 - y) * echelle + decalage_y

    WRITE(texte_complet,"(A,X,I1,X,A,X,I5,X,I5,X,A,A)") "4 0",c,"50 -1 0 8 0.0000 4 90 75",i,j,TRIM(ADJUSTL(texte)),"\\001"
    WRITE(unit,*) TRIM(texte_complet)

  END SUBROUTINE positionner_texte


!----------------------------------------------------------------------!
  !****************** positionner_point *******************************!

  SUBROUTINE positionner_point(unit,x,y,couleur)
  !** Positionne dans le fichier .fig un point aux positions x et y ***!
  ! Exemple de ligne:
	! 1 3 0 1 0 0 50 -1 20 0.000 1 0.0000 2070 2340 30 30 2070 2340 0 0
  ! Explications:
  !     13e et 17e positions (le 2070): position x
  !     14e et 18e positions (2340): positions y
  !     15e et 16e : rayon du point
  !
  ! couleur = le code xfig est 0 pour noir, 1 pour bleu, 4 pour rouge
  !*******************************************************************!

    INTEGER :: unit
    CHARACTER(LEN=*) :: couleur
    DOUBLE PRECISION :: x, y

    DOUBLE PRECISION :: decalage_local
    CHARACTER(LEN=100) :: texte_complet
    INTEGER :: i, j, c

    IF (couleur == "bleu") THEN
      c = 1
    ELSEIF (couleur == "rouge") THEN
      c = 4
    ELSE
      c = 0
    ENDIF

    decalage_local = -20D0

    i = x * echelle + decalage_x + decalage_local
    j = (3 - y) * echelle + decalage_y
    
    WRITE(texte_complet,"(A,X,I5,X,I5,X,A,X,I5,X,I5,X,A)") "1 3 0 1 0 0 50 -1 20 0.000 1 0.0000",i,j,"8 8",i,j,"0 0"
    WRITE(unit,*) TRIM(texte_complet)

  END SUBROUTINE positionner_point


!-----------------------------------------------------------------------!
  !***************** polyl *********************************************!

  SUBROUTINE polyl(unit,nbp,xss,yss,val)
  !******** Code xfig pour colorier une maille selon la valeur de la solution *****
  !** Sens des variables : unit = numero du fichier dans lequel ecrire
  !                       nbp = nombre de sommets de la maille
  !                       xss et yss = coordonnees des sommets de la maille
  !                       val = valeur de la solution sur la maille,
  !                           ramenee entre 0 et 1 via ses min et max globaux
  !**************************************************************************

    INTEGER :: unit, nbp
    DOUBLE PRECISION, DIMENSION(nbp) :: xss,yss
    DOUBLE PRECISION :: val

    INTEGER :: cou, i, k
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ii,jj

    ALLOCATE(ii(nbp),jj(nbp))

    cou = INT(128D0*val)
    IF (cou >= 128) cou = 127
    cou = cou + 32

    DO i=1,nbp
       ii(i) = xss(i) * echelle + decalage_x
       jj(i) = (3 - yss(i)) * echelle + decalage_y
    END DO

    WRITE(unit,*)"2 3 0 1 ",cou, cou," 50 0 20 0.000 0 0 -1 0 0", nbp+1
		! Pour dessiner la grille en noir: commenter la ligne precedente et decommenter la suivante
!    WRITE(unit,*)"2 3 0 1 0",cou," 50 0 20 0.000 0 0 -1 0 0", nbp+1
    WRITE(unit,*) (ii(k),jj(k),k=1,nbp),ii(1),jj(1)

  END SUBROUTINE polyl



!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!------------------------  PARAVIEW ----------------------------------!
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!

  SUBROUTINE dessinergrille_paraview(fichier,detail_grille)
  !******** Creation du fichier paraview qui represente la grille ************
  !** Sens des variables : 
  !       fichier = nom du fichier dans lequel mettre la grille
  !       detail_grille = "voirnumeros" si on veut voir les numeros
  !  Autres: voir mod_grille.f90
  !**************************************

    CHARACTER(LEN=*) :: fichier, detail_grille

    INTEGER :: i, j

    OPEN(unit=21, File=fichier, status="unknown")
    WRITE(21,"(A)") "# vtk DataFile Version 2.0"
    WRITE(21,"(A)") "Grille"
    WRITE(21,"(A)") "ASCII"

    WRITE(21,"(A)") "DATASET POLYDATA"
    WRITE(21,"(A,I10,X,A)") "POINTS",nbsom,"float"
    DO i=1, nbsom
      WRITE(21,"(E12.5,E12.5,E12.5)") xs(i),ys(i),0D0
    ENDDO
    WRITE(21,"(A,I10,X,I10)") "POLYGONS",ncv,ptvois(ncv+1)+ncv
    DO i=1, ncv
      WRITE(21,*) ptvois(i+1)-ptvois(i),(nusom(j)-1,j=ptvois(i)+1,ptvois(i+1))
    ENDDO

    CLOSE(21)

  END SUBROUTINE dessinergrille_paraview


!----------------------------------------------------------------------------!

  SUBROUTINE dessinersolution_paraview(fichier,sol,solsom)
  !*************** Creation du fichier paraview pour voir la solution *******!
  !
  !**************************************************************************!

    CHARACTER(LEN=*) :: fichier
    DOUBLE PRECISION, DIMENSION(:) :: sol, solsom
    
    INTEGER :: i, j

    OPEN(unit=22, File=fichier, status="unknown")
    WRITE(22,"(A)") "# vtk DataFile Version 2.0"
    WRITE(22,"(A)") "Solution"
    WRITE(22,"(A)") "ASCII"

    WRITE(22,"(A)") "DATASET POLYDATA"
    WRITE(22,"(A,I10,X,A)") "POINTS",nbsom,"float"
    DO i=1, nbsom
    	! Selon que l'on veut un dessin plat ou avec relief, decommenter la
    	! premiere ou deuxieme ligne suivante
!      WRITE(22,"(E12.5,E12.5,E12.5)") xs(i),ys(i),0D0
      WRITE(22,"(E12.5,E12.5,E12.5)") xs(i),ys(i),solsom(i)
    ENDDO
    WRITE(22,"(A,I10,X,I10)") "POLYGONS",ncv,ptvois(ncv+1)+ncv
    DO i=1, ncv
      WRITE(22,*) ptvois(i+1)-ptvois(i),(nusom(j)-1,j=ptvois(i)+1,ptvois(i+1))
    ENDDO
    WRITE(22,"(A,I10)") "CELL_DATA",ncv
    WRITE(22,"(A)") "SCALARS sol float 1"
    WRITE(22,"(A)") "LOOKUP_TABLE default"
    DO i=1, ncv
      WRITE(22,"(E12.5)") sol(i)
    ENDDO

    CLOSE(22)



  END SUBROUTINE dessinersolution_paraview
  
  
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
!-------------------- Autres: binaire fortran contenant la solution ---------!

!----------------------------------------------------------------------------!
  !************************* fichier_solution *******************************!

  SUBROUTINE ecrire_fichier_solution(nomfichier,p)
  !***************************************************************************
  ! Ecrit/lit le maillage et la solution dans le fichier "nomfichier"
  !
  !**************************************************************************
    CHARACTER(LEN=*) :: nomfichier
    DOUBLE PRECISION, DIMENSION(:) :: p

    INTEGER :: i


    OPEN(unit=66, file=nomfichier, form="unformatted", status="unknown")
    write(66) ncv, nbsom, nptvois
    write(66) (ptvois(i),i=1,ncv+1)
    write(66) (nusom(i),i=1,nptvois)
    write(66) (nuvois(i),i=1,nptvois)
    write(66) (xs(i),ys(i),i=1,nbsom)
    write(66) (xcv(i),ycv(i),i=1,ncv)
    write(66) (mcv(i),i=1,ncv)
    write(66) (p(i),i=1,ncv)
    CLOSE(66)


!   Les codes pour lire le fichier:
!
!     OPEN(unit=66, file=nomfichier, form="unformatted", status="unknown")
!     read(66) ncv, nbsom, nptvois
!
!     ALLOCATE(ptvois(ncv+1),nusom(nptvois),nuvois(nptvois),xs(nbsom),ys(nbsom),xcv(ncv),ycv(ncv),mcv(ncv),p(ncv))
!
!     read(66) (ptvois(i),i=1,ncv+1)
!     read(66) (nusom(i),i=1,nptvois)
!     read(66) (nuvois(i),i=1,nptvois)
!     read(66) (xs(i),ys(i),i=1,nbsom)
!     read(66) (xcv(i),ycv(i),i=1,ncv)
!     read(66) (mcv(i),i=1,ncv)
!     read(66) (p(i),i=1,ncv)
!
!     CLOSE(66)


  END SUBROUTINE ecrire_fichier_solution


END MODULE mod_dessin

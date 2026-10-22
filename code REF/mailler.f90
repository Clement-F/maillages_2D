 PROGRAM mailler
  !
  ! Mailleur 2D
  !
  !--------------------------------------------------------------------------!
  !     nbsom, ncv, nptvois, ptvois, nusom, nuvois, numar,
  !            nar, xs, ys   --> module mod_maillage.f90
  !
  !     sommets0 --> fonction "divtri/divquad" dans module mod_grille.f90
  !
  !     ndiv  = nombre de divisions des aretes du maillage initial
  !
  !     (xcv(i),ycv(i)) = centre de la maille i
  !     mcv(i) = aire de la maille i
  !     diam(i) = diametre de la maille i
  !     sizem = pas du maillage (sup des diam)
  !
  !--------------------------------------------------------------------------!

	USE mod_variables
  USE mod_maillage
  USE mod_dessin
  
  IMPLICIT NONE

  !------------------ declaration des variables -----------------------------!
  !
  ! Certaines variables "globales" sont definies dans mod_maillage

  ! Pour voir ou non les numeros de maille, aretes...
  CHARACTER(LEN=15) :: detail_grille
	
  ! Mesurer le temps d'execution
  REAL :: Tps1, Tps2, t1, t2

	! A fins de tests
	INTEGER :: ii
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: xcv2, ycv2
	DOUBLE PRECISION :: ecart, ecart_temp, measure, ecart_x, ecart_y

	DOUBLE PRECISION :: trans_a, trans_b

  !-------------------- programme -------------------------------------------!


  CALL CPU_TIME(Tps1)

  !$$$$$$$ GRILLE $$$$$$$$!
  !
  ! Lecture et construction de la grille
  CALL lec_cons_grille

  ! Calcul des pointeurs et numerotation des aretes/mailles
  CALL calcul_pointeurs(1)

	! Re-scale a mesh on the unit square
	trans_a=0.6D0
	trans_b=0.05D0
	IF ( ((typepoly == -9) .OR. (typepoly == -30)) .AND. nutest == 70) THEN
!		xs(:) = 0*ys(:)/4D0 + (1D0-ys(:)/1.05D0)*xs(:)
!		xcv(:) = 0*ycv(:)/4D0 + (1D0-ycv(:)/1.05D0)*xcv(:)

		ys(:) = ys(:) * (2.5D0-ys(:)) / 1.5D0
		xs(:) = (1D0-ys(:))*trans_a + (ys(:)*trans_b + (1D0-ys(:))*(1D0-trans_a))*xs(:)
!		xcv(:) = (1D0-ycv(:))*trans_a + (ycv(:)*trans_b + (1D0-ycv(:))*(1D0-trans_a))*xcv(:)

		CALL calcul_centres(ncv,ptvois,nusom,xs,ys,xcv,ycv,0)

	ENDIF


  ! On dessine la grille, avec les numeros des mailles et aretes s'il y
  ! a moins de 50 mailles
  IF (ncv <= 10) THEN
    detail_grille = "voirnumeros"
  ELSEIF (ncv <= 20) THEN
  	detail_grille = "voircentres"
  ELSE
    detail_grille = "pasnumeros"
  ENDIF
  CALL dessinergrille(type_sorties,"grid",detail_grille)

  ! Calcul de l'aire des mailles et du facteur de regularite
  CALL cal_mcv_reg

	! Si nutest=-1, on cree juste le fichier maillage au format "typ1"
	IF (nutest==-1) THEN
		CALL acces_maillage_typ1("ecrire","mesh.typ1")
		print *,"Fichier maillage typ1 cree"
		STOP
	ENDIF

	! Si nutest=-2, on cree juste le fichier maillage au format "typ2"
	IF (nutest==-2) THEN
		CALL acces_maillage_typ2("ecrire","mesh.typ2")
		print *,"Fichier maillage typ2 cree"
		STOP
	ENDIF

	! Si nutest=-3, on cree juste le fichier maillage au format "dt"
	IF (nutest==-3) THEN
		CALL acces_fichier_maillage("ecrire","mesh.dt",nbsom,ncv,xs,ys,ptvois,nusom,nuvois)
		print *,"Fichier maillage dt cree"
		STOP
	ENDIF

  CALL CPU_TIME(Tps2)
  print *,'\ntemps execution',Tps2-Tps1

  !--------------------------------------------------------------------------!


END PROGRAM mailler

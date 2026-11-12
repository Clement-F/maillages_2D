!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULE THAT READ THE MESH AND INTIALIZE THE GRID ARRAYS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE read_mesh
  USE mod_int_proj
  USE mesh_cart
  IMPLICIT NONE

  INTEGER,PRIVATE :: nb_bdr_edges=0
  INTEGER,DIMENSION(:,:),ALLOCATABLE,PRIVATE :: bdr_edge_nodes
  INTEGER,DIMENSION(:),ALLOCATABLE,PRIVATE :: bdr_edge_ref,bdr_edge_used

CONTAINS

  !------------------------------------------------------------!
  ! READ OF THE MESH (.msh, .amd, polygonal .plt OR cartesian) !
  !------------------------------------------------------------!
  SUBROUTINE init_mesh
    IMPLICIT NONE

    INTEGER :: ncell,npt,pt,nb_cell_sedov,ref
    INTEGER,DIMENSION(40) :: tri_sedov_temp
    REAL(prec),DIMENSION(2) :: XXp

    nb_face_bdr_gmsh=0
    IF (r_name_mesh(1:8)=="tri_cart") THEN
       IF (TRIM(r_name_mesh)=="tri_cart1") THEN
          CALL init_mesh_cart_tri1
       ELSE IF (TRIM(r_name_mesh)=="tri_cart2") THEN
          CALL init_mesh_cart_tri2
       ELSE
          CALL init_mesh_cart_tri3
       END IF

       IF (TRIM(test_case)=="ODD1") CALL odd_even_mesh_cart_upadte
    ELSE IF (r_name_mesh(1:9)=="quad_cart") THEN
       CALL init_mesh_cart_quad
       
       IF (TRIM(test_case)=="ODD1") CALL odd_even_mesh_cart_upadte
    ELSE
       IF (mesh_type=="msh" .OR. mesh_type=="amd" .OR. mesh_type=="plt") THEN
          CALL read_msh_amd
       ELSE
          WRITE(*,*) "MESH TYPE NOT SUPPORTED",name_mesh
          WRITE(*,*) "MAKE USE OF .msh, .amd OR TAGGED POLYGONAL .plt MESHES"
          STOP
       END IF
    END IF

    DO ncell=1,nb_cell
       DO npt=1,mesh(ncell)%nb_face
          IF (mesh(ncell)%neigh(npt).GE.0) CYCLE
          ref=-mesh(ncell)%neigh(npt)
          IF (ref.LT.1 .OR. ref.GT.nbr_bdr) THEN
             WRITE(*,*) "INVALID BOUNDARY REFERENCE",ncell,npt,ref,nbr_bdr
             STOP
          END IF
          IF (TRIM(bdr_cond(ref))=="peri" .OR. TRIM(bdr_cond(ref))=="per") THEN
             WRITE(*,*) "UNMATCHED PERIODIC FACE",ncell,npt,ref
             STOP
          END IF
       END DO
    END DO

    !! Reorganize the node neighboring cell set !!
    CALL reorg_node_neigh
    
    decal_sol=10._prec**prec
    nb_cell_sedov=0
    tri_sedov_temp=0
    vol_sedov=0._prec
    DO ncell=1,nb_cell
       !CALL reorg_mesh(ncell)
       CALL init_neighborhood_ext(ncell)
       
       mesh(ncell)%Barycenter=0._prec
       DO npt=1,mesh(ncell)%nb_face
          pt=mesh(ncell)%node(npt)
          decal_sol=MIN(decal_sol,DOT_PRODUCT(Vit_adv,point(pt)%Coord))
          mesh(ncell)%Barycenter=mesh(ncell)%Barycenter+point(pt)%Coord
       END DO
       mesh(ncell)%Barycenter=mesh(ncell)%Barycenter/REAL(mesh(ncell)%nb_face,prec)

       IF (test_case(1:5)=="SEDOV") THEN
          DO npt=1,mesh(ncell)%nb_face
             pt=mesh(ncell)%node(npt)
             XXp=point(pt)%Coord
             IF (DOT_PRODUCT(XXp,XXp).LE.eps0) THEN
                nb_cell_sedov=nb_cell_sedov+1
                tri_sedov_temp(nb_cell_sedov)=ncell
                vol_sedov=vol_sedov+mesh(ncell)%volume
                length_sedov=mesh(ncell)%length(npt)*SQRT(1._prec+SQRT(2._prec)/2._prec)
             END IF
          END DO
       END IF
    END DO

    IF (test_case(1:5)=="SEDOV") THEN
       ALLOCATE(tri_sedov(nb_cell_sedov))
       tri_sedov=tri_sedov_temp(1:nb_cell_sedov)
       
       ene_sedov=ene_sedov_0
       IF (xdeb.LT.-eps0) ene_sedov=ene_sedov*2._prec
       IF (ydeb.LT.-eps0) ene_sedov=ene_sedov*2._prec
    END IF

    WRITE(*,*) " "
    WRITE(*,*) "MESHS AND NODES ALLOCATION DONE"
    WRITE(*,*) "THE TOTAL NUMBER OF CELLS IS",nb_cell

  END SUBROUTINE init_mesh



  !--------------------------------!
  ! ROTATE MESH OF ANGLE rot_angle !
  !--------------------------------!
  SUBROUTINE rot_mesh
    IMPLICIT NONE

    INTEGER :: pt,ncell,nbase,nbb,npt
    REAL(prec),DIMENSION(2) :: XXp,XX0
    REAL(prec),DIMENSION(2,2) :: Mat_modif

    XX0=(/xdeb,ydeb/)

    Mat_modif=0._prec
    Mat_modif(1,1)=COS(rot_angle)
    Mat_modif(1,2)=-SIN(rot_angle)
    Mat_modif(2,1)=SIN(rot_angle)
    Mat_modif(2,2)=COS(rot_angle)

    DO pt=1,nb_pts
       XXp=point(pt)%Coord
       point(pt)%Coord=XX0+MATMUL(Mat_modif,XXp-XX0)
    END DO

    decal_sol=10._prec**prec
    DO ncell=1,nb_cell
       XXp=mesh(ncell)%Barycenter
       mesh(ncell)%Barycenter=XX0+MATMUL(Mat_modif,XXp-XX0)
       DO npt=1,mesh(ncell)%nb_face
          mesh(ncell)%Normal(npt,:)=MATMUL(Mat_modif,mesh(ncell)%Normal(npt,:))
          pt=mesh(ncell)%node(npt)
          decal_sol=MIN(decal_sol,DOT_PRODUCT(Vit_adv,point(pt)%Coord))
       END DO
    END DO

    WRITE(*,*) "ROTATION OF THE MESH DONE"
    WRITE(*,*) " "

  END SUBROUTINE rot_mesh

  

  !-----------------------------------------!
  ! COMMON GEOMETRY FOR .msh, .amd AND .plt !
  !-----------------------------------------!
  SUBROUTINE read_msh_amd
    IMPLICIT NONE

    INTEGER :: ncell,nvois,vois
    INTEGER :: pt1,pt2,pt3,npt,npt2
    INTEGER,DIMENSION(:),POINTER ::cell_pts
    REAL(prec) :: long,vol,xmin,xmax,ymin,ymax,llg,decalx,ll
    REAL(prec),DIMENSION(2) :: XX1,XX2,XXm,nn,tt
    REAL(prec),DIMENSION(:,:),ALLOCATABLE :: XX_pts

    !! read_mesh and read_amd to allocate and init  !!
    !! the point coordinates and the mesh nodes !!
    IF (mesh_type=="plt") THEN
       CALL read_poly_mesh
    ELSE
       CALL read_msh
    END IF

    decalx=ABS(xdeb-xdeb0)+ABS(xfin-xfin0)+ABS(ydeb-ydeb0)+ABS(yfin-yfin0)

    !! Allocate the node neighboring cells !!
    !!       and mesh deformation          !!
    DO pt1=1,nb_pts
       nvois=point(pt1)%nb_neigh
       ALLOCATE(point(pt1)%neigh(nvois))
       point(pt1)%neigh=0
       point(pt1)%nb_per_neigh=0
       IF (mesh_type/="plt" .AND. point(pt1)%ref/=0) ALLOCATE(point(pt1)%per_neigh(1))

       IF (opt_decal==1 .AND. decalx.GT.epsl) CALL decal_mesh(pt1)
    END DO

    min_size=10._prec**prec
    xmin=10._prec**prec
    ymin=10._prec**prec
    xmax=-10._prec**prec
    ymax=-10._prec**prec

    !! Init the cell face normals and lenght, the !!
    !!       cell volume as well as the nodes     !!
    !!             neighboring cells              !!
    DO ncell=1,nb_cell
       mesh(ncell)%neigh=0

       ALLOCATE(cell_pts(mesh(ncell)%nb_face),XX_pts(mesh(ncell)%nb_face,2))
       DO npt=1,mesh(ncell)%nb_face
          cell_pts(npt)=mesh(ncell)%node(npt)
          XX_pts(npt,:)=point(cell_pts(npt))%Coord
       END DO

       vol=0._prec
       DO npt=1,mesh(ncell)%nb_face
          pt1=mesh(ncell)%node(npt)
          IF (npt==mesh(ncell)%nb_face) THEN
             pt2=mesh(ncell)%node(1)
          ELSE
             pt2=mesh(ncell)%node(npt+1)
          END IF
          XX1=point(pt1)%Coord
          XX2=point(pt2)%Coord
          
          XXm=0.5_prec*(XX1+XX2)
          tt=XX2-XX1          
          nn=(/tt(2),-tt(1)/)
          
          vol=vol+DOT_PRODUCT(XXm,nn)
       END DO
       vol=0.5_prec*vol

       mesh(ncell)%volume=ABS(vol)

       IF (.NOT.(ABS(vol).GT.0._prec .AND. ABS(vol).LE.HUGE(vol))) THEN
          WRITE(*,*) "ZERO OR INVALID CELL AREA",ncell,vol
          STOP
       END IF

       ! Sort the nodes in the counter-clockwise !
       IF (vol.LE.0._prec) THEN
          DO npt=1,mesh(ncell)%nb_face
             mesh(ncell)%node(npt)=cell_pts(mesh(ncell)%nb_face-npt+1)
          END DO
       END IF
       
       DEALLOCATE(cell_pts,XX_pts)
       
       ! Compute the normals and init the node neighboring cells!
       llg=0._prec
       DO npt=1,mesh(ncell)%nb_face
          pt1=mesh(ncell)%node(npt)
          IF (npt==mesh(ncell)%nb_face) THEN
             pt2=mesh(ncell)%node(1)
          ELSE
             pt2=mesh(ncell)%node(npt+1)
          END IF
          XX1=point(pt1)%Coord
          XX2=point(pt2)%Coord

          long=SQRT((XX2(1)-XX1(1))**2+(XX2(2)-XX1(2))**2)
!          min_size=MIN(min_size,long)
!          min_size=MIN(min_size,ABS(vol)/long)
          llg=llg+long

          xmin=MIN(xmin,XX1(1));xmax=MAX(xmax,XX1(1))
          ymin=MIN(ymin,XX1(2));ymax=MAX(ymax,XX1(2))

          mesh(ncell)%length(npt)=long
          mesh(ncell)%Normal(npt,1)=(XX2(2)-XX1(2))/long
          mesh(ncell)%Normal(npt,2)=(XX1(1)-XX2(1))/long
          
          nvois=point(pt1)%nb_neigh-COUNT(point(pt1)%neigh==0)
          point(pt1)%neigh(nvois+1)=ncell
       END DO

       min_size=MIN(min_size,2._prec*ABS(vol)/llg)

    END DO
    length_mesh=xmax-xmin
    height_mesh=ymax-ymin
       
    !! Neighbor allocation for periodic boundary counditions !!
    IF (mesh_type=="msh" .AND. ANY(bdr_cond=="peri")) CALL neigh_periodic

    !! Allocate the cell neighboring cells !!
    DO ncell=1,nb_cell

       bcl_neigh:DO npt=1,mesh(ncell)%nb_face

          IF (mesh(ncell)%neigh(npt).GT.0) CYCLE bcl_neigh
          
          pt1=mesh(ncell)%node(npt)
          IF (npt==mesh(ncell)%nb_face) THEN
             pt2=mesh(ncell)%node(1)
          ELSE
             pt2=mesh(ncell)%node(npt+1)
          END IF
          
          IF (mesh_type=="plt") THEN
             ! Two common vertices need not form a common polygon edge
             vois=find_polygon_neighbor(ncell,pt1,pt2)
          ELSE
             vois=match_array(ncell,point(pt1)%neigh,point(pt2)%neigh)
          END IF
          mesh(ncell)%neigh(npt)=vois

       END DO bcl_neigh

       IF (ANY(mesh(ncell)%neigh.LE.0)) THEN
          DO npt=1,mesh(ncell)%nb_face
             IF (mesh(ncell)%neigh(npt).LE.0) THEN
                pt1=mesh(ncell)%node(npt)
                pt2=mesh(ncell)%node(next_npt(ncell,npt))
                mesh(ncell)%neigh(npt)=-find_boundary_edge_ref(pt1,pt2)
             END IF
          END DO
       END IF

    END DO

    IF (mesh_type=="plt") THEN
       IF (ANY(bdr_edge_used/=1)) THEN
          WRITE(*,*) "PLT BOUNDARY TABLE CONTAINS AN UNUSED OR NON-BOUNDARY EDGE"
          STOP
       END IF
    END IF
    DEALLOCATE(bdr_edge_nodes,bdr_edge_ref,bdr_edge_used)
    nb_bdr_edges=0

  END SUBROUTINE read_msh_amd


  !---------------------!
  ! READ OF A gmsh MESH !
  !---------------------!
  SUBROUTINE read_msh
    IMPLICIT NONE

    INTEGER :: nelem,pt,npt,ncell,nface,ios_back,ntype
    INTEGER,DIMENSION(9) :: temp
    REAL(prec) :: tempr
    REAL(prec),DIMENSION(2) :: XX1

    OPEN(11,file=TRIM(name_mesh),status="old")
    CALL JUMP(4,11)

    READ(11,*) nb_pts
    ALLOCATE(point(nb_pts))
    
    xdeb0=10._prec**prec
    ydeb0=10._prec**prec
    xfin0=-10._prec**prec
    yfin0=-10._prec**prec

    !! Grid points allocation !!
    DO pt=1,nb_pts
        READ(11,*) temp(1),XX1,tempr

        point(pt)%Coord=XX1

        point(pt)%ref=0
        point(pt)%nb_neigh=0

        xdeb0=MIN(xdeb0,XX1(1))
        ydeb0=MIN(ydeb0,XX1(2))
        xfin0=MAX(xfin0,XX1(1))
        yfin0=MAX(yfin0,XX1(2))
    END DO
    CALL JUMP(2,11)


    !! Determine the number of tri !!
    READ(11,*) nelem

    ALLOCATE(bdr_edge_nodes(nelem,2),bdr_edge_ref(nelem),bdr_edge_used(nelem))
    nb_bdr_edges=0
    bdr_edge_used=0

    bcl_face:DO nface=1,nelem
       READ(11,*) temp(1:7)

       IF (temp(2)/=1) EXIT bcl_face

       nb_bdr_edges=nb_bdr_edges+1
       bdr_edge_nodes(nb_bdr_edges,:)=(/temp(6),temp(7)/)
       bdr_edge_ref(nb_bdr_edges)=temp(4)
       
       point(temp(6))%ref=temp(4)
    END DO bcl_face

    nb_face_bdr_gmsh=nb_bdr_edges     ! Number of bdry faces !
    nb_cell=nelem-nb_face_bdr_gmsh    ! Number of cells      !
    ALLOCATE(mesh(nb_cell))

    CLOSE(11)
    OPEN(11,file=TRIM(name_mesh),status="old")
    CALL JUMP(8+nb_pts+nb_face_bdr_gmsh,11)

    !! Mesh nodes and number of neighboring cells of the node !!
    DO ncell=1,nb_cell
       READ(11,*) temp(1:8)
       BACKSPACE(11, IOSTAT=ios_back)
       IF (ios_back /= 0) THEN
          WRITE(*,*) "ERROR IN BACKSPACE IN  SUBROUTINE read_msh"
          STOP
       END IF
       
       ntype=temp(2)+1
       mesh(ncell)%nb_face=ntype
       ALLOCATE(mesh(ncell)%node(ntype),mesh(ncell)%neigh(ntype))
       ALLOCATE(mesh(ncell)%length(ntype),mesh(ncell)%Normal(ntype,2))
       
       READ(11,*) temp(1:5+ntype)
       mesh(ncell)%ref=temp(4)
       mesh(ncell)%node(1:ntype)=temp(6:5+ntype)

       DO npt=1,ntype
          pt=mesh(ncell)%node(npt)
          point(pt)%nb_neigh=point(pt)%nb_neigh+1
       END DO
    END DO

    CLOSE(11)

  END SUBROUTINE read_msh



  !----------------------------!
  ! POLYGONAL ASCII PLT MESHES !
  !----------------------------!
  SUBROUTINE read_poly_mesh
    IMPLICIT NONE

    INTEGER :: ios,pt,ncell,npt,nv,ib,jb,a,b,ref
    CHARACTER(LEN=65536) :: line
    REAL(prec),DIMENSION(2) :: XX

    OPEN(11,FILE=TRIM(name_mesh),STATUS="old",ACTION="read",IOSTAT=ios)
    IF (ios/=0) THEN
       WRITE(*,*) "CANNOT OPEN POLYGONAL PLT MESH",TRIM(name_mesh)
       STOP 1
    END IF

    CALL read_plt_record(11,line,1)
    nb_pts=0
    READ(line,*,IOSTAT=ios) nb_pts
    IF (ios/=0 .OR. nb_pts.LT.3) THEN
       WRITE(*,*) "INVALID PLT POINT COUNT"
       STOP 1
    END IF
    ALLOCATE(point(nb_pts))
    xdeb0=HUGE(1._prec); ydeb0=HUGE(1._prec)
    xfin0=-HUGE(1._prec); yfin0=-HUGE(1._prec)

    DO pt=1,nb_pts
       CALL read_plt_record(11,line,2)
       XX=0._prec
       READ(line,*,IOSTAT=ios) XX
       IF (ios/=0 .OR. .NOT.ALL(ABS(XX).LE.HUGE(1._prec))) THEN
          WRITE(*,*) "INVALID PLT POINT COORDINATES",pt
          STOP 1
       END IF
       point(pt)%Coord=XX
       point(pt)%ref=0
       point(pt)%nb_neigh=0
       point(pt)%nb_per_neigh=0
       NULLIFY(point(pt)%neigh,point(pt)%per_neigh)
       xdeb0=MIN(xdeb0,XX(1)); xfin0=MAX(xfin0,XX(1))
       ydeb0=MIN(ydeb0,XX(2)); yfin0=MAX(yfin0,XX(2))
    END DO
    IF (xfin0.LE.xdeb0 .OR. yfin0.LE.ydeb0) THEN
       WRITE(*,*) "DEGENERATE TWO-DIMENSIONAL PLT DOMAIN"
       STOP 1
    END IF

    CALL read_plt_record(11,line,1)
    nb_cell=0
    READ(line,*,IOSTAT=ios) nb_cell
    IF (ios/=0 .OR. nb_cell.LT.1) THEN
       WRITE(*,*) "INVALID PLT CELL COUNT"
       STOP 1
    END IF
    ALLOCATE(mesh(nb_cell))

    DO ncell=1,nb_cell
       CALL read_plt_record(11,line)
       nv=0
       READ(line,*,IOSTAT=ios) nv
       IF (ios/=0 .OR. nv.LT.3 .OR. nv.GT.nb_pts) THEN
          WRITE(*,*) "INVALID PLT POLYGON SIZE",ncell,nv
          STOP 1
       END IF
       IF (plt_token_count(line)/=nv+1) THEN
          WRITE(*,*) "PLT CELL MUST CONTAIN nv DISTINCT POINT INDICES",ncell,nv
          STOP 1
       END IF
       mesh(ncell)%nb_face=nv
       mesh(ncell)%ref=0
       mesh(ncell)%nb_neigh_ext=0
       NULLIFY(mesh(ncell)%neigh_ext)
       ALLOCATE(mesh(ncell)%node(nv),mesh(ncell)%neigh(nv))
       ALLOCATE(mesh(ncell)%length(nv),mesh(ncell)%Normal(nv,2))
       mesh(ncell)%node=0
       READ(line,*,IOSTAT=ios) nv,mesh(ncell)%node
       IF (ios/=0 .OR. ANY(mesh(ncell)%node.LT.1) .OR. ANY(mesh(ncell)%node.GT.nb_pts)) THEN
          WRITE(*,*) "INVALID PLT CELL CONNECTIVITY",ncell
          STOP 1
       END IF
       DO npt=1,nv
          pt=mesh(ncell)%node(npt)
          IF (COUNT(mesh(ncell)%node==pt)/=1) THEN
             WRITE(*,*) "REPEATED VERTEX IN PLT CELL",ncell,pt
             STOP 1
          END IF
          point(pt)%nb_neigh=point(pt)%nb_neigh+1
       END DO
    END DO
    DO pt=1,nb_pts
       IF (point(pt)%nb_neigh==0) THEN
          WRITE(*,*) "UNUSED PLT POINT",pt
          STOP 1
       END IF
    END DO

    CALL read_plt_record(11,line,1)
    IF (TRIM(line)/="BOUNDARY_EDGES") THEN
       WRITE(*,*) "PLT MESH REQUIRES AN EXPLICIT BOUNDARY_EDGES SECTION"
       STOP 1
    END IF
    CALL read_plt_record(11,line,1)
    nb_bdr_edges=0
    READ(line,*,IOSTAT=ios) nb_bdr_edges
    IF (ios/=0 .OR. nb_bdr_edges.LT.1) THEN
       WRITE(*,*) "INVALID PLT BOUNDARY EDGE COUNT"
       STOP 1
    END IF
    ALLOCATE(bdr_edge_nodes(nb_bdr_edges,2),bdr_edge_ref(nb_bdr_edges),bdr_edge_used(nb_bdr_edges))
    bdr_edge_used=0
    ! Keep the Gmsh-specific output counter separate from the PLT table.
    nb_face_bdr_gmsh=0

    DO ib=1,nb_bdr_edges
       CALL read_plt_record(11,line,3)
       a=0; b=0; ref=0
       READ(line,*,IOSTAT=ios) a,b,ref
       IF (ios/=0 .OR. MIN(a,b).LT.1 .OR. MAX(a,b).GT.nb_pts .OR. a==b) THEN
          WRITE(*,*) "INVALID PLT BOUNDARY EDGE",ib,a,b
          STOP 1
       END IF
       IF (ref.LT.1 .OR. ref.GT.nbr_bdr) THEN
          WRITE(*,*) "PLT REFERENCE OUTSIDE 1:nbr_bdr",ib,ref,nbr_bdr
          STOP 1
       END IF
       IF (TRIM(bdr_cond(ref))=="peri" .OR. TRIM(bdr_cond(ref))=="per") THEN
          WRITE(*,*) "PERIODIC PAIRING IS NOT DEFINED IN THIS PLT FORMAT",ib,ref
          STOP 1
       END IF
       bdr_edge_nodes(ib,:)=(/MIN(a,b),MAX(a,b)/)
       bdr_edge_ref(ib)=ref
       DO jb=1,ib-1
          IF (ALL(bdr_edge_nodes(ib,:)==bdr_edge_nodes(jb,:))) THEN
             WRITE(*,*) "DUPLICATE PLT BOUNDARY EDGE",a,b
             STOP 1
          END IF
       END DO
       ! A point flag is only a boundary marker; the face owns its CDL.
       IF (point(a)%ref==0) point(a)%ref=ref
       IF (point(b)%ref==0) point(b)%ref=ref
    END DO
    CALL read_plt_record(11,line,1)
    IF (TRIM(line)/="END_BOUNDARY_EDGES") THEN
       WRITE(*,*) "MISSING END_BOUNDARY_EDGES IN PLT MESH"
       STOP 1
    END IF
    CLOSE(11)

    WRITE(*,*) "POLYGONAL PLT: POINTS, CELLS, BOUNDARY EDGES",nb_pts,nb_cell,nb_bdr_edges
  END SUBROUTINE read_poly_mesh


  ! Read one nonempty record, stripping comments. One record per line.
  SUBROUTINE read_plt_record(unit,line,ntokens)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: unit
    CHARACTER(LEN=*),INTENT(OUT) :: line
    INTEGER,OPTIONAL,INTENT(IN) :: ntokens
    INTEGER :: ios,i,n

    DO
       READ(unit,'(A)',IOSTAT=ios) line
       IF (ios/=0) THEN
          WRITE(*,*) "UNEXPECTED END OR READ ERROR IN PLT MESH"
          WRITE(*,*) "EXPECTED POINTS, CELLS AND A COMPLETE BOUNDARY_EDGES SECTION"
          STOP 1
       END IF
       IF (LEN_TRIM(line)==LEN(line)) THEN
          WRITE(*,*) "PLT RECORD TOO LONG; INCREASE THE LINE BUFFER"
          STOP 1
       END IF
       i=INDEX(line,"#")
       IF (i.GT.0) line(i:)=" "
       DO i=1,LEN_TRIM(line)
          IF (line(i:i)==ACHAR(9)) line(i:i)=" "
       END DO
       line=ADJUSTL(line)
       IF (LEN_TRIM(line)==0) CYCLE
       ! Reject list-directed null/repetition/termination syntax.
       IF (SCAN(line,",/*")/=0) THEN
          WRITE(*,*) "PLT DATA MUST USE EXPLICIT SPACE-SEPARATED VALUES"
          STOP 1
       END IF
       EXIT
    END DO
    IF (PRESENT(ntokens)) THEN
       n=plt_token_count(line)
       IF (n/=ntokens) THEN
          WRITE(*,*) "WRONG NUMBER OF VALUES IN PLT RECORD",n,ntokens
          WRITE(*,*) TRIM(line)
          STOP 1
       END IF
    END IF
  END SUBROUTINE read_plt_record


  INTEGER FUNCTION plt_token_count(line)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: line
    INTEGER :: i
    LOGICAL :: in_token
    plt_token_count=0
    in_token=.FALSE.
    DO i=1,LEN_TRIM(line)
       IF (line(i:i)==" ") THEN
          in_token=.FALSE.
       ELSE IF (.NOT.in_token) THEN
          plt_token_count=plt_token_count+1
          in_token=.TRUE.
       END IF
    END DO
  END FUNCTION plt_token_count


  ! Only consecutive polygon vertices define a face. Reject non-manifold
  ! edges and equal orientations of the shared edge after CCW ordering.
  INTEGER FUNCTION find_polygon_neighbor(ncell,pt1,pt2)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ncell,pt1,pt2
    INTEGER :: iv,vois,k,a,b
    find_polygon_neighbor=-1
    DO iv=1,point(pt1)%nb_neigh
       vois=point(pt1)%neigh(iv)
       IF (vois==ncell) CYCLE
       DO k=1,mesh(vois)%nb_face
          a=mesh(vois)%node(k)
          b=mesh(vois)%node(next_npt(vois,k))
          IF (.NOT.((a==pt1 .AND. b==pt2) .OR. (a==pt2 .AND. b==pt1))) CYCLE
          IF (find_polygon_neighbor.GT.0) THEN
             WRITE(*,*) "NON-MANIFOLD PLT EDGE",pt1,pt2
             STOP 1
          END IF
          IF (a==pt1) THEN
             WRITE(*,*) "INCONSISTENT ORIENTATION OF SHARED PLT EDGE",pt1,pt2
             STOP 1
          END IF
          find_polygon_neighbor=vois
       END DO
    END DO
  END FUNCTION find_polygon_neighbor


  INTEGER FUNCTION find_boundary_edge_ref(pt1,pt2)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: pt1,pt2
    INTEGER :: ib
    find_boundary_edge_ref=0
    DO ib=1,nb_bdr_edges
       IF ((bdr_edge_nodes(ib,1)==pt1 .AND. bdr_edge_nodes(ib,2)==pt2) .OR. &
            & (bdr_edge_nodes(ib,1)==pt2 .AND. bdr_edge_nodes(ib,2)==pt1)) THEN
          find_boundary_edge_ref=bdr_edge_ref(ib)
          bdr_edge_used(ib)=bdr_edge_used(ib)+1
          RETURN
       END IF
    END DO
    WRITE(*,*) "PHYSICAL BOUNDARY EDGE MISSING FROM MESH TABLE",pt1,pt2
    STOP 1
  END FUNCTION find_boundary_edge_ref


  !--------------------------------------------------!
  ! DECAL MESH TO BE BETWEEN xdeb,xfin,ydeb and yfin !
  !--------------------------------------------------!
  SUBROUTINE decal_mesh(pt)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: pt
    REAL(prec),DIMENSION(2) :: XXp,XX0,XX1
    REAL(prec),DIMENSION(2,2) :: Mat_modif

    XXp=point(pt)%Coord
    XX0=(/xdeb0,ydeb0/)
    XX1=(/xdeb,ydeb/)

    Mat_modif=0._prec
    Mat_modif(1,1)=(xfin-xdeb)/(xfin0-xdeb0)
    Mat_modif(2,2)=(yfin-ydeb)/(yfin0-ydeb0)

    point(pt)%Coord=XX1+MATMUL(Mat_modif,XXp-XX0)

  END SUBROUTINE decal_mesh



  !--------------------------------------------!
  ! INIT THE BIGGER CELL NEIGHBORING CELLS SET !
  !--------------------------------------------!
  SUBROUTINE init_neighborhood_ext(ncell)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ncell
    INTEGER :: nneigh,npt,pt,vpt,vois,nmax
    INTEGER,DIMENSION(:),ALLOCATABLE :: list_neigh

    nmax=0
    DO npt=1,mesh(ncell)%nb_face
       pt=mesh(ncell)%node(npt)
       nmax=nmax+point(pt)%nb_neigh
    END DO
    ALLOCATE(list_neigh(nmax))
    nneigh=0
    list_neigh=0
    DO npt=1,mesh(ncell)%nb_face
       pt=mesh(ncell)%node(npt)
       DO vpt=1,point(pt)%nb_neigh
          vois=point(pt)%neigh(vpt)
          IF (vois/=ncell) THEN
             IF (vois.GT.0 .AND. .NOT.ANY(list_neigh==vois)) THEN
                nneigh=nneigh+1
                list_neigh(nneigh)=vois
             END IF
          END IF
       END DO
    END DO
    mesh(ncell)%nb_neigh_ext=nneigh
    ALLOCATE(mesh(ncell)%neigh_ext(nneigh))
    mesh(ncell)%neigh_ext=list_neigh(1:nneigh)
    DEALLOCATE(list_neigh)

  END SUBROUTINE init_neighborhood_ext

  

  !-----------------------------------------------------!
  ! REORGANIZE CELL NODES STARTING FROM THE WIDER ANGLE !
  !-----------------------------------------------------!
  SUBROUTINE reorg_mesh(ncell)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ncell
    INTEGER :: npt,nptp,nptm,nptmax,nn_face
    INTEGER,DIMENSION(mesh(ncell)%nb_face) :: cell_neigh,cell_node
    REAL(prec) :: angle_max
    REAL(prec),DIMENSION(2) :: Np,Npm
    REAL(prec),DIMENSION(mesh(ncell)%nb_face) :: cell_length,cell_angle
    REAL(prec),DIMENSION(mesh(ncell)%nb_face,2) :: cell_Normal

    nn_face=mesh(ncell)%nb_face
    
    cell_neigh=mesh(ncell)%neigh
    cell_node=mesh(ncell)%node
    cell_length=mesh(ncell)%length
    cell_normal=mesh(ncell)%Normal

    DO npt=1,nn_face
       nptm=npt-1
       IF (nptm==0) nptm=nn_face

       Np=cell_normal(npt,:)
       Npm=cell_normal(nptm,:)

       cell_angle(npt)=ACOS(-DOT_PRODUCT(Np,Npm))
    END DO
    angle_max=MAXVAL(cell_angle)

    bcl_angle:DO nptmax=1,nn_face
       IF (cell_angle(nptmax).GE.angle_max-eps0) EXIT bcl_angle
    END DO bcl_angle

    IF (nptmax.GT.1) THEN
       DO npt=1,nn_face-nptmax+1
          mesh(ncell)%neigh(npt)=cell_neigh(nptmax+npt-1)
          mesh(ncell)%node(npt)=cell_node(nptmax+npt-1)
          mesh(ncell)%length(npt)=cell_length(nptmax+npt-1)
          mesh(ncell)%Normal(npt,:)=cell_Normal(nptmax+npt-1,:)
       END DO
       DO npt=nn_face-nptmax+2,nn_face
          mesh(ncell)%neigh(npt)=cell_neigh(npt-nn_face+nptmax-1)
          mesh(ncell)%node(npt)=cell_node(npt-nn_face+nptmax-1)
          mesh(ncell)%length(npt)=cell_length(npt-nn_face+nptmax-1)
          mesh(ncell)%Normal(npt,:)=cell_Normal(npt-nn_face+nptmax-1,:)
       END DO
    END IF

  END SUBROUTINE reorg_mesh


  !----------------------------------------------------------!
  ! FIND THE NODES NEIGHBORING CELLS FOR PERIODIC BOUNDARIES !
  !----------------------------------------------------------!
  SUBROUTINE neigh_periodic
    IMPLICIT NONE

    INTEGER :: pt,pt2,pt3,pref,pref2
    INTEGER :: nb_neigh1,nb_neigh2,nb_neigh3,is_corner,indv
    INTEGER :: nb_corner1,nb_corner2,nb_corner3,nb_corner4
    INTEGER,DIMENSION(5) :: nneigh_c1,nneigh_c2,nneigh_c3,nneigh_c4
    INTEGER,DIMENSION(20) :: nneigh1,nneigh2
    INTEGER,DIMENSION(nb_pts) :: p_ind
    REAL(prec) :: eps
    REAL(prec),DIMENSION(2) :: XXp,XXp2,XXp3,mdf

    eps=min_size/100._prec

    p_ind=0
    bcl01:DO pt=1,nb_pts
       
       pref=point(pt)%ref
       IF (pref==0) CYCLE bcl01
       IF (p_ind(pt)==1) CYCLE bcl01

       is_corner=0
       XXp=point(pt)%Coord
       mdf=0._prec
       IF (pref==1 .AND. ABS(XXp(2)-ydeb).LE.eps) THEN
          mdf(2)=height_mesh
          IF (ABS(XXp(1)-xdeb).LE.eps) is_corner=1
       ELSE IF (pref==2 .AND. ABS(XXp(1)-xfin).LE.eps) THEN
          mdf(1)=-length_mesh
       ELSE IF (pref==3 .AND. ABS(XXp(2)-yfin).LE.eps) THEN
          mdf(2)=-height_mesh
       ELSE IF (pref==4 .AND. ABS(XXp(1)-xdeb).LE.eps) THEN
          mdf(1)=length_mesh
       ELSE
          WRITE(*,*) "PROBLEM 1 IN SUBROUTINE neigh_periodic"
          STOP
       END IF
          
       IF (is_corner==0) THEN
          IF (TRIM(bdr_cond(pref))/="peri") CYCLE bcl01
          XXp2=XXp+mdf
          bcl02:DO pt2=1,nb_pts
             pref2=point(pt2)%ref
             IF (pref2==0) CYCLE bcl02
             IF (TRIM(bdr_cond(pref2))/="peri") CYCLE bcl02
             IF (p_ind(pt2)==1) CYCLE bcl02
             
             XXp3=point(pt2)%Coord
             IF (SQRT(DOT_PRODUCT(XXp3-XXp2,XXp3-XXp2)).LE.eps) THEN
                nb_neigh1=point(pt)%nb_neigh
                nb_neigh2=point(pt2)%nb_neigh
                IF (nb_neigh1.GT.20 .OR. nb_neigh2.GT.20) THEN
                   WRITE(*,*) "PROBLEM 2 IN SUBROUTINE neigh_periodic"
                   STOP
                END IF
                nneigh1(1:nb_neigh1)=point(pt)%neigh
                nneigh2(1:nb_neigh2)=point(pt2)%neigh

                point(pt)%nb_per_neigh=1
                point(pt2)%nb_per_neigh=1
                point(pt)%per_neigh(1)=pt2
                point(pt2)%per_neigh(1)=pt
                
                DEALLOCATE(point(pt)%neigh,point(pt2)%neigh)
                nb_neigh3=nb_neigh1+nb_neigh2
                point(pt)%nb_neigh=nb_neigh3
                point(pt2)%nb_neigh=nb_neigh3
                ALLOCATE(point(pt)%neigh(nb_neigh3),point(pt2)%neigh(nb_neigh3))
                point(pt)%neigh(1:nb_neigh1)=nneigh1(1:nb_neigh1)
                point(pt)%neigh(nb_neigh1+1:nb_neigh3)=nneigh2(1:nb_neigh2)
                point(pt2)%neigh(1:nb_neigh2)=nneigh2(1:nb_neigh2)
                point(pt2)%neigh(nb_neigh2+1:nb_neigh3)=nneigh1(1:nb_neigh1)
                
                p_ind(pt)=1
                p_ind(pt2)=1
             END IF
          END DO bcl02

       ELSE

          nb_corner1=pt
          nb_corner2=0
          nb_corner3=0
          nb_corner4=0

          pt2=nb_corner1
          bcl03:DO pt3=1,nb_pts
             pref2=point(pt3)%ref
             IF (pref2/=4) CYCLE bcl03
             IF (p_ind(pt3)==1) CYCLE bcl03
                
             XXp2=point(pt3)%Coord
             IF (ABS(XXp2(2)-yfin).LE.eps) THEN
                nb_corner4=pt3
                IF (bdr_cond(1)=="peri") THEN
                   nb_neigh1=point(pt2)%nb_neigh
                   nb_neigh2=point(pt3)%nb_neigh
                   IF (nb_neigh1.GT.20 .OR. nb_neigh2.GT.20) THEN
                      WRITE(*,*) "PROBLEM 3 IN SUBROUTINE neigh_periodic"
                      STOP
                   END IF
                   nneigh1(1:nb_neigh1)=point(pt2)%neigh
                   nneigh2(1:nb_neigh2)=point(pt3)%neigh
                   
                   DEALLOCATE(point(pt2)%neigh,point(pt3)%neigh)
                   nb_neigh3=nb_neigh1+nb_neigh2
                   point(pt2)%nb_neigh=nb_neigh3
                   point(pt3)%nb_neigh=nb_neigh3
                   ALLOCATE(point(pt2)%neigh(nb_neigh3),point(pt3)%neigh(nb_neigh3))
!!$                   point(pt2)%neigh(1:nb_neigh1)=nneigh1(1:nb_neigh1)
!!$                   point(pt2)%neigh(nb_neigh1+1:nb_neigh3)=nneigh2(1:nb_neigh2)
                   point(pt3)%neigh(1:nb_neigh2)=nneigh2(1:nb_neigh2)
                   point(pt3)%neigh(nb_neigh2+1:nb_neigh3)=nneigh1(1:nb_neigh1)
                   point(pt2)%neigh(1:nb_neigh3)=point(pt3)%neigh(1:nb_neigh3)
                END IF
                EXIT bcl03
             END IF
          END DO bcl03

          bcl04:DO pt2=1,nb_pts
             pref2=point(pt2)%ref
             IF (pref2/=2) CYCLE bcl04
             
             XXp2=point(pt2)%Coord
             IF (ABS(XXp2(2)-ydeb).LE.eps) EXIT bcl04
          END DO bcl04
          nb_corner2=pt2
          
          bcl05:DO pt3=1,nb_pts
             pref2=point(pt3)%ref
             IF (pref2/=3) CYCLE bcl05
             IF (p_ind(pt3)==1) CYCLE bcl05
             
             XXp2=point(pt3)%Coord
             IF (ABS(XXp2(1)-xfin).LE.eps) THEN
                nb_corner3=pt3
                IF (bdr_cond(1)=="peri") THEN
                   nb_neigh1=point(pt2)%nb_neigh
                   nb_neigh2=point(pt3)%nb_neigh
                   IF (nb_neigh1.GT.20 .OR. nb_neigh2.GT.20) THEN
                      WRITE(*,*) "PROBLEM 4 IN SUBROUTINE neigh_periodic"
                      STOP
                   END IF
                   nneigh1(1:nb_neigh1)=point(pt2)%neigh
                   nneigh2(1:nb_neigh2)=point(pt3)%neigh
                   
                   DEALLOCATE(point(pt2)%neigh,point(pt3)%neigh)
                   nb_neigh3=nb_neigh1+nb_neigh2
                   point(pt2)%nb_neigh=nb_neigh3
                   point(pt3)%nb_neigh=nb_neigh3
                   ALLOCATE(point(pt2)%neigh(nb_neigh3),point(pt3)%neigh(nb_neigh3))
                   point(pt2)%neigh(1:nb_neigh1)=nneigh1(1:nb_neigh1)
                   point(pt2)%neigh(nb_neigh1+1:nb_neigh3)=nneigh2(1:nb_neigh2)
!!$                   point(pt3)%neigh(1:nb_neigh2)=nneigh2(1:nb_neigh2)
!!$                   point(pt3)%neigh(nb_neigh2+1:nb_neigh3)=nneigh1(1:nb_neigh1)
                   point(pt3)%neigh(1:nb_neigh3)=point(pt2)%neigh(1:nb_neigh3)
                END IF
                EXIT bcl05
             END IF
          END DO bcl05

          IF (bdr_cond(2)=="peri") THEN
             pt2=nb_corner1
             pt3=nb_corner2
             nb_neigh1=point(pt2)%nb_neigh
             nb_neigh2=point(pt3)%nb_neigh
             IF (nb_neigh1.GT.20 .OR. nb_neigh2.GT.20) THEN
                WRITE(*,*) "PROBLEM 5 IN SUBROUTINE neigh_periodic"
                STOP
             END IF
             nneigh1(1:nb_neigh1)=point(pt2)%neigh
             nneigh2(1:nb_neigh2)=point(pt3)%neigh
                   
             DEALLOCATE(point(pt2)%neigh,point(pt3)%neigh)
             nb_neigh3=nb_neigh1+nb_neigh2
             point(pt2)%nb_neigh=nb_neigh3
             point(pt3)%nb_neigh=nb_neigh3
             ALLOCATE(point(pt2)%neigh(nb_neigh3),point(pt3)%neigh(nb_neigh3))
             point(pt2)%neigh(1:nb_neigh1)=nneigh1(1:nb_neigh1)
             point(pt2)%neigh(nb_neigh1+1:nb_neigh3)=nneigh2(1:nb_neigh2)
!!$             point(pt3)%neigh(1:nb_neigh2)=nneigh2(1:nb_neigh2)
!!$             point(pt3)%neigh(nb_neigh2+1:nb_neigh3)=nneigh1(1:nb_neigh1)
             point(pt3)%neigh(1:nb_neigh3)=point(pt2)%neigh(1:nb_neigh3)

             pt2=nb_corner3
             pt3=nb_corner4
             nb_neigh1=point(pt2)%nb_neigh
             nb_neigh2=point(pt3)%nb_neigh
             IF (nb_neigh1.GT.20 .OR. nb_neigh2.GT.20) THEN
                WRITE(*,*) "PROBLEM 6 IN SUBROUTINE neigh_periodic"
                STOP
             END IF
             nneigh1(1:nb_neigh1)=point(pt2)%neigh
             nneigh2(1:nb_neigh2)=point(pt3)%neigh
             
             DEALLOCATE(point(pt2)%neigh,point(pt3)%neigh)
             nb_neigh3=nb_neigh1+nb_neigh2
             point(pt2)%nb_neigh=nb_neigh3
             point(pt3)%nb_neigh=nb_neigh3
             ALLOCATE(point(pt2)%neigh(nb_neigh3),point(pt3)%neigh(nb_neigh3))
             point(pt2)%neigh(1:nb_neigh1)=nneigh1(1:nb_neigh1)
             point(pt2)%neigh(nb_neigh1+1:nb_neigh3)=nneigh2(1:nb_neigh2)
!!$             point(pt3)%neigh(1:nb_neigh2)=nneigh2(1:nb_neigh2)
!!$             point(pt3)%neigh(nb_neigh2+1:nb_neigh3)=nneigh1(1:nb_neigh1)
             point(pt3)%neigh(1:nb_neigh3)=point(pt2)%neigh(1:nb_neigh3)
          END IF

          p_ind(nb_corner1)=1
          p_ind(nb_corner2)=1
          p_ind(nb_corner3)=1
          p_ind(nb_corner4)=1

          DEALLOCATE(point(nb_corner1)%per_neigh,point(nb_corner2)%per_neigh)
          DEALLOCATE(point(nb_corner3)%per_neigh,point(nb_corner4)%per_neigh)
          ALLOCATE(point(nb_corner1)%per_neigh(3),point(nb_corner2)%per_neigh(3))
          ALLOCATE(point(nb_corner3)%per_neigh(3),point(nb_corner4)%per_neigh(3))
          point(nb_corner1)%nb_per_neigh=3
          point(nb_corner2)%nb_per_neigh=3
          point(nb_corner3)%nb_per_neigh=3
          point(nb_corner4)%nb_per_neigh=3
          point(nb_corner1)%per_neigh=(/nb_corner2,nb_corner3,nb_corner4/)
          point(nb_corner2)%per_neigh=(/nb_corner3,nb_corner4,nb_corner1/)
          point(nb_corner3)%per_neigh=(/nb_corner4,nb_corner1,nb_corner2/)
          point(nb_corner4)%per_neigh=(/nb_corner1,nb_corner2,nb_corner3/)

       END IF

    END DO bcl01

  END SUBROUTINE neigh_periodic

  

  !---------------------------------------------------------!
  ! REORGANIZE NODES NEIGHBORINGS CELL SET IN DIRECT MANNER !
  !---------------------------------------------------------!
  SUBROUTINE reorg_node_neigh
    IMPLICIT NONE

    INTEGER :: pt,npt,nptm,dim,ncell,ncell2,ii
    LOGICAL :: has_physical_boundary

    !$OMP PARALLEL DO PRIVATE(dim,ii,ncell,npt,nptm,ncell2,has_physical_boundary) &
    !$OMP & SCHEDULE(DYNAMIC) NUM_THREADS(num_thrd)
    DO pt=1,nb_pts

       dim=point(pt)%nb_neigh
       
       ! An open fan starts at its outgoing physical boundary edge.
       ! The two boundary edges may carry different references.
       has_physical_boundary=.FALSE.
       BCL_neigh1:DO ii=1,dim
          ncell=point(pt)%neigh(ii)
          npt=find_npt(ncell,pt)
          IF (mesh(ncell)%neigh(npt).LT.0) THEN
             has_physical_boundary=.TRUE.
             EXIT BCL_neigh1
          END IF
       END DO BCL_neigh1
       IF (has_physical_boundary) CALL reorg_array(point(pt)%neigh,ncell)

       DO ii=1,dim-2
          
          ncell=point(pt)%neigh(ii)
          npt=find_npt(ncell,pt)
          nptm=prev_npt(ncell,npt)
          ncell2=mesh(ncell)%neigh(nptm)
          CALL reorg_array(point(pt)%neigh(ii+1:dim),ncell2)
       END DO
       
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE reorg_node_neigh

END MODULE read_mesh
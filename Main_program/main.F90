PROGRAM main
    use mod_mesh
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_dual

    Msh%DIM = 2
    Msh%max_neigh_Cell=3

    print *, "Reading mesh file..."
    CALL ReadMesh("Meshes/carre_base.mesh",Msh)
    print *, "Mesh file read successfully."
    
    CALL IsDelauney(Msh)

    ! CALL OUTMESH(Msh)
    ! CALL Correct_Mesh(Msh)
    ! CALL OUTMESH(Msh)

    CALL Calc_Surface(Msh)
    CALL Calc_Centroid(Msh)

    print *,"Dual mesh writing"
    CALL Voronoi(Msh,Msh_dual)
    print *,"wrote dual mesh"

    print *, "Writing mesh file..."
    CALL Write_vtk("carre_out.vtk",Msh)
    print *, "Mesh file written successfully."

    print *, "Writing mesh dual file..."
    CALL Write_vtk("carre_dual_out.vtk",Msh_dual)
    print *, "Mesh file written successfully."

END PROGRAM main
PROGRAM main
    use mod_mesh
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_dual

    Msh%DIM = 2
    Msh%maxneigh=3

    print *, "Reading mesh file..."
    CALL ReadMesh("carre_test.mesh",Msh)
    print *, "Mesh file read successfully."

    ! CALL OUTMESH(Msh)
    ! CALL Correct_Mesh(Msh)

    print *,"Dual mesh writing"
    CALL Voronoi(Msh,Msh_dual)
    print *,"wrote dual mesh"

    CALL OUTMESH(Msh_dual)
    ! CALL Correct_Mesh(Msh_Dual)
    ! CALL OUTMESH(Msh_dual)


    print *, "Writing mesh file..."
    CALL Write_vtk("carre_out.vtk",Msh)
    print *, "Mesh file written successfully."


    print *, "Writing mesh dual file..."
    CALL Write_vtk("carre_dual_out.vtk",Msh_dual)
    print *, "Mesh file written successfully."

END PROGRAM main
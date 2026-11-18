Program main
    use mod_mesh
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_del

    Msh%DIM = 2
    Msh%max_neigh_Cell=3

    print *, "Reading mesh file..."
    CALL ReadMesh("Meshes/carre_base.mesh",Msh)
    print *, "Mesh file read successfully."
 


    ! CALL REALLOC_Mesh(Msh, 16, 13, 29)

    CALL Mesh_Creation(Msh,Msh_del)

    print *, "Writing mesh file..."
    CALL Write_vtk("carre_out.vtk",Msh_del)
    print *, "Mesh file written successfully."
END program
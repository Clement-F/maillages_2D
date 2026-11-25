Program main
    use mod_Delaunay
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_del

    Msh%DIM = 2
    Msh%max_neigh_Cell=3

    print *, "Reading mesh file..."
    CALL ReadMesh("Meshes/ls89.mesh",Msh)
    print *, "Mesh file read successfully."
 
    Call WriteMesh("carre_comp.mesh", Msh)
    CALL Write_vtk("carre_comp.vtk",Msh)


    ! CALL REALLOC_Mesh(Msh, 16, 13, 29)
    ! CALL add_Point(Msh, (/0.1_prec,0.1_prec/))
    ! CALL add_Point(Msh, (/0.1_prec,0.2_prec/))
    CALL Mesh_Creation(Msh,Msh_del)

    ! CALL Copy_Mesh(Msh_del,Msh)

    print *, "Writing mesh file..."
    CALL WriteMesh("carre_out.mesh",Msh_del)
    CALL Write_vtk("carre_out.vtk",Msh_del)
    print *, "Mesh file written successfully."
END program
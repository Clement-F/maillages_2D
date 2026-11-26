Program main
    use mod_Delaunay
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_Vor, Msh_Ll

    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    Msh%DIM = 2
    Msh%max_neigh_Cell=3

    print *, "Reading mesh file..."
    CALL ReadMesh("Meshes/carre_test.mesh",Msh)
    print *, "Mesh file read successfully."
 

    CALL Voronoi(Msh, Msh_Vor)

    CALL Write_vtk("carre_comp.vtk",Msh_Vor)

    CALL Lloyd(Msh_Vor, Msh_Ll,Msh)

    print *, "Writing mesh file..."
    CALL WriteMesh("carre_out.mesh",Msh_Ll)
    CALL Write_vtk("carre_out.vtk",Msh_Ll)
    print *, "Mesh file written successfully."

END program
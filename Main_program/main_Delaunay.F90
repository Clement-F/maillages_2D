Program Main_Delaunay
    use mod_mesh
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_dual

    Msh%DIM = 2
    Msh%max_neigh_Cell=3

    print *, "Reading mesh file..."
    CALL ReadMesh("Meshes/carre_base.mesh",Msh)
    print *, "Mesh file read successfully."
 


    ! CALL REALLOC_Mesh(Msh, 16, 13, 29)

    CALL add_Point(Msh,(/-0.12746987_Prec,0.1587634987_prec/))
    CALL add_Point(Msh,(/-0.12746987_Prec,0.1887634987_prec/))
    CALL add_Point(Msh,(/-0.12746987_Prec,0.2387634987_prec/))
    CALL add_Point(Msh,(/-0.12746987_Prec,0.2587634987_prec/))

    print *, "Writing mesh file..."
    CALL Write_vtk("carre_out.vtk",Msh)
    print *, "Mesh file written successfully."
END program
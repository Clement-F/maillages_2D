PROGRAM main
    use mod_mesh
    IMPLICIT NONE

    TYPE(MESH) :: Msh,Msh_dual

    Msh%DIM = 2
    Msh%maxneigh=3

    print *, "Reading mesh file..."
    CALL ReadMesh("carre_base.mesh",Msh)
    print *, "Mesh file read successfully."


    ! print *, "Reading mesh file..."
    ! CALL Readvtk("carre_base.vtk",Msh)
    ! print *, "Mesh file read successfully."


    print *, "Creating Dual mesh..."
    CALL Voronoi(Msh_P=Msh,Msh=Msh_dual)
    print *, "Dual mesh created successfully"

    ! print *, "Writing mesh file..."
    ! CALL WriteMesh("carre_dual_out.mesh",nb_cell, cref, Sommets, nb_Ver, Vref, Coord)
    ! print *, "Mesh file written successfully."

    print *, "Creating Paraview file..."
    CALL Write_vtk("carre_dual_out.vtk", Msh_dual)
    print *, "Paraview file created successfully."

    ! print *, "Writing mesh file..."
    ! CALL WriteMesh("carre_2h_out.mesh",nb_cell0, cref0, Sommets0, nb_Ver0, Vref0, Coord0)
    ! print *, "Mesh file written successfully."

    ! print *, "Creating Paraview file..."
    ! CALL Write_vtk("carre_2h_out.vtk", nb_Ver0, Coord0(:,1), Coord0(:,2), nb_cell0, Sommets0, nb_EBdr0, Bdr0)
    ! print *, "Paraview file created successfully."

    ! print *, "Reading mesh file..."
    ! CALL Readvtk("carre_2h_out.vtk")
    ! print *, "Mesh file read successfully."

END PROGRAM main
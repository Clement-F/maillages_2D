PROGRAM main
    use mod_mesh
    IMPLICIT NONE

    DIM = 2
    maxneigh=3
    ! print *, "Reading mesh file..."
    ! CALL ReadMesh("carre_debug.mesh")
    ! print *, "Mesh file read successfully."

    ! print *, "Creating Paraview file..."
    ! CALL Write_vtk("carre_debug_out.vtk", nb_Ver0, Coord0(:,1), Coord0(:,2), nb_cell0, Sommets0, nb_Bdr0, Bdr0)
    ! print *, "Paraview file created successfully."

    print *, "Reading mesh file..."
    CALL ReadMesh("carre_2h.mesh")
    print *, "Mesh file read successfully."

    print *, "Writing mesh file..."
    CALL WriteMesh("carre_2h_out.mesh",nb_cell0, cref0, Sommets0, nb_Ver0, Vref0, Coord0)
    print *, "Mesh file written successfully."

    print *, "Creating Paraview file..."
    CALL Write_vtk("carre_2h_out.vtk", nb_Ver0, Coord0(:,1), Coord0(:,2), nb_cell0, Sommets0, nb_Bdr0, Bdr0)
    print *, "Paraview file created successfully."

    print *, "Reading mesh file..."
    CALL Readvtk("carre_2h_out.vtk")
    print *, "Mesh file read successfully."

END PROGRAM main
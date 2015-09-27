using FemMesh
using FactCheck

verbose = true

facts("\nTesting generate_mesh function on solids") do
    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI3)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 121

    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI6)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 441

    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD4)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 121

    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 341

    bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=HEX8)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 1331

    bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=HEX20)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 4961

    rm("out.vtk")
end

facts("\nTesting generate_mesh function on trusses") do

    coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
    conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]
    bl   = BlockTruss(coord, conn)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 6

    coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]  
    conn  = [ 1 3; 1 2; 2 3] 
    bl   = BlockTruss(coord, conn)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 3

    rm("out.vtk")
end

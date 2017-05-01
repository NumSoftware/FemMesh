using FemMesh
using FactCheck

verbose = true

facts("\nMesh generation on solids") do
    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI3)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 121

    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI6)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 441

    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD4)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 121

    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 341

    bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=HEX8)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 1331

    bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=HEX20)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 4961

    rm("out.vtk")
end

facts("\nMesh generation on trusses") do
    coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
    conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]
    bl   = BlockTruss(coord, conn)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 6

    coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]  
    conn  = [ 1 3; 1 2; 2 3] 
    bl   = BlockTruss(coord, conn)
    mesh = Mesh(bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 3

    rm("out.vtk")
end

facts("\nMesh operators") do
    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI3)
    mesh = Mesh(bl, verbose=verbose)
    move(bl, x=1)
    mesh = Mesh(mesh, bl, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.points) --> 231

    rm("out.vtk")
end

facts("\nMesh smoothing") do
    bl   = Block2D( [0 0; 1 1], nx=2, ny=2, shape=QUAD4)
    mesh = Mesh(bl, verbose=verbose)
    smooth!(mesh, verbose=verbose)
    @fact mesh.quality --> roughly(1.0, atol=1e-2)

    bl   = Block2D( [0 0; 1 1], nx=2, ny=2, shape=TRI3)
    mesh = Mesh(bl, verbose=verbose)
    smooth!(mesh, verbose=verbose)
    @fact mesh.quality --> roughly(0.97, atol=1e-2)

    bl = Block3D( [ 0 0 0; 2 0 0; 2 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1 ], nx=3, ny=3, nz=3, shape=HEX8)
    mesh = Mesh(bl, verbose=verbose)
    smooth!(mesh, verbose=verbose)
    @fact mesh.quality --> roughly(0.98, atol=1e-2)

    bl = Block2D( [ 1.0 0.0; 2.0 0.0; 0.0 2.0; 0.0 1.0; 1.5 0.0; 1.5 1.5; 0.0 1.5; 0.7 0.7 ], nx=3, ny=6, shape=QUAD8)
    mesh = Mesh(bl, verbose=verbose)
    smooth!(mesh, verbose=verbose)
    @fact mesh.quality --> roughly(0.99, atol=1e-2)

    bl = Block2D( [ 1.0 0.0; 2.0 0.0; 0.0 2.0; 0.0 1.0; 1.5 0.0; 1.5 1.5; 0.0 1.5; 0.7 0.7 ], nx=3, ny=6, shape=TRI6)
    mesh = Mesh(bl, verbose=verbose)
    smooth!(mesh, eps=1e-3, epsmin=1e-2, verbose=verbose)
    @fact mesh.quality --> roughly(0.94, atol=1e-2)
end

facts("\nMesh with embedded cells") do
    bl = Block3D( [0 0 0; 1 1 1], nx=8, ny=8, nz=8, shape=HEX8)
    bli = BlockInset( [0 0 0; 1 1 1] )
    mesh = Mesh(bl, bli, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.cells[:lines]) --> 8

    bl = Block3D( [0 0 0; 1 1 1], nx=8, ny=8, nz=8, shape=HEX20)
    bli = BlockInset( [0 0 0; 1 1 1] )
    mesh = Mesh(bl, bli, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.cells[:lines]) --> 8

    rm("out.vtk")
end

facts("\nMesh extrude") do
    bl = Block2D( [0 0; 1 1], nx=3, ny=3, shape=QUAD4)
    mesh = Mesh(bl, verbose=verbose)
    mesh = extrude(mesh, axis=[0,0,1], len=4, n=10)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.cells) --> 90

    bl = Block2D( [0 0; 1 1], nx=3, ny=3, shape=QUAD4)
    ble = extrude(bl, axis=[0,0,1], len=4, n=10)
    mesh = Mesh(ble, verbose=verbose)
    save(mesh, "out.vtk", verbose=verbose)
    @fact length(mesh.cells) --> 90

    rm("out.vtk")
end

facts("\nMesh rotate") do
    bl = Block2D( [0 0; 1 1], nx=4, ny=4, shape=QUAD4)
    rotate(bl, base = [0.5, 0.5], axis=[1,1], angle=45)
    mesh = Mesh(bl)
    save(mesh, "out.vtk")
    @fact length(mesh.cells) --> 16

    bl = Block2D( [0 0; 1 1], nx=4, ny=4, shape=QUAD4)
    bls = polar(bl, base = [0, 0],  n=4)
    mesh = Mesh(bls)
    save(mesh, "out.vtk")
    @fact length(mesh.cells) --> 64

    rm("out.vtk")
end

facts("\nMesh split using joints") do
    verbose = true
    bl  = Block2D( [0 0; 1 1], nx=4, ny=4, shape=TRI3)
    bli = BlockInset( [ 0 0; 1 1] )
    mesh = Mesh(bl, bli, verbose=verbose)
    mesh = split!(mesh)
    #save(mesh, "out1.vtk", verbose=verbose)
    @fact length(mesh.cells) --> 80

    bl  = Block3D( [0 0 0; 1.0 2.0 1.0], nx=2, ny=4, nz=2, shape=TET4)
    bli = BlockInset( [0 0 0; 1.0 2.0 1.0] )
    mesh = Mesh(bl, bli, verbose=verbose)
    mesh = split!(mesh)
    save(mesh, "out.vtk", verbose=verbose)
    mesh = Mesh("out.vtk")
    save(mesh, "out2.vtk")
    @fact length(mesh.cells) --> 264

    rm("out.vtk")
end

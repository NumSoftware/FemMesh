using FemMesh
using FactCheck

verbose = true

facts("\nMesh generation on solids") do
    bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI3)
    mesh = generate_mesh(bl, verbose=verbose)
    save(mesh, "out0.vtk", verbose=verbose)

    mesh = disjoin!(mesh)
    save(mesh, "out1.vtk", verbose=verbose)

    @fact length(mesh.points) --> 121
end


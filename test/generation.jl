using FemMesh
using Test

printstyled("\nMesh generation on solids\n", color=:cyan)
println("\nMesh using TRI3")
bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI3)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 121

println("\nMesh using TRI6")
bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=TRI6)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 441

println("\nMesh using QUAD4")
bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD4)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 121

println("\nMesh using QUAD8")
bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 341

println("\nMesh using QUAD9")
bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD9)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 441

println("\nMesh using HEX8")
bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=HEX8)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 1331

println("\nMesh using HEX20")
bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=HEX20)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 4961

println("\nMesh using TET4")
bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=TET4)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 1331

println("\nMesh using TET10")
bl = Block3D( [0 0 0; 1 1 1], nx=10, ny=10, nz=10, shape=TET10)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 9261

println("\nMesh using HEX8 in BlockCylinder")
bl = BlockCylinder( [0 0 0; 5 5 5], r=2.0, nr=6, n=4, shape=HEX8)
mesh = Mesh(bl)
@test length(mesh.points) == 445

println("\nMesh using HEX20 in BlockCylinder")
bl = BlockCylinder( [0 0 0; 5 5 5], r=2.0, nr=6, n=4, shape=HEX20)
mesh = Mesh(bl)
@test length(mesh.points) == 1641

rm("out.vtk")


printstyled("\nMesh generation on trusses\n", color=:cyan)
coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]
bl   = BlockTruss(coord, conn)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 6

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]  
conn  = [ 1 3; 1 2; 2 3] 
bl   = BlockTruss(coord, conn)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 3

rm("out.vtk")


printstyled("\nMesh with embedded cells\n", color=:cyan)
bl = Block3D( [0 0 0; 1 1 1], nx=8, ny=8, nz=8, shape=HEX8)
bli = BlockInset( [0 0 0; 1 1 1] )
mesh = Mesh(bl, bli)
save(mesh, "out.vtk")
@test length(mesh.cells[:lines]) == 8

bl = Block3D( [0 0 0; 1 1 1], nx=8, ny=8, nz=8, shape=HEX20)
bli = BlockInset( [0 0 0; 1 1 1] )
mesh = Mesh(bl, bli)
save(mesh, "out.vtk")
@test length(mesh.cells[:lines]) == 8

rm("out.vtk")

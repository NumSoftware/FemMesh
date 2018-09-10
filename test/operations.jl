
using FemMesh
using Test

printstyled("\nMesh move\n", color=:cyan)
bl = Block2D( [0 0; 1 1], nx=10, ny=10, cellshape=TRI3)
mesh = Mesh(bl)
move!(bl, dx=1)
mesh = Mesh(mesh, bl)
save(mesh, "out.vtk")
@test length(mesh.points) == 231

printstyled("\nMesh extrude\n", color=:cyan)
bl = Block2D( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
mesh = Mesh(bl)
mesh = extrude(mesh, len=4, n=10)
save(mesh, "out.vtk")
@test length(mesh.cells) == 90

bl = Block2D( [0 0; 1 1], nx=3, ny=3, cellshape=QUAD4)
ble = extrude(bl, len=4, n=10)
mesh = Mesh(ble)
save(mesh, "out.vtk")
@test length(mesh.cells) == 90

printstyled("\nMesh rotate\n", color=:cyan)
bl = Block2D( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
rotate!(bl, base = [0.5, 0.5], axis=[1,1], angle=45)
mesh = Mesh(bl)
save(mesh, "out.vtk")
@test length(mesh.cells) == 16

bl = Block2D( [0 0; 1 1], nx=4, ny=4, cellshape=QUAD4)
bls = polar(bl, base = [0, 0],  n=4)
mesh = Mesh(bls)
save(mesh, "out.vtk")
@test length(mesh.cells) == 64

printstyled("\nMesh split using joints\n", color=:cyan)
bl  = Block2D( [0 0; 1 1], nx=4, ny=4, cellshape=TRI3)
bli = BlockInset( [ 0 0; 1 1] )
mesh = Mesh(bl, bli)
mesh = generate_joints!(mesh)
@test length(mesh.cells) == 80

bl  = Block3D( [0 0 0; 1.0 2.0 1.0], nx=2, ny=4, nz=2, cellshape=TET4)
bli = BlockInset( [0 0 0; 1.0 2.0 1.0] )
mesh = Mesh(bl, bli)
mesh = generate_joints!(mesh)
save(mesh, "out.vtk")
mesh = Mesh("out.vtk")
save(mesh, "out.vtk")
@test length(mesh.cells) == 264

rm("out.vtk")

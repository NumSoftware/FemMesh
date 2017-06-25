using FemMesh
using Base.Test

print_with_color(:cyan, "\nMesh smoothing\n")
bl   = Block2D( [0 0; 1 1], nx=2, ny=2, shape=QUAD4)
mesh = Mesh(bl)
smooth!(mesh)
@test mesh.quality ≈ 1.0 atol=1e-2

bl   = Block2D( [0 0; 1 1], nx=2, ny=2, shape=TRI3)
mesh = Mesh(bl)
smooth!(mesh)
@test mesh.quality ≈ 0.97 atol=1e-2

bl = Block3D( [ 0 0 0; 2 0 0; 2 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1 ], nx=3, ny=3, nz=3, shape=HEX8)
mesh = Mesh(bl)
smooth!(mesh)
@test mesh.quality ≈ 0.98 atol=1e-2
laplacian_smooth!(mesh)
@test mesh.quality ≈ 0.98 atol=1e-2

bl = Block2D( [ 1.0 0.0; 2.0 0.0; 0.0 2.0; 0.0 1.0; 1.5 0.0; 1.5 1.5; 0.0 1.5; 0.7 0.7 ], nx=3, ny=6, shape=QUAD8)
mesh = Mesh(bl)
smooth!(mesh)
@test mesh.quality ≈ 0.99 atol=1e-2

bl = Block2D( [ 1.0 0.0; 2.0 0.0; 0.0 2.0; 0.0 1.0; 1.5 0.0; 1.5 1.5; 0.0 1.5; 0.7 0.7 ], nx=3, ny=6, shape=TRI6)
mesh = Mesh(bl)
smooth!(mesh, eps=1e-3, epsmin=1e-2)
@test mesh.quality ≈ 0.94 atol=1e-2


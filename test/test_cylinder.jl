using FemMesh

verbose = true

bl = BlockCylinder( [0 0 0; 1 1 1], r=2.0, nr=6, n=4, shape=HEX20)

mesh = Mesh(bl)

#smooth!(mesh)
#laplacian_smooth!(mesh)

save(mesh, "out0.vtk")
save(mesh, "out1.vtk")

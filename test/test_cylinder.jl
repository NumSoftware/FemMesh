using FemMesh

verbose = true

bl = BlockCylinder( [0 0 0; 0 0 10], r=2.0, nr=12, n=10, shape=HEX8)

mesh = generate_mesh(bl)

#smooth!(mesh)
laplacian_smooth!(mesh)

save(mesh, "out.vtk")

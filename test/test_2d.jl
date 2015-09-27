using FemMesh
using FactCheck

bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD8)
#bl = Block2D( [0 0; 1 1], nx=10, ny=10, shape=QUAD4)
mesh = generate_mesh(bl, verbose=false)

save(mesh, "out.vtk")


#facts("\nTest 2D:") do
    #@fact dom.nodes[end].dofdict[:uy].U --> roughly(-9.55, atol=6e-3)
#end


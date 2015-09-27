using FemLab

coord = [ 0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0]   # matriz de coordenadas
conn = [ 1 3; 1 2;2 3]  # matriz de conectividades

blt = BlockTruss(coord, conn)
mesh = generate_mesh(blt, verbose=false)

dom = Domain(mesh)

set_mat(dom.elems, Truss(E=2.0e8, A=0.002) )

set_bc( dom.nodes[:(x==0 && y==0 && z==0)] , ux=0)
set_bc( dom.nodes[:(x==0 && y==1 && z==0)] , ux=0, uy=0, uz=0)
set_bc( dom.nodes[:(x==0 && y==1 && z==1)] , ux=0, uy=0)
set_bc( dom.nodes[:(x==0 && y==0)] , fz=-50.)

solve!(dom, verbose=false)

facts("\nTest Truss 3D") do
    @fact 1 --> 1
end

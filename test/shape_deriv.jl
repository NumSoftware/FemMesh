using FemMesh
using Base.Test

print_with_color(:cyan, "\nShape functions\n")

for shape in ALL_SHAPES
    print("shape : ", shape.name)
    n  = shape.npoints
    ndim = shape.ndim

    # Check at nodes
    RR = [ shape.nat_coords[i,:] for i=1:n ]
    NN = shape.func.(RR)

    I = hcat(NN...) # should provida an identity matrix
    @test I ≈ eye(n) atol=1e-10

    # Check at default set of integration points
    Q  = shape.quadrature[0]
    nip, _ = size(Q)
    #@show Q
    RR = [ Q[i,:] for i=1:nip ]
    NN = shape.func.(RR)
    @test sum(sum(NN)) ≈ nip atol=1e-10
    println("  ok")
end

print_with_color(:cyan, "\nShape functions derivatives\n")

for shape in ALL_SHAPES
    print("shape : ", shape.name)
    n    = shape.npoints
    ndim = shape.ndim
    RR = [ shape.nat_coords[i,:] for i=1:n ]
    f  = shape.func

    # numerical derivative
    δ  = 1e-8
    for R in RR
        RI = R .+ eye(ndim)*δ
        fR = f(R)
        D  = zeros(ndim, n)
        for i=1:ndim
            Di     = 1/δ*(f(RI[:,i]) - fR)
            D[i,:] = Di
        end
        @test D ≈ shape.deriv(R) atol=1e-7
    end
    println("  ok")
end

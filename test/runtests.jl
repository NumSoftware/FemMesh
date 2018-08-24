using FemMesh
using Test

FILES = [
    "shape_deriv.jl"
    "generation.jl"
    "operations.jl"
    "extrapolation.jl"
    "smoothing.jl"
]

@testset begin
    for f in FILES
        printstyled( "\nRunning file ", f,"...\n", color=:white)
        include(f)
    end
end





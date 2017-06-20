using FemMesh
using Base.Test

FILES = [
    "shape_deriv.jl"
    "generation.jl"
    "operations.jl"
    "extrapolation.jl"
    "smoothing.jl"
]

@testset begin
    for f in FILES
        print_with_color(:white, "\nRunning file ", f,"...\n")
        include(f)
    end
end





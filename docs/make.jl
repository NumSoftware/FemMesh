using Documenter, FemMesh

#makedocs()
makedocs(
    modules = [FemMesh],
    format = Documenter.Formats.HTML,
    sitename = "FemMesh",
    pages = Any[
        "Home" => "index.md",
        "Index" => "lib/public.md"
    ],
    doctest = false
)

deploydocs(
    repo = "github.com/NumSoftware/FemMesh.git",
    target = "build",
    julia  = "1.0",
    deps = nothing,
    make = nothing,
)

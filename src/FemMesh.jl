# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

__precompile__()

"""
**FemMesh.jl**

FemMesh module implements functions and types related to mesh generation for
finite element analyses. 

**Important data types**

Block2D, Block3D, BlockTruss, BlockInset, Point, Cell, Face, Edge, Mesh.

**Important functions** 

copy, move!, array, rotate!, polar, extrude.

"""
module FemMesh
using Printf, Statistics, LinearAlgebra, SparseArrays
using JSON, DataStructures

# Mesh module
include("tools/linalg.jl")
include("tools/expr.jl")
include("tools/show.jl")
include("tools/iteration.jl")
include("tools/table.jl")
export unfold

# Generic exports
export getindex

include("shape.jl")
export ShapeType, TagType, ALL_SHAPES, ShapeFamily
export get_ip_coords, get_shape_from_vtk
export inverse_map, extrapolator

include("entities.jl")
export Point, Cell, hash, get_x, get_y, get_z
export get_coords, get_point, get_points, get_faces, cell_extent, cell_quality
export tag!, iptag!, update!

include("ugrid.jl")
export UnstructuredGrid, save_vtk, read_ugrid_vtk

include("mesh.jl")
export Mesh, update!, quality!, reorder!, save, get_surface, get_neighbors

include("block.jl")
export Block2D, Block3D, BlockTruss, BlockCoords, BlockCylinder

include("operators.jl")
export move!, array, copy, mirror, rotate!, polar, roll_axes!

include("extrude.jl") 
export extrude

include("smooth.jl") 
include("split.jl") 

include("embedded.jl") 
export generate_embedded_cells!

include("mplot.jl") 
export mplot


# show functions for common structures and arrays
@show_function ShapeType
@show_array_function ShapeType
@show_function Point
@show_array_function Point
@show_function Cell
@show_array_function Cell
@show_function Block
@show_array_function Block

@show_function Mesh
@show_function UnstructuredGrid

# precompile hint
#bl = Block2D( [0 0; 1 1], nx=2, ny=2, shape=TRI3)
#mesh = Mesh(bl, verbose=false)
#bl = Block2D( [0 0; 1 1], nx=2, ny=2, shape=QUAD4)
#mesh = Mesh(bl, verbose=false)
#bl = Block3D( [0 0 0; 1 1 1], nx=2, ny=2, nz=2, shape=HEX8)
#mesh = Mesh(bl, verbose=false)

end#module

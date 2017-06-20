##############################################################################
#    FemMesh - Finite Element Library                                        #
#    Copyright (C) 2014 Raul Durand <raul.durand at gmail.com>               #
#                                                                            #
#    This file is part of FemMesh.                                           #
#                                                                            #
#    FemMesh is free software: you can redistribute it and/or modify         #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    any later version.                                                      #
#                                                                            #
#    FemMesh is distributed in the hope that it will be useful,              #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with FemMesh.  If not, see <http://www.gnu.org/licenses/>.        #
##############################################################################

#VERSION >= v"0.4.0-dev+6521" && __precompile__()
__precompile__()

"""
**FemMesh.jl**

FemMesh module implements functions and types related to mesh generation for
finite element analyses. 

**Important data types**

Block2D, Block3D, BlockTruss, BlockInset, Point, Cell, Face, Edge, Mesh.

**Important functions** 

generate_mesh, copy, move, array, rotate, polar, extrude.

"""
module FemMesh
using JSON

# Mesh module
include("tools/linalg.jl")

# Generic exports
export getindex

include("shape.jl")
export ShapeType, ALL_SHAPES, ShapeClass
export get_ip_coords, get_shape_from_vtk
export inverse_map, extrapolator

include("entities.jl")
export Point, Cell, hash, get_coords, get_point, get_faces, cell_extent, cell_quality
export update!

include("mesh.jl")
export Mesh, update!, quality!, reorder!, save, get_surface, get_neighbors

include("block.jl")
export Block2D, Block3D, BlockTruss, BlockCoords, BlockCylinder

include("operators.jl")
export move, array, copy, rotate, polar

include("extrude.jl") 
export extrude

include("smooth.jl") 
include("split.jl") 

include("embedded.jl") 
export generate_embedded_cells!

include("draw.jl") 
export draw

end#module

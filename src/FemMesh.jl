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

# Constants
if is_linux()
    const RED     = "\x1b[31m"
    const GREEN   = "\x1b[32m"
    const YELLOW  = "\x1b[33m"
    const BLUE    = "\x1b[34m"
    const MAGENTA = "\x1b[35m"
    const CYAN    = "\x1b[36m"
    const WHITE   = "\x1b[37m"
    const BOLD    = "\x1b[1m"
    const DEFAULT = "\x1b[0m"
else
    const RED     = ""
    const GREEN   = "" 
    const BLUE    = ""
    const MAGENTA = ""
    const CYAN    = ""
    const WHITE   = ""
    const BOLD    = "" 
    const DEFAULT = "" 
end

# Types
typealias Vect Array{Float64, 1}
typealias Matx Array{Float64, 2}

# Mesh module
include("tools/linalg.jl")
include("mesh.jl")

end#module

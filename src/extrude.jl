
##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) 2014 Raul Durand <raul.durand at gmail.com>               #
#                                                                            #
#    This file is part of FemLab.                                            #
#                                                                            #
#    FemLab is free software: you can redistribute it and/or modify          #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    any later version.                                                      #
#                                                                            #
#    FemLab is distributed in the hope that it will be useful,               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with FemLab.  If not, see <http://www.gnu.org/licenses/>.         #
##############################################################################

export extrude

# Generates a new mesh obtained by extrusion of a 2D mesh
function extrude{T}(mesh::Mesh, dir::Array{T,1}, len::Number, ndiv::Int=10, verbose::Bool=true)
    verbose && println(BOLD, CYAN, "Mesh extrude:", DEFAULT)

    V = dir/norm(dir)
    δ = len/ndiv
    inipoints = mesh.points
    inicells  = mesh.cells
    newmesh   = Mesh()

    length(inicells)>0 || error("Extrude: Cannot extrude mesh with no cells.")
    
    # check if all cells are the same
    shape = inicells[1].shape
    any(Bool[ c.shape!=shape for c in inicells ]) && error("Extrude: Input mesh shoud have same type of cells.")

    # check if cells are QUAD4 or QUAD8
    (shape == QUAD4 || shape == QUAD8) || error("Error: Only can extrude meshes of QUAD4 and QUAD8 cells.")

    # generate extra nodes and cells
    if shape == QUAD4
        # generate new points
        p_arr = [ Array(Point,ndiv+1) for p in inipoints ]
        for (i,p) in enumerate(inipoints)
            pp = p_arr[i]
            X = [ p.x, p.y, p.z ]
            for j=1:length(pp)
                C = X + V*δ*(j-1)
                newp  = Point(C)
                pp[j] = newp
                push!(newmesh.points, newp)
            end
        end

        # generate new cells
        c_arr = [ Array(Cell, ndiv) for c in inicells ]
        for (i,c) in enumerate(inicells)
            cc = c_arr[i]
            for j=1:ndiv
                cpoints = Array(Point, 8)
                cpoints[1] = p_arr[ c.points[1].id ][j]
                cpoints[2] = p_arr[ c.points[2].id ][j]
                cpoints[3] = p_arr[ c.points[3].id ][j]
                cpoints[4] = p_arr[ c.points[4].id ][j]
                cpoints[5] = p_arr[ c.points[1].id ][j+1]
                cpoints[6] = p_arr[ c.points[2].id ][j+1]
                cpoints[7] = p_arr[ c.points[3].id ][j+1]
                cpoints[8] = p_arr[ c.points[4].id ][j+1]

                cpoints[1] = p_arr[ c.points[1].id ][j+1]
                cpoints[2] = p_arr[ c.points[2].id ][j+1]
                cpoints[3] = p_arr[ c.points[3].id ][j+1]
                cpoints[4] = p_arr[ c.points[4].id ][j+1]
                cpoints[5] = p_arr[ c.points[1].id ][j]
                cpoints[6] = p_arr[ c.points[2].id ][j]
                cpoints[7] = p_arr[ c.points[3].id ][j]
                cpoints[8] = p_arr[ c.points[4].id ][j]

                cell = Cell(HEX8, cpoints, c.tag)
                push!(newmesh.cells, cell)
            end
        end
    end

    # generate extra nodes and cells
    if shape == QUAD8

        # flags for middle points
        middle = trues(length(inipoints))
        for c in inicells
            for i=1:4
                middle[ c.points[i].id ] = false
            end
        end

        # generate new points
        p_arr = [ Array(Point,2*ndiv+1) for p in inipoints ]
        for (i,p) in enumerate(inipoints)
            pp = p_arr[i]
            X = [ p.x, p.y, p.z ]
            for j=1:length(pp)
                (middle[i] && iseven(j)) && continue
                C = X + V*δ*(j-1)
                newp  = Point(C)
                pp[j] = newp
                push!(newmesh.points, newp)
            end
        end

        # generate new cells
        c_arr = [ Array(Cell, ndiv) for c in inicells ]
        for (i,c) in enumerate(inicells)
            cc = c_arr[i]
            for j=1:2:2*ndiv
                cpoints = Array(Point, 20)
                cpoints[1 ] = p_arr[ c.points[1].id ][j]
                cpoints[2 ] = p_arr[ c.points[2].id ][j]
                cpoints[3 ] = p_arr[ c.points[3].id ][j]
                cpoints[4 ] = p_arr[ c.points[4].id ][j]
                cpoints[5 ] = p_arr[ c.points[1].id ][j+2]
                cpoints[6 ] = p_arr[ c.points[2].id ][j+2]
                cpoints[7 ] = p_arr[ c.points[3].id ][j+2]
                cpoints[8 ] = p_arr[ c.points[4].id ][j+2]
                cpoints[9 ] = p_arr[ c.points[5].id ][j]
                cpoints[10] = p_arr[ c.points[6].id ][j]
                cpoints[11] = p_arr[ c.points[7].id ][j]
                cpoints[12] = p_arr[ c.points[8].id ][j]
                cpoints[13] = p_arr[ c.points[5].id ][j+2]
                cpoints[14] = p_arr[ c.points[6].id ][j+2]
                cpoints[15] = p_arr[ c.points[7].id ][j+2]
                cpoints[16] = p_arr[ c.points[8].id ][j+2]
                cpoints[17] = p_arr[ c.points[1].id ][j+1]
                cpoints[18] = p_arr[ c.points[2].id ][j+1]
                cpoints[19] = p_arr[ c.points[3].id ][j+1]
                cpoints[20] = p_arr[ c.points[4].id ][j+1]
                cell = Cell(HEX20, cpoints, c.tag)
                push!(newmesh.cells, cell)
            end
        end
    end

    update!(newmesh)

    if verbose
        @printf "  %5d points obtained\n" length(newmesh.points)
        @printf "  %5d cells obtained\n" length(newmesh.cells)
        #if genfacets
            #@printf "  %5d faces obtained\n" nfaces
        #end
        #if genedges
            #@printf "  %5d surface edges obtained\n" nedges
        #end
        println("  done.")
    end

    return newmesh
end


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

# This file includes the code for adding joints between cells

#include("mesh.jl")
export disjoin!

function joint_shape(shape::ShapeType)
    if shape == LIN2  ; return JLIN2  end
    if shape == LIN3  ; return JLIN3  end
    if shape == LIN4  ; return JLIN4  end
    if shape == TRI3  ; return JTRI3  end
    if shape == TRI6  ; return JTRI6  end
    if shape == QUAD4 ; return JQUAD4 end
    if shape == QUAD8 ; return JQUAD8 end
    error("No joint for shape $shape")
end

# Adds joint cells over all shared faces
function disjoin!(mesh::Mesh)
    cells  = mesh.cells
    solids = filter( c -> is_solid(c.shape), cells)

    newpoints = Point[]

    # Splitting
    for c in solids
        is_solid(c.shape) || continue
        for (i,p) in enumerate(c.points)
            newp = Point(p.x, p.y, p.z)
            p.id = 0 # set id=0 as a flag for points to be removed
            newp.extra = c.id
            push!(newpoints, newp)
            c.points[i] = newp
        end
    end

    # List all repeated faces
    face_pairs = Tuple{Cell, Cell}[]

    # Joints generation
    facedict = Dict{UInt64, Cell}()

    # Get paired faces 
    for cell in solids
        for face in get_faces(cell)
            hs = hash(face)
            f  = get(facedict, hs, nothing)
            if f==nothing
                facedict[hs] = face
            else
                push!(face_pairs, (face, f))
                delete!(facedict, hs)
            end
        end
    end

    # Generate joint elements
    jcells = Cell[]
    for (f1, f2) in face_pairs
        n   = length(f1.points)
        con = Array(Point, 2*n)
        for (i,p1) in enumerate(f1.points)
            for p2 in f2.points
                if hash(p1)==hash(p2)
                    con[i]   = p1
                    con[n+i] = p2
                    break
                end
            end
        end

        jshape = joint_shape(f1.shape)
        cell = Cell(jshape, con, "")
        push!(jcells, cell)
    end

    # Fix 1d joint elements
    cdict = Dict{UInt64, Cell}()
    for c in solids
        if c.crossed
            cdict[ hash(c) ] = c
        end
    end

    for c in cells
        if c.shape in (LINK2, LINK3)
            njpoints = elem.shape==LINK2? 2 : 3
            ncpoints = length(c.points)  # number of points of crossed cell
            hs    = sum([ hash(p) for p in c.points[1:end-njpoints]])
            ccell = cdict[hs]  # crossed cell
            c.points[1:ncpoints] = ccell.points[:]
        end
    end

    # Numbering new joint cells
    #n = length(cells)
    #for (i,c) in enumerate(jcells)
        #c.id = n+i
    #end

    # Updating old points numbers
    #points = filter( p -> p.id!=0, mesh.points)
    #for (i,p) in enumerate(points)
        #p.id = i
    #end

    # Updating new points numbers
    #n = length(points)
    #for (i,p) in enumerate(newpoints)
        #p.id = n+i
    #end

    # Updating mesh
    points      = filter( p -> p.id!=0, mesh.points)
    mesh.points = vcat(points, newpoints)
    mesh.cells  = vcat(cells, jcells)
    update!(mesh)

    return mesh
    
end





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

"""
`copy(block)` 

Creates a copy of a `block` object.
"""
copy(bl::BlockTruss) = BlockTruss(copy(bl.coords), bl.conns, shape=bl.shape, tag=bl.tag)

copy(bl::Block2D) = Block2D(copy(bl.coords), nx=bl.nx, ny=bl.ny, shape=bl.shape, tag=bl.tag)

copy(bl::Block3D) = Block3D(copy(bl.coords), nx=bl.nx, ny=bl.ny, nz=bl.nz, shape=bl.shape, tag=bl.tag)

"""
`move(block, [x=0.0,] [y=0.0,] [z=0.0])` 

Changes the coordinates of a `block`. Also returns a reference.
"""
function move!(bl::Block;x=0.0, y=0.0, z=0.0, dx=0.0, dy=0.0, dz=0.0)
    bl.coords[:, 1] += x
    bl.coords[:, 2] += y
    bl.coords[:, 3] += z
    #return bl
end


"""
`move(blocks, [x=0.0,] [y=0.0,] [z=0.0])` 

Changes the coordinates of an array of blocks. Also returns a reference.
"""
function move!(blocks::Array; x=0.0, y=0.0, z=0.0, dx=0.0, dy=0.0, dz=0.0)
    for bl in blocks
        bl.coords[:, 1] += x
        bl.coords[:, 2] += y
        bl.coords[:, 3] += z
    end
    #return blocks 
end


function scale!(bl::Block; factor=1.0, base=[0.,0,0])
    bl.coords = ( base .+ factor*(bl.coords' .- base) )'
end


function scale!(blocks::Array{Block,1}; factor=1.0, base=[0.,0,0])
    for bl in blocks
        bl.coords = ( base .+ factor*(bl.coords' .- base) )'
    end
end

function mirror(block::Block; face=[0. 0 0; 0 1 0; 0 0 1])
    nr, nc = size(face)
    if nc==2
        face = [ face zeros(nr) ]
    end
    if nr==2
        face = [ face; [face[1,1] face[1,2] 1.0] ]
    end

    p1 = face[1,:]
    p2 = face[2,:]
    p3 = face[3,:]
    normal = cross(p2-p1, p3-p1)
    normal = normal/norm(normal)

    bl = copy(block)

    distances    = (bl.coords .- p1')*normal          # d = n^.(xi - xp)
    bl.coords = bl.coords .- 2*distances.*normal'  # xi = xi - 2*d*n^

    # fix coordinates in bl to keep anti-clockwise numbering
    npts = size(bl.coords)[1]
    ndim = typeof(bl)==Block2D ? 2 : 3

    if npts==8 && ndim==2
        idxs = [ 4, 3, 2, 1, 7, 6, 5, 8 ]
    elseif npts==8 && ndim==3
        idxs = [ 5:8; 1:4 ]
    elseif npts==20 && ndim==3
        idxs = [ 5:8; 1:4; 13:16; 9:12; 17:20 ]
    else
        idxs = [ npts:-1:1; ]  # reverse
    end
    bl.coords = bl.coords[idxs,:]

    return bl
end

function mirror(blocks::Array{Block,1}; face=[0. 0 0; 0 1 0; 0 0 1])
    for bl in newblocks
        mirror(bl, face=face)
    end
end


"""
`array(block, [n=2,] [x=0.0,] [y=0.0,] [z=0.0])` 

Creates `n-1` copies of a `block` separated by distances `x`, `y` and/or `z` along respective axes.
Returns n objects.
"""
function array(bl::Block; n=2, x=0.0, y=0.0, z=0.0)
    blocks = [ bl ]
    for i=1:n-1
        dx = i*x
        dy = i*y
        dz = i*z
        cp = copy(bl)
        move(cp, x=dx, y=dy, z=dz)
        push!(blocks, cp)
    end
    return blocks
end


"""
`rotate!(block, [base=[0,0,0],] [axis=[0,0,1],] [angle=90.0])`

Rotates a `block` according to a `base` point, an `axis` vector and an `angle`.
"""
function rotate!(bl::Block; base=[0.,0,0], axis=[0.,0,1], angle=90.0 )

    length(axis)==2 && ( axis=vcat(axis, 0.0) )
    length(base)==2 && ( base=vcat(base, 0.0) )

    # unit vector
    axis = axis/norm(axis)
    a, b, c = axis
    d = sqrt(b^2+c^2)

    # unit vector for rotation
    l = cos(angle*pi/180)
    m = sin(angle*pi/180)

    # Rotation matrices
    if d != 0.0
        Rx  = [  1.    0.    0.
                 0.   c/d  -b/d 
                 0.   b/d   c/d ]

        Rxi = [  1.    0.    0.
                 0.   c/d   b/d 
                 0.  -b/d   c/d ]
    end

    Ry  = [   d    0.  -a
             0.    1.  0.
              a    0.   d ]
           
    Ryi = [   d    0.   a
             0.    1.  0.
             -a    0.   d ]

    Rz  = [   l   -m   0.
              m    l   0.
             0.   0.   1. ]

    # all rotations matrix
    if d != 0.0
        R = Rxi*Ryi*Rz*Ry*Rx
    else
        R = Ryi*Rz*Ry
    end

    # equation: p2 = base + R*(p-base)
    bl.coords = ( base .+ R*(bl.coords' .- base) )'
    #return bl
end

function rotate!{T<:Block}(blocks::Array{T,1}; base=[0.,0,0], axis=[0.,0,1], angle=90.0 )
    for bl in blocks
        rotate!(bl, base=base, axis=axis, angle=angle)
    end
end


"""
`polar(block, [base=[0,0,0],] [axis=[0,0,1],] [angle=360.0,] [n=2])`

Creates `n-1` copies of a `block` and places them using polar distribution based on 
a `base` point, an `axis` vector, a total `angle`.
"""
function polar{T<:Block}(bl::T; base=[0.,0,0], axis=[0.,0,1], angle=360, n=2 )::Array{T,1}
    blocks::Array{T,1} = [ bl ]
    angle = angle/n
    for i=1:n-1
        bli = copy(bl)
        rotate!(bli, base=base, axis=axis, angle=angle*i)
        push!(blocks, bli)
    end
    return blocks
end

function polar{T<:Block}(blocks::Array{T,1}; base=[0.,0,0], axis=[0.,0,1], angle=360, n=2 )::Array{T,1}
    rblocks::Array{T,1} = []

    for bl in blocks
        bls = polar(bl, base=base, axis=axis, angle=angle, n=n)
        append!(rblocks, bls)
    end

    return rblocks
end


"""
`move(mesh, [dx=0.0,] [dy=0.0,] [dz=0.0])` 

Moves a Mesh object `mesh`. Also returns a reference.
"""
function move!(mesh::Mesh; dx=0.0, dy=0.0, dz=0.0)
    for p in mesh.points
        p.x += dx
        p.y += dy
        p.z += dz
    end
    #return mesh
end



function scale!(msh::Mesh; factor=1.0, base=[0.,0,0])
    for pt in mesh.points
        p.x, p.y, p.z = base + ([p.x, p.y, p.z] - base)*factor
    end
    #return mesh
end


"""
`rotate(mesh; base=[0,0,0], axis=[0,0,1], angle=90.0)`

Rotates a Mesh object `mesh` according to a `base` point, an `axis` vector and an `angle`.
"""
function rotate!(mesh::Mesh; base=[0.,0,0], axis=[0.,0,1], angle=90.0 )

    length(axis)==2 && ( axis=vcat(axis, 0.0) )
    length(base)==2 && ( base=vcat(base, 0.0) )
    norm(axis) < 1e-10 && error("rotate: Invalid axis $axis")

    # unit vector
    axis = axis/norm(axis)
    a, b, c = axis
    d = sqrt(b^2+c^2)
    d==0.0 && ( d=1.0 )

    # unit vector for rotation
    l = cos(angle*pi/180)
    m = sin(angle*pi/180)

    # Rotation matrices
    Rx  = [  1.    0.    0.
             0.   c/d  -b/d 
             0.   b/d   c/d ]

    Rxi = [  1.    0.    0.
             0.   c/d   b/d 
             0.  -b/d   c/d ]

    Ry  = [   d    0.  -a
             0.    1.  0.
              a    0.   d ]
           
    Ryi = [   d    0.   a
             0.    1.  0.
             -a    0.   d ]

    Rz  = [   l   -m   0.
              m    l   0.
             0.   0.   1. ]

    # all rotations matrix
    R = Rxi*Ryi*Rz*Ry*Rx

    for p in mesh.points
        p.x, p.y, p.z = base + R*([p.x, p.y, p.z] - base)
    end

    #return mesh
end


# TESTING
function get_surface_alt(cells::Array{Cell,1})
    # Actually slower....
    # Get all points
    pointsd = Dict{UInt64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    points = values(pointsd)

    # Get incidence matrix (shares) (fast)
    np = length(points)
    N = [ Cell[] for i=1:np]
    for cell in cells
        for pt in cell.points
            push!(N[pt.id], cell)
        end
    end

    @show 1
    # Get matrix of cells faces
    F = [ get_faces(cell) for cell in cells]
    nc = length(cells)
    #CF = Array(Array{Array{Int64,1},1}, nc)
    CF = Array(Array{UInt64,1}, nc)
    for cell in cells # fast
        #CF[cell.id] = [ sort([pt.id for pt in face.points]) for face in F[cell.id]]
        CF[cell.id] = [ hash(face) for face in F[cell.id]]
    end

    @show 2
    # Get cells boundary flag matrix
    CB = [ trues(length(CF[cell.id])) for cell in cells]
    for cell in cells
        for (i,fcon) in enumerate(CF[cell.id])
            for pid in fcon
                for cl in N[pid]
                    if cl.id == cell.id; continue end
                    if fcon in CF[cl.id]
                        CB[cell.id][i] = false
                    end
                end
            end
        end
    end

    @show 3
    # Get list of boundary faces (almost fast)
    facets = Cell[]
    for cell in cells
        for (i,face) in enumerate(F[cell.id])
            if CB[cell.id][i]
                push!(facets, face)
            end
        end
    end
    #return facets
end

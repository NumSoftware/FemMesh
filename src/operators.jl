
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

export move, array, copy, rotate, polar

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
function move(bl::Block;x=0.0, y=0.0, z=0.0, dx=0.0, dy=0.0, dz=0.0)
    bl.coords[:, 1] += x
    bl.coords[:, 2] += y
    bl.coords[:, 3] += z
    return bl
end


"""
`move(blocks, [x=0.0,] [y=0.0,] [z=0.0])` 

Changes the coordinates of an array of blocks. Also returns a reference.
"""
function move(blocks::Array; x=0.0, y=0.0, z=0.0, dx=0.0, dy=0.0, dz=0.0)
    for bl in blocks
        bl.coords[:, 1] += x
        bl.coords[:, 2] += y
        bl.coords[:, 3] += z
    end
    return blocks 
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
`rotate(block, [base=[0,0,0],] [axis=[0,0,1],] [angle=90.0])`

Rotates a `block` according to a `base` point, an `axis` vector and an `angle`.
"""
function rotate(bl::Block; base=[0.,0,0], axis=[0.,0,1], angle=90.0 )

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
    return bl
end

"""
`polar(block, [base=[0,0,0],] [axis=[0,0,1],] [angle=360.0,] [n=2])`

Creates `n-1` copies of a `block` and places them using polar distribution based on 
a `base` point, an `axis` vector, a total `angle`.
"""
function polar(bl::Block; base=[0.,0,0], axis=[0.,0,1], angle=360, n=2 )
    blocks = [ bl ]
    angle = angle/n
    for i=1:n-1
        bli = copy(bl)
        rotate(bli, base=base, axis=axis, angle=angle*i)
        push!(blocks, bli)
    end
    return blocks
end

function polar(blocks::Array; base=[0.,0,0], axis=[0.,0,1], angle=360, n=2 )
    blocks3D = []

    for bl in blocks
        bls = polar(bl, base=base, axis=axis, angle=angle, n=n)
        append!(blocks3D, bls)
    end

    return blocks3D
end



function scale(msh::Mesh, value; base=[0.,0,0])
    return msh
end

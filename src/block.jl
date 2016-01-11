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

export Block2D, Block3D, BlockTruss, BlockCoords
export move, copy, rotate, polar
import Base.copy

### Type Block
abstract Block

### Block types:

include("block_inset.jl")

"""
`BlockTruss(coords, conns, [shape=LIN2,] [tag="",])`

Generates a block object for the mesh generation of trusses:
"""
type BlockTruss <: Block
    coords::Array{Float64,2}
    conns ::Array{Int64,2}
    shape ::ShapeType
    tag::AbstractString
    id::Int64

    function BlockTruss(coords::Array{Float64,2}, conns::Array{Int64,2}; shape=LIN2, tag="", id=-1)
        ncols = size(coords,2)
        if !(ncols in (2,3)); error("Invalid coordinates matrix for BlockTruss") end
        if ncols==2
            C = [coords  zeros(size(coords,1)) ]
        else
            C = coords
        end
        this = new(C, conns, shape, tag, id)
        return this
    end
end

"""
`copy(block)` 

Creates a copy of a `block` object.
"""
copy(bl::BlockTruss) = BlockTruss(copy(bl.coords), bl.conns, shape=bl.shape, tag=bl.tag)


"""
`BlockCoords(coords, conns, [shape=LIN2,] [tag="",])`

Generates a 2D or 3D block object based on a matrix of coordinates and a matrix of connectivities.
"""
type BlockCoords <: Block
    coords::Array{Float64,2}
    conns ::Array{Int64,2}
    tag::AbstractString
    id::Int64

    function BlockCoords(coords::Array{Float64,2}, conns::Array{Int64,2}; tag="", id=-1)
        ncols = size(coords,2)
        if !(ncols in (2,3)); error("Invalid coordinates matrix for BlockCoords") end
        if ncols==2
            C = [coords  zeros(size(coords,1)) ]
        else
            C = coords
        end
        this = new(C, conns, tag, id)
        return this
    end
end

copy(bl::BlockTruss) = BlockTruss(copy(bl.coords), bl.conns, shape=bl.shape, tag=bl.tag)


"""
`Block2D(coords, [nx=1,] [ny=1,] [shape=QUAD4,] [tag=""] )`

Generates a block object for the mesh generation of 2D meshes.
`shape` can be TRI3, TRI6, QUAD4, QUAD8.
"""
type Block2D <: Block
    coords::Array{Float64,2}
    nx::Int64
    ny::Int64
    shape::ShapeType
    tag::AbstractString
    id::Int64
    function Block2D(coords; nx=1, ny=1, shape=QUAD4, tag="", id=-1)
        if size(coords,1)==2
            C = box_coords(vec(coords[1,:]), vec(coords[2,:]))
        elseif size(coords,2)==2
            C = [coords zeros(size(coords,1)) ]
        else
            C = coords
        end
        this = new(C, nx, ny, shape, tag, id)
        return this
    end
end

copy(bl::Block2D) = Block2D(copy(bl.coords), nx=bl.nx, ny=bl.ny, shape=bl.shape, tag=bl.tag)


"""
`Block3D(coords, [nx=1,] [ny=1,] [nz=1,] [shape=HEX8,] [tag=""] )`

Generates a block object for the mesh generation of 3D meshes.
"""
type Block3D <: Block
    coords::Array{Float64,2}
    nx::Int64
    ny::Int64
    nz::Int64
    shape::ShapeType
    tag::AbstractString
    id::Int64
    function Block3D(coords; nx=1, ny=1, nz=1, shape=HEX8, tag="", id=-1)
        C    = size(coords,1)==2? box_coords(vec(coords[1,:]), vec(coords[2,:])): coords
        this = new(C, nx, ny, nz, shape, tag, id)
        return this
    end
end

copy(bl::Block3D) = Block3D(copy(bl.coords), nx=bl.nx, ny=bl.ny, nz=bl.nz, shape=bl.shape, tag=bl.tag)



# Functions for blocks

"""
`move(block, [x=0.0,] [y=0.0,] [z=0.0])` 

Changes de coordinates of a `block`. Also returns a reference.
"""
function move(bl::Block;x=0.0, y=0.0, z=0.0)
    n = size(bl.coords, 1)
    bl.coords[1:n, 1] += x
    bl.coords[1:n, 2] += y
    bl.coords[1:n, 3] += z
    return bl
end


"""
`array(block, [n=1,] [x=0.0,] [y=0.0,] [z=0.0])` 

Creates `n-1` copies of a `block` separated by distances `x`, `y` and/or `z` along respective axes.
"""
function array(bl::Block; n=1, x=0.0, y=0.0, z=0.0)
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

Rotates a `block` according to a `base` point, a `axis` vector and an `angle`.
"""
function rotate(bl::Block; base=[0.,0,0], axis=[0.,0,1], angle=90.0 )

    length(axis)==2 && ( axis=vcat(axis, 0.0) )
    length(base)==2 && ( base=vcat(base, 0.0) )

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


function box_coords{T1<:Number, T2<:Number}(C1::Array{T1,1}, C2::Array{T2,1})
    C = Array(Float64, 8, 3)
    x1 = C1[1]
    y1 = C1[2]
    lx = C2[1] - C1[1]
    ly = C2[2] - C1[2]

    if length(C1)==2
        return [
            x1     y1     0.0
            x1+lx  y1     0.0
            x1+lx  y1+ly  0.0
            x1     y1+ly  0.0 ]
    else
        z1 = C1[3]
        lz = C2[3] - C1[3]
        return [
            x1      y1      z1 
            x1+lx   y1      z1 
            x1+lx   y1+ly   z1 
            x1      y1+ly   z1 
            x1      y1      z1+lz 
            x1+lx   y1      z1+lz 
            x1+lx   y1+ly   z1+lz 
            x1      y1+ly   z1+lz ]
    end
end

# Splits a 2D block
# TODO: replace msh::Mesh by points, bpoints and cells
# TODO: optimize matrix products
function split_block(bl::Block2D, msh::Mesh)
    nx, ny = bl.nx, bl.ny
    shape  = bl.shape # cell shape
    bshape = size(bl.coords,1)==4? QUAD4:QUAD8 # block shape

    if shape==QUAD4
        p_arr = Array(Point, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = shape_func(bshape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); 
                    push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:ny
            for i = 1:nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]
                p4 = p_arr[i  , j+1]

                cell = Cell(shape, [p1, p2, p3, p4], bl.tag)
                push!(msh.cells, cell)
            end
        end
        return
    end

    if shape == QUAD8 || shape == QUAD9
        p_arr = Array(Point, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                if shape==QUAD8 && iseven(i) && iseven(j) continue end

                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = shape_func(bshape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny
            for i = 1:2:2*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                if shape==QUAD8
                    cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8], bl.tag)
                else
                    p9   = p_arr[i+1, j+1]
                    cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9], bl.tag)
                end
                push!(msh.cells, cell)
            end
        end
        return
    end

    if shape == QUAD12
        p_arr = Array(Point, 3*nx+1, 3*ny+1)
        for j = 1:3*ny+1
            for i = 1:3*nx+1
                if shape==QUAD12 && (i-1)%3>0 && (j-1)%3>0 continue end

                r = ((2/3)/nx)*(i-1) - 1.0
                s = ((2/3)/ny)*(j-1) - 1.0
                N = shape_func(bshape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, 3*nx+1) || j in (1, 3*ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:3:3*ny
            for i = 1:3:3*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+3, j  ]
                p3 = p_arr[i+3, j+3]
                p4 = p_arr[i  , j+3]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+3, j+1]
                p7 = p_arr[i+2, j+3]
                p8 = p_arr[i  , j+2]

                p9  = p_arr[i+2, j  ]
                p10 = p_arr[i+3, j+2]
                p11 = p_arr[i+1, j+3]
                p12 = p_arr[i  , j+1]

                cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], bl.tag)
                push!(msh.cells, cell)
            end
        end
        return
    end

    if shape == TRI3
        p_arr = Array(Point, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = shape_func(bshape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:ny
            for i = 1:nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]
                p4 = p_arr[i  , j+1]

                cell1 = Cell(shape, [p1, p2, p3], bl.tag)
                cell2 = Cell(shape, [p4, p1, p3], bl.tag)
                push!(msh.cells, cell1)
                push!(msh.cells, cell2)
            end
        end
        return
    end

    if shape == TRI6

        #=   4       7       3
               @-----@-----@
               |         / |
               |       /   |
             8 @     @     @ 6
               |   /  9    |
               | /         |
               @-----@-----@
             1       5       2     =#

        p_arr = Array(Point, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = shape_func(bshape, [r, s])
                C = round(N'*bl.coords, 8)
                C = reshape(C, 3)
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    p = get_point(msh.bpoints, C)
                    if p==nothing
                        p = Point(C); push!(msh.points, p)
                        msh.bpoints[hash(p)] = p
                    end
                else
                    p = Point(C); push!(msh.points, p)
                end
                p_arr[i,j] = p
            end
        end

        for j = 1:2:2*ny
            for i = 1:2:2*nx
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p4 = p_arr[i  , j+2]

                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p7 = p_arr[i+1, j+2]
                p8 = p_arr[i  , j+1]

                p9   = p_arr[i+1, j+1]

                cell1 = Cell(shape, [p1, p2, p3, p5, p6, p9], bl.tag)
                cell2 = Cell(shape, [p4, p1, p3, p8, p9, p7], bl.tag)
                push!(msh.cells, cell1)
                push!(msh.cells, cell2)
            end
        end
        return
    end

    error("block: Can not discretize using shape $shape")
end


function split_block(bl::Block3D, msh::Mesh)
    nx, ny, nz = bl.nx, bl.ny, bl.nz
    shape  = bl.shape
    bshape = size(bl.coords,1)==8? HEX8:HEX20 # block shape

    if shape==HEX8
        p_arr = Array(Point, nx+1, ny+1, nz+1)
        for k = 1:nz+1
            for j = 1:ny+1
                for i = 1:nx+1
                    r = (2.0/nx)*(i-1) - 1.0
                    s = (2.0/ny)*(j-1) - 1.0
                    t = (2.0/nz)*(k-1) - 1.0
                    N = shape_func(bshape, [r, s, t])
                    C = round(N'*bl.coords, 8)
                    C = reshape(C, 3)
                    p::Any = nothing
                    if i in (1, nx+1) || j in (1, ny+1) || k in (1, nz+1)
                        p = get_point(msh.bpoints, C)
                        if p==nothing
                            p = Point(C); push!(msh.points, p)
                            msh.bpoints[hash(p)] = p
                        end
                    else
                        p = Point(C); push!(msh.points, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:nz
            for j = 1:ny
                for i = 1:nx
                    conn = [
                        p_arr[i  , j  , k  ],
                        p_arr[i+1, j  , k  ],
                        p_arr[i+1, j+1, k  ],
                        p_arr[i  , j+1, k  ],
                        p_arr[i  , j  , k+1],
                        p_arr[i+1, j  , k+1],
                        p_arr[i+1, j+1, k+1],
                        p_arr[i  , j+1, k+1]]

                    cell = Cell(shape, conn, bl.tag)
                    push!(msh.cells, cell)
                end
            end
        end
    elseif shape == HEX20
        p_arr = Array(Point, 2*nx+1, 2*ny+1, 2*nz+1)
        for k = 1:2*nz+1
            for j = 1:2*ny+1
                for i = 1:2*nx+1
                    if iseven(i) && iseven(j) continue end
                    if iseven(j) && iseven(k) continue end
                    if iseven(k) && iseven(i) continue end

                    r = (1.0/nx)*(i-1) - 1.0
                    s = (1.0/ny)*(j-1) - 1.0
                    t = (1.0/nz)*(k-1) - 1.0
                    N = shape_func(bshape, [r, s, t])
                    C = round(N'*bl.coords, 8)
                    C = reshape(C, 3)
                    p::Any = nothing
                    if i in (1, 2*nx+1) || j in (1, 2*ny+1) || k in (1, 2*nz+1)
                        p = get_point(msh.bpoints, C)
                        if p==nothing
                            p = Point(C); push!(msh.points, p)
                            msh.bpoints[hash(p)] = p
                        end
                    else
                        p = Point(C); push!(msh.points, p)
                    end
                    p_arr[i,j,k] = p
                end
            end
        end

        for k = 1:2:2*nz
            for j = 1:2:2*ny
                for i = 1:2:2*nx
                    conn = [
                        p_arr[i  , j  , k  ],
                        p_arr[i+2, j  , k  ],
                        p_arr[i+2, j+2, k  ],
                        p_arr[i  , j+2, k  ],
                        p_arr[i  , j  , k+2],
                        p_arr[i+2, j  , k+2],
                        p_arr[i+2, j+2, k+2],
                        p_arr[i  , j+2, k+2],
                                            
                        p_arr[i+1, j  , k  ],
                        p_arr[i+2, j+1, k  ],
                        p_arr[i+1, j+2, k  ],
                        p_arr[i  , j+1, k  ],
                        p_arr[i+1, j  , k+2],
                        p_arr[i+2, j+1, k+2],
                        p_arr[i+1, j+2, k+2],
                        p_arr[i  , j+1, k+2],
                                           
                        p_arr[i  , j  , k+1],
                        p_arr[i+2, j  , k+1],
                        p_arr[i+2, j+2, k+1],
                        p_arr[i  , j+2, k+1]]

                    cell = Cell(shape, conn, bl.tag)
                    push!(msh.cells, cell)
                end
            end
        end

    end
end


function split_block(bl::BlockTruss, msh::Mesh)
    n = size(bl.coords, 1) # number of points
    m = size(bl.conns , 1) # number of truss cells
    p_arr = Array(Point, n)
    for i=1:n
        C = reshape(bl.coords[i,:], 3)
        p = get_point(msh.bpoints, C)
        if p==nothing; 
            p = Point(C) 
            msh.bpoints[hash(p)] = p
            push!(msh.points, p)
        end
        p_arr[i] = p
    end
    for i=1:m
        p1 = p_arr[bl.conns[i, 1]]
        p2 = p_arr[bl.conns[i, 2]]
        cell = Cell(bl.shape, [p1, p2], bl.tag)
        push!(msh.cells, cell)
    end
end

function split_block(bl::BlockCoords, msh::Mesh)
    n = size(bl.coords, 1) # number of points
    m = size(bl.conns , 1) # number of cells
    p_arr = Array(Point, n)
    for i=1:n
        C = vec(bl.coords[i,:])
        p = get_point(msh.bpoints, C)
        if p==nothing; 
            p = Point(C) 
            msh.bpoints[hash(p)] = p
            push!(msh.points, p)
        end
        p_arr[i] = p
    end
    
    for i=1:m
        points = [ p_arr[j] for j in bl.conns[i,:] ] 
        #TODO: update shape calculation
        shape = [nothing, LIN2, TRI3, QUAD4, nothing, nothing, nothing ][length(points)]
        cell = Cell(shape, points, bl.tag)
        push!(msh.cells, cell)
    end
end

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

import Base.copy

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

    function BlockTruss(coords::Array{<:Real}, conns::Array{Int64,2}; shape=LIN2, tag="", id=-1)
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
`BlockCoords(coords, conns, [shape=LIN2,] [tag="",])`

Generates a 2D or 3D block object based on a matrix of coordinates and a matrix of connectivities.
"""
type BlockCoords <: Block
    coords::Array{Float64,2}
    conns ::Array{Int64,2}
    tag::AbstractString
    id::Int64

    function BlockCoords(coords::Array{<:Real}, conns::Array{Int64,2}; tag="", id=-1)
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

    function Block2D(coords::Array{<:Real}; nx=1, ny=1, shape=QUAD4, tag="", id=-1)
        shape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9, QUAD12) || error("Block2D: shape must be TRI3, TRI6, QUAD4, QUAD8, QUAD9 or QUAD12")
        if size(coords,1)==2
            C = box_coords(coords[1,:], coords[2,:])
        elseif size(coords,2)==2
            C = [coords zeros(size(coords,1)) ]
        else
            C = coords
        end
        this = new(C, nx, ny, shape, tag, id)
        return this
    end
end


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

    function Block3D(coords::Array{<:Real}; nx=1, ny=1, nz=1, shape=HEX8, tag="", id=-1)
        shape in (TET4, TET10, HEX8, HEX20) || error("Block3D: shape must be TET4, TET10, HEX8 or HEX20")
        C    = size(coords,1)==2? box_coords(coords[1,:], coords[2,:]): coords
        this = new(C, nx, ny, nz, shape, tag, id)
        return this
    end
end

type BlockCylinder <: Block
    coords::Array{Float64,2} # two end points
    r::Float64
    nr::Int64
    n::Int64
    shape::ShapeType # HEX8, HEX20
    tag::AbstractString
    id::Int64

    function BlockCylinder(coords::Array{<:Real}; r=1.0, nr=3, n=2, shape=HEX8, tag="", id=-1)
        size(coords,1) != 2 && error("Invalid coordinates matrix for BlockCylinder")
        nr<2 && error("Invalid nr=$nr value for BlockCylinder")
        shape in (HEX8, HEX20) || error("BlockCylinder: shape must be HEX8 or HEX20")
        this = new(coords, r, nr, n, shape, tag, id)
        return this
    end
end

function box_coords(C1::Array{<:Real,1}, C2::Array{<:Real,1})
    C = Array{Float64}(8, 3)
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
        p_arr = Array{Point}(nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = bshape.func([r, s])
                C = N'*bl.coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, 8)
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
        p_arr = Array{Point}(2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                if shape==QUAD8 && iseven(i) && iseven(j) continue end

                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bshape.func([r, s])
                C = N'*bl.coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, 8)
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
        p_arr = Array{Point}(3*nx+1, 3*ny+1)
        for j = 1:3*ny+1
            for i = 1:3*nx+1
                if shape==QUAD12 && (i-1)%3>0 && (j-1)%3>0 continue end

                r = ((2/3)/nx)*(i-1) - 1.0
                s = ((2/3)/ny)*(j-1) - 1.0
                N = bshape.func([r, s])
                C = N'*bl.coords
                p::Any = nothing
                if i in (1, 3*nx+1) || j in (1, 3*ny+1)
                    C = round.(C, 8)
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
        p_arr = Array{Point}(nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = bshape.func([r, s])
                C = N'*bl.coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, 8)
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

        p_arr = Array{Point}(2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bshape.func([r, s])
                C = N'*bl.coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, 8)
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

    error("block: Can not discretize using shape $(shape.name)")
end

# Split for 3D meshes

function split_block(bl::Block3D, msh::Mesh)
    nx, ny, nz = bl.nx, bl.ny, bl.nz
    shape  = bl.shape 
    #shape == TET10 && error("block: Cannot generate mesh with TET10 cells jet")
    bshape = size(bl.coords,1)==8? HEX8:HEX20 # block shape (not cell shape)

    if shape==HEX8 || shape==TET4
        p_arr = Array{Point}(nx+1, ny+1, nz+1)
        for k = 1:nz+1
            for j = 1:ny+1
                for i = 1:nx+1
                    r = (2.0/nx)*(i-1) - 1.0
                    s = (2.0/ny)*(j-1) - 1.0
                    t = (2.0/nz)*(k-1) - 1.0
                    N = bshape.func([r, s, t])
                    C = N'*bl.coords
                    p::Any = nothing
                    if i in (1, nx+1) || j in (1, ny+1) || k in (1, nz+1)
                        C = round.(C, 8)
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
                    p1 = p_arr[i  , j  , k  ]
                    p2 = p_arr[i+1, j  , k  ]
                    p3 = p_arr[i+1, j+1, k  ]
                    p4 = p_arr[i  , j+1, k  ]
                    p5 = p_arr[i  , j  , k+1]
                    p6 = p_arr[i+1, j  , k+1]
                    p7 = p_arr[i+1, j+1, k+1]
                    p8 = p_arr[i  , j+1, k+1]

                    if shape==HEX8
                        cell = Cell(shape, [p1, p2, p3, p4, p5, p6, p7, p8], bl.tag)
                        push!(msh.cells, cell)
                    end
                    if shape==TET4
                        push!( msh.cells, Cell(shape, [p2, p4, p1, p8], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p1, p5, p8], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p5, p6, p8], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p6, p7, p8], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p3, p4, p8], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p7, p3, p8], bl.tag) )
                    end
                end
            end
        end
        return
    end

    if shape == HEX20 || shape == TET10
        p_arr = Array{Point}(2*nx+1, 2*ny+1, 2*nz+1)
        for k = 1:2*nz+1
            for j = 1:2*ny+1
                for i = 1:2*nx+1
                    if shape==HEX20
                        if iseven(i) && iseven(j) continue end
                        if iseven(j) && iseven(k) continue end
                        if iseven(k) && iseven(i) continue end
                    end

                    r = (1.0/nx)*(i-1) - 1.0
                    s = (1.0/ny)*(j-1) - 1.0
                    t = (1.0/nz)*(k-1) - 1.0
                    N = bshape.func([r, s, t])
                    C = N'*bl.coords
                    p::Any = nothing
                    if i in (1, 2*nx+1) || j in (1, 2*ny+1) || k in (1, 2*nz+1)
                        C = round.(C, 8)
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

                    if shape == HEX20
                        cell = Cell(shape, conn, bl.tag)
                        push!(msh.cells, cell)
                    end
                    if shape == TET10

                        p1  = p_arr[i  , j  , k  ]
                        p2  = p_arr[i+2, j  , k  ]
                        p3  = p_arr[i+2, j+2, k  ]
                        p4  = p_arr[i  , j+2, k  ]
                        p5  = p_arr[i  , j  , k+2]
                        p6  = p_arr[i+2, j  , k+2]
                        p7  = p_arr[i+2, j+2, k+2]
                        p8  = p_arr[i  , j+2, k+2]

                        p9  = p_arr[i+1, j  , k  ]
                        p10 = p_arr[i+2, j+1, k  ]
                        p11 = p_arr[i+1, j+2, k  ]
                        p12 = p_arr[i  , j+1, k  ]
                        p13 = p_arr[i+1, j  , k+2]
                        p14 = p_arr[i+2, j+1, k+2]
                        p15 = p_arr[i+1, j+2, k+2]
                        p16 = p_arr[i  , j+1, k+2]

                        p17 = p_arr[i  , j  , k+1]
                        p18 = p_arr[i+2, j  , k+1]
                        p19 = p_arr[i+2, j+2, k+1]
                        p20 = p_arr[i  , j+2, k+1]

                        p21 = p_arr[i+1, j+1, k  ]
                        p22 = p_arr[i+1, j+1, k+2]
                        p23 = p_arr[i+1, j  , k+1]
                        p24 = p_arr[i+2, j+1, k+1]
                        p25 = p_arr[i+1, j+2, k+1]
                        p26 = p_arr[i  , j+1, k+1]
                        p27 = p_arr[i+1, j+1, k+1]

                        push!( msh.cells, Cell(shape, [p2, p4, p1, p8, p21, p12, p9, p27, p20, p26], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p1, p5, p8, p9, p17, p23, p27, p26, p16], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p5, p6, p8, p23, p13, p18, p27, p16, p22], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p6, p7, p8, p18, p14, p24, p27, p22, p15], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p3, p4, p8, p10, p11, p21, p27, p25, p20], bl.tag) )
                        push!( msh.cells, Cell(shape, [p2, p7, p3, p8, p24, p19, p10, p27, p15, p25], bl.tag) )
                    end
                end
            end
        end
        return
    end
    error("block: Can not discretize using shape $(shape.name)")
end


function split_block(bl::BlockTruss, msh::Mesh)
    n = size(bl.coords, 1) # number of points
    m = size(bl.conns , 1) # number of truss cells
    p_arr = Array{Point}(n)
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
    p_arr = Array{Point}(n)
    for i=1:n
        C = bl.coords[i,:]
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



function split_block(bl::BlockCylinder, msh::Mesh)

    nx1 = round(Int, bl.nr/3)
    nx2 = bl.nr - nx1
    shape2D = bl.shape==HEX8 ? QUAD4 : QUAD8

    coords = bl.r*[ 0 0; 1/3 0; 1/3 1/3; 0 1/3; 1/6 0; 1/3 1/6; 1/6 1/3; 0 1/6 ]
    bl1 = Block2D(coords, nx=nx1, ny=nx1, shape= shape2D, tag=bl.tag)

    s45  = sin(45*pi/180)
    c45  = s45
    s225 = sin(22.5*pi/180)
    c225 = cos(22.5*pi/180)

    coords = bl.r*[ 1/3 0; 1 0; c45 s45; 1/3 1/3; 2/3 0; c225 s225; (c45+1/3)/2 (s45+1/3)/2; 1/3 1/6 ]
    bl2    = Block2D(coords, nx=nx2, ny=nx1, shape= shape2D, tag=bl.tag)

    coords = bl.r*[ 0 1/3; 1/3 1/3; c45 s45; 0 1; 1/6 1/3; (c45+1/3)/2 (s45+1/3)/2; s225 c225; 0 2/3 ]
    bl3    = Block2D(coords, nx=nx1, ny=nx2, shape= shape2D, tag=bl.tag)

    blocks = [bl1, bl2, bl3 ]

    # polar and move
    blocks = polar(blocks, n=4)
    move!(blocks, x=bl.coords[1,1], y=bl.coords[1,2], z=bl.coords[1,3])

    # extrude
    len      = norm(bl.coords[1,:] - bl.coords[2,:])
    blocks3D = extrude(blocks, len=len, n=bl.n)

    # rotation
    zv    = [0.0, 0.0, 1.0]
    axis  = bl.coords[2,:] - bl.coords[1,:]
    angle = acos( dot(zv, axis)/(norm(zv)*norm(axis)) )*180/pi
    raxis = cross(zv, axis)
    norm(raxis)>1e-10 && rotate!(blocks3D, base=bl.coords[1,:], axis=raxis, angle=angle)

    # split
    for bl in blocks3D
        split_block(bl, msh)
    end

end

# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh


include("block_inset.jl")

"""
`BlockTruss(coords, conns, [shape=LIN2,] [tag="",])`

Generates a block object for the mesh generation of trusses:
"""
mutable struct BlockTruss <: Block
    points::Array{Point,1}
    conns ::Array{Int64,2}
    shape ::ShapeType
    cellshape::ShapeType
    tag::String
    id::Int64

    function BlockTruss(coords::Array{<:Real}, conns::Array{Int64,2}; cellshape=LIN2, tag="", id=-1)
        ncols = size(coords,2)
        size(conns,2)==2 || error("BlockTruss: Invalid size for connectivities matrix")

        local points

        if ncols==2
            points = [ Point(coords[i,1], coords[i,2], 0.0) for i=1:size(coords,1) ]
        elseif ncols==3
            points = [ Point(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]
        else
            error("BlockTruss: Invalid coordinates matrix")
        end

        for (i,p) in enumerate(points); p.id=i end

        return new(points, conns, POLYV, cellshape, tag, id)
    end
end


function box_points(C1::Array{<:Real,1}, C2::Array{<:Real,1})
    x1 = C1[1]
    y1 = C1[2]
    lx = C2[1] - C1[1]
    ly = C2[2] - C1[2]

    if length(C1)==2
        return [
                 Point(x1   , y1   , 0.0),
                 Point(x1+lx, y1   , 0.0),
                 Point(x1+lx, y1+ly, 0.0),
                 Point(x1   , y1+ly, 0.0),
                ]
    else
        z1 = C1[3]
        lz = C2[3] - C1[3]
        return [
                 Point(x1   , y1   , z1 ),
                 Point(x1+lx, y1   , z1 ),
                 Point(x1+lx, y1+ly, z1 ),
                 Point(x1   , y1+ly, z1 ),
                 Point(x1   , y1   , z1+lz ),
                 Point(x1+lx, y1   , z1+lz ),
                 Point(x1+lx, y1+ly, z1+lz ),
                 Point(x1   , y1+ly, z1+lz ),
                ]
    end
end


"""
`BlockCoords(coords, conns, [cellshape=LIN2,] [tag="",])`

Generates a 2D or 3D block object based on a matrix of coordinates and a matrix of connectivities.
"""
mutable struct BlockCoords <: Block
    coords::Array{Float64,2}
    conns ::Array{Int64,2}
    tag::String
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
`Block2D(coords, [nx=1,] [ny=1,] [cellshape=QUAD4,] [tag=""] )`

Generates a block object for the mesh generation of 2D meshes.
`shape` can be TRI3, TRI6, QUAD4, QUAD8.
"""
mutable struct Block2D <: Block
    points::Array{Point,1}
    shape::ShapeType
    cellshape::ShapeType
    nx::Int64
    ny::Int64
    tag::String
    id::Int64

    function Block2D(coords::Array{<:Real}; nx::Int=1, ny::Int=1, cellshape::ShapeType=QUAD4, tag="", id=-1, shape=nothing)
        if shape != nothing
            @warn "Block2D: argument shape was deprecated. Please use cellshape instead"
            cellshape = shape
        end
        cellshape in (TRI3, TRI6, QUAD4, QUAD8, QUAD9, QUAD12) || error("Block2D: shape must be TRI3, TRI6, QUAD4, QUAD8, QUAD9 or QUAD12")

        if size(coords,1)==2
            points = box_points(coords[1,:], coords[2,:])
        elseif size(coords,2)==2
            points = [ Point(coords[i,1], coords[i,2], 0.0) for i=1:size(coords,1) ]
        else
            points = [ Point(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]
        end

        npoints = length(points)
        shape = npoints==4 ? QUAD4 : QUAD8
        for (i,p) in enumerate(points); p.id=i end

        return new(points, shape, cellshape, nx, ny, tag, id)
    end
end


"""
`Block3D(coords, [nx=1,] [ny=1,] [nz=1,] [cellshape=HEX8,] [tag=""] )`

Generates a block object for the mesh generation of 3D meshes.
"""
mutable struct Block3D <: Block
    points::Array{Point,1}
    shape::ShapeType
    cellshape::ShapeType
    nx::Int64
    ny::Int64
    nz::Int64
    tag::String
    id::Int64

    function Block3D(coords::Array{<:Real}; nx=1, ny=1, nz=1, cellshape=HEX8, tag="", id=-1, shape=nothing)
        if shape != nothing
            @warn "Block3D: argument `shape` was deprecated. Please use cellshape instead"
            cellshape = shape
        end
        cellshape in (TET4, TET10, HEX8, HEX20) || error("Block3D: shape must be TET4, TET10, HEX8 or HEX20")

        if size(coords,1)==2
            points = box_points(coords[1,:], coords[2,:])
        else
            points = [ Point(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]
        end

        npoints = length(points)
        shape = npoints==8 ? HEX8 : HEX20
        for (i,p) in enumerate(points); p.id=i end

        return new(points, shape, cellshape, nx, ny, nz, tag, id)
    end
end

#type BlockCylinder <: Block
    #points::Array{Point,1}
    #shape::ShapeType # LIN2
    #cellshape::ShapeType # HEX8, HEX20
    #function Block3D(coords::Array{<:Real}; nx=1, ny=1, nz=1, shape=HEX8, tag="", id=-1)
        #shape in (TET4, TET10, HEX8, HEX20) || error("Block3D: shape must be TET4, TET10, HEX8 or HEX20")
        #C    = size(coords,1)==2 ? box_coords(coords[1,:], coords[2,:]) : coords
        #this = new(C, nx, ny, nz, shape, tag, id)
        #return this
    #end
#end

#=
mutable struct BlockArc <: Block
    coords::Array{Float64,2}
    r1::Float64
    th::Float64
    n1::Int64
    n2::Int64
    shape::ShapeType
    tag::String
    id::Int64

    function BlockArc(coords::Array{<:Real}; r1=1.0, th=45, nr=3, n=2, shape=QUAD4, tag="", id=-1)
        size(coords,1) == 2 || error("Invalid coordinates matrix for BlockArc")
        shape.ndim==2 || error("BlockArc: shape must be a 2d shape")
        #r2 = norm(coords[2,:]-coords[1,:])
        this = new(coords, r1, th, nr, n, shape, tag, id)
        return this
    end

    function BlockArc(center::Array{<:Real}=[0.0,0.0]; r1=1.0, r2=2.0, th1=45, th2=90, nr=3, n=2, shape=QUAD4, tag="", id=-1)
        size(coords,1) == 1 || error("Invalid coordinates matrix for BlockArc")
        shape.ndim==2 || error("BlockArc: shape must be a 2d shape")
        th1r = deg2rad(th1)
        V = [cos(th1r),sin(th1r)]
        th2r = deg2rad(th2)
        X1 = center + r1*V
        X2 = center + r2*V
        coords = [X1'; X2']
        this = new(coords, r1, th, nr, n, shape, tag, id)
        return this
    end
end
=#

mutable struct BlockCylinder <: Block
    points::Array{Point,1}
    shape::ShapeType # LIN2
    cellshape::ShapeType # HEX8, HEX20
    #coords::Array{Float64,2} # two end points
    r::Float64
    nr::Int64
    n::Int64
    tag::String
    id::Int64

    function BlockCylinder(coords::Array{<:Real}; r=1.0, nr=3, n=2, cellshape=HEX8, tag="", id=-1)
        size(coords,1) != 2 && error("Invalid coordinates matrix for BlockCylinder")
        nr<2 && error("Invalid nr=$nr value for BlockCylinder")
        cellshape in (HEX8, HEX20) || error("BlockCylinder: cellshape must be HEX8 or HEX20")
        points = [ Point(coords[i,1], coords[i,2], coords[i,3]) for i=1:size(coords,1) ]
        return new(points, LIN2, cellshape, r, nr, n, tag, id)
    end
end

# Splits a 2D block
# TODO: replace msh::Mesh by points, bpoints and cells
# TODO: optimize matrix products
function split_block(bl::Block2D, msh::Mesh)
    nx, ny = bl.nx, bl.ny
    shape  = bl.shape # cell shape
    coords = getcoords(bl.points)
    cellshape = bl.cellshape
    quadrilat_bshape = true

    if cellshape==QUAD4
        p_arr = Array{Point}(undef, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, digits=8)
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

                cell = Cell(cellshape, [p1, p2, p3, p4], tag=bl.tag)
                push!(msh.cells, cell)
            end
        end
        return
    end

    if cellshape == QUAD8 || cellshape == QUAD9
        p_arr = Array{Point}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                if cellshape==QUAD8 && iseven(i) && iseven(j) continue end

                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, digits=8)
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

                if cellshape==QUAD8
                    cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=bl.tag)
                else
                    p9   = p_arr[i+1, j+1]
                    cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8, p9], tag=bl.tag)
                end
                push!(msh.cells, cell)
            end
        end
        return
    end

    if cellshape == QUAD12
        p_arr = Array{Point}(undef, 3*nx+1, 3*ny+1)
        for j = 1:3*ny+1
            for i = 1:3*nx+1
                if cellshape==QUAD12 && (i-1)%3>0 && (j-1)%3>0 continue end

                r = ((2/3)/nx)*(i-1) - 1.0
                s = ((2/3)/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 3*nx+1) || j in (1, 3*ny+1)
                    C = round.(C, digits=8)
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

                cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], tag=bl.tag)
                push!(msh.cells, cell)
            end
        end
        return
    end

    if cellshape == TRI3 && quadrilat_bshape
        p_arr = Array{Point}(undef, nx+1, ny+1)
        for j = 1:ny+1
            for i = 1:nx+1
                r = (2.0/nx)*(i-1) - 1.0
                s = (2.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1)
                    C = round.(C, digits=8)
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

                cell1 = Cell(cellshape, [p1, p2, p3], tag=bl.tag)
                cell2 = Cell(cellshape, [p4, p1, p3], tag=bl.tag)
                push!(msh.cells, cell1)
                push!(msh.cells, cell2)
            end
        end
        return
    end

    if cellshape == TRI3 && !quadrilat_bshape
        p_arr = Array{Point}(undef, nx+1, nx+1)
        for j = 1:ny+1
            for i = 1:nx+1
                j>i && continue
                r = 1 - (1.0/nx)*(i-1)
                s =  (1.0/ny)*(j-1)
                #r = (2.0/nx)*(i-1) - 1.0
                #s = (2.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, nx+1) || j in (1, ny+1) || i==j
                    C = round.(C, digits=8)
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
                j>i && continue
                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+1, j  ]
                p3 = p_arr[i+1, j+1]

                cell1 = Cell(cellshape, [p1, p2, p3], tag=bl.tag)
                push!(msh.cells, cell1)
                
                if i>j
                    p4 = p_arr[i  , j+1]
                    cell2 = Cell(cellshape, [p4, p1, p3], tag=bl.tag)
                    push!(msh.cells, cell2)
                end
            end
        end
        return
    end

    if cellshape == TRI6 && quadrilat_bshape

        #=   4       7       3
               @-----@-----@
               |         / |
               |       /   |
             8 @     @     @ 6
               |   /  9    |
               | /         |
               @-----@-----@
             1       5       2     =#

        p_arr = Array{Point}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                r = (1.0/nx)*(i-1) - 1.0
                s = (1.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1)
                    C = round.(C, digits=8)
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

                cell1 = Cell(cellshape, [p1, p2, p3, p5, p6, p9], tag=bl.tag)
                cell2 = Cell(cellshape, [p4, p1, p3, p8, p9, p7], tag=bl.tag)
                push!(msh.cells, cell1)
                push!(msh.cells, cell2)
            end
        end
        return
    end

    if cellshape == TRI6 && !quadrilat_bshape

        #=   4       7       3
               @-----@-----@
               |         / |
               |       /   |
             8 @     @     @ 6
               |   /  9    |
               | /         |
               @-----@-----@
             1       5       2     =#

        p_arr = Array{Point}(undef, 2*nx+1, 2*ny+1)
        for j = 1:2*ny+1
            for i = 1:2*nx+1
                j>i && continue
                r = (1.0/nx)*(i-1) # TODO
                s = (1.0/ny)*(j-1)
                #r = (1.0/nx)*(i-1) - 1.0
                #s = (1.0/ny)*(j-1) - 1.0
                N = bl.shape.func([r, s])
                C = N'*coords
                p::Any = nothing
                if i in (1, 2*nx+1) || j in (1, 2*ny+1) || i==j
                    C = round.(C, digits=8)
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
                j>i && continue

                p1 = p_arr[i  , j  ]
                p2 = p_arr[i+2, j  ]
                p3 = p_arr[i+2, j+2]
                p5 = p_arr[i+1, j  ]
                p6 = p_arr[i+2, j+1]
                p9 = p_arr[i+1, j+1]

                cell1 = Cell(cellshape, [p1, p2, p3, p5, p6, p9], tag=bl.tag)
                push!(msh.cells, cell1)

                if i>j
                    p4 = p_arr[i  , j+2]
                    p7 = p_arr[i+1, j+2]
                    p8 = p_arr[i  , j+1]
                    cell2 = Cell(cellshape, [p4, p1, p3, p8, p9, p7], tag=bl.tag)
                    push!(msh.cells, cell2)
                end
            end
        end
        return
    end

    error("block: Can not discretize using shape $(cellshape.name)")
end

# Split for 3D meshes

function split_block(bl::Block3D, msh::Mesh)
    nx, ny, nz = bl.nx, bl.ny, bl.nz
    shape  = bl.shape 
    coords = getcoords(bl.points)
    cellshape = bl.cellshape

    if cellshape==HEX8 || cellshape==TET4
        p_arr = Array{Point}(undef, nx+1, ny+1, nz+1)
        for k = 1:nz+1
            for j = 1:ny+1
                for i = 1:nx+1
                    r = (2.0/nx)*(i-1) - 1.0
                    s = (2.0/ny)*(j-1) - 1.0
                    t = (2.0/nz)*(k-1) - 1.0
                    N = bl.shape.func([r, s, t])
                    C = N'*coords
                    p::Any = nothing
                    if i in (1, nx+1) || j in (1, ny+1) || k in (1, nz+1)
                        C = round.(C, digits=8)
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

                    if cellshape==HEX8
                        cell = Cell(cellshape, [p1, p2, p3, p4, p5, p6, p7, p8], tag=bl.tag)
                        push!(msh.cells, cell)
                    end
                    if cellshape==TET4
                        push!( msh.cells, Cell(cellshape, [p2, p4, p1, p8], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p1, p5, p8], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p5, p6, p8], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p6, p7, p8], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p3, p4, p8], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p7, p3, p8], tag=bl.tag) )
                    end
                end
            end
        end
        return
    end

    if cellshape == HEX20 || cellshape == TET10
        p_arr = Array{Point}(undef, 2*nx+1, 2*ny+1, 2*nz+1)
        for k = 1:2*nz+1
            for j = 1:2*ny+1
                for i = 1:2*nx+1
                    if cellshape==HEX20
                        if iseven(i) && iseven(j) continue end
                        if iseven(j) && iseven(k) continue end
                        if iseven(k) && iseven(i) continue end
                    end

                    r = (1.0/nx)*(i-1) - 1.0
                    s = (1.0/ny)*(j-1) - 1.0
                    t = (1.0/nz)*(k-1) - 1.0
                    N = bl.shape.func([r, s, t])
                    C = N'*coords
                    p::Any = nothing
                    if i in (1, 2*nx+1) || j in (1, 2*ny+1) || k in (1, 2*nz+1)
                        C = round.(C, digits=8)
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

                    if cellshape == HEX20
                        cell = Cell(cellshape, conn, tag=bl.tag)
                        push!(msh.cells, cell)
                    end
                    if cellshape == TET10

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

                        push!( msh.cells, Cell(cellshape, [p2, p4, p1, p8, p21, p12, p9, p27, p20, p26], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p1, p5, p8, p9, p17, p23, p27, p26, p16], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p5, p6, p8, p23, p13, p18, p27, p16, p22], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p6, p7, p8, p18, p14, p24, p27, p22, p15], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p3, p4, p8, p10, p11, p21, p27, p25, p20], tag=bl.tag) )
                        push!( msh.cells, Cell(cellshape, [p2, p7, p3, p8, p24, p19, p10, p27, p15, p25], tag=bl.tag) )
                    end
                end
            end
        end
        return
    end
    error("block: Can not discretize using shape $(cellshape.name)")
end


function split_block(bl::BlockTruss, msh::Mesh)
    n = length(bl.points) # number of points
    m = size(bl.conns,1)  # number of truss cells
    p_arr = Array{Point}(undef, n)

    for (i,point) in enumerate(bl.points)
        C  = [ point.x, point.y, point.z ]
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
        cell = Cell(bl.cellshape, [p1, p2], tag=bl.tag)
        push!(msh.cells, cell)
    end
end


function split_block(bl::BlockCoords, msh::Mesh)
    n = size(coords, 1) # number of points
    m = size(bl.conns , 1) # number of cells
    p_arr = Array{Point}(undef, n)
    for i=1:n
        C = coords[i,:]
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
        cellshape = [nothing, LIN2, TRI3, QUAD4, nothing, nothing, nothing ][length(points)]
        cell = Cell(cellshape, points, tag=bl.tag)
        push!(msh.cells, cell)
    end
end


function split_block(bl::BlockCylinder, msh::Mesh)

    nx1 = round(Int, bl.nr/3)
    nx2 = bl.nr - nx1
    shape2D = bl.cellshape==HEX8 ? QUAD4 : QUAD8

    # constructing quadratic blocks
    coords = bl.r*[ 0 0; 1/3 0; 1/3 1/3; 0 1/3; 1/6 0; 1/3 1/6; 1/6 1/3; 0 1/6 ]
    bl1 = Block2D(coords, nx=nx1, ny=nx1, cellshape= shape2D, tag=bl.tag)

    s45  = sin(45*pi/180)
    c45  = s45
    s225 = sin(22.5*pi/180)
    c225 = cos(22.5*pi/180)

    coords = bl.r*[ 1/3 0; 1 0; c45 s45; 1/3 1/3; 2/3 0; c225 s225; (c45+1/3)/2 (s45+1/3)/2; 1/3 1/6 ]
    bl2    = Block2D(coords, nx=nx2, ny=nx1, cellshape= shape2D, tag=bl.tag)

    coords = bl.r*[ 0 1/3; 1/3 1/3; c45 s45; 0 1; 1/6 1/3; (c45+1/3)/2 (s45+1/3)/2; s225 c225; 0 2/3 ]
    bl3    = Block2D(coords, nx=nx1, ny=nx2, cellshape= shape2D, tag=bl.tag)

    blocks = [bl1, bl2, bl3 ]

    # polar and move
    blocks = polar(blocks, n=4)
    coords = getcoords(bl.points)
    move!(blocks, dx=coords[1,1], dy=coords[1,2], dz=coords[1,3])

    # extrude
    len      = norm(coords[1,:] - coords[2,:])
    blocks3D = extrude(blocks, len=len, n=bl.n)

    # rotation
    zv    = [0.0, 0.0, 1.0]
    axis  = coords[2,:] - coords[1,:]
    angle = acos( dot(zv, axis)/(norm(zv)*norm(axis)) )*180/pi
    raxis = cross(zv, axis)
    norm(raxis)>1e-10 && rotate!(blocks3D, base=coords[1,:], axis=raxis, angle=angle)

    # split
    for bl in blocks3D
        split_block(bl, msh)
    end

end

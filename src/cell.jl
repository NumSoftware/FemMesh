# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

abstract type AbstractCell end
abstract type Block<:AbstractCell end

# Cell
# ====

"""
A geometry type that represent a definite line, area or volume.
"""
mutable struct Cell<:AbstractCell
    shape  ::ShapeType
    points ::Array{Point,1}
    tag    ::TagType
    id     ::Integer
    ndim   ::Integer
    quality::Float64              # quality index: surf/(reg_surf) 
    embedded::Bool                # flag for embedded cells
    crossed::Bool                 # flag if cell crossed by linear inclusion
    ocell  ::Union{AbstractCell,Nothing}     # owner cell if this cell is a face/edge
    nips   ::Int                  # number of integration points if required
    iptag  ::TagType # tag for integration points if required
    linked_cells::Array{Cell,1}   # neighbor cells in case of joint cell
    function Cell(shape::ShapeType, points::Array{Point,1}, tag::TagType=0, ocell=nothing)
        this = new(shape, points, tag, -1)
        this.ndim = 0
        this.quality = 0.0
        this.embedded= false
        this.crossed = false
        this.ocell   = ocell
        this.nips  = 0
        this.iptag = 0
        this.linked_cells = []
        return this
    end
end

const Face=Cell


### Cell methods

Base.hash(c::Cell) = sum(hash(p) for p in c.points)

function getcoords(c::Cell, ndim=3)::Array{Float64,2}
    n = length(c.points)
    C = Array{Float64}(undef, n, ndim)
    for (i,p) in enumerate(c.points)
        C[i,1] = p.x
        if ndim>1 C[i,2] = p.y end
        if ndim>2 C[i,3] = p.z end
    end
    return C
end

# Return all points in cells
function get_points(cells::Array{Cell,1})::Array{Point,1}
    #pointsd = Dict{UInt64, Point}()
    #for cell in cells
        #for point in cell.points
            #pointsd[hash(point)] = point
        #end
    #end
    #points = [ values(pointsd)... ]

    points = Set{Point}()
    for cell in cells
        for point in cell.points
            push!(points, point)
        end
    end
    return Point[point for point in points]
end

# Updates cell internal data
function update!(c::Cell)::Nothing
    c.quality = cell_quality(c)
    return nothing
end


# Index operator for a collection of elements
# This function is not type stable
function Base.getindex(cells::Array{Cell,1}, s::Symbol)
    s == :all      && return cells
    s == :solids   && return filter(cell -> cell.shape.family==SOLID_SHAPE, cells)
    s == :lines    && return filter(cell -> cell.shape.family==LINE_SHAPE, cells)
    s == :joints   && return filter(cell -> cell.shape.family==JOINT_SHAPE, cells)
    s == :joints1D && return filter(cell -> cell.shape.family==JOINT1D_SHAPE, cells)
    s == :points   && return get_points(cells)
    error("Cell getindex: Invalid filter symbol $s")
end


function Base.getindex(cells::Array{Cell,1}, cond::Expr)
    condm = fix_comparison_arrays(cond)
    funex = :( (x,y,z,id,tag) -> false )
    funex.args[2].args[2] = condm
    fun = nothing
    try
        fun = eval(funex)
    catch
        error("Cell getindex: Invalid filter expression ", cond)
    end

    result = Cell[]
    for cell in cells
        x = [ p.x for p in cell.points ]
        y = [ p.y for p in cell.points ]
        z = [ p.z for p in cell.points ]
        if Base.invokelatest(fun, x, y, z, cell.id, cell.tag)
            push!(result, cell) 
        end
    end
    return result
end


function tag!(object::Union{Point, Cell}, tag::TagType)
    object.tag = tag
end

function tag!(arr::Array{T,1}, tag::TagType) where T<:Union{Point,Cell}
    for obj in arr
        obj.tag = tag
    end
end

function iptag!(object::Cell, tag::TagType)
    object.iptag = tag
end

function iptag!(arr::Array{Cell,1}, tag::TagType)
    for obj in arr
        obj.iptag = tag
    end
end

# Gets the coordinates of a bounding box for an array of points
function bounding_box(points::Array{Point,1})
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for point in points
        point.x < minx && (minx = point.x)
        point.y < miny && (miny = point.y)
        point.z < minz && (minz = point.z)
        point.x > maxx && (maxx = point.x)
        point.y > maxy && (maxy = point.y)
        point.z > maxz && (maxz = point.z)
    end
    return [ minx miny minz; maxx maxy maxz ]
end

# Gets the coordinates of a bounding box for a cell
function bounding_box(cell::Cell)
    return bounding_box(cell.points)
end

# Gets the coordinates of a bounding box for an array of cells
function bounding_box(cells::Array{Cell,1})
    points = unique( Point[ p for c in cells for p in c.points ] )
    return bounding_box(points)
end


# gets all facets of a cell
function get_faces(cell::AbstractCell)
    faces  = Cell[]

    all_faces_idxs = cell.shape.facet_idxs
    facet_shape    = cell.shape.facet_shape
    
    if facet_shape==() return faces end

    sameshape  = typeof(facet_shape) == ShapeType # check if all facets have the same shape

    # Iteration for each facet
    for (i, face_idxs) in enumerate(all_faces_idxs)
        points = cell.points[face_idxs]
        shape  = sameshape ? facet_shape : facet_shape[i]
        face   = Cell(shape, points, cell.tag, cell)
        push!(faces, face)
    end

    return faces
end

# gets all edges of a cell
function get_edges(cell::Cell)
    if cell.shape.ndim==2 return get_faces(cell) end

    edges  = Cell[]
    all_edge_idxs = cell.shape.edge_idxs

    for edge_idx in all_edge_idxs
        points = cell.points[edge_idx]
        shape  = (LIN2, LIN3, LIN4)[length(points)-1]
        edge   = Cell(shape, points, cell.tag, cell)
        push!(edges, edge)
    end

    return edges
end


# Pseudo determinant of non-square matrices
function norm2(J::Array{Float64,2})

    if ndims(J)==1; return norm(J) end

    r, c = size(J)
    if r==1; return norm(J) end
    if r==2 && c==3
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
        j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian determinant
    end
    if r==c; return det(J) end
    error("No rule to calculate norm2 of a $r x $c matrix")
end

# Returns the volume/area/length of a cell, returns zero if jacobian<=0 at any ip
function cell_extent(c::Cell)
    IP = get_ip_coords(c.shape)
    nip = size(IP,1)
    nldim = c.shape.ndim # cell basic dimension


    # get coordinates matrix
    C = getcoords(c)
    J = Array{Float64}(undef, nldim, size(C,2))

    # calc metric
    vol = 0.0
    for i=1:nip
        R    = vec(IP[i,1:3])
        dNdR = c.shape.deriv(R)

        @gemm J = dNdR*C
        w    = IP[i,4]
        normJ = norm2(J)
        normJ <= 0.0 && return 0.0
        vol += normJ*w
    end
    return vol
end

# Returns the surface/perimeter of a regular element given the volume/area of a cell
function regular_surface(metric::Float64, shape::ShapeType)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ] 
        A = metric
        a = 2.0*√(A/√3.0)
        return 3*a
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ] 
        A = metric
        a = √A
        return 4*a
    end
    if shape in [ PYR5 ] 
        V = metric
        a = ( 3.0*√2.0*V )^(1.0/3.0)
        return (1 + √3.0)*a^2
    end
    if shape in [ TET4, TET10 ] 
        V = metric
        a = ( 6.0*√2.0*V )^(1.0/3.0)
        return √3.0*a^2
    end
    if shape in [ HEX8, HEX20 ] 
        V = metric
        a = V^(1.0/3.0)
        return 6.0*a^2.0
    end
    if shape in [ WED6, WED15 ]
        V = metric
        a2 = (16.0/3.0*V^2)^(1.0/3.0) 
        return (3.0 + √3.0/2.0)*a2
    end
    error("No regular surface/perimeter value for shape $(get_name(shape))")
end

# Returns the area/volume of a regular element given the perimeter/surface
function regular_vol(metric::Float64, shape::ShapeType)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ] 
        p = metric
        a = p/3
        return a^2/4*√3.0    
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ] 
        p = metric
        a = p/4
        return a^2
    end
    #if shape in [ PYR5 ] 
        #V = metric
        #a = ( 3.0*√2.0*V )^(1.0/3.0)
        #return (1 + √3.0)*a^2
    #end
    if shape in [ TET4, TET10 ] 
        s = metric
        A = s/4
        a = 2.0*√(A/√3.0)
        return a^3/(6*√2.0)
    end
    if shape in [ HEX8, HEX20 ] 
        s = metric
        A = s/6
        a = √A
        return a^3
    end
    #if shape in [ WED6, WED15 ]
        #s = metric
        #return (3.0 + √3.0/2.0)*a2
    #end
    error("No regular area/volume value for shape $(get_name(shape))")
end


#= Returns the cell aspect ratio
function cell_aspect_ratio(c::Cell)
    # get faces
    faces = get_faces(c)
    if length(faces)==0 
        return 1.0
    end

    # cell surface
    fmetrics = [ cell_extent(f) for f in faces ]
    surf = sum(fmetrics)
    ar = minimum(fmetrics)/maximum(fmetrics)

    # quality calculation
    metric = cell_extent(c) # volume or area
    rsurf  = regular_surface(metric, c.shape)
    return rsurf/surf
end =#



# Returns the cell quality ratio as vol/reg_vol
function cell_quality_2(c::Cell)::Float64
    # get faces
    faces = get_faces(c)
    length(faces)==0 && return 1.0

    # cell surface
    surf = sum( cell_extent(f) for f in faces )

    # quality calculation
    vol = cell_extent(c) # volume or area
    rvol = regular_vol(surf, c.shape)
    c.quality = min(vol/rvol, 1.0) # << updates cell property!!!
    return c.quality
end



# Returns the cell quality ratio as reg_surf/surf
function cell_quality(c::Cell)::Float64
    # get faces
    faces = get_faces(c)
    if length(faces)==0 
        return 1.0
    end

    # cell surface
    surf = 0.0
    for f in faces
        surf += cell_extent(f)
    end

    # quality calculation
    metric = cell_extent(c) # volume or area
    rsurf  = regular_surface(metric, c.shape)

    k = 2
    q = (rsurf/surf)^k # << updates cell property!!!
    if q>1
        q = 2 - q
    end
    c.quality = q
    return q
end

function cell_aspect_ratio(c::Cell)::Float64
    faces = get_faces(c)
    len = [ cell_extent(f) for f in faces ]
    c.quality = minimum(len)/maximum(len)
    return c.quality 
end

export getproperty
function Base.getproperty(c::Cell, s::Symbol)
    if s == :coords
        return getcoords(c)
    elseif s == :faces
        return get_faces(c)
    elseif s == :edges
        return get_edges(c)
    elseif s == :extent
        return cell_extent(c)
    else
        return getfield(c, s)
    end
end

# Get an array with shares for all points
function get_patches(cells::Array{Cell,1})
    # get all points from cells if needed
    pointsd = Dict{UInt64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    
    points = collect(values(pointsd))
    np     = length(points)

    # backup points ids
    bk_pt_id = [ pt.id for pt in points ]
    for i=1:np
        points[i].id = i
    end

    # get incidence array
    patches  = [ Cell[] for i=1:np ]
    for cell in cells
        for pt in cell.points
            push!(patches[pt.id], cell)
        end
    end

    # restore points ids
    for i=1:np
        points[i].id = bk_pt_id[i]
    end

    return points, patches
end


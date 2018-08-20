# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

abstract type Block end

# Point
# =====

mutable struct Point
    x    ::Float64
    y    ::Float64
    z    ::Float64
    tag  ::TagType
    id   ::Int64
    extra::Int64 # TODO: check where is used
    function Point(x::Real, y::Real, z::Real=0.0, tag::TagType=0)
        const NDIG = 14
        # zero is added to avoid negative bit sign for zero signed values
        x = round(x, NDIG) + 0.0
        y = round(y, NDIG) + 0.0
        z = round(z, NDIG) + 0.0
        return new(x, y, z, tag,-1,-1)
    end
    function Point(C::AbstractArray{<:Real}; tag=0)
        # zero is added to avoid negative bit sign for zero signed values
        if length(C)==2
            return new(C[1]+0.0, C[2]+0.0, 0.0, tag, -1, -1)
        else
            return new(C[1]+0.0, C[2]+0.0, C[3]+0.0, tag, -1, -1)
        end
    end
end


### Point methods

# The functions below can be used in conjuntion with sort
get_x(point::Point) = point.x
get_y(point::Point) = point.y
get_z(point::Point) = point.z

#Base.hash(p::Point) = round(UInt, 1000000 + p.x*1001 + p.y*10000001 + p.z*100000000001)
#Base.hash(p::Point) = hash( (p.x==0.0?0.0:p.x, p.y==0.0?0.0:p.y, p.z==0.0?0.0:p.z) ) # comparisons used to avoid signed zero
Base.hash(p::Point) = hash( (p.x, p.y, p.z) ) 
Base.hash(points::Array{Point,1}) = sum(hash(p) for p in points)

Base.:(==)(p1::Point, p2::Point) = hash(p1)==hash(p2)

getcoords(p::Point) = [p.x, p.y, p.z]

# Sentinel instance for not found point
const NO_POINT = Point(1e+308, 1e+308, 1e+308)

function get_point(points::Dict{UInt64,Point}, C::AbstractArray{<:Real})
    hs = hash(Point(C))
    return get(points, hs, nothing)
end

function get_coords(points::Array{Point,1}, ndim::Int64=3)
    np = length(points)
    C  = Array{Float64}(np, ndim)
    for i=1:np
        C[i,1] = points[i].x
        C[i,2] = points[i].y
        ndim>2 && (C[i,3] = points[i].z)
    end

    return C
end



# Cell
# ====

type Cell
    shape  ::ShapeType
    points ::Array{Point,1}
    tag    ::TagType
    id     ::Integer
    ndim   ::Integer
    quality::Float64              # quality index: surf/(reg_surf) 
    embedded::Bool                # flag for embedded cells
    crossed::Bool                 # flag if cell crossed by linear inclusion
    ocell  ::Union{Cell,Void}     # owner cell if this cell is a face/edge
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

function get_coords(c::Cell, ndim=3)::Array{Float64,2}
    n = length(c.points)
    C = Array{Float64, 2}(n, ndim)
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
function update!(c::Cell)::Void
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

    result = Array{Cell}(0)
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


type Bins
    bins::Array{Array{Cell,1},3}
    bbox::Array{Float64,2}
    lbin::Float64
    function Bins(nx=0, ny=0, nz=0, bbox=nothing)
        this = new()
        this.bins = Array{Array{Cell,1}}(nx, nx, nx)
        this.bbox = zeros(0,0)
        # this.lbin = ..
        return this
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

    #=
    # Get all points
    #pointsd = Dict{UInt64, Point}()
    #for cell in cells
        #for point in cell.points
            #pointsd[hash(point)] = point
        #end
    #end
    #points = values(pointsd)

    # Get bounding box
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for point in points
        if point.x<minx; minx = point.x end
        if point.y<miny; miny = point.y end
        if point.z<minz; minz = point.z end
        if point.x>maxx; maxx = point.x end
        if point.y>maxy; maxy = point.y end
        if point.z>maxz; maxz = point.z end
    end
    return [ minx miny minz; maxx maxy maxz ]
end
=#

function get_bin_cells(bins::Bins, p::Point) # returns an array of cells: TODO: untested function
    minx, miny, minz = bins.bbox[1,:]
    lbin = bins.lbin
    ix = floor(Int, (p.x - minx)/lbin) + 1
    iy = floor(Int, (p.y - miny)/lbin) + 1
    iz = floor(Int, (p.z - minz)/lbin) + 1
    return bins.bins[ix, iy, iz]
end

function add_cell(bins::Bins, cell::Cell) # TODO: untested function
    minx, miny, minz = bins.bbox[1,:]
    lbin = bins.lbin
    bbox = bounding_box(cell)
    X, Y, Z  = bbox[:,1], bbox[:,2], bbox[:,3]
    verts    = [ [x y z] for z in Z, y in Y, x in X ]
    cell_pos = Set()
    for V in verts
        x, y, z = V
        ix = floor(Int, (x - minx)/lbin) + 1
        iy = floor(Int, (y - miny)/lbin) + 1
        iz = floor(Int, (z - minz)/lbin) + 1
        push!(cell_pos, (ix, iy, iz))
    end

    for (ix, iy, iz) in cell_pos
        push!(bins.bins[ix, iy, iz], cell)
    end
end

# Grepares a group of bins that contain given cells
function build_bins(cells::Array{Cell,1}, bins::Bins)
    # Get all points
    bins.bbox = bounding_box(cells)
    minx, miny, minz = bins.bbox[1,:]

    # Get max global lengths
    Lx, Ly, Lz = bins.bbox[2,:] - bins.bbox[1,:]
    max_L = max(Lx, Ly, Lz)

    # Get max cell lengths
    max_l = 0.0
    for cell in cells
        if cell.shape.family!=SOLID_SHAPE continue end
        bbox = bounding_box(cell)
        l    = maximum(bbox[2,:] - bbox[1,:])
        if l>max_l; max_l = l end
    end

    # Get number of divisions
    ndiv = min(50, 1*floor(Int, max_L/max_l)) # calibrate for bins efficiency

    lbin = max_L/ndiv     # Get bin length
    bins.lbin = lbin

    nx = floor(Int, Lx/lbin) + 1
    ny = floor(Int, Ly/lbin) + 1
    nz = floor(Int, Lz/lbin) + 1

    # Allocate bins
    bins.bins = Array{Array{Cell,1}}(nx, ny, nz)
    for k=1:nz, j=1:ny, i=1:nx
        bins.bins[i,j,k] = Array{Cell}(0)
    end

    # Fill bins
    for cell in cells
        bbox = bounding_box(cell)
        X, Y, Z  = bbox[:,1], bbox[:,2], bbox[:,3]
        verts    = [ (x, y, z) for z in Z, y in Y, x in X ]
        cell_loc = Set() # cells can be in more than one bin
        for (x, y, z) in verts
            ix = floor(Int, (x - minx)/lbin) + 1
            iy = floor(Int, (y - miny)/lbin) + 1
            iz = floor(Int, (z - minz)/lbin) + 1
            push!(cell_loc, (ix, iy, iz))
        end

        for (ix, iy, iz) in cell_loc
            push!(bins.bins[ix, iy, iz], cell)
        end
    end
end

# Find the cell that contains a given point
function find_cell(X::Array{Float64,1}, cells::Array{Cell,1}, bins::Bins, tol::Float64, exc_cells::Array{Cell,1})
    # Point coordinates
    x, y, z = vcat(X, 0)[1:3]
    lbin = bins.lbin

    # Build bins if empty
    length(bins.bins) == 0 && build_bins(cells, bins)

    for attempt=1:2 
        Cmin = reshape(bins.bbox[1,:],3)
        Cmax = reshape(bins.bbox[2,:],3)
        lbin = bins.lbin

        if any(X.<Cmin-tol) || any(X.>Cmax+tol)
            error("find_cell: point outside bounding box")
        end

        # Find bin index
        ix = floor(Int, (x - Cmin[1])/lbin) + 1
        iy = floor(Int, (y - Cmin[2])/lbin) + 1
        iz = floor(Int, (z - Cmin[3])/lbin) + 1

        # Search cell in bin
        bin = bins.bins[ix, iy, iz]
        for cell in bin
            coords = get_coords(cell)
            if is_inside(cell.shape, coords, X, tol) && !(cell in exc_cells)
                return cell
            end
        end

        # If not found in the first attempt try then rebuild bins
        if attempt==1
            warn("find_cell: Bin search failed. Rebuilding bins...")
            build_bins(cells, bins) 
        end
    end

    warn("find_cell: Bin search failed")
    return nothing
end


# gets all facets of a cell
function get_faces(cell::Cell)
    faces  = Array{Cell}(0)

    #all_faces_idxs = FACETS_IDXS[cell.shape]
    all_faces_idxs = cell.shape.facet_idxs
    #facet_shape = FACETS_SHAPE[cell.shape]         # facet_shape could be of type ShapeType or Tuple
    facet_shape = cell.shape.facet_shape
    
    if facet_shape==() return faces end

    sameshape  = typeof(facet_shape) == ShapeType # check if all facets have the same shape

    # Iteration for each facet
    for (i, face_idxs) in enumerate(all_faces_idxs)
        points = cell.points[face_idxs]
        shape  = sameshape? facet_shape : facet_shape[i]
        face   = Cell(shape, points, cell.tag, cell)
        push!(faces, face)
    end

    return faces
end

# gets all edges of a cell
function get_edges(cell::Cell)
    if cell.shape.ndim==2 return get_faces(cell) end

    edges  = Array{Cell}(0)
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
    C = get_coords(c)
    J = Array{Float64}(nldim, size(C,2))

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

# Returns the desirable surface/perimeter given the volume/area of a cell
function regular_surface(metric::Float64, shape::ShapeType)
    if shape in [ TRI3, TRI6, TRI9, TRI10 ] 
        A = metric
        a = 2.*√( A/√3.)
        return 3*a
    end
    if shape in [ QUAD4, QUAD8, QUAD9, QUAD12, QUAD16 ] 
        A = metric
        a = √A
        return 4*a
    end
    if shape in [ PYR5 ] 
        V = metric
        a = ( 3.*√2.*V )^(1./3.)
        return (1 + √3.)*a^2
    end
    if shape in [ TET4, TET10 ] 
        V = metric
        a = ( 6.*√2.*V )^(1./3.)
        return √3.*a^2
    end
    if shape in [ HEX8, HEX20 ] 
        V = metric
        a = V^(1./3.)
        return 6.*a^2.
    end
    if shape in [ WED6, WED15 ]
        V = metric
        a2 = (16./3.*V^2)^(1./3.) 
        return (3. + √3./2.)*a2
    end
    error("No regular surface/perimeter value for shape $(get_name(shape))")
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

    c.quality = rsurf/surf # << updates cell property!!!
    return min(c.quality, 1.0)
end

function cell_aspect_ratio(c::Cell)::Float64
    faces = get_faces(c)
    len = [ cell_extent(f) for f in faces ]
    c.quality = minimum(len)/maximum(len)
    return c.quality
end

#end #module

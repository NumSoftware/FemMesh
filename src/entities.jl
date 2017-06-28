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

import Base.getindex
export get_points

abstract type Block end

# Point
# =====

type Point
    x    ::Float64
    y    ::Float64
    z    ::Float64
    tag  ::AbstractString
    id   ::Int64
    extra::Int64
    function Point(x::Real, y::Real, z::Real=0.0; tag::String="")
        const NDIG = 14
        x = round(x, NDIG)
        y = round(y, NDIG)
        z = round(z, NDIG)
        return new(x, y, z, tag,-1,-1)
    end
    function Point(C::AbstractArray{<:Real}; tag::String="")
        if length(C)==2
            return new(C[1], C[2], 0.0, tag, -1, -1)
        else
            return new(C[1], C[2], C[3], tag, -1, -1)
        end
    end
end


### Point methods

import Base.hash
#hash(p::Point) = round(UInt, 1000000 + p.x*1001 + p.y*10000001 + p.z*100000000001)
hash(p::Point) = hash( (p.x==0.0?0.0:p.x, p.y==0.0?0.0:p.y, p.z==0.0?0.0:p.z) ) # comparisons used to avoid signed zero

function hash(points::Array{Point,1})::UInt64
    hs = 0x000000000::UInt64
    for p in points
        hs += hash(p)
    end
    return hs
end

getcoords(p::Point) = [p.x, p.y, p.z]

# Sentinel instance for not found point
const NO_POINT = Point(1e+308, 1e+308, 1e+308)

function get_point(points::Dict{UInt64,Point}, C::AbstractArray{<:Real})
    hs = hash(Point(C))
    return get(points, hs, nothing)
end

function get_coords(points::Array{Point,1}, ndim::Int64=3)::Array{Float64,2}
    np = length(points)
    C  = Array{Float64}(np, ndim)
    for i=1:np
        C[i,1] = points[i].x
        C[i,2] = points[i].y
        if ndim>2
            C[i,3] = points[i].z
        end
    end

    return C
end



# Cell
# ====

type Cell
    shape  ::ShapeType
    points ::Array{Point,1}
    tag    ::AbstractString
    id     ::Integer
    ndim   ::Integer
    quality::Float64              # quality index: surf/(reg_surf) 
    crossed::Bool                 # flag if cell crossed by linear inclusion
    ocell  ::Union{Cell,Void}     # owner cell in the case of a face
    linked_cells::Array{Cell,1}   # neighbor cells in case of joint cell
    function Cell(shape::ShapeType, points::Array{Point,1}, tag::AbstractString="", ocell=nothing)
        this = new(shape, points, tag, -1)
        this.ndim = 0
        this.quality = 0.0
        this.crossed = false
        this.ocell   = ocell
        this.linked_cells = []
        return this
    end
end

Face=Cell


### Cell methods


hash(c::Cell) = sum([ hash(p) for p in c.points])

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
function getindex(cells::Array{Cell,1}, s::Symbol)::Array{Cell,1}
    if s == :solids
        return filter(cell -> cell.shape.class==SOLID_SHAPE, cells)
    end
    if s == :lines
        return filter(cell -> cell.shape.class==LINE_SHAPE, cells)
    end
    if s == :joints
        return filter(cell -> cell.shape.class==JOINT_SHAPE, cells)
    end
    if s == :joints1D
        return filter(cell -> cell.shape.class==JOINT1D_SHAPE, cells)
    end
    if s == :points
        return get_points(cells)
    end
    error("Cell getindex: Invalid symbol $s")
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

# Gets the coordinates of a bounding box for an array of cells
function bounding_box(cells::Array{Cell,1})
    # Get all points
    pointsd = Dict{UInt64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    points = values(pointsd)

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

# Gets the coordinates of a bounding box for a cell
function bounding_box(cell::Cell)
    minx = miny = minz =  Inf
    maxx = maxy = maxz = -Inf
    for point in cell.points
        if point.x<minx; minx = point.x end
        if point.y<miny; miny = point.y end
        if point.z<minz; minz = point.z end
        if point.x>maxx; maxx = point.x end
        if point.y>maxy; maxy = point.y end
        if point.z>maxz; maxz = point.z end
    end
    return [ minx miny minz; maxx maxy maxz ]
end

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
        if cell.shape.class!=SOLID_SHAPE continue end
        bbox = bounding_box(cell)
        max  = maximum(bbox[2,:] - bbox[1,:])
        if max>max_l; max_l = max end
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
end

# Find the cell that contains a given point
function find_cell(X::Array{Float64,1}, cells::Array{Cell,1}, bins::Bins, Tol::Float64, exc_cells::Array{Cell,1})
    # Point coordinates
    x, y, z = vcat(X, 0)[1:3]
    lbin = bins.lbin

    # Build bins if empty
    if length(bins.bins) == 0
        build_bins(cells, bins)
    end

    for attempt=1:2 
        Cmin = reshape(bins.bbox[1,:],3)
        Cmax = reshape(bins.bbox[2,:],3)
        lbin = bins.lbin

        if any(X.<Cmin-Tol) || any(X.>Cmax+Tol)
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
            if is_inside(cell.shape, coords, X, Tol) && !(cell in exc_cells)
                return cell
            end
        end

        # If not found in the first try then rebuild bins
        if attempt==1
            warn("Bin search failed. Rebuilding bins...")
            build_bins(cells, bins) 
        end
    end

    warn("Bin search failed")
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

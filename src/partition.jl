# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

mutable struct SpacePartition{T}
    bins::Array{Array{T,1},3} # list of bins
    bbox::Array{Float64,2}    # bounding box
    lbin::Float64             # length of a bin

    function SpacePartition{T}() where T
        this = new{T}()
        this.bins = Array{Array{T,1}}(undef, 0, 0, 0)
        this.bbox = zeros(0,0)
        this.lbin = 0.0
        return this
    end

end


function SpacePartition{Point}(n::Int64, bbox::Array{Float64,2})
    # n: number of points hint

    par = SpacePartition{Point}()
    par.bbox = bbox

    Lx, Ly, Lz = bbox[2,:] - bbox[1,:]
    ndim = 0
    V    = 1.0
    for L in (Lx, Ly, Lz)
        L == 0 && continue
        V *= L
        ndim += 1
    end

    mpoints = 5
    lbin = (V/n*mpoints)^(1/ndim) # estimated bin length
    par.lbin = lbin

    nx = floor(Int, Lx/lbin) + 1
    ny = floor(Int, Ly/lbin) + 1
    nz = floor(Int, Lz/lbin) + 1

    par.bins = Array{Array{Point,1}}(undef, nx, ny, nz)

    return par
end


function Base.push!(par::SpacePartition{Point}, point::Point)
    @assert(par.lbin>0)

    lbin = par.lbin
    x, y, z, = point.x, point.y, point.z
    ix = floor(Int, (x - minx)/lbin) + 1
    iy = floor(Int, (y - miny)/lbin) + 1
    iz = floor(Int, (z - minz)/lbin) + 1
    push!(par.bins[ix, iy, iz], point)
end


function Base.get(par::SpacePartition{Point}, point::Point, default=nothing)
    @assert(par.lbin>0)

    lbin = par.lbin
    x, y, z, = point.x, point.y, point.z
    ix = floor(Int, (x - minx)/lbin) + 1
    iy = floor(Int, (y - miny)/lbin) + 1
    iz = floor(Int, (z - minz)/lbin) + 1

    for p in par.bins[ix, iy, iz]
        point == p && return p
    end

    return nothing
end


function build_partition!(points::Array{Point,1})
    # Get all points
    npoints = length(points)
    bbox = bounding_box(points)

    par = SpacePartition{Point}(npoints, bbox)

    # Fill bins
    for point in points
        push!(par, point)
    end

end



mutable struct Bins
    bins::Array{Array{Cell,1},3}
    bbox::Array{Float64,2}
    lbin::Float64
    function Bins(nx=0, ny=0, nz=0, bbox=nothing)
        this = new()
        this.bins = Array{Array{Cell,1}}(undef, nx, nx, nx)
        this.bbox = zeros(0,0)
        # this.lbin = ..
        return this
    end
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
    bins.bins = Array{Array{Cell,1}}(undef, nx, ny, nz)
    for k=1:nz, j=1:ny, i=1:nx
        bins.bins[i,j,k] = Cell[]
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

        if any(X .< Cmin .- tol) || any(X .> Cmax .+ tol)
            error("find_cell: point outside bounding box")
        end

        # Find bin index
        ix = floor(Int, (x - Cmin[1])/lbin) + 1
        iy = floor(Int, (y - Cmin[2])/lbin) + 1
        iz = floor(Int, (z - Cmin[3])/lbin) + 1

        # Search cell in bin
        bin = bins.bins[ix, iy, iz]
        for cell in bin
            coords = getcoords(cell)
            if is_inside(cell.shape, coords, X, tol) && !(cell in exc_cells)
                return cell
            end
        end

        # If not found in the first attempt try then rebuild bins
        if attempt==1
            @warn "find_cell: Bin search failed. Rebuilding bins..."
            build_bins(cells, bins) 
        end
    end

    @warn "find_cell: Bin search failed"
    return nothing
end

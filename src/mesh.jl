# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

#include("delaunay.jl")

### Type Mesh
"""
    Mesh()

constructs an unitialized mesh object to be used in finite element analyses.
It contains geometric fields as: points, cells, faces, edges, ndim, quality, etc.
"""
mutable struct Mesh
    ndim   ::Int
    points ::Array{Point,1}
    cells  ::Array{Cell,1}

    faces  ::Array{Cell,1}
    edges  ::Array{Cell,1}
    bpoints::Dict{UInt64,Point}
    bins   ::Bins
    quality::Float64
    qmin   ::Float64 

    # Data
    point_scalar_data::Dict{String,Array}
    cell_scalar_data ::Dict{String,Array}
    point_vector_data::Dict{String,Array}

    function Mesh()
        this = new()
        this.points = []
        this.bpoints = Dict{UInt64, Point}()
        this.cells  = []
        this.faces  = []
        this.edges  = []
        this.ndim   = 0
        this.bins   = Bins()
        this.quality = 0.0
        this.qmin    = 0.0

        this.point_scalar_data = Dict()
        this.cell_scalar_data  = Dict()
        this.point_vector_data = Dict()
        return this
    end
end

function get_surface(cells::Array{Cell,1})::Array{Cell,1}
    surf_dict = Dict{UInt64, Cell}()

    # Get only unique faces. If dup, original and dup are deleted
    for cell in cells
        for face in get_faces(cell)
            #conns = [ pt.id for pt in face.points ]
            #hs = hash(conns)
            hs = hash(face)
            if haskey(surf_dict, hs)
                delete!(surf_dict, hs)
            else
                surf_dict[hs] = face
            end
        end
    end

    return [ face for face in values(surf_dict) ]
end

function get_edges(surf_cells::Array{Cell,1})::Array{Cell,1}
    edges_dict = Dict{UInt64, Cell}()

    # Get only unique edges
    for cell in surf_cells
        for edge in get_faces(cell)
            edge.ocell = cell.ocell
            edges_dict[hash(edge)] = edge
        end
    end

    return [ edge for edge in values(edges_dict) ] 
end

# Return a list of neighbors for each cell
function get_neighbors(cells::Array{Cell,1})::Array{Cell,1}
    faces_dict = Dict{UInt64, Cell}()
    neighbors = [ Cell[] for i=1:length(cells) ]

    # Get cell faces. If dup, original and dup are deleted but neigh info saved
    for cell in cells
        for face in get_faces(cell)
            hs = hash(face)
            other = get(faces_dict, hs, nothing)
            if other == nothing
                faces_dict[hs] = face
            else
                push!(neighbors[face.ocell.id], other.ocell)
                push!(neighbors[other.ocell.id], face.ocell)
                delete!(faces_dict, hs)
            end
        end
    end

    return neighbors

end

# Return the cell patch for each point
function get_patches(mesh::Mesh)
    patches = [ Cell[] for i=1:length(mesh.points) ]
    for cell in mesh.cells
        for pt in cell.points
            push!(patches[pt.id], cell)
        end
    end

    return mesh.points, patches
end

# Reverse Cuthill–McKee algorithm (RCM) 
function reorder!(mesh::Mesh; sort_degrees=true, reversed=false)

    # Get all mesh edges
    all_edges = Dict{UInt64, Cell}()
    for cell in mesh.cells

        # adding cell edges
        if cell.shape.family == SOLID_SHAPE #is_solid(cell.shape)
            for edge in get_edges(cell)
                hs = hash(edge)
                all_edges[hs] = edge
            end

            #check for lagrangian elements
            if cell.shape==QUAD9
                edge = Cell(POLYV, [ cell.points[1], cell.points[end] ] )
                hs   = hash(edge)
                all_edges[hs] = edge
            end
            continue
        end

        # joint1D cells (semi-embedded approach)
        if cell.shape in (JLINK2, JLINK3)
            edge = Cell(POLYV, [ cell.points[1], cell.points[end-1] ])
            hs = hash(edge)
            all_edges[hs] = edge
            continue
        end

        # embedded line cells
        if cell.shape.family == LINE_SHAPE && length(cell.linked_cells)>0
            edge1 = Cell(cell.shape, cell.points)
            edge2 = Cell(LIN2, [ cell.points[1], cell.linked_cells[1].points[1] ])
            all_edges[hash(edge1)] = edge1
            all_edges[hash(edge2)] = edge2
            continue
        end

        # all other cells
        edge = Cell(cell.shape, cell.points)
        hs = hash(edge)
        all_edges[hs] = edge

    end

    # Get neighbors ids
    npoints = length(mesh.points)
    neighs_ids = Array{Int64}[ [] for i in 1:npoints ]

    for edge in values(all_edges)
        points = edge.points
        np = length(points)
        for i=1:np-1
            for j=i+1:np
                push!(neighs_ids[points[i].id], points[j].id)
                push!(neighs_ids[points[j].id], points[i].id)
            end
        end
    end

    # removing duplicates and get neighbors
    neighs = Array{Point}[ mesh.points[unique(list)] for list in neighs_ids ]
    #neighs = Array{Point}[ mesh.points[unique(list)] for list in neighs_ids ]

    # list of degrees per point
    degrees = Int64[ length(list) for list in neighs]
    mindeg, idx  = findmin(degrees)

    if mindeg == 0
        # Case of overlapping elements where edges have at least one point with the same coordinates
        @warn "reorder!: Reordering nodes failed! Check for overlapping cells or non used points."
        return
    end

    N = [ mesh.points[idx] ] # new list of points
    L = Dict{Int64,Point}()  # last levelset. Use ids as keys instead of hash to avoid collisions of points with same coordinates
    L[idx] = mesh.points[idx]
    LL = Dict{Int64,Point}()  # levelset before the last one

    while length(N) < npoints
        # Generating current levelset A
        A = Dict{Int64,Point}() 

        for p in values(L)
            for q in neighs[p.id]
                (haskey(L, q.id) || haskey(LL, q.id)) && continue
                A[q.id] = q
            end
        end
        if length(A)==0
            @error "reorder!: Reordering nodes failed! Possible error with cell connectivities."
            return
        end
        
        # Convert A into an array RA
        RA = collect(values(A))
        if sort_degrees
            D  = [ degrees[point.id] for point in RA ]
            RA = RA[sortperm(D)]
        end

        append!(N, RA)
        LL, L = L, A
    end

    # Reverse list of new nodes
    if reversed
        N = reverse(N)
    end

    # Renumerating
    for (i,p) in enumerate(N)
        p.id = i
    end

    mesh.points = N

    return nothing

end


# Updates numbering, faces and edges in a Mesh object
function update!(mesh::Mesh; verbose::Bool=false, genfacets::Bool=true, genedges::Bool=true)

    # Get ndim
    ndim = 1
    for point in mesh.points
        point.y != 0.0 && (ndim=2)
        point.z != 0.0 && (ndim=3; break)
    end
    mesh.ndim = ndim
    
    # Numberig nodes
    for (i,p) in enumerate(mesh.points) p.id = i end

    # Numberig cells 
    for (i,c) in enumerate(mesh.cells ) 
        c.id = i; 
        c.ndim=ndim; 
    end

    # Facets
    if genfacets
        verbose && print("  finding facets...   \r")
        mesh.faces = get_surface(mesh.cells)
    end
    ndim==2 && (mesh.edges=mesh.faces)

    # Edges
    if genedges && ndim==3
        verbose && print("  finding edges...   \r")
        mesh.edges = get_edges(mesh.faces)
    end

    # Quality
    Q = Float64[]
    for c in mesh.cells
        c.quality = cell_quality(c)
        push!(Q, c.quality)
    end

    # Reset data
    mesh.point_scalar_data = Dict()
    mesh.cell_scalar_data  = Dict()
    mesh.point_vector_data = Dict()
    mesh.cell_scalar_data["quality"] = Q

    return nothing
end

# Mesh quality
function quality!(mesh::Mesh)
    for c in mesh.cells
        c.quality = cell_quality(c)
    end
    Q = Float64[ c.quality for c in mesh.cells]
    mesh.quality = sum(Q)/length(mesh.cells)
    mesh.qmin    = minimum(Q)
    return nothing
end

# Adds m2 to m1
function join_mesh!(m1::Mesh, m2::Mesh)
    #m = Mesh()

    # copy m1 to m
    #m.points = copy(m1.points)
    #m.cells  = copy(m1.cells)
    #m.bpoints = copy(m1.bpoints)

    pid = length(m1.points)
    cid = length(m1.cells)

    #@show length(m1.points)
    #@show length(m2.points)

    #for (h,p) in m1.bpoints
        #@show (p.id, p.x, p.y, p.z)
    #end

    # Add new points from m2
    for p in m2.points
        hs = hash(p)
        if !haskey(m1.bpoints, hs)
            #@show (p.x, p.y, p.z)
            pid += 1
            p.id = pid
            push!(m1.points, p)
        end
    end

    #@show length(m1.points)

    # Fix m2 cells connectivities for points at m1 border
    for c in m2.cells
        for (i,p) in enumerate(c.points)
            hs = hash(p)
            if haskey(m1.bpoints, hs)
                pp = m1.bpoints[hs]
                c.points[i] = pp
            else
                # update bpoints dict
                if haskey(m2.bpoints, hs)
                    m1.bpoints[hs] = p
                end
            end
        end

        cid += 1
        c.id = cid
        push!(m1.cells, c)
    end

    return nothing
end


function Mesh(coords::Array{Float64}, conns::Array{Array{Int64,1},1}, cellshapes::Array{ShapeType,1}=ShapeType[])
    n = size(coords, 1) # number of points
    m = size(conns , 1) # number of cells

    points = Point[]
    for i=1:n
        C = coords[i,:]
        push!(points, Point(C))
    end

    # Get ndim
    ndim = size(coords,2)
    
    cells = Cell[]
    for i=1:m
        pts = points[conns[i]]
        if length(cellshapes)>0
            shape = cellshapes[i]
        else
            shape = get_shape_from_ndim_npoints(length(pts), ndim)
        end
        cell = Cell(shape, points)
        push!(cells, cell)
    end

    mesh = Mesh()
    mesh.points = points
    mesh.cells  = cells
    update!(mesh) # no node ordering

    return mesh
end


# flattens a list of nested lists
function flatten(x, y)
    ty = typeof(x)
    if ty <: Tuple || ty <: Array
        for item in x
            flatten(item, y)
        end
    else
        push!(y, x)
    end
    return y
end
flatten(x)=flatten(x, [])


"""
    Mesh(blocks, [verbose=true], [genfacets=true], [genedges=false], [reorder=true])

Generates a mesh based on an array of geometry blocks
"""
function Mesh(items::Union{Mesh, Block, Array{<:Union{Block, Array},1}}...; 
              verbose::Bool=true, genfacets::Bool=true, genedges::Bool=true, reorder=true)

    # Flatten items list
    fitems = flatten(items)

    # Get list of blocks and update mesh objects
    blocks = []
    meshes = []
    for item in fitems
        if isa(item, Block)
            push!(blocks, item)
        elseif isa(item, Mesh)
            push!(meshes, copy(item))
        else
            error("Mesh: item of type $(typeof(item)) cannot be used as Block or Mesh")
        end
    end

    nmeshes = length(meshes)
    nblocks = length(blocks)
    if verbose
        printstyled("Mesh generation:\n", bold=true, color=:cyan)
        nmeshes>0 && @printf "  %5d meshes\n" nmeshes
        @printf "  %5d blocks\n" nblocks
    end

    # New mesh object
    mesh = Mesh()

    # Join meshes
    for m in meshes
        join_mesh!(mesh, m)
    end

    # Split blocks: generates points and cells
    for (i,b) in enumerate(blocks)
        b.id = i
        split_block(b, mesh)
        verbose && print("  spliting block ", i, "...    \r")
    end

    # Updates numbering, quality, facets and edges
    update!(mesh, verbose=verbose, genfacets=genfacets, genedges=genedges)

    # Reorder nodal numbering
    if reorder
        verbose && print("  reordering points...             \r")
        reorder!(mesh)
    end

    if verbose
        npoints = length(mesh.points)
        ncells  = length(mesh.cells)
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        @printf "  %4dd mesh                             \n" mesh.ndim
        @printf "  %5d points\n" npoints
        @printf "  %5d cells\n" ncells
        if genfacets
            @printf "  %5d faces\n" nfaces
        end
        if genedges
            @printf "  %5d surface edges\n" nedges
        end
        icount = 0
        for block in blocks
            if typeof(block) == BlockInset; icount += block.icount end
        end
        if icount>0
            @printf "  %5d intersections\n" icount
        end
        println("  done.")
    end

    return mesh
end


function Mesh(cells::Array{Cell,1})
    # New mesh object
    mesh = Mesh()
    mesh.points = cells[:points]
    mesh.cells  = cells
    update!(mesh)
    reorder!(mesh)

    return mesh
end

function Base.convert(::Type{UnstructuredGrid}, mesh::Mesh)
    npoints = length(mesh.points)
    ncells  = length(mesh.cells)

    # points
    points  = zeros(npoints, 3)
    for (i,P) in enumerate(mesh.points)
        points[i, 1] = P.x
        points[i, 2] = P.y
        points[i, 3] = P.z
    end
    cells   = [ [point.id for point in C.points] for C in mesh.cells ] # TODO: make it not dependent of point.id
    cell_tys= [ Int64(C.shape.vtk_type) for C in mesh.cells ]
    
    # point id
    point_ids = [ P.id for P in mesh.points ]

    # cell data
    cell_scalar_data = copy(mesh.cell_scalar_data)
    # cell id
    cell_scalar_data["cell-id"] = [ C.id for C in mesh.cells ]
    # cell quality
    cell_scalar_data["quality"] = [ C.quality for C in mesh.cells ]
    # cell crossed
    cell_crs  = [ C.crossed for C in mesh.cells ]
    if any(cell_crs)
        cell_scalar_data["crossed"] = Int.(cell_crs)
    end

    point_scalar_data = copy(mesh.point_scalar_data)
    point_scalar_data["point-id"] = point_ids
    point_vector_data = copy(mesh.point_vector_data)

    # title
    title = "FemMesh output"

    return UnstructuredGrid(title, points, cells, cell_tys,
                                 point_scalar_data = point_scalar_data,
                                 cell_scalar_data  = cell_scalar_data,
                                 point_vector_data = point_vector_data)
end



function Mesh(ugrid::UnstructuredGrid)

    mesh = Mesh()

    # Setting points
    npoints = size(ugrid.points,1)
    ndim = 2
    for i=1:npoints
        X = ugrid.points[i,:]
        point = Point(X[1], X[2], X[3])
        point.id = i
        push!(mesh.points, point)
        if X[3] != 0.0
            ndim = 3
        end
    end

    # Set ndim
    mesh.ndim = ndim

    # Setting cells
    ncells = length(ugrid.cells)
    has_polyvertex = false

    for i=1:ncells
        conn = mesh.points[ ugrid.cells[i] ]
        vtk_shape = VTKCellType(ugrid.cell_types[i])
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYV
            has_polyvertex = true
        else
            shape = get_shape_from_vtk( vtk_shape, length(conn), ndim )
        end
        cell  = Cell(shape, conn)
        push!(mesh.cells, cell)
    end

    # update mesh and get faces and edges
    update!(mesh)

    # Setting data
    mesh.point_scalar_data = ugrid.point_scalar_data
    mesh.point_vector_data = ugrid.point_vector_data
    mesh.cell_scalar_data  = ugrid.cell_scalar_data

    # Fix shape for polyvertex cells
    if has_polyvertex
        # mount dictionary of cells
        cdict = Dict{UInt64, Cell}() 
        for cell in mesh.cells
            hs = hash(cell)
            cdict[hs] = cell
        end

        # check cells
        for cell in mesh.cells
            if cell.shape == POLYV 
                n = length(cell.points)
                # look for joints1D and fix shape
                if n>=5
                    # check if cell is related to a JLINK2
                    hss = hash( [ cell.points[i] for i=1:n-2] )
                    if haskey(cdict, hss) 
                        cell.shape = JLINK2
                        hs0   = hash( [ cell.points[i] for i=n-1:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_cells = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end

                    # check if cell is related to a JLINK3
                    hss = hash( [ cell.points[i] for i=1:n-3] )
                    if haskey(cdict, hss) 
                        cell.shape = JLINK3
                        hs0   = hash( [ cell.points[i] for i=n-2:n] )
                        hcell = cdict[hss]
                        lcell = cdict[hs0]
                        cell.linked_cells = [hcell, lcell]
                        hcell.crossed = true
                        continue
                    end
                end

                # look for conventional joints and fix shape
                if n%2==0 && n>=4
                    isjoint = true
                    for i=1:div(n,2)
                        p1 = cell.points[i]
                        p2 = cell.points[div(n,2)+i]
                        if hash(p1) != hash(p2) 
                            isjoint = false
                            break
                        end
                    end
                    if isjoint
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim, 2)
                        continue
                    end
                end

                # look for joint elements with 3 layers and fix shape
                if n%3==0 && n>=6
                    stride= div(n,3)
                    delta = 0.0
                    for i=1:stride
                        p1 = cell.points[i]
                        p2 = cell.points[stride+i]
                        if hash(p1) != hash(p2) 
                            isjoint = false
                            break
                        end
                    end
                    if isjoint
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim, 3)
                        continue
                    end
                end

                # remaining polyvertex cells
                cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)
            end
        end

        # Update linked cells in joints

        # check if there are joints
        has_joints = any( C -> C.shape.family==JOINT_SHAPE, mesh.cells )

        if has_joints
            # generate dict of faces
            facedict = Dict{UInt64, Cell}()
            for cell in mesh.cells
                for face in get_faces(cell)
                    hs = hash(face)
                    f  = get(facedict, hs, nothing)
                    facedict[hs] = face
                end
            end

            for cell in mesh.cells
                if cell.shape.family == JOINT_SHAPE
                    n = length(cell.points)
                    hs1 = hash( [ cell.points[i] for i=1:div(n,2)] )
                    hs2 = hash( [ cell.points[i] for i=div(n,2)+1:n] )
                    cell1 = facedict[hs1].ocell
                    cell2 = facedict[hs2].ocell
                    cell.linked_cells = [ cell1, cell2 ]
                end
            end
        end
    end

    # Fix linked_cells for embedded elements (line cells not connected to other elements)
    has_line = any(C->C.shape.family==LINE_SHAPE, mesh.cells)
    if has_line
        # TODO: find the owner of orphan line cells OR use parent id information
    end

    # Update boundary points: bpoints
    for f in mesh.faces
        for p in f.points
            mesh.bpoints[hash(p)] = p
        end
    end

    return mesh
end




"""
    save(mesh, filename, [verbose=true])

Saves a mesh object into a file in VTK legacy format
"""
function save(mesh::Mesh, filename::String; verbose::Bool=true)
    ugrid = convert(UnstructuredGrid, mesh)

    # Warning for lost of embedded flag
    has_embedded = any( C -> C.embedded, mesh.cells )
    if has_embedded
        # the function read_mesh_vtk cannot setup jet embedded cells
        @warn "save: the flag for embedded cells will be lost after saving mesh to file $filename"
    end

    save_vtk(ugrid, filename)
    verbose && printstyled( "  file $filename written (Mesh)\n", color=:cyan)
end


function load_mesh_json(filename; format="json")
    mesh = Mesh()

    data = JSON.parsefile(filename)
    points_key = haskey(data, "points") ? "points" : "verts"
    shape_key  = haskey(data["cells"][1], "shape") ? "shape" : "geo"
    verts = data[points_key]
    cells = data["cells"]

    # Loading points
    for vert in verts
        C = vert["c"]
        push!(C, 0.0)
        P = Point(C[1], C[2], C[3])
        if points_key=="points"
            P.id  = vert["id"]
        else
            P.id  = vert["id"]+1
        end
        P.tag = string(vert["tag"])
        push!(mesh.points, P)

    end

    # Get ndim
    ndim = 2
    for point in mesh.points
        if point.z != 0.0; ndim = 3; break end
    end
    mesh.ndim = ndim

    # Loading cells
    for cell in cells
        conn  = cell[points_key]
        id =  cell["id"]

        if shape_key=="shape"
            shapetype = eval(parse(cell[shape_key]))
        else
            conn  = [ i+1 for i in conn]
            id    = id+1
            shapetype = get_shape_from_geo(cell[shape_key])
        end
        points = Point[ mesh.points[i] for i in conn ]

        C = Cell(shapetype, points, string(cell["tag"]))

        # check for embedded joint
        if haskey(cell, "jlinid") # is joint element
            lincell = mesh.cells[cell["jlinid"]+1] # TODO
            sldcell = mesh.cells[cell["jsldid"]+1]
            C.linked_cells = [sldcell, lincell]
            if length(lincell.points)==2
                C.shape = JLINK2
            else
                C.shape = JLINK3
            end
        end

        C.id  = id
        C.ndim=ndim; 
        push!(mesh.cells, C)
    end

    # Mesh quality
    for c in mesh.cells
        update!(c)
    end

    Q = [ c.quality for c in mesh.cells]
    mesh.quality = sum(Q)/length(mesh.cells)
    mesh.qmin    = minimum(Q)

    # Generating faces
    mesh.faces = get_surface(mesh.cells)
    surf_dict = Dict(hash(F) => F for F in mesh.faces)

    all_faces = Face[]
    for (i,C) in enumerate(mesh.cells)
        if C.shape.family == SOLID_SHAPE #is_solid(C.shape)
            faces = get_faces(C)
            ftags = get(cells[i],"ftags", [])
            if length(ftags) > 0
                for (j,F) in enumerate(faces)
                    tag = string(ftags[j])
                    if ftags != "0"
                        F.tag = tag
                        surf_dict[hash(F)] = F
                    end
                end
            end
        end
    end

    mesh.faces = Cell[ face for face in values(surf_dict) ]

    return mesh
end


"""
`Mesh(filename)` constructs a mesh object based on a file in VTK legacy format or JSON format.
"""
function Mesh(filename::String; verbose::Bool=true, reorder::Bool=false)
    if verbose
        printstyled("Mesh loading:\n", bold=true, color=:cyan)
    end

    basename, ext = splitext(filename)
    if ext==".vtk"
        verbose && print("  Reading VTK legacy format...\n")
        mesh = Mesh(read_ugrid_vtk(filename))
    elseif ext==".json"
        verbose && print("  Reading JSON format...\n")
        mesh = load_mesh_json(filename)
    else
        verbose && print("  Reading tetgen output files...\n")
        mesh = Mesh(read_ugrid_tetgen(filename))
    end

    # Reorder nodal numbering
    if reorder
        verbose && print("  reordering points...\r")
        reorder!(mesh)
    end

    if verbose
        npoints = length(mesh.points)
        ncells  = length(mesh.cells)
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        println("  ", mesh.ndim, "d found             ")
        @printf "  %5d points read\n" npoints
        @printf "  %5d cells read\n" ncells
        @printf "  %5d faces obtained\n" nfaces
        @printf "  %5d surface edges obtained\n" nedges
    end

    return mesh
end

# Read a mesh from TetGen output files
# TODO: deprecate this function. Consider read_ugrid_tetgen function
export tetreader
function tetreader(filekey::String)
    # reading .node file
    f = open(filekey*".node")
    data = readlines(f)
    npoints, ndim, att, hasmarker = parse.(split(data[1]))
    points = Array{Point,1}(npoints)

    for i=1:npoints
        pdata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? pdata[end] : ""
        point = Point(pdata[2], pdata[3], pdata[4], tag)
        point.id = i
        points[i] = point
    end

    # reading .ele file
    f = open(filekey*".ele")
    data = readlines(f)
    ncells, ntetpts, hasatt = parse.(split(data[1]))
    celltype = ntetpts==4 ? TET4 : TET10
    cells = Array{Cell,1}(ncells)

    for i=1:ncells
        cdata = parse.(split(data[i+1]))
        tag   = hasatt==1 ? cdata[end] : 0
        cellpoints = points[cdata[2:ntetpts+1]]
        cell = Cell(celltype, cellpoints, tag)
        cell.id = i
        cells[i] = cell
    end

    # reading .face file
    f = open(filekey*".face")
    data = readlines(f)
    nfaces, hasmarker = parse.(split(data[1]))
    nfacepts = ntetpts==4 ? 3 : 6
    celltype = ntetpts==4 ? TRI3 : TRI6
    faces = Array{Cell,1}(nfaces)

    for i=1:nfaces
        fdata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? fdata[end] : 0
        facepts = points[fdata[2:nfacepts+1]]
        nfacepts==6 && (facepts = facepts[[1,2,3,6,4,5]])
        faces[i] = Cell(celltype, facepts, tag)
    end

    # reading .edge file
    f    = open(filekey*".edge")
    data = readlines(f)
    nedges, hasmarker = parse.(split(data[1]))
    nedgepts = ntetpts==4 ? 2 : 3
    celltype = ntetpts==4 ? LIN2 : LIN3
    edges = Array{Cell,1}(nedges)

    for i=1:nedges
        edata = parse.(split(data[i+1]))
        tag   = hasmarker==1 ? edata[end] : 0
        edgepts = points[edata[2:nedgepts+1]]
        edges[i] = Cell(celltype, edgepts, tag)
    end

    mesh = Mesh()
    mesh.points = points
    mesh.cells  = cells
    mesh.faces  = faces
    mesh.edges  = edges

    return mesh
end

# Gets a part of a mesh filtering elements
function Base.getindex(mesh::Mesh, filter_ex::Expr)
    # filter cells
    cells  = mesh.cells[filter_ex]
    # get points
    points = get_points(cells)

    # ids from selected cells and poitns
    cids = [ c.id for c in cells ]
    pids = [ p.id for p in points ]

    # new mesh object
    new_mesh = Mesh()
    new_mesh.points = points
    new_mesh.cells = cells

    # select relevant data
    for (key,vals) in mesh.point_scalar_data
        new_mesh.point_scalar_data[key] = vals[pids]
    end

    for (key,vals) in mesh.cell_scalar_data
        new_mesh.cell_scalar_data[key] = vals[cids]
    end

    for (key,vals) in mesh.point_vector_data
        new_mesh.point_vector_data[key] = vals[pids,:]
    end

    # update node numbering, facets and edges
    update!(new_mesh)

    return new_mesh

end


function threshold!(mesh::Mesh, field::Union{Symbol,String}, minval::Float64, maxval::Float64)
    field = string(field)

    # check if field exists
    found = haskey(mesh.cell_scalar_data, field)
    if found
        vals = mesh.cell_scalar_data[field]
    else
        found = haskey(mesh.point_scalar_data, field)
        found || error("threshold: field $field not found")
        data  = mesh.point_scalar_data[field]
        vals = [ mean(data[i]) for i=1:ncells ]
    end

    # filter cells
    cells = Cell[]
    for (cell, val) in zip(mesh.cells, vals)
        if minval <= val <= maxval
            push!(cells, cell)
        end
    end

    # get points
    points = get_points(cells)

    # ids from selected cells and points
    cids = [ c.id for c in cells ]
    pids = [ p.id for p in points ]

    # new mesh object
    new_mesh = Mesh()
    new_mesh.points = points
    new_mesh.cells = cells

    # select relevant data
    for (key,vals) in mesh.point_scalar_data
        new_mesh.point_scalar_data[key] = vals[pids]
    end

    for (key,vals) in mesh.cell_scalar_data
        new_mesh.cell_scalar_data[key] = vals[cids]
    end

    for (key,vals) in mesh.point_vector_data
        new_mesh.point_vector_data[key] = vals[pids,:]
    end

    # update node numbering, facets and edges
    update!(new_mesh)

    return new_mesh

end

export randmesh
function randmesh(l::Real...)
    ndim = length(l)
    if ndim==2
        lx, ly = l
        nx, ny = rand(4:7, 2)
        cellshape = rand((TRI3, TRI6, QUAD4, QUAD8))
        m = Mesh(Block2D([0.0 0.0; lx ly], nx=nx, ny=ny, cellshape=cellshape), verbose=false)
    else
        lx, ly, lz = l
        nx, ny, nz = rand(4:7, 3)
        cellshape = rand((TET4, TET10, HEX8, HEX20))
        m = Mesh(Block3D([0.0 0.0 0.0; lx ly lz], nx=nx, ny=ny, nz=nz, cellshape=cellshape), verbose=false)
    end
end

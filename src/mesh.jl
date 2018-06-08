# This file is part of FemMesh package. See copyright license in https://github.com/NumSoftware/FemMesh

#include("delaunay.jl")

### Type Mesh
"""
    Mesh()

constructs an unitialized mesh object to be used in finite element analyses.
It contains geometric fields as: points, cells, faces, edges, ndim, quality, etc.
"""
type Mesh
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
    #point_scalar_data::Dict{String,Array}
    #point_vector_data::Dict{String,Array}
    #cell_scalar_data ::Dict{String,Array}

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
function get_patches(mesh::Mesh)::Array{Array{Cell,1},1}
    patches = [ Cell[] for i=1:length(mesh.points) ]
    for cell in mesh.cells
        for pt in cell.points
            push!(patches[pt.id], cell)
        end
    end

    return patches
end

# Reverse Cuthill–McKee algorithm (RCM) 
function reorder!(mesh::Mesh; sort_degrees=true, reversed=false)::Void

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
            edge = Cell(LIN2, [ cell.points[1], cell.points[end-1] ])
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

    # Get points neighbors
    npoints = length(mesh.points)
    neighs  = Array{Point}[ [] for i in 1:npoints  ]

    for edge in values(all_edges)
        points = edge.points
        np = length(points)
        for i=1:np-1
            for j=i+1:np
                push!(neighs[points[i].id], points[j])
                push!(neighs[points[j].id], points[i])
            end
        end
    end

    # removing duplicates
    neighs  = Array{Point}[ unique(list) for list in neighs  ]

    # list of degrees per point
    degrees = Int64[ length(list) for list in neighs]
    mindeg, idx  = findmin(degrees)

    if mindeg == 0
        # Case of overlapping elements where edges have at least one point with the same coordinates
        warn("reorder!: Reordering nodes failed! Check for overlapping cells.")
        return
    end

    N = [ mesh.points[idx] ] # new list of cells
    R = [ mesh.points[idx] ] # first levelset

    R0 = [ ] # to be the last levelset

    while length(N) < npoints
        # Generating current levelset A
        A = Set{Point}([ ])

        for p in R
            for q in neighs[p.id]
                if q in R0 || q in R continue end
                push!(A, q)
            end
        end
        
        # Convert set A into an array RA
        RA = collect(A)
        if sort_degrees
            D  = [ degrees[point.id] for point in RA ]
            RA = RA[sortperm(D)]
        end

        append!(N, RA)
        R0 = R
        R  = A
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
function update!(mesh::Mesh; verbose::Bool=false, genfacets::Bool=true, genedges::Bool=true)::Void

    # Get ndim
    ndim = 2
    for point in mesh.points
        if point.z != 0.0; ndim = 3; break end
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
    #if 0<length(mesh.cells)<=1000
        quality!(mesh)
    #end

    return nothing
end

# Mesh quality
function quality!(mesh::Mesh)::Void
    Q = Float64[ cell_quality(c) for c in mesh.cells]
    mesh.quality = sum(Q)/length(mesh.cells)
    mesh.qmin    = minimum(Q)
    return nothing
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
function Mesh(items::Union{Block, Mesh, Array}...; verbose::Bool=true, genfacets::Bool=true, genedges::Bool=true, reorder=true)

    # New mesh object
    mesh = Mesh()

    # Flatten items list
    fitems = flatten(items)

    # Get list of blocks and update mesh objects
    blocks = []
    for item in fitems
        if isa(item, Block)
            push!(blocks, item)
        elseif isa(item, Mesh)
            append!(mesh.points, item.points)
            append!(mesh.cells, item.cells)
            mesh.bpoints = merge(mesh.bpoints, item.bpoints)
        else
            error("Mesh: item of type $(typeof(item)) cannot be used as Block or Mesh")
        end
    end

    nblocks = length(blocks)
    if verbose
        print_with_color(:cyan, "Mesh generation:\n", bold=true)
        #println("  analyzing $nblocks block(s)") 
    end

    # Split blocks: generates points and cells
    @printf "  %5d blocks\n" nblocks
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

import Base.convert
function convert(::Type{VTK_unstructured_grid}, mesh::Mesh)
    npoints = length(mesh.points)
    ncells  = length(mesh.cells)

    # points
    points  = zeros(npoints, 3)
    for (i,P) in enumerate(mesh.points)
        points[i, 1] = P.x
        points[i, 2] = P.y
        points[i, 3] = P.z
    end
    cells   = [ [point.id for point in C.points] for C in mesh.cells ]
    cell_tys= [ C.shape.vtk_type for C in mesh.cells ]
    
    # point id
    point_ids = [ P.id for P in mesh.points ]

    # cell data
    cell_scalar_data = Dict()
    # cell id
    cell_scalar_data["cell-id"] = [ C.id for C in mesh.cells ]
    # cell quality
    cell_scalar_data["quality"] = [ C.quality for C in mesh.cells ]
    # cell crossed
    cell_crs  = [ C.crossed for C in mesh.cells ]
    if any(cell_crs)
        cell_scalar_data["crossed"] = Int.(cell_crs)
    end

    # title
    title = "FemMesh output"

    return VTK_unstructured_grid(title, points, cells, cell_tys,
                                 point_scalar_data = Dict("point-id"=>point_ids),
                                 cell_scalar_data  = cell_scalar_data)

end



"""
    save(mesh, filename, [verbose=true])

Saves a mesh object into a file in VTK legacy format
"""
function save(mesh::Mesh, filename::String; verbose::Bool=true)
    vtk_data = convert(VTK_unstructured_grid, mesh)

    # Warning for lost of embedded flag
    has_embedded = any( C -> C.embedded, mesh.cells )
    if has_embedded
        # the function read_mesh_vtk cannot setup jet embedded cells
        warn("save: the flag for embedded cells will be lost after saving mesh to file $filename")
    end

    save_vtk(vtk_data, filename)
    verbose && print_with_color(:green, "  file $filename written (Mesh)\n")
end


function read_mesh_vtk(filename::String, verbose::Bool=false)

    vtk_data = read_VTK_unstructured_grid(filename)
    mesh = Mesh()

    # Setting points
    npoints = size(vtk_data.points,1)
    ndim = 2
    for i=1:npoints
        X = vtk_data.points[i,:]
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
    ncells = length(vtk_data.cells)
    has_polyvertex = false

    for i=1:ncells
        conn = mesh.points[ vtk_data.cells[i] ]
        vtk_shape = VTKCellType(vtk_data.cell_types[i])
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYV
            has_polyvertex = true
        else
            shape = get_shape_from_vtk( vtk_shape, length(conn), ndim )
        end
        cell  = Cell(shape, conn)
        push!(mesh.cells, cell)
    end

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
                    delta = 0.0
                    for i=1:div(n,2)
                        p1 = cell.points[i]
                        p2 = cell.points[div(n,2)+i]
                        delta += abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z)
                    end
                    if delta<1e-10
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)
                        continue
                    end
                end

                # look for joint elements (with 2 or 3 layers)
                # TODO
                #=
                if n%3==0 && n>=6
                    string = div(n,3)
                    delta = 0.0
                    for i=1:stride
                        p1 = cell.points[i]
                        p2 = cell.points[stride+i]
                        delta += abs(p1.x - p2.x) + abs(p1.y - p2.y) + abs(p1.z - p2.z)
                    end
                    if delta<1e-10
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, 2*stride, ndim, pointlayers=3)
                        continue
                    end
                end
                =#

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

    # update mesh and get faces and edges
    update!(mesh)

    # Update boundary points: bpoints
    for f in mesh.faces
        for p in f.points
            mesh.bpoints[hash(p)] = p
        end
    end

    return mesh
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
            #@show shapetype.name
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
            #@show " ", shapetype.name
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
function Mesh(filename::String; verbose::Bool=true, reorder::Bool=false)::Mesh
    if verbose
        print_with_color(:cyan, "Mesh loading:\n", bold=true)
    end

    basename, ext = splitext(filename)
    if ext==".vtk"
        verbose && print("  Reading VTK legacy format...\n")
        mesh = read_mesh_vtk(filename)
    else # try JSON
        verbose && print("  Reading JSON format...\n")
        mesh = load_mesh_json(filename)
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


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

include("delaunay.jl")

### Type Mesh
"""
`Mesh()` constructs an unitialized mesh object to be used in finite element analyses.
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
    point_scalars::Dict{Symbol,Array{Float64}}
    point_vectors::Dict{Symbol,Array{Array{Float64}}}
    cell_scalars ::Dict{Symbol,Array{Float64}}

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
        if cell.shape.class == SOLID_SHAPE #is_solid(cell.shape)
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

        # joint1D cells
        if cell.shape in (JLINK2, JLINK3)
            edge = Cell(LIN2, [ cell.points[1], cell.points[end-1] ])
            hs = hash(edge)
            all_edges[hs] = edge
            continue
        end

        # all other cells
        edge = Cell(cell.shape, cell.points)
        hs = hash(edge)
        all_edges[hs] = edge

    end

    #myedges = [ [ed.points[1].id, ed.points[2].id] for ed in values(all_edges)  ]
    #@show myedges

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
        # Case of overlapping elements where edges have at least a point with same coordinates
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
        verbose && print("  finding facets...\r")
        mesh.faces = get_surface(mesh.cells)
    end
    ndim==2 && (mesh.edges=mesh.faces)

    # Edges
    if genedges && ndim==3
        verbose && print("  finding edges...\r")
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
Generates a mesh based on an array of geometry blocks:
```
Mesh(blocks, [verbose=true,] [genfacets=true,] [genedges=false,] [initial_mesh=nothing]) -> mesh_object
```
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
        end
    end

    nblocks = length(blocks)
    if verbose
        print_with_color(:cyan, "Mesh generation:\n", bold=true)
        println("  analyzing $nblocks block(s)") 
    end

    # Split blocks: generates points and cells
    for (i,b) in enumerate(blocks)
        b.id = i
        split_block(b, mesh)
        verbose && print("  spliting block ",i,"...\r")
    end

    # Updates numbering, quality, facets and edges
    update!(mesh, verbose=verbose, genfacets=genfacets, genedges=genedges)

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
        @printf "  %5d points obtained\n" npoints
        @printf "  %5d cells obtained\n" ncells
        if genfacets
            @printf "  %5d faces obtained\n" nfaces
        end
        if genedges
            @printf "  %5d surface edges obtained\n" nedges
        end
        icount = 0
        for block in blocks
            if typeof(block) == BlockInset; icount += block.icount end
        end
        if icount>0
            @printf "  %5d intersections found\n" icount
        end
    end

    verbose && println("  done.")

    return mesh
end

#function generate_mesh(items::Union{Block, Array}...; verbose::Bool=true, genfacets::Bool=true, genedges::Bool=false)
    #println("generate_mesh function deprecated. Use Mesh instead")
    #return Mesh(items..., verbose=verbose, genfacets=genfacets, genedges=genedges)
#end


"""
Saves a mesh object into a file in VTK legacy format:
```
save(mesh, filename, [verbose=true])
```
"""
function save(mesh::Mesh, filename::AbstractString; verbose::Bool=true)::Void
    # Saves the mesh information in vtk format

    npoints = length(mesh.points)
    ncells  = length(mesh.cells)

    ndata = 0
    for cell in mesh.cells
        ndata += 1 + length(cell.points)
    end

    has_crossed = any(Bool[cell.crossed for cell in mesh.cells])

    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, "pyfem output ")
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", npoints, " float64")

    # Write points
    for point in mesh.points
        @printf f "%23.15e %23.15e %23.15e \n" point.x point.y point.z
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ",ncells, " ", ndata)
    for cell in mesh.cells
        print(f, length(cell.points), " ")
        for point in cell.points
            print(f, point.id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write cell types
    println(f, "CELL_TYPES ", ncells)
    for cell in mesh.cells
        println(f, cell.shape.vtk_type)
    end
    println(f)

    # Write point data
    println(f, "POINT_DATA ", npoints)

    # Write point numbers
    println(f, "SCALARS ", "Point-ID", " int 1")
    println(f, "LOOKUP_TABLE default")
    for point in mesh.points
        @printf f "%5d " point.id
    end
    println(f, )

    println(f, "CELL_DATA ", ncells)

    # Write cell numbers
    println(f, "SCALARS ", "Cell-ID", " int 1")
    println(f, "LOOKUP_TABLE default")
    for cell in mesh.cells  #naelems
        @printf f "%5d " cell.id
    end
    println(f, )

    # Write cells quality
    println(f, "SCALARS quality float 1")
    println(f, "LOOKUP_TABLE default")
    for cell in mesh.cells
        println(f, cell.quality)
    end
    println(f, )

    # Write cell type as cell data
    #println(f, "SCALARS cell_type int 1")
    #println(f, "LOOKUP_TABLE default")
    #for cell in mesh.cells
        #println(f, Int64(cell.shape))
    #end
    #println(f, )

    # Write flag for crossed cells
    if has_crossed
        println(f, "SCALARS crossed int 1")
        println(f, "LOOKUP_TABLE default")
        for cell in mesh.cells
            println(f, Int(cell.crossed))
        end
        println(f, )
    end

    if verbose
        print_with_color(:green, "  file $filename written (Mesh)\n")
    end

    close(f)

    return nothing

end

function load_mesh_vtk(filename::String)::Mesh

    file = open(filename)
    mesh = Mesh()

    # read nodal information
    alltext = readstring(filename)
    data    = split(alltext)

    # read header
    idx  = 1
    while data[idx] != "DATASET"
        idx += 1
    end
    idx += 1

    gridtype = data[idx]
    gridtype == "UNSTRUCTURED_GRID" || error("This reader only support VTK UNSTRUCTURED_GRID.")

    # read number of points
    while data[idx] != "POINTS"
        idx += 1
    end
    npoints = parse(data[idx+1])
    idx += 3

    # read points
    for i=1:npoints
        x = parse(Float64, data[idx])
        y = parse(Float64, data[idx+1])
        z = parse(Float64, data[idx+2])
        point = Point(x,y,z)
        point.id = i
        push!(mesh.points, point)
        idx += 3
    end

    # get ndim
    ndim = 2
    for point in mesh.points
        if point.z != 0.0; ndim = 3; break end
    end
    mesh.ndim = ndim

    # read number of cells
    while data[idx] != "CELLS"
        idx += 1
    end
    ncells = parse(data[idx+1])
    ncdata = parse(data[idx+2])
    idx += 3

    # read cells connectivities
    conns = Array{Point}[]
    for i=1:ncells
        npoints = parse(data[idx])
        idx += 1
        conn = Point[]
        for j=1:npoints
            point = mesh.points[ parse(data[idx]) + 1 ]
            idx  += 1
            push!(conn, point)
        end
        push!(conns, conn)
    end

    # read type of cells
    types = Int[]
    while data[idx] != "CELL_TYPES"
        idx += 1
    end
    idx += 2

    # read type of cells (cont.)
    has_polyvertex = false
    for i=1:ncells
        vtk_shape = VTKCellType(parse(data[idx]))
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYV
        else
            shape = get_shape_from_vtk( vtk_shape, length(conns[i]), ndim )
        end
        cell  = Cell(shape, conns[i])
        push!(mesh.cells, cell)
        idx  += 1
        if shape == POLYV   has_polyvertex = true end
    end

    # end of reading
    close(file)

    # Fixing shape field for polyvertex cells
    if has_polyvertex
        # mount dictionary of cells
        cdict = Dict{UInt64, Cell}() 
        for cell in mesh.cells
            #conn = [ point.id for point in cell.points]
            #hs = hash(conn)
            hs = hash(cell)
            cdict[hs] = cell
        end

        # check cells
        for cell in mesh.cells
            if cell.shape == POLYV 
                n = length(cell.points)
                # look for joints1D and fix shape
                if n>=5
                    #conn = [ point.id for point in cell.points]
                    #hss = hash(conn[1:end-2])
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
                    #hss = hash(conn[1:end-3])
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

                # look for joints and fix shape
                if n%2==0 && n>=4
                    delta = 0.0
                    for i=1:div(n,2)
                        p1 = cell.points[i]
                        p2 = cell.points[div(n,2)+i]
                        delta += abs(p1.x - p2.x)
                        delta += abs(p1.y - p2.y)
                        delta += abs(p1.z - p2.z)
                    end
                    if delta<1e-10
                        cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)
                        continue
                    end
                end

                # remaining polyvertex cells
                cell.shape = get_shape_from_vtk(VTK_POLY_VERTEX, n, ndim)
            end
        end

        # Update linked cells in joints

        # check if there are joints
        has_joints = false
        for cell in mesh.cells
            if cell.shape.class == JOINT_SHAPE
                has_joints = true
                break
            end
        end

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
                if cell.shape.class == JOINT_SHAPE
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
        if C.shape.class == SOLID_SHAPE #is_solid(C.shape)
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
        mesh = load_mesh_vtk(filename)
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


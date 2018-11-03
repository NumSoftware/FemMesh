#using FemMesh
#export mplot

const MOVETO = 1
const LINETO = 2
const CURVE3 = 3
const CURVE4 = 4
const CLOSEPOLY = 79

function plot_data_for_cell2d(points::Array{Array{Float64,1},1}, ctype::Int)
    npoints = length(points)
    shape = get_shape_from_vtk(VTKCellType(ctype), npoints, 2)

    if shape==LIN2
        verts = points
        codes = [ MOVETO, LINETO ]
    elseif shape == LIN3
        p1, p2, p3 = points
        cp    = 2*p3 - 0.5*p1 - 0.5*p2
        verts = [ p1, cp, p2 ]
        codes = [ MOVETO, CURVE3, CURVE3]
    elseif shape in (TRI3, QUAD4)
        n = shape==TRI3 ? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i=1:n
            p2 = i<n ? points[i+1] : points[1]
            push!(verts, p2)
            push!(codes, LINETO)
        end
    elseif shape in (TRI6, QUAD8, QUAD9)
        n = shape==TRI6 ? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i=1:n
            p1 = points[i]
            p2 = i<n ? points[i+1] : points[1]
            p3 = points[i+n]
            cp = 2*p3 - 0.5*p1 - 0.5*p2
            append!(verts, [cp, p2])
            append!(codes, [CURVE3, CURVE3])
        end
    elseif shape in (QUAD12, QUAD16)
        n = 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i=1:n
            p1 = points[i]
            p2 = i<n ? points[i+1] : points[1]
            p3 = points[2*i+3]
            p4 = points[2*i+4]
            cp2 = 1/6*(-5*p1+18*p2-9*p3+2*p4)
            cp3 = 1/6*( 2*p1-9*p2+18*p3-5*p4)
            append!(verts, [cp2, cp3, p2])
            append!(codes, [CURVE4, CURVE4, CURVE4])
        end
    else
        error("plot_data_for_cell2d: Not implemented for ", shape)
    end

    return verts, codes
end

function plot_data_for_cell3d(points::Array{Array{Float64,1},1}, ctype::Int)
    npoints = length(points)
    shape = get_shape_from_vtk(VTKCellType(ctype), npoints, 2)
    if shape == LIN2
        verts = points
    elseif shape in (TRI3, QUAD4)
        verts = points
    elseif shape == TRI6
        verts = points[[1,4,2,5,3,6]]
    elseif shape == QUAD8
        verts = points[[1,5,2,6,3,7,4,8]]
    end
    return verts
end


function mplot(items::Union{Block, Array}...; args...)
    # Get list of blocks and check type
    blocks = unfold(items)

    for item in blocks
        isa(item, Block) || error("mplot: Block object expected")
    end

    # Using Point and Cell types
    points = Array{Point,1}()
    cells  = Array{AbstractCell,1}()

    for bl in blocks
        append!(points, bl.points)

        if bl.shape.family==SOLID_SHAPE
            if bl.shape.ndim==2
                push!(cells, bl)
            else
                faces = get_faces(bl)
                append!(cells, faces)
            end
        elseif bl.shape.family==LINE_SHAPE
            lines = [ Cell(LIN2, bl.points[i-1:i]) for i=2:length(bl.points)]
            append!(cells, lines)
        else
            continue
        end

    end

    # points, conns and types
    #points  = [ p for bl in blocks for p in bl.points]
    for (i,p) in enumerate(points); p.id = i end
    pts_arr = [ [p.x, p.y, p.z] for p in points ]
    coords  = [ pts_arr[i][j] for i=1:length(pts_arr), j=1:3]
    conns   = [ [ p.id for p in c.points ] for c in cells]
    stypes  = [ Int(c.shape.vtk_type) for c in cells]
    #colors  = Dict("family" => [ ty==LIN2 ? 0.0 : 1.0 for ty in stypes])

    ugrid = UnstructuredGrid("Blocks", coords, conns, stypes)
    mplot(ugrid; args...)
end


function mplot(mesh::Mesh, filename::String=""; args...)
    ugrid = convert(UnstructuredGrid, mesh)
    if mesh.ndim==3
        # get surface cells and update ugrid
        scells = get_surface(mesh.cells)
        spoints = [ p for c in scells for p in c.points ]
        pt_ids = [ p.id for p in spoints ]
        oc_ids = [ c.ocell.id for c in scells ]

        # update data
        for (field, data) in ugrid.point_scalar_data
            ugrid.point_scalar_data[field] = data[pt_ids]
        end
        for (field, data) in ugrid.cell_scalar_data
            ugrid.cell_scalar_data[field] = data[oc_ids]
        end
        for (field, data) in ugrid.point_vector_data
            ugrid.point_vector_data[field] = data[pt_ids, :]
        end

        # renumerate points
        id_dict = Dict{Int, Int}( p.id => i for (i,p) in enumerate(spoints) )

        # points and cells
        pts_arr      = [ [p.x, p.y, p.z] for p in spoints ]
        ugrid.points = [ pts_arr[i][j] for i=1:length(pts_arr), j=1:3]
        ugrid.cells  = [ Int[ id_dict[p.id] for p in c.points ] for c in scells ]
        ugrid.cell_types = [ Int(c.shape.vtk_type) for c in scells ]
    end
    mplot(ugrid, filename; args...)
end


using PyCall # required

function mplot(ugrid::UnstructuredGrid, filename::String=""; axis=true, 
               pointmarkers=false, pointlabels=false, vectorfield=nothing, arrowscale=0.0,
               celllabels=false, fieldlims=nothing, cmap=nothing, colorbarscale=0.9, field=nothing,
               alpha=1.0, warpscale=0.0, highlightcell=0, elev=30.0, azim=45.0, dist=10.0,
               figsize=(4,2.5), leaveopen=false)

    @eval import PyPlot:plt, matplotlib, figure, art3D, Axes3D


    plt[:close]("all")

    plt[:rc]("font", family="serif", size=7)
    plt[:rc]("lines", lw=0.5)
    plt[:rc]("legend", fontsize=7)
    #plt[:rc]("figure", figsize=(4.5, 3))
    plt[:rc]("figure", figsize=figsize)

    ndim = all(p->p==0.0, ugrid.points[:,3]) ? 2 : 3
    ncells = length(ugrid.cells)
    npoints = size(ugrid.points,1)

    # All points coordinates
    XYZ = ugrid.points
    if warpscale>0
        found = haskey(ugrid.point_vector_data, "U")
        found || error("mplot: vector field U not found for warp")
        XYZ = ugrid.points .+ warpscale.*ugrid.point_vector_data["U"]
    end
    X = XYZ[:,1]
    Y = XYZ[:,2]
    Z = XYZ[:,3]

    limX = collect(extrema(X))
    limY = collect(extrema(Y))
    limZ = collect(extrema(Z))
    limX = limX + 0.05*[-1, 1]*norm(limX)
    limY = limY + 0.05*[-1, 1]*norm(limY)
    limZ = limZ + 0.05*[-1, 1]*norm(limZ)
    L = max(norm(limX), norm(limY), norm(limZ))

    # Configure plot
    if ndim==3
        ax = @eval Axes3D(figure())
        ax[:set_aspect]("equal")
        
        # Set limits
        meanX = mean(limX)
        meanY = mean(limY)
        meanZ = mean(limZ)
        limX = [meanX-L/2, meanX+L/2]
        limY = [meanY-L/2, meanY+L/2]
        limZ = [meanZ-L/2, meanZ+L/2]
        ax[:set_xlim]( meanX-L/2, meanX+L/2)
        ax[:set_ylim]( meanY-L/2, meanY+L/2)
        ax[:set_zlim]( meanZ-L/2, meanZ+L/2)
        #ax[:scatter](limX, limY, limZ, color="w", marker="o", alpha=0.0)

        # Labels
        ax[:set_xlabel]("x")
        ax[:set_ylabel]("y")
        ax[:set_zlabel]("z")

        if axis == false
            ax[:set_axis_off]()
        end
    else
        ax = plt[:gca]()
        plt[:axes]()[:set_aspect]("equal", "datalim")

        # Set limits
        ax[:set_xlim](limX[1], limX[2])
        ax[:set_ylim](limY[1], limY[2])

        # Labels
        ax[:set_xlabel]("x")
        ax[:set_ylabel]("y")
        if axis == false
            plt[:axis]("off")
        end
    end

    #cm = colors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],256)
    #colors =  matplotlib[:colors]
    cm = matplotlib[:colors][:ListedColormap]([(1,0,0),(0,1,0),(0,0,1)],256)

    cdict = Dict("red"   => [(0.0,  0.8, 0.8), (0.5, 0.7, 0.7), (1.0, 0.0, 0.0)],
                 "green" => [(0.0,  0.2, 0.2), (0.5, 0.7, 0.7), (1.0, 0.2, 0.2)],
                 "blue"  => [(0.0,  0.0, 0.0), (0.5, 0.7, 0.7), (1.0, 0.6, 0.6)])

    cm = matplotlib[:colors][:LinearSegmentedColormap]("my_colormap",cdict,256)

    if cmap==nothing # cmap may be "bone", "plasma", "inferno", etc.
        cmap = cm
    end

    has_field = field != nothing
    if has_field
        field = string(field)
        found = haskey(ugrid.cell_scalar_data, field)
        if found
            fvals = ugrid.cell_scalar_data[field]
        else
            found = haskey(ugrid.point_scalar_data, field)
            found || error("mplot: field $field not found")
            data  = ugrid.point_scalar_data[field]
            fvals = [ mean(data[ugrid.cells[i]]) for i=1:ncells ]
        end
        fieldlims==nothing && (fieldlims = extrema(fvals))
    end

    # Plot cells
    if ndim==3
        # Plot line cells
        for i=1:ncells 
            ctype = ugrid.cell_types[i]
            ctype in (3,21) || continue
            con = ugrid.cells[i]
            X = XYZ[con, 1]
            Y = XYZ[con, 2]
            Z = XYZ[con, 3]
            plt[:plot](X, Y, Z, color="red", lw=1.0)
        end

        # Plot surface cells
        all_verts  = []
        for i=1:ncells 
            ctype = ugrid.cell_types[i]
            ctype in (10,12,13,14,24,25) && @warn "mplot: cannot directly plot 3d cell" vtk_type=ctype
            ctype in (5,22,9,23,28) || continue
            con = ugrid.cells[i]
            points = [ XYZ[i,1:3] for i in con ]
            verts = plot_data_for_cell3d(points, ctype)
            push!(all_verts, verts)
        end

        #cltn = art3D[:Poly3DCollection](all_verts, cmap=cmap, facecolor="aliceblue", edgecolor="black", lw=2)
        cltn = @eval art3D[:Poly3DCollection]($all_verts, cmap=$cmap, facecolor="aliceblue", edgecolor="black",
                                              lw=0.5, alpha=$alpha)

        if has_field
            cltn[:set_array](fvals)
            cltn[:set_clim](fieldlims)
            #cbar = plt[:colorbar](cltn, label=field, shrink=0.9)
            cbar = plt[:colorbar](cltn, label=field, shrink=colorbarscale, aspect=20*colorbarscale)
            cbar[:ax][:tick_params](labelsize=7)
        end
        @eval $ax[:add_collection3d]($cltn)

    else
        all_patches = []
        for i=1:ncells
            ctype = ugrid.cell_types[i]
            ctype in (3,21,5,22,9,23,28) || continue
            #lw = ugrid.cell_types[i] in (3,21) ? 1.0 : 0.5
            
            con = ugrid.cells[i]
            points = [ XYZ[i,1:2] for i in con ]
            verts, codes = plot_data_for_cell2d(points, ctype)
            path  = matplotlib[:path][:Path](verts, codes)
            patch = matplotlib[:patches][:PathPatch](path)
            push!(all_patches, patch)

            if highlightcell==i
                patch = matplotlib[:patches][:PathPatch](path, facecolor="cadetblue", edgecolor="black")
                ax[:add_patch](patch)
            end

        end
        cltn = matplotlib[:collections][:PatchCollection](all_patches, cmap=cmap, edgecolor="black", 
                                                          facecolor="aliceblue", lw=0.5)
        if has_field
            cltn[:set_array](fvals)
            cltn[:set_clim](fieldlims)
            #cbar = plt[:colorbar](cltn, label=field, shrink=0.9)
            cbar = plt[:colorbar](cltn, label=field, shrink=colorbarscale, aspect=0.8*20*colorbarscale)
            cbar[:ax][:tick_params](labelsize=7)
            cbar[:outline][:set_linewidth](0.5)
        end
        ax[:add_collection](cltn)
    end

    # Draw points
    if pointmarkers
        if ndim==3
            ax[:scatter](X, Y, Z, color="k", marker="o", s=1)
        else
            plt[:plot](X, Y, color="black", marker="o", markersize=3, lw=0)
        end
    end

    # Draw arrows
    if vectorfield!=nothing && ndim==2
        data = ugrid.point_vector_data[vectorfield]
        color = "blue"
        if arrowscale==0
            plt[:quiver](X, Y, data[:,1], data[:,2], color=color)
        else
            plt[:quiver](X, Y, data[:,1], data[:,2], color=color, scale=1.0/arrowscale)
        end
    end

    # Draw point numbers
    if pointlabels
        npoints = length(X)
        for i=1:npoints
            x = X[i] + 0.01*L
            y = Y[i] - 0.01*L
            z = Z[i] - 0.01*L
            if ndim==3
                ax[:text](x, y, z, i, va="center", ha="center", backgroundcolor="none")
            else
                ax[:text](x, y, i, va="top", ha="left", backgroundcolor="none")
            end
        end
    end

    # Draw cell numbers
    if celllabels && ndim==2
        for i=1:ncells
            coo = ugrid.points[ ugrid.cells[i], : ]
            x = mean(coo[:,1])
            y = mean(coo[:,2])
            ax[:text](x, y, i, va="top", ha="left", color="blue", backgroundcolor="none", size=8)
        end
    end

    if ndim==3 
        ax[:view_init](elev=elev, azim=azim)
        ax[:dist] = dist
    end

    if filename==""
        plt[:show]()
    else
        plt[:savefig](filename, bbox_inches="tight", pad_inches=0.00, format="pdf")
    end

    # close the figure
    #@show plt[:isinteractive]()
    #@show isinteractive()

    #plt[:isinteractive]() || plt[:close]("all")
    #isinteractive() || plt[:close]("all")

    # Do not close if in IJulia
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return
    end
    if !leaveopen
        plt[:close]("all")
    end

    return
end

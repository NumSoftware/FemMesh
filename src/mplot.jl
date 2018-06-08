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
        n = shape==TRI3? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i=1:n
            p2 = i<n? points[i+1] : points[1]
            push!(verts, p2)
            push!(codes, LINETO)
        end
    elseif shape in (TRI6, QUAD8, QUAD9)
        n = shape==TRI6? 3 : 4
        codes = [ MOVETO ]
        verts = [ points[1] ]
        for i=1:n
            p1 = points[i]
            p2 = i<n? points[i+1] : points[1]
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
            p2 = i<n? points[i+1] : points[1]
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
    cells  = Array{Cell,1}()

    for bl in blocks
        npts = size(bl.coords,1)
        pts = [ Point(bl.coords[i,:]) for i=1:npts ]
        append!(points, pts)
        bl_type = typeof(bl)
        ndim = 0

        if bl_type==BlockInset
            shapetype = LIN2
        elseif npts==4 && bl_type==Block2D
            shapetype = QUAD4
            ndim = 2
        elseif npts==8 && bl_type==Block2D
            ndim = 2
            shapetype = QUAD8
        elseif npts==8 && bl_type==Block3D
            ndim = 3
            shapetype = HEX8
        elseif npts==20 && bl_type==Block3D
            ndim = 3
            shapetype = HEX20
        end

        if ndim==2
            cell = Cell(shapetype, pts)
            push!(cells, cell)
            append!(points, pts)
        else
            if shapetype == LIN2
                lines = [ Cell(shapetype, [pts[i-1], pts[i]]) for i=2:npts ]
                append!(cells, lines)
            else
                cell  = Cell(shapetype, pts)
                faces = get_faces(cell)
                append!(cells, faces)
            end
        end
    end

    # points, cells and types
    for (i,p) in enumerate(points); p.id = i end
    pts_arr = [ [p.x, p.y, p.z] for p in points ]
    coords  = [ pts_arr[i][j] for i=1:length(pts_arr), j=1:3]
    conns   = [ [ p.id for p in c.points ] for c in cells ]
    stypes  = [ c.shape.vtk_type for c in cells ]
    colors  = Dict("family" => [ ty==LIN2? 0.0 : 1.0 for ty in stypes])

    ugrid = VTK_unstructured_grid("Blocks", coords, conns, stypes, cell_scalar_data=colors)
    mplot(ugrid; alpha=0.4, args...)
end


function mplot(mesh::Mesh; args...)
    ugrid = convert(VTK_unstructured_grid, mesh)
    if mesh.ndim==3
        # get surface cells and update ugrid
        scells = get_surface(mesh.cells)
        spoints = get_points(scells)
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
        ugrid.cell_types = [ c.shape.vtk_type for c in scells ]
    end
    mplot(ugrid; args...)
end


using PyCall

function mplot(ugrid::VTK_unstructured_grid, filename::String=""; axis=true, 
               pointmarkers=false, pointlabels=false, celllabels=false, elev=30, azim=45,
               fieldlims=nothing, cmap=nothing, field=nothing, alpha=1.0, warpscale=0.0)

    @eval import PyPlot:plt, matplotlib, figure, art3D, Axes3D

    plt[:rc]("font", family="serif", size=10)

    ndim = all(p->p==0.0, ugrid.points[:,3]) ? 2: 3
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
        ax[:scatter](limX, limY, limZ, color="w", marker="o", alpha=0.0)

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

    if cmap==nothing
        cmap = cm
    end

    has_field = field != nothing
    if has_field
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
            ctype in (5,22,9,23,28) || continue
            con = ugrid.cells[i]
            points = [ XYZ[i,1:3] for i in con ]
            verts = plot_data_for_cell3d(points, ctype)
            push!(all_verts, verts)
        end

        #cltn = art3D[:Poly3DCollection](all_verts, cmap=cmap, facecolor="aliceblue", edgecolor="black", lw=2)
        cltn = @eval art3D[:Poly3DCollection]($all_verts, cmap=$cmap, facecolor="aliceblue", edgecolor="black",
                                              lw=1., alpha=$alpha)

        if has_field
            cltn[:set_array](fvals)
            cltn[:set_clim](fieldlims)
            plt[:colorbar](cltn, label=field, shrink=0.7)
        end
        @eval $ax[:add_collection3d]($cltn)

    else
        all_patches = []
        for i=1:ncells
            ctype = ugrid.cell_types[i]
            ctype in (3,21,5,22,9,23,28) || continue
            lw = ugrid.cell_types[i] in (3,21)? 1 : 0.5
            
            #verts, codes = plot_data_for_cell2d(ugrid, i)
            con = ugrid.cells[i]
            #coords = XYZ[con, 1:2]
            points = [ XYZ[i,1:2] for i in con ]
            verts, codes = plot_data_for_cell2d(points, ctype)
            path  = matplotlib[:path][:Path](verts, codes)
            patch = matplotlib[:patches][:PathPatch](path, lw=lw)
            push!(all_patches, patch)

        end
        cltn = matplotlib[:collections][:PatchCollection](all_patches, cmap=cmap, edgecolor="black", 
                                                          facecolor="aliceblue", lw=1)
        if has_field
            cltn[:set_array](fvals)
            cltn[:set_clim](fieldlims)
            plt[:colorbar](cltn, label=field, shrink=0.7)
        end
        ax[:add_collection](cltn)
    end

    # Draw points
    if pointmarkers
        if ndim==3
            ax[:scatter](X, Y, Z, color="k", marker="o")
        else
            plt[:plot](X, Y, color="black", marker="o", markersize=4, lw=0)
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
            ax[:text](x, y, i, va="top", ha="left", color="blue", backgroundcolor="none")
        end
    end

    ndim==3 && ax[:view_init](elev=elev, azim=azim)

    if filename==""
        plt[:show]()
    else
        plt[:savefig](filename, bbox_inches="tight", pad_inches=0.25, format="png")
    end

    # close de figure
    isinteractive() || plt[:close]("all")
    return nothing
end


#coord = [ 0 0; 9 0; 18 0; 0 9; 9 9; 18 9.]
#conn  = [ 1 2; 1 5; 2 3; 2 6; 2 5; 2 4; 3 6; 3 5; 4 5; 5 6]
#blk  = BlockTruss(coord, conn)
#msh = Mesh(blk, reorder=false)

#blk = Block2D([0 0; 8 4], nx=8, ny=4)
#blb = BlockInset([0.5 0.5; 3 1; 3.5 2.5; 7.5 3.5 ], curvetype=3)
#msh = Mesh(blk, blb)
#split!(msh)

#blk = Block2D( [ 0.0 0.0 ; 1.0 0.0 ; 0.707 0.707; 0.0 1.0 ; 0.5 0.0 ; 0.923 0.382; 0.382 0.923 ; 0.0 0.5 ], shape=TRI6, nx=6, ny=6)
#msh = Mesh(blk)
#smooth!(msh, maxit=4)

#blk = Block2D( [ 0.0 0.0 ; 1.0 0.0 ; 0.707 0.707; 0.0 1.0 ; 0.5 0.0 ; 0.923 0.382; 0.382 0.923 ; 0.0 0.5 ], shape=TET4, nx=6, ny=6)
#blk = extrude(blk, len=0.5, n=3)
#blk.shape = TET4
#blk = rotate(blk, base=[0,0,0], axis=[1,0,0], angle=90)
#msh = Mesh(blk, reorder=false)
#smooth!(msh, maxit=8)

#blk = Block3D( [0 0 0 ; 1 1 1], nx=3, ny=3, shape=TET10)
#msh = Mesh(blk)
#save(msh, "msh.vtk")
#smooth!(msh)

#draw(msh, axis=false, pointlabels=true, quality=true)
#draw(msh, points=true, pointlabels=true, quality=true)

#blk = Block2D( [ 0.0 0.0 ; 1.0 0.0 ; 0.707 0.707; 0.0 1.0 ; 0.5 0.0 ; 0.923 0.382; 0.382 0.923 ; 0.0 0.5 ], shape=TRI6, 
              #nx=6, ny=6)
#blk = Block3D( [0 0 0 ; 1 1 1], nx=6, ny=6, nz=6, shape=TET10)

#msh = Mesh(blk)
#mplot(convert(VTK_unstructured_grid,msh), axis=false, cellfield="quality")
#mplot(msh, axis=false, field="point-id", pointmarkers=true, pointlabels=true, celllabels=true)
#mplot(msh, axis=false, cellfield="quality")
#mplot(msh, axis=false, cellfield="quality")

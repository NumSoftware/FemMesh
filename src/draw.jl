#using FemMesh

using PyCall

const MOVETO = 1
const LINETO = 2
const CURVE3 = 3
const CURVE4 = 4
const CLOSEPOLY = 79

function draw_cell2d(c::Cell)

    if c.shape == LIN2
        verts = [ [p.x, p.y] for p in c.points ]
        codes = [ MOVETO, LINETO ]
    elseif c.shape == LIN3
        verts = [ [p.x, p.y] for p in c.points ]
        pc    = 2*verts[3] - 0.5*verts[1] - 0.5*verts[2]
        verts[2], verts[3] = pc, verts[2]
        codes = [ MOVETO, CURVE3, CURVE3]
    elseif c.shape in (TRI3, QUAD4)
        codes = []
        verts = []
        for edge in get_faces(c)
            everts = [ [p.x, p.y] for p in edge.points ]
            if length(codes)==0
                push!(verts, everts[1])
                push!(codes, MOVETO)
            end
            push!(verts, everts[2])
            push!(codes, LINETO)
        end
    elseif c.shape in (TRI6, QUAD8, QUAD9)
        codes = []
        verts = []
        for edge in get_faces(c)
            everts = [ [p.x, p.y] for p in edge.points ]
            if length(codes)==0
                push!(verts, everts[1])
                push!(codes, MOVETO)
            end
            pc    = 2*everts[3] - 0.5*everts[1] - 0.5*everts[2]
            push!(verts, pc)
            push!(verts, everts[2])
            append!(codes, [CURVE3, CURVE3])
        end
    elseif c.shape in (QUAD12, QUAD16)
        codes = []
        verts = []
        for edge in get_faces(c)
            everts = [ [p.x, p.y] for p in edge.points ]
            if length(codes)==0
                push!(verts, everts[1])
                push!(codes, MOVETO)
            end
            # Calculate Bezier control points
            p1 = everts[1]; p2 = everts[3]; p3 = everts[4]; p4 = everts[2]
            pc2 = 1/6*(-5*p1+18*p2-9*p3+2*p4  )
            pc3 = 1/6*( 2*p1-9*p2+18*p3-5*p4  )
            push!(verts, pc2)
            push!(verts, pc3)
            push!(verts, p4)
            append!(codes, [CURVE4, CURVE4, CURVE4])
        end
    else
        error("draw_cell: Not implemented for ", c.shape)
    end

    return verts, codes
end

function draw_cell3d(c::Cell)
    if c.shape in (TRI3, QUAD4)
        verts = [ Float64[p.x, p.y, p.z] for p in c.points ]
        #verts = push!([], verts)
    elseif c.shape == TRI6
        verts = [ Float64[p.x, p.y, p.z] for p in c.points[ [1,4,2,5,3,6] ] ]
        #verts = push!([], verts)
    elseif c.shape == QUAD8
        verts = [ Float64[p.x, p.y, p.z] for p in c.points[ [1,5,2,6,3,7,4,8] ] ]
        #verts = push!([], verts)
    end
    return verts
end

function draw(mesh::Mesh, filename::String=""; axis=true, points=false, pointlabels=false, celllabels=false, quality=false, lims=nothing,
    cmap=nothing)
    if is_windows() eval(Expr(:using, :PyPlot)) end

    @pyimport matplotlib.pyplot as plt
    @pyimport matplotlib.path as path
    #Path = path.pymember("Path")
    Path = path.Path
    @pyimport matplotlib.patches as patches
    @pyimport matplotlib.collections as collections

    @pyimport matplotlib.colors as colors
    @pyimport matplotlib.colorbar as colorbar


    plt.rc("font", family="serif",size=15)


    # All points coordinates
    coords = [ getfield(point, field) for point in mesh.points, field in (:x, :y, :z)]
    X = coords[:,1]
    Y = coords[:,2]
    Z = coords[:,3]
    ndim = mesh.ndim

    limX = collect(extrema(X))
    limY = collect(extrema(Y))
    limZ = collect(extrema(Z))
    limX = limX + 0.05*[-1, 1]*norm(limX)
    limY = limY + 0.05*[-1, 1]*norm(limY)
    limZ = limZ + 0.05*[-1, 1]*norm(limZ)
    L = max(norm(limX), norm(limY), norm(limZ))

    # Plot
    if ndim==3
        @pyimport mpl_toolkits.mplot3d as mplot3d
        @pyimport matplotlib.collections as collections
        @pyimport mpl_toolkits.mplot3d.art3d as art3d

        ax = plt.gca(projection="3d")
        #Axes3D = mplot3d.Axes3D
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
        fig = plt.figure()
        #ax = fig[:add_axes]()
        ax = plt.gca()
        plt.axes()[:set_aspect]("equal", "datalim")

        # Set limits
        ax[:set_xlim](limX[1], limX[2])
        ax[:set_ylim](limY[1], limY[2])

        # Labels
        #plt.plot([1,2]
        plt.xlabel("x")
        plt.ylabel("y")
        if axis == false
            plt.axis("off")
        end
    end

    cm = colors.ListedColormap([(1,0,0),(0,1,0),(0,0,1)],256)

    cdict = Dict("red"   => [(0.0,  0.8, 0.8), (0.5, 0.7, 0.7), (1.0, 0.0, 0.0)],
                 "green" => [(0.0,  0.2, 0.2), (0.5, 0.7, 0.7), (1.0, 0.2, 0.2)],
                 "blue"  => [(0.0,  0.0, 0.0), (0.5, 0.7, 0.7), (1.0, 0.6, 0.6)])

    cm = colors.LinearSegmentedColormap("my_colormap",cdict,256)

    if cmap==nothing
        cmap = cm
    end

    # Draw cells
    facecolor = "aliceblue"
    if lims==nothing
        lims = [mesh.qmin, 1.0]
    end
    cbnorm = colors.Normalize(vmin=lims[1], vmax=lims[2])
    if ndim==3
        all_verts = []
        surface = get_surface(mesh.cells)
        for cell in surface
            verts = draw_cell3d(cell)
            #if quality facecolor = cm(cell.ocell.quality) end
            if cell.shape.class==LINE_SHAPE facecolor="none" end

            push!(all_verts, verts)
            #poly = art3d.Poly3DCollection(verts, facecolor=facecolor)
            #ax[:add_collection3d]( poly )
        end
        colors = [cell.ocell.quality for cell in surface ]
        cltn   = art3d.Poly3DCollection(all_verts, cmap=cmap)
        cltn[:set_array](colors)
        cltn[:set_clim](lims)
        ax[:add_collection3d](cltn)
        if quality
            plt.colorbar(cltn)
        end
    else
        all = []
        for cell in mesh.cells
            if Int(cell.shape.class)>2 continue end
            verts, codes = draw_cell2d(cell)
            path  = Path(verts, codes)
            #if quality facecolor = cm(cbnorm(cell.quality)) end
            lw = 0.5
            if cell.shape.class == LINE_SHAPE
                facecolor="none" 
                lw = 3
            end
            patch = patches.PathPatch(path, facecolor=facecolor, lw=lw)
            #ax[:add_patch](patch)
            push!(all, patch)
        end
        #colors = [cbnorm(cell.quality) for cell in mesh.cells ]
        colors = [cell.quality for cell in mesh.cells ]
        cltn   = collections.PatchCollection(all, cmap=cmap)
        cltn[:set_array](colors)
        cltn[:set_clim](lims)
        ax[:add_collection](cltn)
        if quality
            plt.colorbar(cltn)
        end

    end

    # Colorbar
    if false # quality
        #fig = plt.figure()
        ax1 = fig[:add_axes]([0.9, 0.1, 0.05, 0.8])
        ax1[:locator_params](tight=true, nbins=5)
        #cax = plt.inset_axes(ax, width="8%", height="70%", loc=4)
        cb = colorbar.ColorbarBase(ax1, cmap=cm, norm=cbnorm, orientation="vertical")
        cb[:set_label]("Cell quality")

    end

    # Draw points
    if points
        if ndim==3
            ax[:scatter](X, Y, Z, color="k", marker="o")
        else
            plt.plot(X, Y, color="black", marker="o", markersize=4)
        end
    end

    # Draw point numbers
    if pointlabels
        for point in mesh.points
            x = point.x + 0.01*L
            y = point.y - 0.01*L
            z = point.z - 0.01*L
            if ndim==3
                ax[:text](x, y, z, point.id, va="center", ha="center", backgroundcolor="none")
            else
                ax[:text](x, y, point.id, va="top", ha="left", backgroundcolor="none")
            end
        end
    end

    # Draw cell numbers
    if celllabels
        for cell in mesh.cells
            C = get_coords(cell)
            x = mean(C[:,1])
            y = mean(C[:,2])
            z = mean(C[:,3])
            ax[:text](x, y, cell.id, va="top", ha="left", color="blue", backgroundcolor="none")
        end
    end

    #ax[:view_init](elev=20, azim=-49)

    if filename==""
        plt.show()
    else
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.25, format="pdf")
    end

    # close de figure
    plt.close("all")
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


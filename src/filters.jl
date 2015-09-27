
# Get an array with shares for all points
function get_shares(cells::Array{Cell,1}, points::Point=Array(Point,0))
    # points should be numbered

    # get all points from cells if needed
    if length(points==0)
        pointsd = Dict{Uint64, Point}()
        for cell in cells
            for point in cell.points
                pointsd[hash(point)] = point
            end
        end
        points = values(pointsd)
    end

    # get incidence array (shares) (fast)
    np = length(points)
    N  = [ Cell[] for i=1:np]
    for cell in cells
        for pt in cell.points
            push!(N[pt.id], cell)
        end
    end

    return N
end

# Get an array with shares for all points
function get_shares(mesh::Mesh)
    # get incidence array (shares) (fast)
    np = length(mesh.points)
    C  = [ Cell[] for i=1:np]
    for cell in mesh.cells
        for pt in cell.points
            push!(C[pt.id], cell)
        end
    end

    return C
end

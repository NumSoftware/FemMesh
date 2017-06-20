

# Get an array with shares for all points
function get_patches(cells::Array{Cell,1})
    # get all points from cells if needed
    pointsd = Dict{UInt64, Point}()
    for cell in cells
        for point in cell.points
            pointsd[hash(point)] = point
        end
    end
    
    points = values(pointsd)
    np     = length(points)

    # backup points ids
    bk_pt_id = [ pt.id for pt in points ]
    for i=1:np
        points[i].id = i
    end

    # get incidence array
    P  = [ Cell[] for i=1:np ]
    for cell in cells
        for pt in cell.points
            push!(P[pt.id], cell)
        end
    end

    # restore points ids
    for i=1:np
        points[i].id = bk_pt_id[i]
    end

    return P
end

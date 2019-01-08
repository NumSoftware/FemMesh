
mutable struct UnstructuredGrid
    title::String
    points::Array{Float64,2}
    cells ::Array{Array{Int64,1},1}
    cell_types::Array{Int64,1}
    point_scalar_data::Dict{String,Array}
    cell_scalar_data ::Dict{String,Array}
    point_vector_data::Dict{String,Array}
    version::String
    function UnstructuredGrid(title, points, cells, cell_types; point_scalar_data=Dict(), cell_scalar_data=Dict(), point_vector_data=Dict())
        return new(title, points, cells, cell_types, point_scalar_data, cell_scalar_data, point_vector_data, "3.0")
    end
end

function save_vtk(vtk_data::UnstructuredGrid, filename::String)
    # Saves a UnstructuredGrid
    npoints = size(vtk_data.points, 1)
    ncells = length(vtk_data.cells)

    # Number of total connectivities
    nconns = 0
    for cell in vtk_data.cells
        nconns += 1 + length(cell)
    end


    # Open filename
    f = open(filename, "w")

    println(f, "# vtk DataFile Version $(vtk_data.version)")
    println(f, vtk_data.title)
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", npoints, " float64")

    # Write nodes
    for i=1:npoints
        @printf f "%23.15e %23.15e %23.15e \n" vtk_data.points[i,1] vtk_data.points[i,2] vtk_data.points[i,3]
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ", ncells, " ", nconns)
    for cell in vtk_data.cells
        print(f, length(cell), " ")
        for id in cell
            print(f, id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write elem types
    println(f, "CELL_TYPES ", ncells)
    for ty in vtk_data.cell_types
        println(f, ty)
    end
    println(f)

    has_point_scalar_data = !isempty(vtk_data.point_scalar_data)
    has_point_vector_data = !isempty(vtk_data.point_vector_data)
    has_point_data = has_point_scalar_data || has_point_vector_data
    has_cell_data  = !isempty(vtk_data.cell_scalar_data)

    # Write point data
    if has_point_data
        println(f, "POINT_DATA ", npoints)
        # Write scalar data
        if has_point_vector_data
            for (field,D) in vtk_data.point_vector_data
                isempty(D) && continue
                dtype = eltype(D)<:Integer ? "int" : "float64"
                println(f, "VECTORS ", "$field $dtype")
                for i=1:npoints
                    @printf f "%23.15e %23.15e %23.15e \n" D[i,1] D[i,2] D[i,3]
                end
            end
        end
        # Write vector data
        if has_point_scalar_data
            for (field,D) in vtk_data.point_scalar_data
                isempty(D) && continue
                dtype = eltype(D)<:Integer ? "int" : "float64"
                println(f, "SCALARS $field $dtype 1")
                println(f, "LOOKUP_TABLE default")
                if dtype=="float64"
                    for i=1:npoints
                        @printf f "%23.10e" D[i]
                    end
                else
                    for i=1:npoints
                        @printf f "%10d" D[i]
                    end
                end
                println(f)
            end
        end
    end

    # Write cell data
    if has_cell_data
        println(f, "CELL_DATA ", ncells)
        for (field,D) in vtk_data.cell_scalar_data
            isempty(D) && continue
            dtype = eltype(D)<:Integer ? "int" : "float64"
            println(f, "SCALARS $field $dtype 1")
            println(f, "LOOKUP_TABLE default")
            if dtype=="float64"
                for i=1:ncells
                    @printf f "%23.10e" D[i]
                end
            else
                for i=1:ncells
                    @printf f "%10d" D[i]
                end
            end
            println(f)
        end
    end

    close(f) 

    return nothing
end


function read_ugrid_vtk(filename::String)
    #file = open(filename)

    # read nodal information
    alltext = read(filename, String)
    data    = split(alltext)
    #close(file)

    local points, cells, cell_types, npoints, ncells
    point_scalar_data = Dict{String,Array}()
    point_vector_data = Dict{String,Array}()
    cell_scalar_data  = Dict{String,Array}()
    reading_point_data = false
    reading_cell_data  = false

    TYPES = Dict("float32"=>Float32, "float64"=>Float64, "int"=>Int64)

    idx = 1
    while idx<=length(data)
        if data[idx] == "DATASET"
            gridtype = data[idx+1]
            gridtype == "UNSTRUCTURED_GRID" || error("load_VTK_unstructured_grid: this reader only support files of VTK UNSTRUCTURED_GRID")
        end

        # read points
        if data[idx] == "POINTS"
            npoints = parse(Int64, data[idx+1]) # read number of points
            points  = zeros(npoints,3)
            idx += 2
            for i=1:npoints
                points[i,1] = parse(Float64, data[idx+1])
                points[i,2] = parse(Float64, data[idx+2])
                points[i,3] = parse(Float64, data[idx+3])
                idx += 3
            end
        end

        # read cells connectivities
        if data[idx] == "CELLS"
            ncells = parse(Int64, data[idx+1])
            ncdata = parse(Int64, data[idx+2])
            idx += 2

            cells = Array{Int,1}[]
            for i=1:ncells
                npts = parse(Int64, data[idx+1])
                idx += 1
                conn = Int[]
                for j=1:npts
                    idx += 1
                    id = parse(Int64, data[idx]) + 1
                    push!(conn, id)
                end
                push!(cells, conn)
            end
        end

        # read type of cells
        if data[idx] == "CELL_TYPES"
            idx += 1
            cell_types = Int[]
            for i=1:ncells
                idx += 1
                vtk_shape = parse(Int64, data[idx])
                push!(cell_types, vtk_shape)
            end
        end

        if data[idx] == "POINT_DATA"
            idx += 1
            reading_point_data = true
            reading_cell_data  = false
        end

        if data[idx] == "CELL_DATA"
            idx += 1
            reading_cell_data  = true
            reading_point_data = false
        end

        if data[idx] == "VECTORS" && reading_point_data
            label = data[idx+1]
            idx += 2
            vectors = zeros(npoints,3)
            for i=1:npoints
                vectors[i,1] = parse(Float64, data[idx+1])
                vectors[i,2] = parse(Float64, data[idx+2])
                vectors[i,3] = parse(Float64, data[idx+3])
                idx += 3
            end
            point_vector_data[label] = vectors
        end

        if data[idx] == "SCALARS" && reading_point_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], npoints)
            for i=1:npoints
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            point_scalar_data[label] = scalars
        end

        if data[idx] == "SCALARS" && reading_cell_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], ncells)
            for i=1:ncells
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            cell_scalar_data[label] = scalars
        end

        idx += 1

    end

    ugrid = UnstructuredGrid("VTK unstructured grid", points, cells, cell_types, 
                                point_scalar_data=point_scalar_data,
                                point_vector_data=point_vector_data,
                                cell_scalar_data =cell_scalar_data)

    return ugrid
end


function read_ugrid_tetgen(filekey::String)
    # reads files .node and .ele

    alltext  = read(filename*".node", String)
    nodedata = split(alltext, "\n")
    alltext  = read(filename*".ele", String)
    celldata = split(alltext, "\n")

    local points, cells, cell_types, npoints, ncells
    point_scalar_data = Dict{String,Array}()
    point_vector_data = Dict{String,Array}()
    cell_scalar_data  = Dict{String,Array}()

    # read nodal information

    f=open(filekey*".node")
    line = readline(f)
    npoints, ndim, _, _ = parse.(Int, split(line))
    points  = zeros(npoints,3)

    for i=1:npoints
        line = readline(f)
        items = split(line)
        points[i,1] = parse(Float64, items[2])
        points[i,2] = parse(Float64, items[3])
        points[i,3] = parse(Float64, items[4])
    end
    close(f)

    # read cell information

    f=open(filekey*".cell")
    line = readline(f)
    ncells, npts, _ = parse.(Int, split(line))
    cells = Array{Int,1}[]

    for i=1:npoints
        line = readline(f)
        items = split(line)
        pts = parse.(Int, items[2:end])
        push!(cells, pts)
    end
    close(f)

    # cell types
    #
    cell_types = 10*ones(Int, ncells) 


    ugrid = UnstructuredGrid("Unstructured grid from tetgen", points, cells, cell_types)

    return ugrid
end

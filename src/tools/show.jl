
import Base.show
export show

function print_short(io::IO, obj::Any)
    ty = typeof(obj) # update type
    print(io, ty, " ")
    if isbits(ty) || ty in [Symbol, Expr]
        print(io, obj)
    elseif ty<:AbstractString
        print(io, "\"", obj, "\"")
    elseif ty<:AbstractArray
        print_short_array(obj)
    end
end

function print_short_array(io::IO, array::AbstractArray)
    maxn = 8
    if length(array)>maxn
        str = string(array[1:maxn])[1:end-1]*", ...]"
        if str[1] != '['
            str = str[findfirst(str,'['): end]
        end
        print(io, str)
    else
        print(io, array)
    end
end


function print_datatype_fields(io::IO, obj::Any)
    ctype = typeof(obj)
    print(io, ctype, " with:")
    for (field, ty) in zip(fieldnames(ctype), ctype.types )
        print(io, "\n  ", field, ": ")
        if !isdefined(obj, field)
            #print(io, "#undef ", ty, " object")
            print(io, ty, "#undef ")
            continue
        end

        item = getfield(obj, field)
        ty   = typeof(item) # update type
        #if ty<:Real || ty in [Symbol, Expr]
            #print(io, item)
        #elseif ty<:AbstractString
            #print(io, "\"", item, "\"")
        #elseif :name in fieldnames(ty) 
            #if isdefined(item, :name)
                #print(io, ty, " ", getfield(item, :name))
            #end
        #elseif 

        if :name in fieldnames(ty) 
            if isdefined(item, :name)
                print(io, "name field ", getfield(item, :name))
            end
        elseif :id in fieldnames(ty) 
            if isdefined(item, :id )
                print(io, "id field ", getfield(item, :id ))
            end
        elseif ty<:AbstractArray && ndims(ty)==1
            array = getfield(obj, field)
            if eltype(ty)<:Real && length(size(array))==1
                print_short_array(io, array)
                #maxn = 8
                #if length(array)>maxn
                    #str = string(array[1:maxn])[1:end-1]*", ...]"
                    #print(io, str)
                #else
                    #print(io, array)
                #end
            else
                #print(io, length(array), "-element ", ty, " object")
                print(io, ty, " length ", length(array))
            end

            n = length(array)
            if n>0
                maxn = 10
                # add items IDs if available
                has_id = all( :id in fieldnames(typeof(item)) for item in array[1:min(maxn,n)] )
                if has_id
                    ids = [ item.id for item in array ]
                    print(io, ", id fields ")
                    print_short_array(io, ids)
                    #print(io, " with id field ")
                    #if length(ids)>maxn
                        #str = string(ids[1:maxn-1])[1:end-1]*", ...]"
                        #print(io, str)
                    #else
                        #print(io, ids)
                    #end
                end
                # add items names if available
                has_name = all( :name in fieldnames(typeof(item)) for item in array[1:min(maxn,n)] )
                if has_name
                    names = [ item.name for item in array ]
                    print(io, ", name fields ")
                    print_short_array(io, names)
                    #print(io, " with name field ")
                    #if length(names)>maxn
                        #str = string(names[1:maxn-1])[1:end-1]*", ...]"
                    #else
                        #str = string(names)
                    #end
                    #str = str[findfirst(str,'['): end]
                    #print(io, str)
                end
            end
        elseif ty<:AbstractArray
            array = getfield(obj, field)
            str   = replace(string(size(array)), ", " => "×")[2:end-1]
            #print(io, str, " ", ty, " object")
            print(io, ty, " ", str)
        else
            #print(io, ty)
            print_short(io, item)
        end
    end
    return nothing
end

# Reuses show function to display array of objects
function print_datatype_array(io::IO, array::AbstractArray)
    print(io, length(array), "-element ",typeof(array), ":")
    n = length(array)
    maxn = 8
    half = div(maxn,2)
    idx = n<=maxn ? [1:n;] : [1:half; n-half+1:n]
    for i in idx
        str = "\n  "*replace(string(array[i]), "\n" => "\n  ")
        print(io, str)
        if n>maxn && i==half
            print(io, "\n    ⋮")
        end
    end
    return nothing
end


# Generates show functions for data types
macro show_function(datatype)
    return quote
        function $(esc(:show))(io::IO, obj::$datatype)
            print_datatype_fields(io, obj)
        end
    end
end

# Reuses show function to display array of objects
macro show_array_function(datatype)
    return quote
        function $(esc(:show))(io::IO, array::Array{<:$(datatype),1})
            print_datatype_array(io, array)
        end
    end
end


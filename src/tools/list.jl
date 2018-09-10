import Base.show, Base.getindex, Base.push!, Base.delete!, Base.start, Base.next, Base.done, Base.length

abstract type AbstractNode end

type Node{T} where T
    data::T
    next::Node{T}
    prev::Node{T}
    function Node()
        return new()
    end
    function Node(data::T)
        return new(data)
    end
end

function show{T}(io::IO, node::Node{T})
    print(io, "Node(data: $(node.data))")
end

function getindex(node::Node, index::Int)
    nodei = node
    if index>0
        for i=1:index
            nodei = nodei.next
        end
    end

    if index<0
        index = abs(index)
        for i=1:index
            nodei = nodei.prev
        end
    end

    return nodei
end

#function Node{T}(data::T)
    #@show 1
    #node = Node{T}()
    #@show 1
    #node.data = data
#end


type List{T}
    size  ::Int
    first::Node{T}
    last ::Node{T}
    function List()
        return new(0)
    end
end

function length(list::List)
    return list.size
end


#type Iterator{T}
    #container::T
#end

#iterator(container) = Iterator(container)

# DummyNode used as start point for iteration
type DummyNode{T}
    next::Node{T}
end

Base.start{T}(list::List{T}) = DummyNode{T}(list.first)
Base.next{T}(list::List{T}, node) = node.next, node.next
Base.done{T}(list::List{T}, node) = node == list.last


function push!{T}(list::List{T}, node::Node{T})
    if list.size==0
        list.first= node
        list.last = node
        node.prev = node
        node.next = node
    else
        node.prev = list.last
        node.next = list.first # circular
        list.first.prev = node
        list.last.next  = node
        list.last = node
    end
    list.size += 1
end

function push!{T}(list::List{T}, data::T)
    push!(list, Node{T}(data))
end

# insert(list, cnode, new)
function insert!{T}(list::List{T}, nodepos::Node{T}, newnode::Node{T})
    # inserts before nodepos
    #if nodepos==list.first
        #list.first.prev = nodepos
    #end
    @assert list.size>0

    newnode.prev = nodepos.prev
    newnode.prev.next = newnode
    newnode.next = nodepos
    nodepos.prev = newnode

    #if nodepos==list.last
        #push!(list, newnode)
        #return
    #end

    #newnode.prev = nodepos.prev
    #newnode.next = nodepos
    #newnode.prev.next = newnode
    #nodepos.prev = newnode
    list.size += 1

end

function delete!{T}(list::List{T}, node::Node{T})
    @assert list.size > 0
    node.prev.next = node.next
    node.next.prev = node.prev
    list.size -= 1
end


#=

l = List{Int}()
n = Node{Int}(1)
push!(l, n)
n = Node{Int}(2)
push!(l, n)
n = Node{Int}(3)
push!(l, n)
n = Node{Int}(4)
push!(l, n)
ii = n
n = Node{Int}(44)
insert!(l, ii, n)

delete!(l, n)

for node in l
    @show "========="
    @show node.prev.data
    @show node.data
    @show node.next.data
    #@show (node.data, node.prev, node.next)
end

node = l.first
@show node[1]
@show node[2]
@show node[3]
node = l.last
@show node[-1]
@show node[-2]
@show node[-3]

#A*B = Aik*Bkj
#A*B*C = Aik*Bkm*Cmj

=#

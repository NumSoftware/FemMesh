abstract AbstractNode


type List{T}
    len  ::Int
    first::AbstractNode
    last ::AbstractNode
    function List()
        return new(0)
    end
end

# DummyNode used as start point for iteration
type DummyNode <: AbstractNode
    next::AbstractNode
end

Base.start(list::List) = DummyNode(list.first)
Base.next(list::List, node) = node.next, node.next
Base.done(list::List, node) = node == list.last

type Node{T} <: AbstractNode
    data::T
    next::Node{T}
    prev::Node{T}
    list::List{T}
    function Node(data::T)
        return new(data)
    end
end


function push!{T}(list::List, node::Node{T})
    if !isdefined(list, :first)
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

    list.len += 1
end

# insert(list, cnode, new)
function insert!{T}(list::List, nodepos::Node{T}, newnode::Node{T})
    if nodepos==list.last
        push!(list, newnode)
        return
    end

    newnode.prev = nodepos.prev
    newnode.next = nodepos
    newnode.prev.next = newnode
    nodepos.prev = newnode
    list.len += 1

end

function delete!{T}(list::List, node::Node{T})
    node.prev.next = node.next
    node.next.prev = node.prev
    list.len -= 1
end

l = List{Int}()
n = Node{Int}(1)
push!(l, n)
n = Node{Int}(2)
push!(l, n)
n = Node{Int}(3)
push!(l, n)
n = Node{Int}(4)
push!(l, n)

for node in l
    @show "========="
    @show node.prev.data
    @show node.data
    @show node.next.data
    #@show (node.data, node.prev, node.next)
end

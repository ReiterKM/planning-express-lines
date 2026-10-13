using Graphs

ASC = 1
DESC = 2
DIRS = 1:2

EPSILON = 10^(-5)

# node = (t,h)
TIME_IDX = 1
STOP_IDX = 2

# arc = (src_node, dst_node, dir)
DIR_IDX = 3

function safe_div(a, b)
    if b == 0
        return NaN
    else
        return a / b
    end
end
function idx_up(t::Int, tΔ::Int)::Int
    return cld(t - 1, tΔ) + 1
end

function idx_down(t::Int, tΔ::Int)::Int
    return fld(t - 1, tΔ) + 1
end

function reverse_edge(e::Tuple{Int,Int})::Tuple{Int,Int}
    return (e[2], e[1])
end

function index_to_node(index::Int, n::Int)::Tuple{Int,Int}
    row, col = divrem(index - 1, n)
    return (row + 1, col + 1)
end

function node_to_index(row::Int, col::Int, n::Int)::Int
    return (row-1)*n + col
end

function index_to_stop(index::Int, n::Int)::Int
    return index_to_node(index, n)[STOP_IDX]
end

function index_to_time(index::Int, n::Int)::Int
    return index_to_node(index, n)[TIME_IDX]
end

function arc_to_index(arc::Tuple{Int,Int,Int}, Aset::Vector)::Int
    idx = findfirst(==(arc), Aset)
    @assert !isnothing(idx) "No index found for arc $arc."
    @assert length(idx) == 1 "More than one index found for arc $arc."
    return idx
end

function arc_to_index(src::Int, dst::Int, dir::Int, Aset::Vector)::Int
    arc = (src, dst, dir)
    idx = findfirst(==(arc), Aset)
    @assert length(idx) == 1 "No unique index found for arc $arc."
    return idx
end

function edge_to_time_arcs(src::Int, dst::Int, Aset::Vector, Nset::Matrix{Int})::Vector{Int}
    src_nodes = Nset[:, src]
    dst_nodes = Nset[:, dst]

    return [a_idx for (a_idx, a) in enumerate(Aset) if (a[1] in src_nodes && a[2] in dst_nodes) || (a[2] in src_nodes && a[1] in dst_nodes)]
end

function edge_to_time_arcs(src::Int, dst::Int, Aset::Vector, Nset::Matrix{Int}, dir::Int)::Vector{Int}
    src_nodes = Nset[:, src]
    dst_nodes = Nset[:, dst]

    if (src == dst == 1 || src == dst == size(Nset)[2])
        @assert (dir == ASC) "Only ASC waiting edges on nodes 1 and n."
    end

    if dir == ASC
        return [a_idx for (a_idx, a) in enumerate(Aset) if (a[1] in src_nodes && a[2] in dst_nodes) && a[3] == dir]
    else
        return [a_idx for (a_idx, a) in enumerate(Aset) if (a[2] in src_nodes && a[1] in dst_nodes) && a[3] == dir]
    end
end

function time_expanded_arc_to_edge(arc::Tuple{Int,Int,Int}, n::Int)::Tuple{Int,Int}
    return index_to_stop(arc[1], n), index_to_stop(arc[2], n)
end

function get_departure_time_from_arc(arc::Tuple{Int,Int,Int}, n::Int)::Int
    return index_to_node(arc[1], n)[TIME_IDX]
end

function get_arrival_time_from_arc(arc::Tuple{Int,Int,Int}, n::Int)::Int
    return index_to_node(arc[2], n)[TIME_IDX]
end

function get_directed_arcs(arclist::Vector{Tuple{Int,Int,Int}}, dir::Int)::Vector{Int}
    @assert dir in [1,2] "Unknown direction."
    return [idx for (idx,a) in enumerate(arclist) if a[DIR_IDX] == dir]
end

function get_outgoing_arcs_from_station(src::Int, n::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # outgoing from station, i.e., all a = (v,w): v_h = src
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[1], n) == src]
end

function get_outgoing_arcs_from_station(src::Int, n::Int, dir::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    @assert dir ∈ [ASC,DESC] "Unknown direction $dir."
    # outgoing from station, i.e., all a = (v,w): v_h = src and dir == dir
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[1], n) == src && a[DIR_IDX]== dir]
end

function get_outgoing_arcs_from_station(src::Int, n::Int, t::Int, f::Function, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # outgoing from station, i.e., all a = (v,w): v_h = src
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[1], n) == src && f(index_to_time(a[1], n), t)]
end

function get_outgoing_arcs_from_station(src::Int, n::Int, t::Vector{Int}, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # outgoing from station, i.e., all a = (v,w): v_h = src
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[1], n) == src && index_to_time(a[1], n) ≥ t[1] && index_to_time(a[1],n) ≤ t[2]]
end

# note: includes cyclic arcs
function get_incoming_arcs_to_node(ten::TimeExpandedNetwork, node::Int)::Vector{Int}
    return ten.in_arcs[node]
end

function get_incoming_arcs_to_node(ten::TimeExpandedNetwork, node::Int, dir::Int)::Vector{Int}
    return [a for a in ten.in_arcs[node] if ten.arc_dir[a] == dir]
end

function get_outgoing_arcs_from_node(ten::TimeExpandedNetwork, node::Int)::Vector{Int}
    return ten.out_arcs[node]
end

function get_outgoing_arcs_from_node(ten::TimeExpandedNetwork, node::Int, dir::Int)::Vector{Int}
    return [a for a in ten.out_arcs[node] if ten.arc_dir[a] == dir]
end

function get_outgoing_driving_arcs(ten::TimeExpandedNetwork, src::Int)::Vector{Int}
    return [a for a in ten.out_arcs[src] if ten.is_driving[a]]
end

function get_outgoing_driving_arcs(ten::TimeExpandedNetwork, src::Int, dir::Int)::Vector{Int}
    return [a for a in ten.out_arcs[src] if ten.is_driving[a] && ten.arc_dir[a] == dir]
end

function get_outgoing_waiting_arcs(ten::TimeExpandedNetwork, src::Int)::Vector{Int}
    return [a for a in ten.out_arcs[src] if ten.is_waiting[a]]
end

function get_outgoing_waiting_arcs(ten::TimeExpandedNetwork, src::Int, dir::Int)::Vector{Int}
    return [a for a in ten.out_arcs[src] if ten.is_waiting[a] && ten.arc_dir[a] == dir]
end

#=
function get_outgoing_arcs_from_node(src::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # outgoing from node, i.e., all a : a \in \delta^+(src)
    return [idx for (idx,a) in enumerate(arclist) if a[1] == src]
end

function get_outgoing_arcs_from_node(src::Int, arclist::Vector{Tuple{Int,Int,Int}}, dir::Int)::Vector{Int}
    # outgoing from node, i.e., all a : a \in \delta^+(src) in direction dir
    return [idx for (idx,a) in enumerate(arclist) if a[1] == src && a[3] == dir]
end

function get_outgoing_arcs_from_node_with_offset(src::Int, arclist::Vector{Tuple{Int,Int,Int}}, offset::Int)::Vector{Int}
    return [e + offset for e in get_outgoing_arcs_from_node(src, arclist)]
end
=#

function get_outgoing_edges(h::Int, edgelist::Vector{Edge{Int}})::Vector{Int}
    return [i for (i,e) in enumerate(edgelist) if src(e) == h]
end

function get_outgoing_edges(h::Int, edgelist::Vector{Tuple{Int,Int}})::Vector{Int}
    return [i for (i,e) in enumerate(edgelist) if e[1] == h]
end

function get_outgoing_edges(h::Int, edgelist::Vector{Edge{Int}}, dir::Int)::Vector{Int}
    if dir == ASC
        return [i for (i,e) in enumerate(edgelist) if src(e) == h && dst(e) >= h]
    elseif dir == DESC
        return [i for (i,e) in enumerate(edgelist) if src(e) == h && dst(e) <= h]
    else
        @error "Unknown direction $dir."
    end
end

function get_outgoing_edges(h::Int, edgelist::Vector{Tuple{Int,Int}}, dir::Int)::Vector{Int}
    if dir == ASC
        return [i for (i,e) in enumerate(edgelist) if e[1] == h && e[2] >= h]
    elseif dir == DESC
        return [i for (i,e) in enumerate(edgelist) if e[1] == h && e[2] <= h]
    else
        @error "Unknown direction $dir."
    end
end

function get_outgoing_edges(h::Int, max_h::Int, edgelist::Vector{Edge{Int}}, dir::Int)::Vector{Int}
    if dir == ASC
        return [i for (i,e) in enumerate(edgelist) if src(e) == h && dst(e) >= h && dst(e) <= max_h]
    elseif dir == DESC
        return [i for (i,e) in enumerate(edgelist) if src(e) == h && dst(e) <= h && dst(e) >= max_h]
    else
        @error "Unknown direction $dir."
    end
end

function get_incoming_arcs_to_station(dst::Int,  n::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # incoming to station, i.e., all a = (v,w): w_h = dst
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[2],n) == dst]
end

function get_incoming_arcs_to_station(dst::Int,  n::Int, dir::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    @assert dir ∈ [ASC,DESC] "Unknown direction $dir."
    # incoming to station, i.e., all a = (v,w): w_h = dst and dir==dir
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[2],n) == dst && a[DIR_IDX] == dir]
end

function get_incoming_arcs_to_station(dst::Int, n::Int, t::Int, f::Function, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # incoming to station, i.e., all a = (v,w): w_h = dst
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[2], n) == dst && f(index_to_time(a[2], n), t)]
end

function get_incoming_arcs_to_station(dst::Int, n::Int, t::Vector{Int}, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # incoming to station, i.e., all a = (v,w): w_h = dst
    return [idx for (idx,a) in enumerate(arclist) if index_to_stop(a[2], n) == dst && index_to_time(a[2], n) ≥ t[1] && index_to_time(a[2], n) ≤ t[2]]
end

function get_incoming_driving_arcs(ten::TimeExpandedNetwork, dst::Int)::Vector{Int}
    return [a for a in ten.in_arcs[dst] if ten.is_driving[a]]
end

function get_incoming_driving_arcs(ten::TimeExpandedNetwork, dst::Int, dir::Int)::Vector{Int}
    return [a for a in ten.in_arcs[dst] if ten.is_driving[a] && ten.arc_dir[a] == dir]
end

function get_incoming_waiting_arcs(ten::TimeExpandedNetwork, dst::Int)::Vector{Int}
    return [a for a in ten.in_arcs[dst] if ten.is_waiting[a]]
end

function get_incoming_waiting_arcs(ten::TimeExpandedNetwork, dst::Int, dir::Int)::Vector{Int}
    return [a for a in ten.in_arcs[dst] if ten.is_waiting[a] && ten.arc_dir[a] == dir]
end

#=
function get_incoming_arcs_to_node(dst::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Int}
    # incoming to node, i.e., all a: a \in \delta^-(dst)
    return [idx for (idx,a) in enumerate(arclist) if a[2] == dst]
end

function get_incoming_arcs_to_node(dst::Int, arclist::Vector{Tuple{Int,Int, Int}}, dir::Int)::Vector{Int}
    # incoming to node, i.e., all a: a \in \delta^-(dst) in direction dir
    return [idx for (idx,a) in enumerate(arclist) if a[2] == dst && a[3] == dir]
end

function get_incoming_arcs_to_node_with_offset(src::Int, arclist::Vector{Tuple{Int,Int,Int}}, offset::Int)::Vector{Int}
    return [e + offset for e in get_incoming_arcs_to_node(src, arclist)]
end
=#

function get_incoming_edges(h::Int, edgelist::Vector{Edge{Int}})::Vector{Int}
    return [i for (i,e) in enumerate(edgelist) if dst(e) == h]
end

function get_incoming_edges(h::Int, edgelist::Vector{Tuple{Int,Int}})::Vector{Int}
    return [i for (i,e) in enumerate(edgelist) if e[2] == h]
end

function get_incoming_edges(h::Int, edgelist::Vector{Edge{Int}}, dir::Int)::Vector{Int}
    @assert dir in DIRS "Unknown direction $dir."
    if dir == ASC
        return [i for (i,e) in enumerate(edgelist) if dst(e) == h && src(e) <= h]
    elseif dir == DESC
        return [i for (i,e) in enumerate(edgelist) if dst(e) == h && src(e) >= h]
    else
        @error "Unknown direction $dir."
    end
end

function get_incoming_edges(h::Int, edgelist::Vector{Tuple{Int,Int}}, dir::Int)::Vector{Int}
    @assert dir in DIRS "Unknown direction $dir."
    if dir == ASC
        return [i for (i,e) in enumerate(edgelist) if e[2] == h && e[1] <= h]
    elseif dir == DESC
        return [i for (i,e) in enumerate(edgelist) if e[2] == h && e[1] >= h]
    else
        @error "Unknown direction $dir."
    end
end

function get_incoming_edges(h::Int, max_h::Int, edgelist::Vector{Edge{Int}}, dir::Int)::Vector{Int}
    @assert dir in DIRS "Unknown direction $dir."
    if dir == ASC
        return [i for (i,e) in enumerate(edgelist) if dst(e) == h && src(e) <= h && src(e) >= max_h]
    elseif dir == DESC
        return [i for (i,e) in enumerate(edgelist) if dst(e) == h && src(e) >= h && src(e) <= max_h]
    else
        @error "Unknown direction $dir."
    end
end

function get_arcs_between_node_and_station(src_node::Int, dst_station::Int, n::Int, arclist::Vector{Tuple{Int,Int,Int}})
    return [a for a in arclist if a[1] == src_node && index_to_stop(a[2], n) == dst_station]
end

function construct_path(n::Int, initial_node::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Tuple{Int,Int,Int}}
    arcs = Vector{Tuple{Int,Int,Int}}()
    node = initial_node
    for h in 2:n
        outgoing_arcs = get_arcs_between_node_and_station(node, h, n, arclist)
        @assert length(outgoing_arcs) == 1 "More than one arc found between $(h-1) and $h: $(outgoing_arcs)."

        push!(arcs, outgoing_arcs[1])
        node = outgoing_arcs[1][2]
    end
    for h in (n-1):-1:1
        outgoing_arcs = get_arcs_between_node_and_station(node, h, n, arclist)
        @assert length(outgoing_arcs) == 1 "More than one arc found between $(h-1) and $h: $(outgoing_arcs)."

        push!(arcs, outgoing_arcs[1])
        node = outgoing_arcs[1][2]
    end

    @assert index_to_stop(node, n) == 1 "Initial path does not return to station 1."
    @assert length(arcs) == 2*(n-1) "Number of arcs does not match number of stations $n."
    return arcs
end

function extend_path(initial_path::Vector{Tuple{Int,Int,Int}}, f::Int, w::Int, n::Int, T::Int, arclist::Vector{Tuple{Int,Int,Int}})::Vector{Tuple{Int,Int,Int}}
    full_path = initial_path

    node = full_path[end][2]

    for wi in 1:w # waiting arcs
        current_time = index_to_time(node, n)
        next_node = node_to_index(current_time + 1, 1, n)
        push!(full_path, (node, next_node,1))
        node = next_node
    end

    for fi in 2:f
        f_path = construct_path(n, node, arclist)
        push!(full_path, f_path...)
        node = full_path[end][2]

        for wi in 1:w # waiting arcs
            current_time = index_to_time(node, n)
            next_time = current_time + 1 ≤ T ? current_time + 1 : 1
            next_node = node_to_index(next_time, 1, n)
            push!(full_path, (node,next_node,1))
            node = next_node
        end
    end

    @assert length(full_path) == f * (2*(n-1) + w) "Length of path does not match."

    return full_path
end

function find_frequency_and_remainder(T::Int, path_length::Int)::Tuple{Int,Int}
    # find the largest frequency that is at most floor(T/path_length)
    @assert T > 0 && path_length > 0 "Time and path length must be positive integers"
    @assert path_length <= T "requires path_length ≤ T"

    limit = div(T, path_length) # upper bound
    freq = 1

    i = 1
    while i * i <= T
        if T % i == 0
            j = div(T, i)
            if i <= limit && i > freq
                freq = i
            end
            if j <= limit && j > freq
                freq = j
            end
        end
        i += 1
    end

    rem = div(T, freq) - path_length
    return (freq, rem)
end
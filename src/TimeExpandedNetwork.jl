struct TimeExpandedNetwork
    tΔ::Int
    Tset::Vector{Int}
    Tidx::UnitRange{Int}
    num_time_steps::Int
    fmax::Int

    Nset::Matrix{Int}
    Aset::Vector{Tuple{Int,Int,Int}}
    Aidx::UnitRange{Int}

    # arc: (tail, head, dir)
    arc_head::Vector{Int}
    arc_tail::Vector{Int}
    arc_dir::Vector{Int}
    arc_len::Vector{Real}
    arc_dep::Vector{Int}
    arc_arr::Vector{Int}

    is_cyclic::BitVector
    is_driving::BitVector
    is_waiting::BitVector
    
    out_arcs::Vector{Vector{Int}}
    in_arcs::Vector{Vector{Int}}

    node_t::Vector{Int}
    node_h::Vector{Int}

    AD::Vector{Tuple{Int,Int,Int}} # driving arcs
    AW::Vector{Tuple{Int,Int,Int}} # waiting arcs
    AD_ASC::Vector{Int}
    AD_DESC::Vector{Int}
    A_CYCLIC::Vector{Int}

    ## == underlying line graph ==
    Eset::Vector{Tuple{Int,Int}}
    Eidx::UnitRange{Int}
    num_edges::Int
    arc_edge::Vector{Int}
    edge_arcs::Vector{Vector{Int}}
end

function TimeExpandedNetwork(D::InstanceData, tΔ::Int)
    Tset, Tidx, Nset, AD, AW, AD_ASC, AD_DESC, A_CYCLIC = build_time_expanded_network(D,tΔ)

    num_time_steps = length(Tidx)

    roundtrip = max(1, ceil(D.w[1, D.n] / tΔ)) + max(1, ceil(D.w[D.n,1] / tΔ))
    fmax = floor(Int, num_time_steps / roundtrip)
    @assert fmax ≥ 1 "Maximum frequency ($fmax) needs to be at least 1."

    @info "Maximum frequency: $fmax"
    @assert maximum(D.OD) <= D.Q * fmax * length(D.K) "Inputs (K: $(length(D.K)), fmax: $(D.fmax), Q: $(D.Q)) are not sufficient to serve maximum demand of $(maximum(D.OD))."

    Aset = vcat(AD,AW)
    num_arcs = length(Aset)
    Aidx = 1:num_arcs
    
    num_nodes = length(Nset)

    # arc = (tail, head, dir)
    arc_tail = [a[1] for a in Aset]
    arc_head = [a[2] for a in Aset]
    arc_dir = [a[3] for a in Aset]

    is_cyclic = falses(num_arcs)
    is_cyclic[A_CYCLIC] .= true

    is_driving = falses(num_arcs)
    is_driving[1:length(AD)] .= true

    is_waiting = falses(num_arcs)
    is_waiting[length(AD)+1:num_arcs] .= true

    out_arcs = [Vector{Int}() for _ in 1:num_nodes]
    in_arcs = [Vector{Int}() for _ in 1:num_nodes]

    for a in Aidx
        push!(in_arcs[arc_head[a]], a)
        push!(out_arcs[arc_tail[a]], a)
    end

    node_t = zeros(Int, num_nodes)
    node_h = zeros(Int, num_nodes)
    for t in Tidx, h in D.H
        v = Nset[t,h]
        node_t[v] = t
        node_h[v] = h
    end

    ascEset = [(h,g) for h in D.H for g in (h+1):D.n]
    descEset = [reverse_edge(e) for e in ascEset]
    Eset = vcat(ascEset, descEset)
    num_edges = length(Eset)
    Eidx = 1:num_edges

    edge_index = Dict(e => i for (i, e) in enumerate(Eset))

    arc_edge = zeros(Int, num_arcs) # arc => edge
    edge_arcs = [Int[] for _ in Eidx] # edge => [arcs]
    arc_len = zeros(num_arcs)
    for a in Aidx
        if is_driving[a]
            e_o, e_d = node_h[arc_tail[a]], node_h[arc_head[a]]
            e = (e_o, e_d)
            arc_edge[a] = edge_index[e]
            push!(edge_arcs[edge_index[e]], a)
            arc_len[a] = D.w[e_o, e_d]
        end
    end

    # departure times (in real time)
    arc_dep = [Tset[node_t[arc_tail[a]]] for a in Aidx]
    arc_arr = [Tset[node_t[arc_head[a]]] for a in Aidx]

    return TimeExpandedNetwork(
        tΔ, Tset, Tidx, num_time_steps, fmax, 
        Nset, Aset, Aidx, 
        arc_head, arc_tail, arc_dir, arc_len, arc_dep, arc_arr,
        is_cyclic, is_driving, is_waiting, 
        out_arcs, in_arcs, 
        node_t, node_h, 
        AD, AW, AD_ASC, AD_DESC, A_CYCLIC, 
        Eset, Eidx, num_edges, arc_edge, edge_arcs)
end

function validate_network(ten::TimeExpandedNetwork, D::InstanceData)
    fwd_in_time = a -> ten.node_t[ten.arc_head[a]] > ten.node_t[ten.arc_tail[a]]
    same_station = a -> ten.node_h[ten.arc_head[a]] == ten.node_h[ten.arc_tail[a]]

    @assert all(fwd_in_time(a) for a in ten.Aidx if !ten.is_cyclic[a])
    @assert all(same_station(a) for a in ten.Aidx if ten.is_waiting[a])
    @assert all(!same_station(a) for a in ten.Aidx if ten.is_driving[a])

    # model invariant: end stations only have ASC waiting arcs
    @assert all(ten.arc_dir[a] == ASC for a in ten.Aidx if ten.is_waiting[a] && ten.node_h[ten.arc_tail[a]] in (1, D.n))

    @assert D.fmin ≤ length(D.K) * ten.fmax "Minimum frequency $(D.fmin) exceeds the maximum achievable frequency $(length(D.K) * ten.fmax)."
end

function build_time_expanded_network(D::InstanceData, tΔ::Int)::Tuple{Vector{Int}, UnitRange{Int}, Matrix{Int}, Vector{Tuple{Int,Int,Int}}, Vector{Tuple{Int,Int,Int}}, Vector{Int}, Vector{Int}, Vector{Int}}
    Tset = collect(1:tΔ:D.T)
    num_time_steps = length(Tset)
    Tidx = 1:num_time_steps
    Nset = Matrix{Int}([(t_idx-1)*D.n + h for t_idx in Tidx, h in D.H])

    # arcs are of the form (start_node, end_node, DIR)
    # where a node is an index in Nset
    AD = Vector{Tuple{Int,Int,Int}}() # driving
    AW = Vector{Tuple{Int,Int,Int}}() # waiting

    AD_ASC = Vector{Int}()
    AD_DESC = Vector{Int}()
    A_CYCLIC = Vector{Int}()

    for src in Nset
        t, h_o = index_to_node(src, D.n)
        for h_d in D.H
            if h_o != h_d
                t_steps = max(1, ceil(Int, D.w[h_o,h_d] / tΔ)) # travel time in time steps
                DIR = h_o < h_d ? ASC : DESC
                if t + t_steps ≤ num_time_steps
                    dst = node_to_index(t + t_steps, h_d, D.n)
                    push!(AD, (src, dst, DIR))
                else
                    # t = (t + travel_time) MOD (D.T / tΔ)
                    t_cyclic = mod1(t + t_steps, num_time_steps)
                    dst = node_to_index(t_cyclic, h_d, D.n)
                    push!(AD, (src, dst, DIR))
                    push!(A_CYCLIC, length(AD))
                end

                push!(DIR == ASC ? AD_ASC : AD_DESC, length(AD))
            end
        end
    end

    offset = length(AD)
    
    for h in D.H[2:end-1]
        for t in Tidx[1:end-1]
            push!(AW, (Nset[t,h], Nset[t+1,h], ASC))
            push!(AW, (Nset[t,h], Nset[t+1,h], DESC))
        end
        push!(AW, (Nset[Tidx[end],h], Nset[Tidx[1], h], ASC))
        push!(A_CYCLIC, offset + length(AW))
        push!(AW, (Nset[Tidx[end],h], Nset[Tidx[1], h], DESC))
        push!(A_CYCLIC, offset + length(AW))
    end

    for h in [1, D.H[end]]
        for t in Tidx[1:end-1]
            push!(AW, (Nset[t,h], Nset[t+1,h], ASC))
        end
        push!(AW, (Nset[Tidx[end],h], Nset[Tidx[1], h], ASC))
        push!(A_CYCLIC, offset + length(AW))
    end
    return Tset, Tidx, Nset, AD, AW, AD_ASC, AD_DESC, A_CYCLIC
end

function reach_nodes(starts::Vector{Int}, adj::Vector{Vector{Int}}, next_node::Vector{Int}, keep_arc, keep_node)
    # Depth-search:: find all reachable nodes, only using "allowed" arcs, from a list of starts nodes
    seen = falses(length(adj))
    stack = collect(starts)
    foreach(s -> seen[s] = true, stack)
    while !isempty(stack)
        v = pop!(stack)
        for a in adj[v]
            keep_arc(a) || continue
            w = next_node[a]
            if !seen[w] && keep_node(w)
                seen[w] = true
                push!(stack, w)
            end
        end
    end
    return seen
end

function get_arcs_per_request(ten::TimeExpandedNetwork, ω::ScenarioData, r::Int)::Vector{Int}
    o_r, d_r = ω.origins[r], ω.destinations[r]
    consider_arc = a -> ten.arc_dir[a] == ω.R_DIRS[r] && !ten.is_cyclic[a] && !(ten.is_waiting[a] && (ten.node_h[ten.arc_tail[a]] == o_r || ten.node_h[ten.arc_head[a]] == d_r))

    e_or = max(1, idx_up(ω.origin_tws[1,r], ten.tΔ))
    l_or = min(ten.num_time_steps, idx_down(ω.origin_tws[2,r], ten.tΔ))
    @assert e_or <= l_or "Origin time window ($e_or, $l_or) of request $r is infeasible."

    e_dr = max(1, idx_up(ω.destination_tws[1,r], ten.tΔ))
    l_dr = min(ten.num_time_steps, idx_down(ω.destination_tws[2,r], ten.tΔ))
    @assert e_dr <= l_dr "Destination time window ($e_dr, $l_dr) of request $r is infeasible."

    origin_nodes = [ten.Nset[t, o_r] for t in e_or:l_or]
    destination_nodes = [ten.Nset[t, d_r] for t in e_dr:l_dr]

    fwd_arcs = reach_nodes(origin_nodes, ten.out_arcs, ten.arc_head, consider_arc, v -> ten.node_t[v] <= ω.destination_tws[2,r]) # l_dr
    bwd_arcs = reach_nodes(destination_nodes, ten.in_arcs, ten.arc_tail, consider_arc, v -> ten.node_t[v] >= ω.origin_tws[1,r]) # e_or

    return [a for a in ten.Aidx if consider_arc(a) && fwd_arcs[ten.arc_tail[a]] && bwd_arcs[ten.arc_head[a]]]
end

function get_arcs_for_all_requests(ten::TimeExpandedNetwork, ω::ScenarioData)::Vector{Vector{Int}}
    return [get_arcs_per_request(ten, ω, r) for r in ω.R]
end

function get_all_requests_per_arc(ten::TimeExpandedNetwork, request_to_arcs::Vector{Vector{Int}})
    arc_to_requests = [Int[] for _ in ten.Aidx]

    for (r, arcs) in enumerate(request_to_arcs), a in arcs
        push!(arc_to_requests[a], r)
    end

    return arc_to_requests
end

function validate_request_arcs(ten::TimeExpandedNetwork, ω::ScenarioData, arcs_per_request::Vector{Vector{Int}})
    for (r, arcs) in enumerate(arcs_per_request)
        @assert all(ten.arc_dir[a] == ω.R_DIRS[r] for a in arcs)
        @assert all(!ten.is_cyclic[a] for a in arcs)
        @assert !any(ten.is_waiting[a] && ten.node_h[ten.arc_tail[a]] in (ω.origins[r], ω.destinations[r]) for a in arcs)
        @assert all(ten.is_driving[a] for a in arcs if ten.node_h[ten.arc_tail[a]] == ω.origins[r])
        @assert all(ten.is_driving[a] for a in arcs if ten.node_h[ten.arc_head[a]] == ω.destinations[r])
    end
end
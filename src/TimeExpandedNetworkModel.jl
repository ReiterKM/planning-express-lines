using JuMP
using Gurobi

function add_ten_objective_function(model::Model, D::InstanceData, Ω::Scenarios, ten::TimeExpandedNetwork, which_obj_func::String, weights::Vector{Float64})
    if which_obj_func ∈ ["dist-pax", "dist-pax-wt"]
        model[:first_stage_objective] = weights[1] * sum(ten.arc_len[a_idx] * model[:x][k,a_idx] for a_idx in eachindex(ten.AD) for k in D.K)
    else
        model[:first_stage_objective] = 0
    end

    if which_obj_func == "dist-pax"
        model[:second_stage_objective] = weights[2] * sum(Ω.probabilities[ω.idx] * (sum(1 - model[:z][ω.idx][r] for r in ω.R; init=0)) for ω in Ω.scenarios; init = 0)
    elseif which_obj_func == "dist-pax-wt"
        model[:second_stage_objective] = sum(Ω.probabilities[ω.idx] * (
            weights[2] * sum(1 - model[:z][ω.idx][r] for r in ω.R; init=0) +
            weights[3] * sum(model[:c][ω.idx][r] for r in ω.R; init=0)
            ) for ω in Ω.scenarios; init = 0)
    elseif which_obj_func == "pax-wt"
        model[:second_stage_objective] = sum(Ω.probabilities[ω.idx] * (
            weights[1] * sum(1 - model[:z][ω.idx][r] for r in ω.R; init=0) +
            weights[2] * sum(model[:c][ω.idx][r] for r in ω.R; init=0)
            ) for ω in Ω.scenarios; init = 0)
    elseif which_obj_func == "pax"
        model[:second_stage_objective] = sum(Ω.probabilities[ω.idx] * (
            sum(1 - model[:z][ω.idx][r] for r in ω.R; init=0)
            ) for ω in Ω.scenarios; init=0)
    else
        model[:second_stage_objective] = 0
    end

    @objective(model, Min, model[:first_stage_objective] + model[:second_stage_objective])
end

function build_time_expanded_network_model(D::InstanceData, Ω::Scenarios, tΔ::Int, directionality_property::Bool, relax_symmetry::Bool, warm_start::Bool, relax_second_stage::Bool, only_second_stage::Bool, drive_all_stops::Bool, verbose::Bool, which_obj_func::String, weights::Vector{Float64}, solution_arcs::Dict{Int, Vector{Any}})::Tuple{Model, TimeExpandedNetwork}
    model = Model(Gurobi.Optimizer)
    if verbose
        set_string_names_on_creation(model, true)
    else
        set_string_names_on_creation(model, false)
    end

    # == helpers == #
    ten = TimeExpandedNetwork(D, tΔ)
    validate_network(ten, D)

    odPairs = [(o,d) for o in D.H for d in D.H if o < d]

    function _build_first_stage()
        # == variables ==
        @variable(model, x[D.K, ten.Aidx], Bin, base_name = "x")
        if !relax_symmetry
            @variable(model, xbar[D.K,ten.Eidx], Bin, base_name = "xbar")
        end
        @variable(model, p[ten.Eidx,D.H,D.H] ≥ 0, base_name = "p")
        @variable(model, ten.fmax ≥ fod[D.K,1:length(odPairs)] ≥ 0, base_name = "f_od")
        @variable(model, f[ten.Eidx] ≥ 0, base_name = "f")

        if warm_start
            foreach(v -> set_start_value(v,0), x)
            if !relax_symmetry
                foreach(v -> set_start_value(v,1), xbar)
            end

            initial_path = construct_path(D.n, node_to_index(1,1,D.n), ten.AD)
            path_length = index_to_time(initial_path[end][2], D.n) - 1

            path_frequency, wait = find_frequency_and_remainder(D.T, path_length)
            full_path = extend_path(initial_path, path_frequency, wait, D.n, D.T, ten.Aset)
            idxs = [arc_to_index(a, ten.Aset) for a in full_path]
            @assert all(!isnothing, idxs) "Path constructed in warm start arcs are unknown."
            for k in D.K, a in idxs
                set_start_value(x[k,a], 1)
            end
        end

        if directionality_property
            for e_idx in ten.Eidx
                src, dst = ten.Eset[e_idx]
                e_lo, e_hi = minmax(src, dst)

                for o_h in D.H, d_h in D.H
                    if o_h != d_h
                        od_lo, od_hi = minmax(o_h, d_h)
                        is_outside_path = e_lo < od_lo || e_hi > od_hi

                        if is_outside_path
                            fix(p[e_idx, o_h, d_h], 0; force=true)
                        end
                    end
                end
            end
        end

        if drive_all_stops
            for h in 1:(D.n-1)
                e_idx = findfirst(ten.Eset, (h,h+1))
                for k in D.K
                    fix(xbar[k,e_idx], 1; force=true)
                end
            end
        end

        # == constraints ==
        if !relax_symmetry
            for h in D.H
                outgoing_edges = get_outgoing_edges(h, ten.Eset, ASC)
                incoming_edges = get_incoming_edges(h, ten.Eset, ASC)
                for k in D.K
                    if h == 1 
                        @constraint(model, sum(xbar[k,e] for e in outgoing_edges) == 1, base_name="1:enable_outgoing_ASC_edge_$(h)_$k") # (2a)
                    elseif h == D.n 
                        @constraint(model, sum(xbar[k,e] for e in incoming_edges) == 1, base_name="1:enable_incoming_ASC_edge_$(h)_$k") # (2a)
                    else 
                        @constraint(model, sum(xbar[k,e] for e in outgoing_edges) - sum(xbar[k,e] for e in incoming_edges) == 0, base_name="1:flow_ASC_edge_$(h)_$k") # (2b)
                    end
                end
            end
            for k in D.K
                num_dir_edges = div(ten.num_edges,2)
                for e in 1:num_dir_edges 
                    @constraint(model, xbar[k,e] == xbar[k,e + num_dir_edges], base_name="1:symmetric_path_$(k)_$(e)") # (2c)
                end
            end
        end

        for h in D.H
            in_edges = get_incoming_edges(h, ten.Eset)
            out_edges = get_outgoing_edges(h, ten.Eset)

            for o in D.H, d in D.H
                if D.OD[o,d] > 0
                    if h == o
                        C = -D.OD[o,d]
                    elseif h == d
                        C = D.OD[o,d]
                    else
                        C = 0.0
                    end
                    @constraint(model, sum(p[e,o,d] for e in in_edges) - sum(p[e,o,d] for e in out_edges) == C, base_name = "1:od_flow_($h,$o,$d)") # (2d)
                    for e in ten.Eidx
                        set_upper_bound(p[e,o,d], D.OD[o,d])
                    end
                end
            end
        end
        for node in ten.Nset
            h = ten.node_h[node]

            # arc lists for vehicle flows (independent of vehicle k)
            if h == 1
                in_desc_arcs = [get_incoming_driving_arcs(ten, node, DESC); get_incoming_waiting_arcs(ten, node)]
                out_asc_arcs = [get_outgoing_driving_arcs(ten, node, ASC); get_outgoing_waiting_arcs(ten, node)]
                flow_pairs = [(in_desc_arcs, out_asc_arcs, 0)]
            elseif h == D.n
                in_asc_arcs = [get_incoming_driving_arcs(ten, node, ASC); get_incoming_waiting_arcs(ten, node)]
                out_desc_arcs = [get_outgoing_driving_arcs(ten, node, DESC); get_outgoing_waiting_arcs(ten, node)]
                flow_pairs = [(in_asc_arcs, out_desc_arcs, 0)]
            else
                flow_pairs = [(get_incoming_arcs_to_node(ten, node, DIR), get_outgoing_arcs_from_node(ten, node, DIR), DIR) for DIR in DIRS]
            end

            for k in D.K
                for (ins, outs, dir) in flow_pairs 
                    @constraint(model, sum(x[k,a] for a in ins) == sum(x[k,a] for a in outs), base_name = "1:flow_constraints_$(node)_$(k)_$(dir)") # (3b,3c,3d)
                end
            end
        end

        for e in eachindex(ten.Eset)
            @constraint(model, sum(p[e,o,d] for o in D.H for d in D.H) ≤ D.Q * f[e], base_name="1:enable_arc_capacity_$(e)") # (2e)
            @constraint(model, sum(x[k,a] for k in D.K for a in ten.edge_arcs[e]) == f[e], base_name="1:set_edge_freq_$(e)") # (3h)        
            if !relax_symmetry
                for k in D.K
                    for a in ten.edge_arcs[e]
                        @constraint(model, x[k,a] ≤ xbar[k,e], base_name="1:enable_arc_$(a)_$k") # (3a)
                    end
                end
            end
        end

        out_arcs_per_t = [reduce(vcat, [get_outgoing_arcs_from_node(ten, n) for n in ten.Nset[t,:]]) for t in ten.Tidx]

        for k in D.K
            # strengthening
            @constraint(model, sum(x[k,a] for a in ten.A_CYCLIC) ≤ 1, base_name="1:one_cyclic_arc_$(k)") # (3f)
            for t in ten.Tidx # strengthening
                @constraint(model, sum(x[k,a] for a in out_arcs_per_t[t]) ≤ 1, base_name="1:one_outgoing_arc_per_time_$(t)_$(k)") # (3e)
            end

            for k_bar in D.K 
                if k_bar > k
                    @constraint(model, sum(ten.arc_len[a_idx] * (x[k,a_idx] - x[k_bar,a_idx]) for (a_idx,a) in enumerate(ten.AD)) ≥ 0, base_name="1:symmetry_breaking_$(k)_$(k_bar)") # (3g)
                end
            end
        end

        out_driving_arcs_from_1 = get_outgoing_arcs_from_station(1, D.n, ASC, ten.AD)
        for (od_idx, (o,d)) in enumerate(odPairs)
            outgoing_asc_arcs = get_outgoing_arcs_from_station(o, D.n, ASC, ten.AD)
            incoming_asc_arcs = get_incoming_arcs_to_station(d, D.n, ASC, ten.AD)

            for k in D.K
                @constraint(model, sum(x[k,a] for a in outgoing_asc_arcs) ≥ fod[k,od_idx], base_name="2:set_fod_outgoing_$(od_idx)_$(k)") # (3i)
                @constraint(model, sum(x[k,a] for a in incoming_asc_arcs) ≥ fod[k,od_idx], base_name="2:set_fod_incoming_$(od_idx)_$(k)") # (3j)
                @constraint(model, sum(x[k,a] for a in out_driving_arcs_from_1) ≥ fod[k,od_idx], base_name="2:set_fod_max_$(od_idx)_$(k)") # (3k)
            end

            @constraint(model, sum(fod[k,od_idx] for k in D.K) ≥ D.fmin, base_name="2:set_fmin_($(o),$(d))") # (3l)
        end
    end

    function _build_second_stage()
        # == helpers ==
        # vector of feasible arcs for each request (respecting time window!)
        arc_sets = []
        for ω in Ω.scenarios
            arcs_per_requests = get_arcs_for_all_requests(ten, ω)
            validate_request_arcs(ten, ω, arcs_per_requests)
            push!(arc_sets, arcs_per_requests)
        end

        # == variables ==
        y = Any[]; z = Any[]; c = Any[]; Bdep = Any[]; Barr = Any[];

        for ω in Ω.scenarios
            push!(c, @variable(model, [ω.R], lower_bound = 0, upper_bound = D.max_waiting_time, base_name = "c_$(ω.idx)"))
            push!(Bdep, @variable(model, [ω.R], lower_bound = 0, base_name = "B_dep_$(ω.idx)"))
            push!(Barr, @variable(model, [ω.R], lower_bound = 0, base_name = "B_arr_$(ω.idx)"))
            if relax_second_stage
                push!(z, @variable(model, [ω.R], lower_bound = 0, upper_bound = 1, base_name = "z_$(ω.idx)"))
                push!(y, [@variable(model, [D.K, arc_sets[ω.idx][r]], lower_bound = 0, upper_bound = 1, base_name = "y_$(ω.idx)_$r") for r in ω.R])
            else
                push!(z, @variable(model, [ω.R], Bin, base_name = "z_$(ω.idx)"))
                push!(y, [@variable(model, [D.K, arc_sets[ω.idx][r]], Bin, base_name = "y_$(ω.idx)_$r") for r in ω.R])
            end

            if warm_start
                foreach(v -> set_start_value(v, 0), c[ω.idx])
                foreach(v -> set_start_value(v, 0), z[ω.idx])
                for yr in y[ω.idx]
                    foreach(v -> set_start_value(v, 0), yr)
                end
            end
        end

        model[:y] = y; model[:z] = z; model[:c] = c; model[:Bdep] = Bdep; model[:Barr] = Barr;

        # == constraints ==
        for ω in Ω.scenarios
            allowed_requests = get_all_requests_per_arc(ten, arc_sets[ω.idx])

            # capacity
            for a in findall(ten.is_driving)
                if !isempty(allowed_requests[a])
                    for k in D.K
                        @constraint(model, sum(ω.loads[r] * y[ω.idx][r][k,a] for r in allowed_requests[a]) ≤ D.Q * model[:x][k,a], base_name="2:capacity_$(ω.idx)_$(k)_$a") # (4d)
                    end
                end
            end

            for r in ω.R
                A_r = arc_sets[ω.idx][r]
                if isempty(A_r)
                    @warn "Request $r in scenario $(ω.idx) has no feasible arcs, fixing z = 0."
                    fix(z[ω.idx][r], 0; force = true)
                    fix(c[ω.idx][r], 0; force = true)
                    continue
                end

                yr = y[ω.idx][r]
                o_r, d_r = ω.origins[r], ω.destinations[r]
                out_origin = [a for a in A_r if ten.node_h[ten.arc_tail[a]] == o_r]  # pick-up arcs
                into_destination = [a for a in A_r if ten.node_h[ten.arc_head[a]] == d_r]  # drop-off arcs

                # pick-up only once
                @constraint(model, sum(yr[k,a] for k in D.K, a in out_origin) == z[ω.idx][r], base_name="2:pickup_once_$(ω.idx)_$r") # (4a)
                # drop-off only once, link z-variable
                @constraint(model, sum(yr[k,a] for k in D.K, a in into_destination) == z[ω.idx][r], base_name="2:dropoff_once_$(ω.idx)_$r") # (4b)

                # waiting time
                @constraint(model, Bdep[ω.idx][r] == sum(ten.arc_dep[a] * yr[k,a] for k in D.K, a in out_origin), base_name="2:dep_time_$(ω.idx)_$r") # (4e)
                @constraint(model, Barr[ω.idx][r] == sum(ten.arc_arr[a] * yr[k,a] for k in D.K, a in into_destination), base_name="2:arr_time_$(ω.idx)_$r") # (4f)
                # waiting time dependent on departure / arrival request (incoming / outgoing)
                if ω.is_departure_request[r] # given departure time
                    @constraint(model, c[ω.idx][r] ≥ Bdep[ω.idx][r] - ω.requested_times[r] * z[ω.idx][r], base_name="2:waiting_time_$(ω.idx)_$r") # (4g)
                else # given arrival time
                    @constraint(model, c[ω.idx][r] ≥ ω.requested_times[r] * z[ω.idx][r] - Barr[ω.idx][r], base_name="2:waiting_time_$(ω.idx)_$r") # (4h) 
                end              
                @constraint(model, ω.max_ride_times[r] ≥ Barr[ω.idx][r] - Bdep[ω.idx][r], base_name="2:max_ride_time_$(ω.idx)_$r") # (4i)
                
                # flow conservation
                in_at = Dict{Int, Vector{Int}}()
                out_at = Dict{Int, Vector{Int}}()
                for a in A_r
                    push!(get!(out_at, ten.arc_tail[a], Int[]), a)
                    push!(get!(in_at,  ten.arc_head[a], Int[]), a)
                end
                for node in sort!(collect(union(keys(in_at), keys(out_at))))
                    if ten.node_h[node] ∉ (o_r, d_r)
                        ins  = get(in_at,  node, Int[])
                        outs = get(out_at, node, Int[])
                        for k in D.K
                            @constraint(model, sum(yr[k,a] for a in outs) == sum(yr[k,a] for a in ins), base_name = "2:flow_$(ω.idx)_$(r)_$(node)_$k") # (4c)
                        end
                    end
                end
            end
        end
    end
    
    if !only_second_stage
        _build_first_stage()
    else
        model[:x] = zeros(Int, (length(D.K), length(ten.Aidx)))
        for k in D.K
            for a_idx in solution_arcs[k]
                model[:x][k,a_idx] = 1
            end
        end
    end
    _build_second_stage()

    model[:total_distance] = sum(ten.arc_len[a_idx] * model[:x][k,a_idx] for a_idx in eachindex(ten.AD) for k in D.K)

    add_ten_objective_function(model, D, Ω, ten, which_obj_func, weights)

    return model, ten
end
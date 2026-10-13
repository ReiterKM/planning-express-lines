using JuMP
using Gurobi
using Graphs

function build_and_solve(D::InstanceData, Ω::Scenarios, time_limit::Int, which_obj_func::String, weights::Vector{Float64}, verbose::Bool, tΔ::Int, directionality_property::Bool, relax_symmetry::Bool, warm_start::Bool, relax_second_stage::Bool, drive_all_stops::Bool)::Dict{String,Any}
    # === Build ===
    @info "Building model."

    start_build_time = time()
    model, ten = build_time_expanded_network_model(D,Ω,tΔ,directionality_property,relax_symmetry,warm_start,relax_second_stage,false,drive_all_stops,verbose,which_obj_func,weights,Dict{Int, Vector{Any}}())
    end_build_time = time()
    @info "Built model in $(round(end_build_time - start_build_time; digits=3)) s."

    # === Solve ===
    @info "Solving model."
    remaining_time = time_limit - (end_build_time - start_build_time)
    if remaining_time < 0
        return Dict{String, Any}("model_build_time" => (end_build_time - start_build_time), "time_limit" => time_limit)
    end

    set_time_limit_sec(model, remaining_time)
    set_optimizer_attribute(model, "LogFile", D.output_directory * D.instance_name * D.comment * "_log.txt")
    if !verbose
        set_optimizer_attribute(model, "LogToConsole", 0)
        # write_to_file(model, D.output_directory * D.instance_name * "_model.mps")
    end

    optimize!(model)
    end_solve_time = time()
    @info "Model solved in $(round(end_solve_time - end_build_time; digits=3)) s."

    has_solution = verify_solution(D, Ω, ten, model, false, relax_second_stage)

    solution_vals = get_solution_values(D, Ω, ten, model, which_obj_func, weights, false, relax_second_stage)
    solution_vals["has_solution"] = has_solution
    solution_vals["model_build_time"] = end_build_time - start_build_time
    solution_vals["model_solve_time"] = end_solve_time - end_build_time
    solution_vals["total_time"] = end_solve_time - start_build_time
    solution_vals["time_limit"] = time_limit
    solution_vals["directionality_property"] = directionality_property
    solution_vals["warm_start"] = warm_start

    return solution_vals
end

function get_all_od_paths(D::InstanceData,model::Model, index_to_edge::Dict)::Dict
    passenger_paths = Dict{Tuple{Int,Int}, Dict{Tuple{Int,Int}, Int}}()
    for o in D.H
        for d in D.H
            if D.OD[o,d] > 0
                passenger_paths[(o,d)] = get_paths_per_od_pair(model, index_to_edge, o, d)
            end
        end
    end
    return passenger_paths
end

function get_paths_per_od_pair(model::Model, index_to_edge::Dict, o::Int, d::Int)::Dict{Tuple{Int,Int},Int}
    p = value.(model[:p][:,o,d])
    p_path = Dict{Tuple{Int,Int},Int}((src(e),dst(e)) => p[e_idx] for (e_idx,e) in index_to_edge if p[e_idx] > 0)
    return p_path
end

function directed_edges(E::Vector{Graphs.SimpleGraphs.SimpleEdge{Int}}, dir::Int)::Vector{Graphs.SimpleGraphs.SimpleEdge{Int}}
    if dir == ASC
        return [e for e in E if src(e) < dst(e)]
    else
        return [e for e in E if src(e) > dst(e)]
    end
end

function verify_solution(D::InstanceData, Ω::Scenarios, ten::TimeExpandedNetwork, model::Model, is_evaluation::Bool, relax_second_stage::Bool)::Bool
    if !(is_solved_and_feasible(model) || (termination_status(model) == TIME_LIMIT && result_count(model) > 0))
        return false
    end
    verify_time_expanded_network_solution(D, Ω, ten, model, is_evaluation, relax_second_stage)
    return true
end

function verify_time_expanded_network_solution(D::InstanceData, Ω::Scenarios, ten::TimeExpandedNetwork, model::Model, is_evaluation::Bool, relax_second_stage::Bool) 
    for ω in Ω.scenarios
        accepted_requests = Set{Int}()
        allowed_arcs = get_arcs_for_all_requests(ten, ω)
        requests_of_arc = get_all_requests_per_arc(ten, allowed_arcs)

        if !is_evaluation
            for k in D.K, a in ten.Aidx
                if ten.is_driving[a] && value(model[:x][k,a]) ≥ 0.9
                    requests_on_arc = [r for r in requests_of_arc[a] if value.(model[:y][ω.idx][r][k,a]) ≥ 0.9]
                    num_pax = sum(ω.loads[r] for r in requests_on_arc; init=0)
                    @assert num_pax ≤ D.Q "Number of passengers $num_pax on arc $a of vehicle $k in scenario $(ω.idx) exceeds capacity."
                end
            end
        end
        @info "== SCENARIO $(ω.idx) with $(ω.m) requests =="

        if !relax_second_stage
            z = model[:z][ω.idx]
            for r in ω.R
                yr = model[:y][ω.idx][r]
                if sum(value.(yr)) ≥ 0.9
                    push!(accepted_requests, r)

                    vehicles = [k for k in D.K if sum(value(yr[k,a]) for a in allowed_arcs[r]) ≥ 0.9]
                    @assert length(vehicles) == 1 "Request $r can only be assigned to one vehicle."
                    k = vehicles[1]

                    @assert value(z[r]) ≥ 0.9 "Request $r of scenario $(ω.idx) accepted but value of z $(value(z[r])) does not match."

                    used_arcs = [a for a in allowed_arcs[r] if value(yr[k,a]) ≥ 0.9]
                    @assert all(value(model[:x][k,a]) ≥ 0.9 for a in used_arcs if ten.is_driving[a]) "Request $r uses a driving arc that vehicle $k does not drive."
                    @assert length(used_arcs) ≥ 1 "Request must be assigned to at least one arc."

                    pickup_arc  = [a for a in used_arcs if ten.node_h[ten.arc_tail[a]] == ω.origins[r]]
                    dropoff_arc = [a for a in used_arcs if ten.node_h[ten.arc_head[a]] == ω.destinations[r]]
                    @assert length(pickup_arc) == 1 && length(dropoff_arc) == 1 "Request $r needs exactly one pick-up and one drop-off arc."

                    dep_time = ten.arc_dep[pickup_arc[1]]
                    arr_time = ten.arc_arr[dropoff_arc[1]]
                    @assert arr_time - dep_time ≤ ω.max_ride_times[r] "Ride time of request $r exceeds its maximum ride time."

                    if ω.is_departure_request[r]
                        waiting_time = dep_time - ω.requested_times[r]
                    else
                        waiting_time = ω.requested_times[r] - arr_time
                    end
                    @assert value(model[:c][ω.idx][r]) ≥ waiting_time - EPSILON "Request $r: c = $(value(model[:c][ω.idx][r])) is smaller than the deviation $(waiting_time)."
                    @assert waiting_time ≤ D.max_waiting_time + EPSILON "Request $r: deviation $(waiting_time) exceeds the maximum waiting time."

                    @info "($(ω.R_DIRS[r])) request $r picked-up by vehicle $k at $(dep_time) ($(ω.origins[r]) --> $(ω.destinations[r])), arrival at $(arr_time)."
                end
            end
        
            @assert length(accepted_requests) == round(Int, sum(value.(z))) "Number of transported requests does not match: accepted = $(length(accepted_requests)), z = $(sum(value.(z)))."
        end
    end
end

function get_solution_values(D::InstanceData, Ω::Scenarios, ten::TimeExpandedNetwork, model::Model, which_obj_func::String, weights::Vector{Float64}, is_evaluation::Bool, relax_second_stage::Bool; first_stage_solution::Dict=Dict())::Dict
    solution_vals = Dict{String, Any}()

    solution_vals["status"] = termination_status(model)
    solution_vals["num_vars"] = num_variables(model)
    solution_vals["num_constr"] = sum(num_constraints(model, F, K) for (F,K) in list_of_constraint_types(model))
    solution_vals["which_model"] = "Time-Expanded Network Model"
    solution_vals["which_obj_func"] = which_obj_func
    solution_vals["weights"] = weights

    @info "Termination status: $(termination_status(model))"

    if is_solved_and_feasible(model) || (termination_status(model) == TIME_LIMIT && result_count(model) > 0)
        @info "Obj. Value: $(objective_value(model))"
        @info "First stage: $(value(model[:first_stage_objective]))"
        @info "Second stage: $(value(model[:second_stage_objective]))"

        model_solution_vals = get_time_expanded_network_solution(D, Ω, ten, model, is_evaluation, relax_second_stage; first_stage_solution)
        solution_vals["num_TEN_nodes"] = length(ten.Nset)
        solution_vals["num_TEN_arcs"] = length(ten.Aidx)
        merge!(solution_vals, model_solution_vals)
    else
        solution_vals["mip_gap"] = missing
        solution_vals["obj_val"] = missing
        solution_vals["first_stage"] = missing
        solution_vals["second_stage"] = missing
        solution_vals["frequencies"] = missing
        solution_vals["B"] = missing
        solution_vals["B_T"] = missing
        solution_vals["od_freq"] = missing
        solution_vals["vehicle_flows"] = missing
        solution_vals["passenger_flows"] = missing
        solution_vals["accepted_req_per_scenario"] = missing
        solution_vals["rejected_re_per_scenario"] = missing
        solution_vals["waiting_time_per_scenario"] = missing
        solution_vals["avg_waiting_time_per_scenario"] = missing


        diagnose_first_stage(D, ten)

        if termination_status(model) ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            @info "Computing iis..." 
            output_iis(D, model)
        else
            @assert termination_status(model) == MOI.TIME_LIMIT
        end
    end
    return solution_vals
end

function get_time_expanded_network_solution(D::InstanceData, Ω::Scenarios, ten::TimeExpandedNetwork, model::Model, is_evaluation::Bool, relax_second_stage::Bool; first_stage_solution::Dict=Dict())::Dict
    solution_vals = Dict{String, Any}()
    solution_vals["mip_gap"] = relative_gap(model)
    solution_vals["obj_val"] = objective_value(model)
    solution_vals["first_stage"] = value(model[:first_stage_objective])
    solution_vals["second_stage"] = value(model[:second_stage_objective])

    paths = Dict{Int, Vector{Int}}()
    arcs = Dict{Int, Vector{Any}}()
    arcs_by_idx = Dict{Int, Vector{Any}}()
    frequencies = Vector{Int}()

    # Statistics
    total_distance = 0.0

    # first stage
    if is_evaluation
        paths = first_stage_solution["vehicle_flows"]
        arcs = first_stage_solution["vehicle_flows_TEN"]
        arcs_by_idx = first_stage_solution["vehicle_flows_TEN_idx"]
        frequencies = first_stage_solution["frequencies"]
        od_frequencies = first_stage_solution["od_freq"]
        total_distance = first_stage_solution["total_distance"]
        solution_vals["total_distance"] = total_distance
    else
        x = model[:x]
        od_frequencies = zeros(Int, D.n, D.n)

        for k in D.K
            # arcs in time order = walk order
            ids = sort([a for a in ten.Aidx if value(x[k,a]) > 0.9], by = a -> ten.arc_tail[a])
            for (i, a) in enumerate(ids)
                nxt = ids[mod1(i + 1, length(ids))]
                @assert ten.arc_head[a] == ten.arc_tail[nxt] "Vehicle $k: arcs $a and $nxt are not connected (route not closed)."
            end

            arcs_by_idx[k] = ids
            arcs[k] = [((ten.node_t[ten.arc_tail[a]], ten.node_h[ten.arc_tail[a]]),
                        (ten.node_t[ten.arc_head[a]], ten.node_h[ten.arc_head[a]])) for a in ids]
            total_distance += sum(ten.arc_len[a] for a in ids; init=0)

            # station path up to the first arrival at the last station (repeats = waiting)
            if length(ids) > 0
                stops = [ten.node_h[ten.arc_tail[ids[1]]]; [ten.node_h[ten.arc_head[a]] for a in ids]]
                cut = something(findfirst(==(D.n), stops), length(stops))
                paths[k] = stops[1:cut]
            else
                paths[k] = []
            end

            # each round trip starts with exactly one driving arc leaving station 1
            f = count(a -> ten.is_driving[a] && ten.node_h[ten.arc_tail[a]] == 1, ids)
            push!(frequencies, f)
            @info "Vehicle $k drives path $(paths[k]) with frequency $f."

            # OD frequency contribution, same logic as the fod constraints
            out_count = zeros(Int, D.n)
            in_count = zeros(Int, D.n)
            for a in ids
                (ten.is_driving[a] && ten.arc_dir[a] == ASC) || continue
                out_count[ten.node_h[ten.arc_tail[a]]] += 1
                in_count[ten.node_h[ten.arc_head[a]]]  += 1
            end
            for o in D.H, d in (o+1):D.n
                od_frequencies[o, d] += min(out_count[o], in_count[d], out_count[1])
            end
        end

        for o in D.H, d in (o+1):D.n
            od_frequencies[d, o] = od_frequencies[o, d]
            @assert od_frequencies[o, d] ≥ D.fmin "OD-frequency of pair ($o,$d) = $(od_frequencies[o,d]) does not meet the minimum frequency $(D.fmin)."
        end
        display(od_frequencies)

        # passenger statistics
        arc_flow = [sum(value(model[:p][e, o, d]) for o in D.H, d in D.H) for e in ten.Eidx]
        pax_km = sum(D.w[e...] * arc_flow[e] for e in ten.Eidx)
        booked_km = sum(D.w[o, d] * D.OD[o, d] for o in D.H, d in D.H) # all pax transported
        filled_distance = sum(ten.arc_len[a] for k in D.K for a in arcs_by_idx[k] if ten.is_driving[a] && arc_flow[ten.arc_edge[a]] > EPSILON; init = 0.0)
        empty_distance  = total_distance - filled_distance

        @assert empty_distance ≥ -EPSILON
        @assert filled_distance ≥ 0
        @assert total_distance ≥ 0
        @assert abs(total_distance - value(model[:total_distance])) ≤ EPSILON "Total distance $(total_distance) does not match model value $(value(model[:total_distance]))."
        @assert pax_km ≥ booked_km - EPSILON "Passagierkilometer $pax_km must be ≥ Buchungskilometer $booked_km."

        solution_vals["vehicle_flows"] = paths
        solution_vals["vehicle_flows_TEN"] = arcs
        solution_vals["vehicle_flows_TEN_idx"] = arcs_by_idx
        solution_vals["frequencies"] = frequencies
        solution_vals["od_freq"] = od_frequencies
        solution_vals["total_distance"] = total_distance
        solution_vals["1_empty_distance"] = empty_distance
        solution_vals["1_filled_distance"] = filled_distance
        solution_vals["1_share_empty"] = safe_div(empty_distance, total_distance)
        solution_vals["1_pax_km"] = pax_km
        solution_vals["1_booked_km"] = booked_km
        solution_vals["1_detour_factor"] = safe_div(pax_km, booked_km)
        solution_vals["1_system_efficiency"] = safe_div(booked_km, total_distance)
        solution_vals["1_occupancy_rate"] = safe_div(pax_km, filled_distance)
        solution_vals["1_pooling_factor"] = safe_div(pax_km, total_distance)

        if all(isfinite, (solution_vals["1_detour_factor"], solution_vals["1_occupancy_rate"], solution_vals["1_share_empty"]))
            @assert abs(solution_vals["1_system_efficiency"] - (1 / solution_vals["1_detour_factor"]) * solution_vals["1_occupancy_rate"] * (1 - solution_vals["1_share_empty"])) < EPSILON
        end
    end

    # second stage

    if relax_second_stage
        accepted_req_per_scenario = Vector{Vector{Any}}()
        rejected_req_per_scenario = Vector{Vector{Any}}()
    else
        accepted_req_per_scenario = Vector{Int}()
        rejected_req_per_scenario = Vector{Int}()
    end

    waiting_time_per_scenario = Vector{Vector{Float64}}()
    avg_waiting_time_per_scenario = Vector{Float64}()
    ride_time_per_scenario = Vector{Vector{Float64}}()
    avg_ride_time_per_scenario = Vector{Float64}()

    # Statistics
    booked_km_per_scenario = Vector{Float64}()
    pax_km_per_scenario = Vector{Float64}()
    filled_distance_per_scenario = Vector{Float64}()
    empty_distance_per_scenario = Vector{Float64}()

    for ω in Ω.scenarios
        allowed_arcs = get_arcs_for_all_requests(ten, ω)
        requests_on_arc = Dict{Tuple{Int,Int}, Vector{Int}}() # (vehicle, arc) => request
        booked_km = 0.0
        pax_km = 0.0
        filled_distance = 0.0
        empty_distance = 0.0
        waiting_times = zeros(ω.m)
        ride_times = zeros(ω.m)

        if relax_second_stage
            push!(accepted_req_per_scenario, [value(model[:z][ω.idx][r]) for r in ω.R])
            push!(rejected_req_per_scenario, [1 - value(model[:z][ω.idx][r]) for r in ω.R])
        else
            num_accepted = 0

            for r in ω.R
                used_vehicle_arcs = [(k, a) for k in D.K for a in allowed_arcs[r] if value(model[:y][ω.idx][r][k, a]) ≥ 0.9]
                if value(model[:z][ω.idx][r]) < 0.5
                    @assert isempty(used_vehicle_arcs) "Request $r is not accepted but uses arcs."
                    continue
                end
                num_accepted += 1
                @assert length(unique(first.(used_vehicle_arcs))) == 1 "Request $r must be assigned to exactly one vehicle."
                for (k, a) in used_vehicle_arcs
                    push!(get!(requests_on_arc, (k, a), Int[]), r)
                end

                pickup_arcs  = [a for (_, a) in used_vehicle_arcs if ten.node_h[ten.arc_tail[a]] == ω.origins[r]]
                dropoff_arcs = [a for (_, a) in used_vehicle_arcs if ten.node_h[ten.arc_head[a]] == ω.destinations[r]]
                @assert length(pickup_arcs) == 1 && length(dropoff_arcs) == 1 "Request $r needs exactly one pick-up and one drop-off arc."


                @assert ω.origin_tws[1, r] ≤ ten.arc_dep[pickup_arcs[1]] ≤ ω.origin_tws[2, r] "Departure at $(ten.arc_dep[pickup_arcs[1]]) of request $r in scenario $(ω.idx) is outside window $(ω.origin_tws[:, r])."
                @assert ω.destination_tws[1, r] ≤ ten.arc_arr[dropoff_arcs[1]] ≤ ω.destination_tws[2, r] "Arrival at $(ten.arc_arr[dropoff_arcs[1]]) of request $r outside window $(ω.destination_tws[:, r])."

                ride_times[r] = ten.arc_arr[dropoff_arcs[1]] - ten.arc_dep[pickup_arcs[1]]
                @assert 0 < ride_times[r] ≤ ω.max_ride_times[r] "Ride time of request $r is invalid: $(ride_times[r])."
                if ω.is_departure_request[r]
                    waiting_times[r] = ten.arc_dep[pickup_arcs[1]] - ω.requested_times[r]
                else
                    waiting_times[r] =  ω.requested_times[r] - ten.arc_arr[dropoff_arcs[1]]
                end
                @assert waiting_times[r] ≥ 0 "Request $r has negative waiting time."
                @assert value(model[:c][ω.idx][r]) ≥ waiting_times[r] - EPSILON "Request $r: $(value(model[:c][ω.idx][r])) is smaller than the actual waiting time."

                booked_km += D.w[ω.origins[r], ω.destinations[r]] * ω.loads[r]
            end
            @assert num_accepted == round(Int, sum(value.(model[:z][ω.idx]))) "Accepted requests $num_accepted do not match the sum of z."

            push!(accepted_req_per_scenario, num_accepted)
            push!(rejected_req_per_scenario, ω.m - num_accepted)
            push!(waiting_time_per_scenario, waiting_times)
            push!(avg_waiting_time_per_scenario, safe_div(sum(waiting_times), num_accepted))
            push!(ride_time_per_scenario, ride_times)
            push!(avg_ride_time_per_scenario, safe_div(sum(ride_times), num_accepted))
            @info "Scenario $(ω.idx): accepted $num_accepted ($(100 * num_accepted / ω.m)%), average waiting time $(safe_div(sum(waiting_times), num_accepted))."
        end
        
        vehicles_out = Dict{String, Any}()
        for k in D.K
            arcs_out = Dict{String, Any}[]
            for a in arcs_by_idx[k]
                t1, h1 = ten.node_t[ten.arc_tail[a]], ten.node_h[ten.arc_tail[a]]
                t2, h2 = ten.node_t[ten.arc_head[a]], ten.node_h[ten.arc_head[a]]
                requests = get(requests_on_arc, (k, a), Int[])

                if !relax_second_stage && ten.is_driving[a]
                    pax = sum(ω.loads[r] for r in requests; init = 0)
                    pax_km += ten.arc_len[a] * pax
                    if pax > 0
                        filled_distance += ten.arc_len[a]
                    else
                        empty_distance += ten.arc_len[a]
                    end
                end

                push!(arcs_out, Dict(
                    "from" => Dict("t" => t1, "h" => h1),
                    "to" => Dict("t" => t2, "h" => h2),
                    "direction" => ten.arc_dir[a],
                    "requests" => requests))
            end
            sort!(arcs_out, by = a -> a["from"]["t"])
            vehicles_out[string(k)] = Dict("arcs" => arcs_out)
        end

        if !relax_second_stage
            @assert pax_km ≥ booked_km - EPSILON "Passagierkilometer $pax_km must be ≥ Buchungskilometer $booked_km."
            @assert isapprox(filled_distance + empty_distance, total_distance; atol = EPSILON) "Filled $filled_distance + empty $empty_distance ≠ total $total_distance."
        end

        push!(booked_km_per_scenario, booked_km)
        push!(pax_km_per_scenario, pax_km)
        push!(filled_distance_per_scenario, filled_distance)
        push!(empty_distance_per_scenario, empty_distance)

        scenario_result = Dict("instance_name" => D.instance_name, "scenario_idx" => ω.idx, "vehicles" => vehicles_out)
        open(joinpath(D.output_directory, "solution_$(D.comment)_scenario_$(ω.idx).json"), "w") do f
            JSON.print(f, scenario_result, 2)
        end
    end

    solution_vals["accepted_req_per_scenario"] = accepted_req_per_scenario
    solution_vals["avg_accepted_req"] = sum(accepted_req_per_scenario) / Ω.num_scenarios
    solution_vals["rejected_req_per_scenario"] = rejected_req_per_scenario
    solution_vals["avg_rejected_req"] = sum(rejected_req_per_scenario) / Ω.num_scenarios
    solution_vals["waiting_time_per_scenario"] = waiting_time_per_scenario
    solution_vals["avg_waiting_time_per_scenario"] = avg_waiting_time_per_scenario
    solution_vals["avg_waiting_time"] = sum(avg_waiting_time_per_scenario) / Ω.num_scenarios
    solution_vals["ride_time_per_scenario"] = ride_time_per_scenario
    solution_vals["avg_ride_time_per_scenario"] = avg_ride_time_per_scenario
    solution_vals["avg_ride_time"] = sum(avg_ride_time_per_scenario) / Ω.num_scenarios

    solution_vals["2_empty_distance_per_scenario"] = empty_distance_per_scenario
    solution_vals["2_avg_empty_distance"] = sum(empty_distance_per_scenario) / Ω.num_scenarios
    solution_vals["2_filled_distance_per_scenario"] = filled_distance_per_scenario
    solution_vals["2_avg_filled_distance"] = sum(filled_distance_per_scenario) / Ω.num_scenarios
    solution_vals["2_share_empty_per_scenario"] = empty_distance_per_scenario ./ total_distance
    solution_vals["2_avg_share_empty"] = sum(solution_vals["2_share_empty_per_scenario"]) / Ω.num_scenarios
    solution_vals["2_pax_km_per_scenario"] = pax_km_per_scenario
    solution_vals["2_avg_pax_km"] = sum(pax_km_per_scenario) / Ω.num_scenarios
    solution_vals["2_booked_km_per_scenario"] = booked_km_per_scenario
    solution_vals["2_avg_booked_km"] = sum(booked_km_per_scenario) / Ω.num_scenarios
    solution_vals["2_detour_factor_per_scenario"] = pax_km_per_scenario ./ booked_km_per_scenario
    solution_vals["2_avg_detour_factor"] = sum(solution_vals["2_detour_factor_per_scenario"]) / Ω.num_scenarios
    solution_vals["2_system_efficiency_per_scenario"] = booked_km_per_scenario ./ total_distance
    solution_vals["2_avg_system_efficiency"] = sum(solution_vals["2_system_efficiency_per_scenario"]) / Ω.num_scenarios
    solution_vals["2_occupancy_rate_per_scenario"] = pax_km_per_scenario ./ filled_distance_per_scenario
    solution_vals["2_avg_occupancy_rate"] = sum(solution_vals["2_occupancy_rate_per_scenario"]) / Ω.num_scenarios
    solution_vals["2_pooling_factor_per_scenario"] = pax_km_per_scenario ./ total_distance
    solution_vals["2_avg_pooling_factor"] = sum(solution_vals["2_pooling_factor_per_scenario"]) / Ω.num_scenarios

    return solution_vals
end

function output_iis(D::InstanceData, model::Model)
    if !is_solved_and_feasible(model)
        compute_conflict!(model)
        if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
            counts = Dict{String,Int}()
            for (F, S) in list_of_constraint_types(model), con in all_constraints(model, F, S)
                MOI.get(model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT || continue
                nm = name(con)
                key = isempty(nm) ? "(bounds / integrality / unnamed)" : first(split(nm, r"_\d"; limit = 2))
                counts[key] = get(counts, key, 0) + 1
            end
            for (k, v) in sort(collect(counts), by = last, rev = true)
                println(rpad(k, 45), v)
            end
            println()
            iis_model, _ = copy_conflict(model)
            print(iis_model)
            write_to_file(iis_model, D.output_directory * D.instance_name * "_iis_model.lp")
        end
    end
end

function diagnose_first_stage(D, ten)
    L = ten.num_time_steps
    hop(h, g) = max(1, ceil(Int, D.w[h,g] / ten.tΔ))
    rt     = hop(1, D.n) + hop(D.n, 1)
    rt_all = sum(hop(h,h+1) + hop(h+1,h) for h in 1:D.n-1)
    println("L=$L K=$(length(D.K)) Q=$(D.Q) fmin=$(D.fmin) fmax=$(ten.fmax)")
    println("direct round trip $rt steps -> at most $(L ÷ rt) per vehicle")
    println("round trip over all stations $rt_all steps -> at most $(L ÷ rt_all) per vehicle")
    cap = D.Q * length(D.K) * (L ÷ rt)
    for h in 1:D.n-1
        asc  = sum(D.OD[o,d] for o in 1:h, d in (h+1):D.n)
        desc = sum(D.OD[o,d] for o in (h+1):D.n, d in 1:h)
        println("cut $h|$(h+1): demand ASC=$asc DESC=$desc, capacity bound=$cap")
    end
end
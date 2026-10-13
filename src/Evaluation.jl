using JuMP, Gurobi

function evaluate_solution(D::InstanceData, Ω::Scenarios, tΔ::Int, solution::Dict{String,Any}, which_obj_func::String, weights::Vector{Float64}, time_limit::Int, verbose::Bool; first_stage_solution::Dict=Dict())::Dict{String, Any}
    @assert "vehicle_flows_TEN" in keys(solution) "Solution does not contain time-expanded-network arcs."

    # === Build ===
    @info "Building evaluation model."

    start_build_time = time()
    model, ten = build_time_expanded_network_model(D, Ω, tΔ, false, false, false, false, true, false, verbose, which_obj_func, weights, solution["vehicle_flows_TEN_idx"])
    end_build_time = time()
    @info "Built evaluation model in $(round(end_build_time - start_build_time; digits=3)) s."

    # === Solve ===
    @info "Running evaluation."
    remaining_time = time_limit - (end_build_time - start_build_time)
    if remaining_time < 0
        return Dict{String, Any}("model_build_time" => (end_build_time - start_build_time), "time_limit" => time_limit)
    end

    set_time_limit_sec(model, remaining_time)
    set_optimizer_attribute(model, "LogFile", D.output_directory * D.instance_name * D.comment * "_evaluation_log.txt")
    if !verbose
        set_optimizer_attribute(model, "LogToConsole", 0)
        # write_to_file(model, D.output_directory * D.instance_name * "_model.mps")
    end

    optimize!(model)
    end_solve_time = time()
    @info "Evaluation model solved in $(round(end_solve_time - end_build_time; digits=3)) s."

    verify_solution(D, Ω, ten, model, true, false)

    solution_vals = get_solution_values(D, Ω, ten, model, which_obj_func, weights, true, false; first_stage_solution)

    solution_vals["model_build_time"] = end_build_time - start_build_time
    solution_vals["model_solve_time"] = end_solve_time - end_build_time
    solution_vals["total_time"] = end_solve_time - start_build_time
    solution_vals["time_limit"] = time_limit

    return solution_vals
end
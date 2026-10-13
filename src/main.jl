#!/usr/bin/env julia
# main.jl

using Pkg
using ArgParse
using Logging: ConsoleLogger, Warn
using JSON

Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

include("PlanningExpressLines.jl")
using .PlanningExpressLines

function main(path_to_config::String, scenario_dir::String, eval_dir::String, input_base_dir::String, output_base_dir::String;
    verbose::Bool=false,
    which_obj_func::String="pax",
    weights::Vector{Float64}=[1.0, 10.0],
    time_step::Int=1,
    comment::String="",
    solve_time_limit::Int=300 * 60,
    eval_time_limit::Int=120 * 60,
    directionality_property::Bool=false,
    relax_symmetry::Bool=false,
    warm_start::Bool=false,
    relax_second_stage::Bool=false,
    drive_all_stops::Bool=false,
)

    if !verbose
        Base.global_logger(ConsoleLogger(Warn))
    end

    @assert isfile(path_to_config) "Cannot find $path_to_config."
    @assert isdir(scenario_dir) "Scenario directory at $(scenario_dir) does not exist."
    @assert isdir(eval_dir) || eval_dir == "" "Evaluation directory at $(eval_dir) does not exist."
    @assert isdir(input_base_dir) || input_base_dir == "" "Input directory at $(input_base_dir) does not exist."
    @assert isdir(output_base_dir) || output_base_dir == "" "Output directory at $(output_base_dir) does not exist."
    @assert time_step ≥ 1 "Time step needs to be at least 1."
    
    which_obj_func = parse_obj_func(which_obj_func, weights)

    @info "Reading file $path_to_config."

    D = parse_input(path_to_config, input_base_dir, output_base_dir, comment)
    if !isdir(D.output_directory)
        mkpath(D.output_directory)
    end
    Ω_training = parse_scenarios(D, scenario_dir)

    @assert mod(D.T, time_step) == 0 "Maximum time of $(D.T) is not divisible by time step $(time_step)."
    @assert 2 * (D.n - 1) * time_step ≤ D.T "Maximum time of $(D.T) is not enough to visit all stations in both directions once with time step $(time_step)."
    @info "Number of training scenarios: $(Ω_training.num_scenarios)."
    
    solution_vals = build_and_solve(D, Ω_training, solve_time_limit, which_obj_func, weights, verbose, time_step, directionality_property, relax_symmetry, warm_start, relax_second_stage, drive_all_stops)

    @info "Finished solve; saving to file."
    save_output_to_file(D, solution_vals, comment, 
        Dict{String, Any}("time_step" => time_step,
        "directionality_property" => directionality_property,
        "relax_symmetry" => relax_symmetry,
        "warm_start" => warm_start,
        "relax_second_stage" => relax_second_stage,
        "drive_all_stops" => drive_all_stops,
        )
    )

    if eval_dir ≠ "" && solution_vals["has_solution"]
        Ω_evaluation = parse_scenarios(D, eval_dir)

        @info "Number of evaluation scenarios: $(Ω_evaluation.num_scenarios)"

        evaluation_vals = evaluate_solution(D, Ω_evaluation, time_step, solution_vals, which_obj_func, weights, eval_time_limit, verbose; first_stage_solution = solution_vals)

        @info "Finished evaluation; saving to file."
        save_output_to_file(D, evaluation_vals, comment * "_evaluation", 
            Dict{String, Any}("eval_time_limit" => eval_time_limit)
        )
    else
        @info "No evaluation."
    end    
end

function parse_obj_func(which_obj_func::String, weights::Vector{Float64})::String
    which_obj_func = lowercase(which_obj_func)
    @assert which_obj_func ∈ ["dist-pax", "dist-pax-wt", "pax-wt", "pax"] "Unknown objective function $(which_obj_func)"
    if which_obj_func ∈ ["dist-pax", "pax-wt"]
        @assert length(weights) == 2 "Number of weights passed ($(length(weights))) does not match objective function $which_obj_func."
    elseif which_obj_func == "dist-pax-wt"
        @assert length(weights) == 3 "Number of weights passed ($(length(weights))) does not match objective function $which_obj_func."
    end
    return which_obj_func
end

function save_output_to_file(D::InstanceData, solution_vals::Dict, comment::String, main_args::Dict)
    output_name = D.instance_name * comment * "_details.json"

    output_path = joinpath(D.output_directory, output_name)

    if isfile(output_path)
        @warn "Output file for $(output_name) exists and will be over-written."
    end

    instance_settings = get_instance_settings(D)

    output = merge(instance_settings, solution_vals, main_args)

    open(output_path, "w") do file
        JSON.json(file, output; allownan=true)
    end
end

# ensure running from BASH will pass all arguments
if abspath(PROGRAM_FILE) == @__FILE__
    using .Main  # ensure module context

    s = ArgParseSettings()

    @add_arg_table s begin
        "path_to_config"
            help = "Full path to config file."
            arg_type = String
            required = true
        "scenario_dir"
            help = "Directory to scenarios which should be evaluated in the 2nd stage."
            arg_type = String
            required = true
        "eval_dir"
            help = "Directory to scenarios which should be used for evaluating the solution."
            arg_type = String
            required = true
        "input_base_dir"
            help = "Base directory where inputs are located."
            arg_type = String
            required = true
        "output_base_dir"
            help = "Base directory where outputs should be saved."
            arg_type = String
            required = true
        "--verbose"
            help = "Enable verbose output."
            action = :store_true
        "--obj_func"
            help = "Which objective function?"
            arg_type = String
            default = "pax"
        "--weights"
            help = "Objective function weights"
            arg_type = Float64
            nargs = '*'
            default = [1.0, 10.0]
        "--time_step"
            help = "Time step for Time-Expanded Network model."
            arg_type = Int
            default = 1
        "--comment"
            help = "Optional comment for output file name."
            arg_type = String
            default = ""
        "--solve_time_limit"
            help = "Solve time limit, in seconds."
            arg_type = Int
            default = 300 * 60
        "--eval_time_limit"
            help = "Evaluation time limit, in seconds."
            arg_type = Int
            default = 120 * 60
        "--directionality_property"
            help = "Enforce directionality property?"
            action = :store_true
        "--relax_symmetry"
            help = "Relax restriction that paths have to be symmetric."
            action = :store_true
        "--warm_start"
            help = "Use a trivial solution to set an initial path for a warm start."
            action = :store_true
        "--relax_second_stage"
            help = "Relax integrality constraints on the second stage."
            action = :store_true
        "--drive_all_stops"
            help = "Force all vehicles to drive the line, i.e., all stops."
            action = :store_true
    end

    parsed_args = parse_args(ARGS, s)
    
    main(
        parsed_args["path_to_config"],
        parsed_args["scenario_dir"],
        parsed_args["eval_dir"],
        parsed_args["input_base_dir"],
        parsed_args["output_base_dir"];
        verbose = parsed_args["verbose"],
        which_obj_func = parsed_args["obj_func"],
        weights = Vector(parsed_args["weights"]),
        comment = parsed_args["comment"],
        solve_time_limit = parsed_args["solve_time_limit"],
        eval_time_limit = parsed_args["eval_time_limit"],
        directionality_property = parsed_args["directionality_property"],
        relax_symmetry = parsed_args["relax_symmetry"],
        warm_start = parsed_args["warm_start"],
        relax_second_stage = parsed_args["relax_second_stage"],
        time_step = parsed_args["time_step"],
        drive_all_stops = parsed_args["drive_all_stops"],
        )

    flush(stdout)
end
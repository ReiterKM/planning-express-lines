using DelimitedFiles
using JSON

@kwdef mutable struct InstanceData
    instance_name::String = ""
    output_directory::String
    comment::String = ""

    n::Int # number of stations
    H::UnitRange{Int} # set of stations
    w::Matrix{Float64} # distance matrix
    T::Int = 0 # time period in minutes
    max_waiting_time::Int = 15
    max_ride_time_factor::Float64 = 1.0
    fmin::Int = 0 # minimum frequency
    Q::Int = 0 # vehicle capacity
    OD::Matrix{Real} # demand matrix
    K::UnitRange{Int} # vehicles
end

function InstanceData(n::Int; Q::Int = 0, num_vehicles::Int = 0, max_waiting_time::Int = 15, max_ride_time_factor::Float64 = 1.0, T::Int = 0, instance_name::String = "", output_directory::String="", comment::String = "")
    @assert n >= 1 "At least one station is required."
    return InstanceData(
        n = n,
        H = 1:n,
        w = zeros(Float64, n,n),
        OD = zeros(n,n),
        Q = Q,
        max_waiting_time = max_waiting_time,
        max_ride_time_factor = max_ride_time_factor, 
        T = T,
        instance_name = instance_name,
        output_directory = output_directory,
        comment = comment,
        K = 1:num_vehicles,
    )
end

function parse_input(infile::String, indir::String, outdir::String, comment::String)::InstanceData
    input_data = JSON.parsefile(infile)

    D = InstanceData(input_data["number_of_stations"]; 
        Q = input_data["capacity"], 
        instance_name = input_data["instance_name"], 
        num_vehicles = input_data["number_of_vehicles"], 
        max_waiting_time = input_data["max_waiting_time"],
        max_ride_time_factor = input_data["max_ride_time_factor"],
        T = input_data["maximum_time"],
        output_directory = joinpath(outdir, input_data["output_directory"]),
        comment = comment)

    read_od_data!(D, joinpath(indir, input_data["path_to_od_matrix"]))

    if isempty(input_data["path_to_distance_matrix"])
        generate_distances!(D)
    else
        generate_distances!(D, joinpath(indir, input_data["path_to_distance_matrix"]))
    end

    # set_max_frequency!(D)
    set_min_frequency!(D)
    validate_input(D)

    return D
end

function set_min_frequency!(D::InstanceData)
    @assert D.max_waiting_time > 0 "Waiting time must be non-negative."
    D.fmin = ceil(Int, D.T / (2 * D.max_waiting_time))
    @info "Minimum frequency: $(D.fmin)."
end

function validate_input(D::InstanceData)
    @assert all(>=(0), D.w) "Distances cannot be negative."
    @assert all(>=(0), D.OD) "Demand cannot be negative."
    @assert D.T > 0 "Maximum time must be positive."
    @assert D.Q > 0 "Capacity must be positive."
    @assert length(D.K) > 0 "Number of vehicles must be positive."
    @assert D.fmin > 0 "Minimum frequency must be positive."
end

function read_od_data!(D::InstanceData, infile::String)
    @assert isfile(infile) "OD Matrix file $infile does not exist."

    D.OD = readdlm(infile)
    @assert size(D.OD) == (D.n, D.n) "Size of OD-matrix does not match the number of stations $(D.n)."
    @assert all(D.OD[i,i] == 0 for i in D.H) "OD-matrix needs to have a zero-diagonal."
end

function generate_distances!(D::InstanceData)
    for i in D.H
        for j in D.H
            if i != j
                D.w[i,j] = abs(i - j)
            end
        end
    end
end

function generate_distances!(D::InstanceData, infile::String)
    @assert isfile(infile) "Distance matrix $infile does not exist."

    D.w = readdlm(infile)
    @assert size(D.w) == (D.n, D.n) "Size of distance matrix does not match the number of stations $(D.n)."
end

function get_instance_settings(D::InstanceData)::Dict
    return Dict{String, Any}(
        "instance_name" => D.instance_name,
        "vehicle_capacity" => D.Q,
        "num_vehicles" => length(D.K),
        "num_stations" => D.n,
        "num_requests" => sum(D.OD),
        "maximum_time" => D.T,
    )
end
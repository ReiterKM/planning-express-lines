using JSON

@kwdef mutable struct ScenarioData
    idx::Int
    m::Int # number of requests
    R::UnitRange{Int} # set of requests

    origins::Vector{Int}
    destinations::Vector{Int}
    loads::Vector{Int}
    origin_tws::Matrix{Int} # time window at origin
    destination_tws::Matrix{Int} # time window at destination
    max_ride_times::Vector{Int} # maximum ride time (alpha * w[o_r, d_r])
    requested_times::Vector{Int}
    is_departure_request::BitVector # true if earliest departure time given

    R_ASC::Vector{Int} # set of ascending requests
    R_DESC::Vector{Int} # set of descending requests
    R_DIRS::Vector{Int} # asc / desc per request

    Γ_ASC::Set{Set{Int}} # set of overlapping requests in ASC direction
    Γ_DESC::Set{Set{Int}} # set of overlapping requests in DESC direction
end

function ScenarioData(m::Int, idx::Int)
    @assert m >= 1 "At least one request is required."
    return ScenarioData(
        idx = idx,
        m = m,
        R = 1:m,
        origins = zeros(Int, m),
        destinations = zeros(Int, m),
        loads = zeros(Int, m),
        origin_tws = zeros(Int, (2,m)),
        destination_tws = zeros(Int, (2,m)),
        max_ride_times = ones(Int, m),
        requested_times = zeros(Int, m),
        is_departure_request = falses(m),
        R_DIRS = Vector{Int}(),
        R_ASC = Vector{Int}(),
        R_DESC = Vector{Int}(),
        Γ_ASC = Set{Set{Int}}(),
        Γ_DESC = Set{Set{Int}}(),
    )
end

function ScenarioData(idx::Int, m::Int, origins::Vector{Any}, destinations::Vector{Any}, loads::Vector{Any})
    @assert m >= 1 "At least one request is required."
    return ScenarioData(
        idx = idx,
        m = m,
        R = 1:m,
        origins = origins,
        destinations = destinations,
        loads = loads,
        origin_tws = zeros(Int, (2,m)),
        destination_tws = zeros(Int, (2,m)),
        max_ride_times = ones(Int, m),
        requested_times = zeros(Int, m),
        is_departure_request = falses(m),
        R_DIRS = Vector{Int}(),
        R_ASC = Vector{Int}(),
        R_DESC = Vector{Int}(),
        Γ_ASC = Set{Set{Int}}(),
        Γ_DESC = Set{Set{Int}}(),
    )
end

@kwdef mutable struct Scenarios
    num_scenarios::Int
    scenarios::Vector{ScenarioData}
    probabilities::Vector{Real}
end

function Scenarios(num_scenarios::Int)
    @assert num_scenarios >= 0 "Number of scenarios cannot be negative."
    return Scenarios(
        num_scenarios = num_scenarios,
        scenarios = Vector{ScenarioData}(),
        probabilities = (1 / num_scenarios) * ones(num_scenarios)
    )
end

function Scenarios(num_scenarios::Int, scenarios::Vector{ScenarioData}, probabilities::Vector{Any})
    @assert num_scenarios >= 0 "Number of scenarios cannot be negative."
    @assert length(scenarios) == num_scenarios "Number of scenarios does not match."
    @assert length(probabilities) == num_scenarios "Number of probabilities does not match."
    @assert sum(probabilities) == 1 "Probabilities do not sum to one."
    @assert length(Set(ω.idx for ω in scenarios)) == num_scenarios "Indices of scenarios are not unique."
    return Scenarios(
        num_scenarios = num_scenarios,
        scenarios = scenarios,
        probabilities = probabilities
    )
end

function parse_scenarios(D::InstanceData, indir::String)::Scenarios
    scenario_files = [f for f in readdir(indir; join = true) if f[end-3:end] == "json" && !contains(f, "config")]
    sort!(scenario_files, by=f->parse(Int,split(basename(f)[1:end-5],"_")[end]))

    S = Scenarios(length(scenario_files))
    @assert sum(S.probabilities) - 1 ≤ EPSILON "Probabilities do not add up to 1."
    @assert D.max_ride_time_factor ≥ 1 "Maximum ride time factor $(D.max_ride_time_factor) must be ≥ 1."

    for (idx,infile) in enumerate(scenario_files)
        read_scenario_data!(D, S, infile, idx)
    end
    return S
end

function get_overlapping_requests(D::InstanceData, SD::ScenarioData, dir::Int)::Set{Set{Int}}
    overlapping_requests = Set{Set{Int}}()

    for h in 1:D.n
        overlap = Set{Int}()
        if dir == ASC
            for r in SD.R_ASC
                if (SD.origins[r] <= h) && (SD.destinations[r] > h)
                    push!(overlap, r)
                end
            end
        else
            for r in SD.R_DESC
                if (SD.origins[r] >= h) && (SD.destinations[r] < h)
                    push!(overlap, r)
                end
            end
        end
        if length(overlap) > 1
            push!(overlapping_requests, overlap)
        end
    end
    return overlapping_requests
end

function get_overlapping_requests(n::Int, SD::ScenarioData, dir::Int)::Set{Set{Int}}
    overlapping_requests = Set{Set{Int}}()

    for h in 1:n
        overlap = Set{Int}()
        if dir == ASC
            for r in SD.R_ASC
                if (SD.origins[r] <= h) && (SD.destinations[r] > h)
                    push!(overlap, r)
                end
            end
        else
            for r in SD.R_DESC
                if (SD.origins[r] >= h) && (SD.destinations[r] < h)
                    push!(overlap, r)
                end
            end
        end
        if length(overlap) > 1
            push!(overlapping_requests, overlap)
        end
    end
    return overlapping_requests
end

function read_scenario_data!(D::InstanceData, S::Scenarios, infile::String, idx::Int)
    data = JSON.parsefile(infile)

    SD = ScenarioData(idx,
        data["num_requests"],
        data["origins"],
        data["destinations"],
        data["loads"]
    )

    add_time_windows!(D,SD, data["departure_times"], data["arrival_times"])
    validate_scenario_data(D, SD)

    for r in SD.R
        if SD.origins[r] < SD.destinations[r]
            push!(SD.R_ASC, r)
            push!(SD.R_DIRS, ASC)
        elseif SD.origins[r] > SD.destinations[r]
            push!(SD.R_DESC, r)
            push!(SD.R_DIRS, DESC)
        else
            @error("Request $r has the same origin and destination $(S.origins[r]).")
        end
    end

    SD.Γ_ASC = get_overlapping_requests(D, SD, ASC)
    SD.Γ_DESC = get_overlapping_requests(D, SD, DESC)

    push!(S.scenarios, SD)
end

function add_time_windows!(D::InstanceData, SD::ScenarioData, departure_times::Vector{Any}, arrival_times::Vector{Any})
    # tightened TWs according to Cordeau (2006)
    for r in SD.R
        SD.max_ride_times[r] = ceil(Int, D.max_ride_time_factor * D.w[SD.origins[r], SD.destinations[r]])

        if !isnothing(departure_times[r]) && isnothing(arrival_times[r])
            @assert departure_times[r] ≥ 1 "Departure time of request $r is before zero: $(departure_times[r])."
            SD.requested_times[r] = departure_times[r]
            SD.is_departure_request[r] = true
            # e_(o_r) given, calculate l_(o_r)
            SD.origin_tws[1, r] = departure_times[r]
            SD.origin_tws[2, r] = min(D.T, departure_times[r] + D.max_waiting_time)
            # calculate e_(d_r), l_(d_r)
            SD.destination_tws[1, r] = min(D.T, ceil(Int, departure_times[r] + D.w[SD.origins[r], SD.destinations[r]]))
            SD.destination_tws[2, r] = min(D.T, departure_times[r] + D.max_waiting_time + SD.max_ride_times[r])
        elseif isnothing(departure_times[r]) && !isnothing(arrival_times[r])
            @assert arrival_times[r] ≤ D.T "Arrival time of request $r is after $(D.T): $(arrival_times[r])."
            SD.requested_times[r] = arrival_times[r]
            SD.is_departure_request[r] = false
            # calculate e_(o_r), l_(o_r)
            SD.origin_tws[1,r] = max(1, arrival_times[r] - D.max_waiting_time - SD.max_ride_times[r])
            SD.origin_tws[2,r] = max(1, ceil(Int, arrival_times[r] - D.w[SD.origins[r], SD.destinations[r]]))
            # l_(d_r) given, calculate e_(d_r)
            SD.destination_tws[1,r] = max(1, arrival_times[r] - D.max_waiting_time)
            SD.destination_tws[2,r] = arrival_times[r]
        else
            error("Cannot set time window for request $r.")
        end
    end
end

function validate_scenario_data(D::InstanceData, SD::ScenarioData)
    @assert SD.m ≥ 1 "Number of requests needs to be positive."
    @assert all(x -> x ≥(0) && x ≤(D.n), SD.origins) "Origins need to be within the number of stations."
    @assert all(x -> x ≥(0) && x ≤(D.n), SD.destinations) "Destinations need to be within the number of stations."
    @assert all(x -> x ≥(1) && x ≤(D.Q), SD.loads) "Loads need to be within capacity."
    @assert all(≥(1), SD.origin_tws) "Departure times need to be greater-equal 1."
    @assert all(≥(1), SD.destination_tws) "Arrival times need to be greater-equal 1."
    for r in SD.R
        @assert SD.origin_tws[1, r] ≤ SD.origin_tws[2, r] "Departure time window $(SD.origin_tws[:, r]) for request $r is invalid."
        @assert SD.destination_tws[1, r] ≤ SD.destination_tws[2, r] "Arrival time window $(SD.destination_tws[:, r]) for request $r is invalid."
        @assert SD.origin_tws[1, r] ≤ SD.destination_tws[1, r] "Arrival of request $r at $(SD.origin_tws[1, r]) is before departure at $(SD.destination_tws[1, r])."
        @assert SD.origin_tws[1, r] ≥ 0 "Departure of request $r at $(SD.origin_tws[1, r]) is before time 0."
        @assert SD.destination_tws[2, r] ≤ D.T "Arrival of request $r at $(SD.destination_tws[2, r]) is after time T = $(D.T)."
    end
end
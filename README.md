# Planning Express Lines with Aggregated Demand and Stochastic Requests
This repository contains the Julia implementation of the time-expanded network model to solve the single express line planning problem with aggregated demand and stochastic requests, as an integrated problem.

## Prerequisites
Make sure to have Julia (version 1.11.6) and Gurobi (version 10.0.1) installed. A valid Gurobi license is required.

Install the dependencies listed in `Project.toml`:
```
julia --project=. -e using Pkg; Pkg.instantiate()
```

## Quickstart
```
julia --project=. main.jl example_config.json scenarios/training scenarios/testing input/ output/
```

## Usage
The project can be run from the command line or directly from the Julia REPL.

From the command line, run the project using `main.jl` and pass optional arguments using the `--` syntax:

```
julia --project=. src/main.jl <path_to_config> <scenario_dir> <eval_dir> <input_base_dir> <output_base_dir> [--option=value ...]
```

Example:
```
julia --project=. src/main.jl input/example_config.json input/training/ input/testing/ input/ output/ --verbose 
```

The module `PlanningExpressLines.jl` bundles all required files.

### Required arguments

| Argument | Type | Description | 
|----------|------|-------------| 
|`path_to_config` | String | path to configuration (JSON) file |
|`scenario_dir` | String | path to directory with training scenarios |
|`eval_dir` | String | path to directory with testing scenarios |
|`input_base_dir` | String | path to directory with input file(s) |
|`output_base_dir` | String | path to directory for output file(s) |

The `scenario_dir` and `eval_dir` fields may be empty strings (`""`), indicating that no training or testing should be carried out.

#### Optional arguments
| Argument | Type | Default | Description |
|----------|------|---------|-------------|
|`verbose`|Flag|false|verbose output|
|`which_obj_func`|String|"pax"|which objective function to use |
|`weights`| Float64, multiple values | `1.0 10.0` | Weights for the scalarized objective, passed space-separated (e.g. `--weights 1.0 10.0`) |
|`time_step` | Int | 1| time step $t_{\Delta}$ |
|`comment`|String|""|comment to add to output file names|
|`solve_time_limit`|Int|18000|time limit to solve the model, in seconds|
|`eval_time_limit`|Int|7200|time limit for testing the solution, in seconds|
|`directionality_property`|Flag|false|enforce directionality property for aggregated requests|
|`relax_symmetry`|Flag|false|allow non-symmetric express lines|
|`warm_start`|Flag|false|use a trivial solution as a warm start|
|`relax_second_stage`|Flag|false|relax integrality constraints in the second stage (stochastic scenarios)|
|`drive_all_stops`|Flag|false|force express lines to drive all stops, i.e., forbid shortcuts |

Arguments of type Flag take no value: pass `--verbose` to enable them; omit them to keep the default (off).

### Objective functions
The objective function passed via `which_obj_func` must be one of the following:

| Value | Description |
|-------|-------------|
| `"pax"` | maximize number of requests transported |
| `"pax-wt"` | maximize number of requests transported, minimize waiting time |
| `"dist-pax"` | minimize total vehicle distance, maximize number of requests transported |
| `"dist-pax-wt"` | minimize total vehicle distance, maximize number of requests transported, minimize waiting time |

All except for the `pax` objective function are treated as a scalarized sum with an appropriate `weights` vector.

## Inputs
### Configuration File Format
The config file is a JSON with the problem instance settings:

| Field | Type | Description |
|-------|------|-------------|
| `instance_name` | String | name of the instance |
| `output_directory` | String | subdirectory of `output_base_dir` for the output |
| `path_to_od_matrix` | String | path to the aggregated OD matrix |
| `path_to_distance_matrix` | String | (optional) path to the distance matrix of the line |
| `number_of_stations` | Int | number of stations $n$ |
| `number_of_vehicles` | Int | number of vehicles |
| `capacity` | Int | vehicle capacity |
| `max_ride_time_factor` | Float | maximum detour factor $\alpha$ |
| `max_waiting_time` | Int | maximum waiting time $\beta$|
| `maximum_time` | Int | maximum period length $T$|

All fields except the `path_to_distance_matrix` are required. See the `example_config.json` provided for an example.

### Aggregated OD Demand File
OD demand must be provided as a complete $n \times n$ matrix with zeros on the diagonal, formatted as

````
0 C_{1,2} ... C_{1,n}
C_{2,1} 0 ... C_{2,n}
...
C_{n,1} C_{n,2}... 0
````

### Distance File
Distances between stations can optionally be provided as a complete $n \times n$ matrix with zeros on the diagonal, formatted as

````
0 w_{1,2} ... w_{1,n}
w_{2,1} 0 ... w_{2,n}
...
w_{n,1} w_{n,2}... 0
````

If `path_to_distance_matrix` is not given, all distances between adjacent stations will be set to 1.

### Scenarios
Training and testing scenarios are provided via a JSON file per scenario, in the following format:

| Field | Type | Description |
|-------|------|-------------|
| `instance_name` | String | name of the instance |
| `num_request` | Int | number of requests $m$ |
| `origins` | Vector{Int} | origins $o_r$ of the $m$ requests |
| `destinations` | Vector{Int} | destinations $d_r$ of the $m$ requests |
| `departure_times` | Vector{Int, null} | earliest departure times, or `null` |
| `arrival_times` | Vector{Int, null} | latest arrival times, or `null` |
| `loads` | Vector{Int} | number of passengers $q_r$ per request |

For each request, either the departure time or the arrival time must be given; if both are `null`, an error is raised.

## Outputs

Results are written to `output_base_dir/<output_directory>/`.

| File | Description |
|------|-------------|
| `*_details.json` | settings, solution and statistics of the training |
| `*_evaluation_details.json` | settings, solution and statistics of the testing |

In addition, the solution for each scenario (training and testing) is saved.

## Support
E-Mail kendra.reiter@uni-wuerzburg.de. Feedback and improvements are welcome!

### Citation
Cite this work as

````
@software{reiter2026expressline,
    author = {Reiter, Kendra},
    title = {Planning Express Lines with Aggregated Demand and Stochastic Requests},
    url = {https://github.com/ReiterKM/planning-express-lines},
    year = {2026},
    note = {Software}
}
````

### Acknowledgement
This work was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – 566982931.

## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />

This work is licensed under [CC BY-NC-SA 4.0 ](https://creativecommons.org/licenses/by-nc-sa/4.0/).
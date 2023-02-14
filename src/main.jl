include("utils.jl")
include("CreatSample.jl")
include("Crossovers.jl")
include("Mutations.jl")
include("DP.jl")
include("GA.jl")
include("Initial.jl")
include("Neighborhood.jl")
include("Neighborhood_local_optima.jl")
include("read_files.jl")


function solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, T::Matrix{Float64}, D::Matrix{Float64};
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)

    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_not_Eligible)
end


function solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, depot::Vector{Float64}, 
    Customers::Matrix{Float64}, tspeed::Float64, dspeed::Float64;
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)

    T, D = Calculate_duration_matrices(tspeed, dspeed, depot, Customers, problem_type)
    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_not_Eligible)
end


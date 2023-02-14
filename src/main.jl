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
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Int = 0, sL::Int = 0)

    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_not_Eligible)
end


function solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, depot::Vector{Float64}, 
    Customers::Matrix{Float64}, tspeed::Float64, dspeed::Float64;
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Int = 0, sL::Int = 0)

    T, D = Calculate_duration_matrices(tspeed, dspeed, depot, Customers, problem_type)
    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_not_Eligible)
end


# depot, Customers, dEligible = Read_TSPLIB_instance(:berlin52)
# T, D = Calculate_duration_matrices(40.0, 40.0, depot, Customers, FSTSP)
# solve_tspd_by_HGA_TAC(FSTSP, 5, T, D, drone_not_Eligible = dEligible, flying_range = 40.0, sR = 1, sL = 1)

# P = solve_tspd_by_HGA_TAC(FSTSP, 5, depot, Customers, 40.0, 40.0, drone_not_Eligible = dEligible, flying_range = 40.0, sR = 1, sL = 1)
# T, D = Calculate_duration_matrices(40.0, 40.0, depot, Customers, FSTSP)
# T, D = read_data_Agatz("uniform-51-n10")
# solve_tspd_by_HGA_TAC(FSTSP, 5, T, D)

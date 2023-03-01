using Distances
using StatsBase
using LKH
using Shuffle

include("utils.jl")
include("CreatSample.jl")
include("Crossovers.jl")
include("Mutations.jl")
include("DP.jl")
include("GA.jl")
include("Initial.jl")
include("Neighborhood.jl")
include("Neighborhood_local_optima.jl")


# This function is compatible with TSPDrone.solve_tspd
# https://github.com/chkwon/TSPDrone.jl
function solve_tspd(
    truck_cost_mtx::Matrix{Float64}, 
    drone_cost_mtx::Matrix{Float64};
    num_runs::Int64=1,
    flying_range::Float64=Inf,
)
    # the size of truck_cost_mtx and drone_cost_mtx is (n_customers + 1) x (n_customers + 1)
    # the first row and column are the depot

    # What need to be T and D? Shouldn't each of T and D be Ct and Cd?
    Ct = [
        truck_cost_mtx          truck_cost_mtx[:, 1];
        truck_cost_mtx[1, :]'    0.0
    ]

    Cd = [
        drone_cost_mtx          drone_cost_mtx[:, 1];
        drone_cost_mtx[1, :]'    0.0
    ]

    result = solve_tspd_by_HGA_TAC(
        TSPD, 
        num_runs, 
        Ct, 
        Cd, 
        flying_range=flying_range
    )

    @assert length(result) == 1
    
    return result[1]

end



# This function is compatible with TSPDrone.solve_tspd
# https://github.com/chkwon/TSPDrone.jl
function solve_tspd(
    x::Vector{Float64}, 
    y::Vector{Float64}, 
    truck_cost_factor::Float64, 
    drone_cost_factor::Float64;
    num_runs::Int64=1,
    flying_range::Float64=Inf,
)

    # the first element in x and y is the depot
    # the rest of the elements are the customers
    depot = [x[1], y[1]]
    customers = [x[2:end] y[2:end]]

    result = solve_tspd_by_HGA_TAC(
        TSPD, 
        num_runs, 
        depot, 
        customers, 
        1 / truck_cost_factor, 
        1 / drone_cost_factor, 
        flying_range=flying_range
    )

    @assert length(result) == 1
    
    return result[1]

end


function solve_tspd_by_HGA_TAC(problem_type::ProblemType, num_runs::Int64, T::Matrix{Float64}, D::Matrix{Float64};
    drone_not_Eligible::Vector{Int}=Int[], flying_range::Float64=Inf, sR::Float64=0.0, sL::Float64=0.0)

    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_not_Eligible)
end


function solve_tspd_by_HGA_TAC(problem_type::ProblemType, num_runs::Int64, depot::Vector{Float64},
    Customers::Matrix{Float64}, tspeed::Float64, dspeed::Float64;
    drone_not_Eligible::Vector{Int}=Int[], flying_range::Float64=Inf, sR::Float64=0.0, sL::Float64=0.0)

    T, D = Calculate_duration_matrices(tspeed, dspeed, depot, Customers, problem_type)

    @show size(T), size(D)
    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_not_Eligible)
end


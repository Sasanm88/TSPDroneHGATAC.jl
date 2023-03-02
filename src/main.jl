using Distances
using StatsBase
using LKH
using Shuffle

include("utils.jl")
include("create_samples.jl")
include("crossovers.jl")
include("mutations.jl")
include("dynamic_programming.jl")
include("genetic_algorithm.jl")
include("initial.jl")
include("neighborhood.jl")
include("neighborhood_local_optima.jl")


# This function is compatible with TSPDrone.solve_tspd
# https://github.com/chkwon/TSPDrone.jl
function solve_tspd(
    truck_cost_mtx::Matrix{Float64},
    drone_cost_mtx::Matrix{Float64};
    num_runs::Int64=1,
    flying_range::Float64=Inf
)

    @assert size(truck_cost_mtx) == size(drone_cost_mtx)

    # check the diagonal elements are all zero
    for i in 1:size(truck_cost_mtx)[1]
        @assert truck_cost_mtx[i, i] == 0.0
        @assert drone_cost_mtx[i, i] == 0.0
    end


    # the size of truck_cost_mtx and drone_cost_mtx is (n_customers + 1) x (n_customers + 1)
    # the first row and column are the depot
    # the rest of the rows and columns are the customers

    # we need to add the depot to the end of the cost matrix    
    T, D = cost_matrices_with_dummy(truck_cost_mtx, drone_cost_mtx)


    routes = solve_tspd_by_HGA_TAC(
        TSPD,
        num_runs,
        T,
        D,
        flying_range=flying_range
    )

    @assert length(routes) == num_runs

    result = prepare_return_value(routes)

    return result

end



# This function is compatible with TSPDrone.solve_tspd
# https://github.com/chkwon/TSPDrone.jl
function solve_tspd(
    x::Vector{Float64},
    y::Vector{Float64};
    truck_cost_factor::Float64=1.0,
    drone_cost_factor::Float64=0.5,
    num_runs::Int64=1,
    flying_range::Float64=Inf
)

    # the first element in x and y is the depot
    # the rest of the elements are the customers
    depot = [x[1], y[1]]
    customers = [x[2:end] y[2:end]]

    routes = solve_tspd_by_HGA_TAC(
        TSPD,
        num_runs,
        depot,
        customers,
        1 / truck_cost_factor,
        1 / drone_cost_factor,
        flying_range=flying_range
    )

    @assert length(routes) == num_runs

    result = prepare_return_value(routes)

    return result

end



function solve_fstsp(
    truck_cost_mtx::Matrix{Float64},
    drone_cost_mtx::Matrix{Float64};
    num_runs::Int64=1,
    flying_range::Float64=Inf,
    drone_ineligible_nodes::Vector{Int}=Int[],
    sR::Float64=0.0,
    sL::Float64=0.0
)

    @assert size(truck_cost_mtx) == size(drone_cost_mtx)

    # check the diagonal elements are all zero
    for i in 1:size(truck_cost_mtx)[1]
        @assert truck_cost_mtx[i, i] == 0.0
        @assert drone_cost_mtx[i, i] == 0.0
    end

    # the size of truck_cost_mtx and drone_cost_mtx is (n_customers + 1) x (n_customers + 1)
    # the first row and column are the depot
    # the rest of the rows and columns are the customers

    # we need to add the depot to the end of the cost matrix    
    T, D = cost_matrices_with_dummy(truck_cost_mtx, drone_cost_mtx)

    routes = solve_tspd_by_HGA_TAC(
        FSTSP,
        num_runs,
        T,
        D,
        flying_range=flying_range,
        drone_ineligible_nodes=drone_ineligible_nodes,
        sR=sR,
        sL=sL
    )

    @assert length(routes) == num_runs

    result = prepare_return_value(routes)

    return result


end


function solve_fstsp(
    x::Vector{Float64},
    y::Vector{Float64};
    num_runs::Int64=1,
    flying_range::Float64=Inf,
    truck_speed::Float64=1.0,
    drone_speed::Float64=2.0,
    drone_ineligible_nodes::Vector{Int}=Int[],
    sR::Float64=0.0,
    sL::Float64=0.0
)

    # the first element in x and y is the depot
    # the rest of the elements are the customers
    depot = [x[1], y[1]]
    customers = [x[2:end] y[2:end]]

    routes = solve_tspd_by_HGA_TAC(
        FSTSP,
        num_runs,
        depot,
        customers,
        truck_speed,
        drone_speed;
        drone_ineligible_nodes=drone_ineligible_nodes,
        flying_range=flying_range,
        sR=sR,
        sL=sL
    )

    @assert length(routes) == num_runs

    result = prepare_return_value(routes)

    return result

end










# Internal functions to run HGA_TAC

function solve_tspd_by_HGA_TAC(
    problem_type::ProblemType,
    num_runs::Int64,
    T::Matrix{Float64},
    D::Matrix{Float64};
    drone_ineligible_nodes::Vector{Int}=Int[],
    flying_range::Float64=Inf,
    sR::Float64=0.0,
    sL::Float64=0.0
)

    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_ineligible_nodes)
end


function solve_tspd_by_HGA_TAC(
    problem_type::ProblemType,
    num_runs::Int64,
    depot::Vector{Float64},
    Customers::Matrix{Float64},
    tspeed::Float64,
    dspeed::Float64;
    drone_ineligible_nodes::Vector{Int}=Int[],
    flying_range::Float64=Inf,
    sR::Float64=0.0,
    sL::Float64=0.0
)

    T, D = Calculate_duration_matrices(tspeed, dspeed, depot, Customers, problem_type)

    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_ineligible_nodes)
end


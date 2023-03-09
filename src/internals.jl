

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

    T, D = calculate_duration_matrices(tspeed, dspeed, depot, Customers, problem_type)

    return run_GA(problem_type, num_runs, T, D, flying_range, sR, sL, drone_ineligible_nodes)
end

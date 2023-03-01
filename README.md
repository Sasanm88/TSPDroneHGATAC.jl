# TSPDroneHGATAC.jl
This is a package that uses Hybrid Genetic Algorithm to solve any TSPD or FSTSP instance. 
The main function that can be utilized by the user is "solve_tspd_by_HGA_TAC". 
This function can take two different set of arguments. 
1. solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, T::Matrix{Float64}, D::Matrix{Float64};
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)
2. solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, depot::Vector{Float64}, 
    Customers::Matrix{Float64}, tspeed::Float64, dspeed::Float64; drone_not_Eligible::Vector{Int} = Int[],
    flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)
    
problem_type is an argument that takes either TSPD or FSTSP. 


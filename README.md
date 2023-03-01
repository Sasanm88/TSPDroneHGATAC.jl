# TSPDroneHGATAC.jl
This is a package that uses Hybrid Genetic Algorithm to solve any TSPD or FSTSP instance. 
The main function that can be utilized by the user is "solve_tspd_by_HGA_TAC". 
This function can take two different set of arguments. 
1. solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, T::Matrix{Float64}, D::Matrix{Float64};
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)
2. solve_tspd_by_HGA_TAC(problem_type::problem, num_runs::Int64, depot::Vector{Float64}, 
    Customers::Matrix{Float64}, tspeed::Float64, dspeed::Float64; drone_not_Eligible::Vector{Int} = Int[],
    flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)
    
* problem_type is an argument that takes either TSPD or FSTSP. 

* num_runs is the number of times that user wants the algorithm to solve the required problem. The output will be a vector of size num_runs. 

* drone_not_Eligible is a vector of integers containing the list of all customers that are not eligible to be served by drone. The default value is an empty vector, we assume if the user does not pass a value for this argument, then all the customers are eligible to be visited by drone. 

* flying_range is a float value that refers to the drone maximum endurance in terms of seconds. The default value is infinity. 
* sR is the retreival time when the drone is landing on the truck. Default value is assumed to be zero. 
* sL is the time needed for luanching the drone. Default value is assumed to be zero. 




# TSPDroneHGATAC.jl


A package that uses Hybrid Genetic Algorithm to solve any TSPD or FSTSP instance

# Installation

To install:
```julia
import Pkg; Pkg.add(url="https://github.com/Sasanm88/TSPDroneHGATAC.jl")
```
or
```julia
] add https://github.com/Sasanm88/TSPDroneHGATAC.jl
```


To install `TSPDrone.jl`:
```julia
] add https://github.com/chkwon/TSPDrone.jl
```


## Examples 

This works....
```julia
using TSPDroneHGATAC
n = 10 
x = rand(n); y = rand(n);

result1 = TSPDroneHGATAC.solve_tspd(x, y, 1.0, 0.5)

# compare with TSPDrone.jl
using TSPDrone
result2 = TSPDrone.solve_tspd(x, y, 1.0, 0.5)

@show result1[1].total_cost
@show result2.total_cost
```

Another way...
```julia
using TSPDroneHGATAC
n = 10 
dist_mtx = rand(n, n)
dist_mtx = dist_mtx + dist_mtx' # symmetric distance only
for i in 1:n # diagonal needs to be zero
    dist_mtx[i, i] = 0.0
end
truck_cost_mtx = dist_mtx .* 1.0
drone_cost_mtx = truck_cost_mtx .* 0.5 

result1 = TSPDroneHGATAC.solve_tspd(truck_cost_mtx, drone_cost_mtx)

# compare with TSPDrone.jl
using TSPDrone
result2 = TSPDrone.solve_tspd(truck_cost_mtx, drone_cost_mtx)

@show result1[1].total_cost
@show result2.total_cost
```


## Internal Solver Functions

This is a package that uses Hybrid Genetic Algorithm to solve any TSPD or FSTSP instance. 
The main function that can be utilized by the user is "solve_tspd_by_HGA_TAC". 
This function can take two different set of arguments. 

```julia
solve_tspd_by_HGA_TAC(problem_type::ProblemType, num_runs::Int64, T::Matrix{Float64}, D::Matrix{Float64};
    drone_not_Eligible::Vector{Int} = Int[], flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)
```

and

```julia
solve_tspd_by_HGA_TAC(problem_type::ProblemType, num_runs::Int64, depot::Vector{Float64}, 
    Customers::Matrix{Float64}, tspeed::Float64, dspeed::Float64; drone_not_Eligible::Vector{Int} = Int[],
    flying_range::Float64 = Inf, sR::Float64 = 0.0, sL::Float64 = 0.0)
```    

* `problem_type` is an argument that takes either `TSPD` or `FSTSP`. 
* `num_runs` is the number of times that user wants the algorithm to solve the required problem. The output will be a vector of size `num_runs`. 
    
* problem_type is an argument that takes either TSPD or FSTSP. 

* num_runs is the number of times that user wants the algorithm to solve the required problem. The output will be a vector of size num_runs. 

* drone_not_Eligible is a vector of integers containing the list of all customers that are not eligible to be served by drone. The default value is an empty vector, we assume if the user does not pass a value for this argument, then all the customers are eligible to be visited by drone. 

* flying_range is a float value that refers to the drone maximum endurance in terms of seconds. The default value is infinity. 
* sR is the retreival time when the drone is landing on the truck. Default value is assumed to be zero. 
* sL is the time needed for luanching the drone. Default value is assumed to be zero. 




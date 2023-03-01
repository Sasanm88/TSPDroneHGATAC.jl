# TSPDroneHGATAC.jl

A package that uses Hybrid Genetic Algorithm to solve any TSPD or FSTSP instance

# Installation

To install:
```julia
import Pkg; Pkg.add("https://github.com/Sasanm88/TSPDroneHGATAC.jl")
```
or
```julia
] add https://github.com/Sasanm88/TSPDroneHGATAC.jl
```

## Examples 

This works....
```julia
using TSPDroneHGATAC
n = 10 
x = rand(n); y = rand(n);

result1 = TSPDroneHGATAC.solve_tspd(x, y, 1.0, 0.5)

import Pkg; Pkg.add("https://github.com/chkwon/TSPDrone.jl")
using TSPDrone
result2 = TSPDrone.solve_tspd(x, y, 1.0, 0.5)

@show result1[1].total_cost
@show result2.total_cost
```


The below makes an error.... 
```julia
using TSPDroneHGATAC
n = 10 
dist_mtx = rand(n, n)
dist_mtx = dist_mtx + dist_mtx' # symmetric distance only
truck_cost_mtx = dist_mtx .* 1.0
drone_cost_mtx = truck_cost_mtx .* 0.5 

result1 = TSPDroneHGATAC.solve_tspd(truck_cost_mtx, drone_cost_mtx)

import Pkg; Pkg.add("https://github.com/chkwon/TSPDrone.jl")
using TSPDrone
result2 = TSPDrone.solve_tspd(truck_cost_mtx, drone_cost_mtx)

@show result1[1].total_cost
@show result2.total_cost
```
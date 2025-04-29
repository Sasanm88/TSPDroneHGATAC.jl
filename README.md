# TSPDroneHGATAC.jl


A package that uses Hybrid Genetic Algorithm to solve any TSPD or FSTSP instance, based on the following paper:


https://arxiv.org/abs/2303.00614

```
@misc{mahmoudinazlou_hybrid,
  author = {Mahmoudinazlou, Sasan and Kwon, Changhyun},
  title = {A Hybrid Genetic Algorithm with Type-Aware Chromosomes for Traveling Salesman Problems with Drone},
  year = {2023},
  doi = {10.48550/ARXIV.2303.00614},  
  url = {https://arxiv.org/abs/2303.00614},
}
```




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


## TSP-D Examples 

Using `x`-`y` coordinates. 
Based on the coordinates, the Euclidean distances are calculated. 
The Eucliden distances will be divided by `truck_speed` and `drone_speed` to calculate the travel cost for the truck and the drone, respectively.


```julia
using TSPDroneHGATAC
n = 10 
x = rand(n); y = rand(n);

result1 = TSPDroneHGATAC.solve_tspd(x, y, truck_speed=1.0, drone_speed=2.0)

# compare with TSPDrone.jl
using TSPDrone
result2 = TSPDrone.solve_tspd(x, y, 1.0, 0.5)

@show result1.best_total_cost
@show result2.total_cost
```

You can also input the distance cost matrices for the truck and the drone directly.
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

@show result1.best_total_cost
@show result2.total_cost
```


## FS-TSP Examples 


Based on the coordinates:
```julia
using TSPDroneHGATAC
n = 10 
x = rand(n); y = rand(n);

result1 = TSPDroneHGATAC.solve_fstsp(x, y; truck_speed=1.0, drone_speed=2.0, flying_range=0.5, retrieval_time=0.1, launching_time=0.1)

@show result1.best_total_cost
```


Based on the cost matrices:
```julia
using TSPDroneHGATAC
n = 10 
dist_mtx = rand(n, n)
dist_mtx = dist_mtx + dist_mtx' # symmetric distance only
for i in 1:n # diagonal needs to be zero
    dist_mtx[i, i] = 0.0
end
truck_cost_mtx = dist_mtx ./ 1.0
drone_cost_mtx = truck_cost_mtx ./ 2.0 

result1 = TSPDroneHGATAC.solve_fstsp(truck_cost_mtx, drone_cost_mtx; flying_range=0.5, retrieval_time=0.1, launching_time=0.1)

@show result1.best_total_cost
```


## Related Packages

- [TSPDroneLIB](https://github.com/chkwon/TSPDroneLIB)
- [TSPDrone.jl](https://github.com/chkwon/TSPDrone.jl)
- [TSPDroneHGATAC.jl](https://github.com/Sasanm88/TSPDroneHGATAC.jl)
- [TSPDroneICP.jl](https://github.com/0505daniel/TSPDroneICP.jl)


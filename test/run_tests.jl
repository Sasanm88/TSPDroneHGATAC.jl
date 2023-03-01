include("../src/main.jl")

include("read_files.jl")
include("test_instances.jl")


# Solve_Agatz_unrestricted("uniform-alpha_1-75-n50", 10)   # 10 test runs
# Solve_Agatz_restricted("uniform-62-n20-maxradius-20", 10)  # 10 test runs
Solve_Bogyrbayeva("Amsterdam-n_nodes-50", 1, 10)  # instance number 1, 10 test runs
Solve_Murray("20140810T123437v1", 20.0, 10) # flying_range = 20, 10 test runs
Solve_Ha("mbA101", 10) # 10 test runs
Solve_TSPLIB(:berlin52, Int[2, 8, 26, 42], 40.0, 40.0, 40.0, 1.0, 1.0, 10) # nodes 2,8,26 and 42 cannot be served by drone 
# flying_range = 40, truck speed = 40, drone speed = 40, service time = 1, retrieval time = 1, 10 test runs
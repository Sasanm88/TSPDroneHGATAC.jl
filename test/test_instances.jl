
include("../src/main.jl")


function Solve_Agatz_unrestricted(sample::String, num_runs::Int64)
    T, D = read_data_Agatz(sample)
    return solve_tspd_by_HGA_TAC(TSPD, num_runs, T, D)
end

function Solve_Agatz_restricted(sample::String, num_runs::Int64)
    T, D, flying_range = read_data_Agatz_restricted(sample)
    return solve_tspd_by_HGA_TAC(TSPD, num_runs, T, D, flying_range = flying_range)
end

function Solve_Bogyrbayeva(file_name::String, sample_number::Int64, num_runs::Int64)
    depot, customers = read_data_Bogyrbayeva(file_name, sample_number)
    return solve_tspd_by_HGA_TAC(TSPD, num_runs, depot, customers, 1.0, 1.0)
end


Solve_Bogyrbayeva("Amsterdam-n_nodes-10", 1, 5)

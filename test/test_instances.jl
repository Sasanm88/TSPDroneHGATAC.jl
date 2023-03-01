
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

function Solve_Murray(file_name::String, flying_range::Float64, num_runs::Int64)
    T, D, drone_not_eligible = read_Murray(file_name)
    return solve_tspd_by_HGA_TAC(FSTSP, num_runs, T, D, drone_not_Eligible = drone_not_eligible,
     flying_range = flying_range, sR = 1.0, sL = 1.0)
end

function Solve_Ha(file_name::String, num_runs::Int64)
    depot, customers, drone_not_eligible, tspeed, dspeed, flying_range, sL, sR = read_data_Ha(file_name)
    return solve_tspd_by_HGA_TAC(FSTSP, num_runs, depot, customers, tspeed, dspeed, drone_not_Eligible = drone_not_eligible,
    flying_range = flying_range, sR = sR, sL = sL)
end

function Solve_TSPLIB(file_name::Symbol, drone_not_eligible::Vector{Int}, flying_range::Float64, tspeed::Float64, dspeed::Float64, sR::Float64, sL::Float64, num_runs::Int64)
    depot, customers = Read_TSPLIB_instance(file_name)
    return solve_tspd_by_HGA_TAC(FSTSP, num_runs, depot, customers, tspeed, dspeed, drone_not_Eligible = drone_not_eligible,
    flying_range = flying_range, sR = sR, sL = sL)
end




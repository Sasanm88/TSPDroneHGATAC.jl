using LKH
using Shuffle

function convert_TSPtour_to_chromosome(tsp_tour::Vector{Int64}, n_nodes::Int64)
    chrm = Vector{Int}()
    if tsp_tour[2] == n_nodes + 2
        chrm = tsp_tour[3:n_nodes+2]
    else
        chrm = tsp_tour[2:n_nodes+1]
    end
    return chrm
end

function find_initial_TSP_tour(TT::Matrix{Float64}, n_nodes::Int64)
    tsp_tour, _ = find_tsp_tour1(TT)
    return convert_TSPtour_to_chromosome(tsp_tour, n_nodes)
end

function find_tsp_tour1(Ct::Matrix{Float64})
    scale_factor = 1000
    dist_mtx = round.(Int, Ct .* scale_factor)

    tsp_tour, tsp_tour_len = LKH.solve_tsp(dist_mtx)

    @assert tsp_tour[1] == 1

    return tsp_tour, tsp_tour_len
end

# MAX_DRONE_RANGE = Inf

function exact_p(initial_tour::Vector{Int64}, Ct::Matrix{Float64}, Cd::Matrix{Float64}, flying_range::Float64, sR::Float64, sL::Float64)
    n, _ = size(Ct)

    r = initial_tour
    Land = zeros(n)
    T = fill(Inf, n, n)
    M = fill(-99, n, n)
    @inbounds for i in 1:n-1
        @inbounds for j in i+1:n
            if j == i + 1
                T[i, j] = Ct[r[i], r[j]]
                M[r[i], r[j]] = -1
            else
                @inbounds for k in i+1:j-1
                    Tk1 = Cd[r[i], r[k]] + Cd[r[k], r[j]] + sR
                    if Tk1 <= flying_range
                        Tk2 = sum([Ct[r[l], r[l+1]] for l in i:k-2]) +
                              Ct[r[k-1], r[k+1]] +
                              sum([Ct[r[l], r[l+1]] for l in k+1:j-1]) + sR + sL * Land[i]
                        if Tk2 <= flying_range
                            Tk = max(Tk1, Tk2)
                            if Tk < T[i, j]
                                T[i, j] = Tk
                                M[r[i], r[j]] = r[k]
                                Land[j] = 1
                            end
                        end
                    end
                end

            end
        end
    end

    V = zeros(n)
    P = fill(-1, n)

    V[1] = 0.0
    @inbounds for i in 2:n
        VV = [V[k] + T[k, i] for k in 1:i-1]
        V[i] = minimum(VV)
        P[i] = r[argmin(VV)]
    end

    # Retrieving solutions.
    combined_nodes = Int[]
    current_idx = n
    current = r[current_idx]
    while current != -1
        pushfirst!(combined_nodes, current)
        current_idx = findfirst(x -> x == current, r)
        current = P[current_idx]
    end
    drone_only_nodes = Int[]
    drone_route = Int[]
    push!(drone_route, combined_nodes[1])
    @assert combined_nodes[1] == r[1]
    @inbounds for i in 1:length(combined_nodes)-1
        j1 = combined_nodes[i]
        j2 = combined_nodes[i+1]
        if M[j1, j2] != -1
            push!(drone_only_nodes, M[j1, j2])
            push!(drone_route, M[j1, j2])
        end
        push!(drone_route, j2)
    end
    truck_route = setdiff(initial_tour, drone_only_nodes)

    # obj_val = objective_value(truck_route, drone_route, Ct, Cd)
    final_time = V[end]

    #     @assert isapprox(obj_val, final_time)
    c = truck_route .- 1
    d = drone_route .- 1
    don = drone_only_nodes .- 1
    @inbounds for dnode in don
        tnode = d[findfirst(x -> x == dnode, d)-1]
        insert!(c, findfirst(x -> x == tnode, c) + 1, -dnode)
    end
    return final_time, c
end

function Build_Initial_chromosome(TT::Matrix{Float64}, DD::Matrix{Float64}, n_nodes::Int64, flying_range::Float64, sR::Float64, sL::Float64)
    tour = find_initial_TSP_tour(TT, n_nodes)
    pushfirst!(tour, 1)
    push!(tour, n_nodes + 2)
    f, c = exact_p(tour, TT, DD, flying_range, sR, sL)
    deleteat!(c, [1, n_nodes + 2])
    return c
end

function Change_initial(c::Vector{Int64})
    n_nodes = length(c)
    cc = copy(c)
    if rand() < 0.4
        @inbounds for i = 1:length(c)-1
            r1 = rand()
            if r1 < 0.1
                cc[i] = -cc[i]
            elseif r1 < 0.2
                cc[i+1], cc[i] = cc[i], cc[i+1]
            end
        end
    else
        idx1 = rand(1:n_nodes)
        idx2 = rand(1:n_nodes)
        if idx1 > idx2
            temp = idx1
            idx1 = idx2
            idx2 = temp
        end
        r = rand()
        if r < 0.3
            cc[idx1:idx2] = reverse(cc[idx1:idx2])
        elseif r < 0.7
            cc[idx1:idx2] = -cc[idx1:idx2]
        else
            cc[idx1:idx2] = shuffle(cc[idx1:idx2])
        end
    end
    return cc
end

function Generate_initial_population(mu::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64},
    n_nodes::Int64, flying_range::Float64, penaltyR::Float64, penaltyM::Float64, sR::Float64, sL::Float64, 
    initial_chrm::Vector{Int64}, problem_type::problem)
    t1 = time()
    Population = Chromosome[]
    # chrm = Build_Initial_chromosome(TT, DD, n_nodes, flying_range, sR, sL)
    chrm = initial_chrm
    f, LLnodes, Real_LLnodes = find_fitness(chrm, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
    push!(Population, Chromosome(chrm, LLnodes, Real_LLnodes, f, 'F', 0.0))
    count_f = 1
    count_infR = 0
    count_infM = 0
    if flying_range == Inf 
        count_infR = mu
    end
    while true
        if count_infM == mu
            if count_infR == mu
                if count_f == mu
                    break
                else
                    cr = Change_initial(chrm)
                    if !Is_feasibleM(cr)
                        cr = make_feasibleM(cr)
                    end
                    violating_drones = Is_feasibleR(cr, DD, TT, dEligible, flying_range, sR, sL, problem_type)
                    if length(violating_drones) > 0
                        cr = make_feasibleR(cr, violating_drones)
                    end
                    f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                    push!(Population, Chromosome(cr, LLnodes, Real_LLnodes, f, 'F', 0.0))
                    count_f += 1
                end
            else
                cr = Change_initial(chrm)
                if !Is_feasibleM(cr)
                    cr = make_feasibleM(cr)
                end
                violating_drones = Is_feasibleR(cr, DD, TT, dEligible, flying_range, sR, sL, problem_type)
                if length(violating_drones) == 0
                    if count_f < mu
                        f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                        push!(Population, Chromosome(cr, LLnodes, Real_LLnodes, f, 'F', 0.0))
                        count_f += 1
                    end
                else
                    f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
                    push!(Population, Chromosome(cr, LLnodes, Real_LLnodes, f, 'R', 0.0))
                    count_infR += 1
                end
            end
        else
            cr = Change_initial(chrm)
            if Is_feasibleM(cr)
                violating_drones = Is_feasibleR(cr, DD, TT, dEligible, flying_range, sR, sL, problem_type)
                if length(violating_drones) == 0
                    if count_f < mu
                        f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                        push!(Population, Chromosome(cr, LLnodes, Real_LLnodes, f, 'F', 0.0))
                        count_f += 1
                    end
                elseif count_infR < mu
                    f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
                    push!(Population, Chromosome(cr, LLnodes, Real_LLnodes, f, 'R', 0.0))
                    count_infR += 1
                end
            else
                f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'M', problem_type)
                push!(Population, Chromosome(cr, LLnodes, Real_LLnodes, f, 'M', 0.0))
                count_infM += 1
            end
        end
    end
    sort!(Population, by=x -> x.fitness)
    best_f = best_objective(Population)
    t2 = time()
    # println("Generation 1", " the best objective is: ", round(best_f, digits=2), " and it took ", round(t2 - t1, digits=2), " seconds.")
    return Population, best_f
end

function travel_cost(path::Vector{Int}, C::Matrix{T}) where {T}
    # @show path
    sum = zero(T)
    @inbounds for i in 1:length(path)-1
        sum += C[path[i], path[i+1]]
    end
    return sum
end

function objective_value(truck_route::Vector{Int64}, drone_route::Vector{Int64}, Ct::Matrix{Float64}, Cd::Matrix{Float64})
    combined_nodes = intersect(truck_route, drone_route)
    obj_val = 0.0
    for i in 1:length(combined_nodes)-1
        j1 = combined_nodes[i]
        j2 = combined_nodes[i+1]

        t_idx1 = findfirst(x -> x == j1, truck_route)
        t_idx2 = findfirst(x -> x == j2, truck_route)
        t_cost = travel_cost(truck_route[t_idx1:t_idx2], Ct)

        d_idx1 = findfirst(x -> x == j1, drone_route)
        d_idx2 = findfirst(x -> x == j2, drone_route)
        d_cost = travel_cost(drone_route[d_idx1:d_idx2], Cd)

        obj_val += max(t_cost, d_cost)
    end
    return obj_val
end

function Diversify(Population::Vector{Chromosome}, mu::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64},
    n_nodes::Int64, flying_range::Float64, sR::Float64, sL::Float64, penaltyR::Float64, penaltyM::Float64, 
    initial_chrm::Vector{Int64}, problem_type::problem)
    # chrm = Build_Initial_chromosome(TT, DD, n_nodes, flying_range, sR, sL)
    chrm = initial_chrm
    n_best = Int(round(0.3 * mu))
    feas_count = 0
    infeasR_count = 0
    infeasM_count = 0
    if flying_range == Inf
        infeasR_count = n_best
    end
    @inbounds for i = 1:length(Population)
        if Population[i].feasible == 'F'
            if feas_count < n_best
                feas_count += 1
            else
                cr = Change_initial(chrm)
                if !Is_feasibleM(cr)
                    cr = make_feasibleM(cr)
                end
                violating_drones = Is_feasibleR(cr, DD, TT, dEligible, flying_range, sR, sL, problem_type)
                if length(violating_drones) > 0
                    cr = make_feasibleR(cr, violating_drones)
                end
                f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                Population[i].genes = cr
                Population[i].fitness = f
                Population[i].LLnodes = LLnodes
                Population[i].Real_LLnodes = Real_LLnodes
            end
        elseif Population[i].feasible == 'R'
            if infeasR_count < n_best
                infeasR_count += 1
            else
                cr = Change_initial(chrm)
                if !Is_feasibleM(cr)
                    cr = make_feasibleM(cr)
                end
                violating_drones = Is_feasibleR(cr, DD, TT, dEligible, flying_range, sR, sL, problem_type)
                if length(violating_drones) > 0
                    f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
                    Population[i].genes = cr
                    Population[i].fitness = f
                    Population[i].LLnodes = LLnodes
                    Population[i].Real_LLnodes = Real_LLnodes
                end
            end
        else
            if infeasM_count < n_best
                infeasM_count += 1
            else
                cr = Change_initial(chrm)
                if !Is_feasibleM(cr)
                    f, LLnodes, Real_LLnodes = find_fitness(cr, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'M', problem_type)
                    Population[i].genes = cr
                    Population[i].fitness = f
                    Population[i].LLnodes = LLnodes
                    Population[i].Real_LLnodes = Real_LLnodes
                end
            end
        end
    end
end


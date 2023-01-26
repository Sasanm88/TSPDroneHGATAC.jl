using StatsBase

mutable struct Chromosome
    genes::Vector{Int64}      #Represents the sequence and the types of nodes, positive number means the node is visited by truck, while negative is for drone 
    LLnodes::Vector{Int64}    #Represents the location of luanch and land nodes so along with genes it gives us the tour
    fitness::Float64          #This is the Makespan of the given sequence found by JOIN algorithm (smaller, better)
    feasible::Char            #'F': feasible, 'R':infeasible type range, 'M':infeasible type multiple drone nodes   (In TSPD with unlimited flight range, there would be no 'R')
    power::Float64            # This is the fitness after considering the diversity contribution (smaller, better)
end


function Creat_Random_Cromosome(n_nodes::Int64)
    chromosome = shuffle!([i for i = 1:n_nodes])
    @inbounds for i = 1:n_nodes
        if rand() < 0.6
            chromosome[i] = -chromosome[i]
        end
    end
    chromosome
end


function Parent_Selection_TS(Population::Vector{Chromosome}, k::Int, popsize::Int)  #Tournament Selection
    idx = sample(1:popsize, k, replace=false)
    return Population[idx[argmin(idx)]]
end

function Select_parents(Population::Vector{Chromosome}, k_tournament::Int64, popsize::Int64)
    return Parent_Selection_TS(Population, k_tournament, popsize), Parent_Selection_TS(Population, k_tournament, popsize)
end

function Reproduce(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64) 
    r = rand(1:5)

    if r == 1
        return Crossover_OX1(parent1, parent2, n_nodes)  
    elseif r == 2
        return Crossover_OX2(parent1, parent2, n_nodes)  
    elseif r == 3
        return Crossover_DX1(parent1, parent2, n_nodes)  
    elseif r == 4
        return Crossover_DX2(parent1, parent2, n_nodes)  
    else
        return Crossover_DX3(parent1, parent2, n_nodes)   
    end
end

function mutate(c::Vector{Int64})
    r = rand()
    if r < 0.1
        if rand() < 0.5
            return sign_mutation(c)
        else
            return tour_mutation(c)
        end
    else
        return c
    end
end

function find_difference(c1::Vector{Int64}, c2::Vector{Int64})  #range between zero and 1, zero when two chromosomes are exactly the same, 1 when all genes are different
    diff1 = 0
    diff2 = 0
    c3 = reverse(c2)
    @inbounds for i = 1:length(c1)
        if c1[i] != c2[i]
            diff1 += 1
        end
        if c1[i] != c3[i]
            diff2 += 1
        end
    end
    return min(diff1, diff2) / length(c1)
end

function Sort_based_on_power(Population::Vector{Chromosome})
    popsize = length(Population)
    diff1 = 0.0
    diff2 = 0.0
    @inbounds for i = 1:popsize
        if i == 1
            diff1 = find_difference(Population[1].genes, Population[2].genes)
            Population[i].power = Population[i].fitness * 0.8^diff1 #(2-diff1)
        elseif i == popsize
            Population[i].power = Population[i].fitness * 0.8^diff1 #(2-diff1)
        else
            diff2 = find_difference(Population[i].genes, Population[i+1].genes)
            Population[i].power = Population[i].fitness * 0.8^((diff1 + diff2) / 2) #(2-(diff1+diff2)/2)
            diff1 = diff2
        end
    end
    sort!(Population, by=x -> x.power)
end

function process_child(Population::Vector{Chromosome}, child::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64},
    dEligible::Vector{Int64}, flying_range::Float64, penaltyR::Float64, penaltyM::Float64, fractionFeasibleLoad::Float64,
    fractionInFeasibleRLoad::Float64, feas_count::Int64, InfR_count::Int64, InfM_count::Int64, sR::Int64, sL::Int64,
     ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, problem_type::String)
    n_nodes = length(child)
    if Is_feasibleM(child)
        violating_drones = Is_feasibleR(child, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            fractionFeasibleLoad = fractionFeasibleLoad * 0.99 + 0.01
            fractionInFeasibleRLoad = fractionInFeasibleRLoad * 0.99
        else
            fractionFeasibleLoad = fractionFeasibleLoad * 0.99
            fractionInFeasibleRLoad = fractionInFeasibleRLoad * 0.99 + 0.01
        end
        fitness, LLnodes = find_fitness(child, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        offspring = Chromosome(child, LLnodes, fitness, 'R', 0.0)
        offspring = Improve_chromosome(offspring, n_nodes, TT, DD, dEligible, ClosenessT, ClosenessD, flying_range, sR, sL, penaltyR, penaltyM, problem_type)
        violating_drones = Is_feasibleR(offspring.genes, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            offspring.feasible = 'F'
            push!(Population, offspring)
            feas_count += 1
        elseif rand() < 0.5
            child = make_feasibleR(offspring.genes, violating_drones)
            fitness, LLnodes = find_fitness(offspring.genes, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            offspring = Chromosome(child, LLnodes, fitness, 'F', 0.0)
            push!(Population, offspring)
            feas_count += 1
        else
            push!(Population, offspring)
            InfR_count += 1
        end
    elseif rand() < 0.5
        fractionFeasibleLoad = fractionFeasibleLoad * 0.99
        fractionInFeasibleRLoad = fractionInFeasibleRLoad * 0.99
        child = make_feasibleM(child)
        fitness, LLnodes = find_fitness(child, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        offspring = Chromosome(child, LLnodes, fitness, 'R', 0.0)
        offspring = Improve_chromosome(offspring, n_nodes, TT, DD, dEligible, ClosenessT, ClosenessD, flying_range, sR, sL, penaltyR, penaltyM, problem_type)
        violating_drones = Is_feasibleR(offspring.genes, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            offspring.feasible = 'F'
            push!(Population, offspring)
            feas_count += 1
        elseif rand() < 0.5
            child = make_feasibleR(offspring.genes, violating_drones)
            fitness, LLnodes = find_fitness(offspring.genes, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            offspring = Chromosome(child, LLnodes, fitness, 'F', 0.0)
            push!(Population, offspring)
            feas_count += 1
        else
            push!(Population, offspring)
            InfR_count += 1
        end
    else
        fractionFeasibleLoad = fractionFeasibleLoad * 0.99
        fractionInFeasibleRLoad = fractionInFeasibleRLoad * 0.99
        fitness, LLnodes = find_fitness(child, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'M', problem_type)
        offspring = Chromosome(child, LLnodes, fitness, 'M', 0.0)
        push!(Population, offspring)
        if offspring.feasible == 'F'
            feas_count += 1
        elseif offspring.feasible == 'R'
            InfR_count += 1
        else
            InfM_count += 1
        end
    end
    return fractionFeasibleLoad, fractionInFeasibleRLoad, feas_count, InfR_count, InfM_count
end

function Survive(Population::Vector{Chromosome}, sigma::Int64, feasiblity::Char, feas_count::Int64, InfR_count::Int64, InfM_count::Int64)
    del_count = 0
    del_idx = Vector{Int64}()
    @inbounds for i = 1:length(Population)-1
        if del_count == sigma
            break
        end
        if Population[i].feasible == feasiblity
            @inbounds for j = i+1:length(Population)
                if del_count == sigma
                    break
                end
                if Population[j].feasible == feasiblity
                    if Population[i].genes == Population[j].genes
                        if !(j in del_idx)
                            push!(del_idx, j)
                            del_count += 1
                        end
                    end
                end
            end
        end
    end
    deleteat!(Population, sort(del_idx))
    del_idx = Vector{Int64}()
    last_index = length(Population)
    index = 0

    @inbounds while del_count < sigma
        i = last_index - index

        if Population[i].feasible == feasiblity
            push!(del_idx, i)
            del_count += 1
        end
        index += 1
    end
    deleteat!(Population, sort(del_idx))
    if feasiblity == 'F'
        feas_count -= sigma
    elseif feasiblity == 'R'
        InfR_count -= sigma
    else
        InfM_count -= sigma
    end
    return feas_count, InfR_count, InfM_count
end

function Perform_Survival_plan(Population::Vector{Chromosome}, mu::Int64, sigma::Int64, feas_count::Int64, InfR_count::Int64, InfM_count::Int64)
    if InfR_count >= mu + sigma
        feas_count, InfR_count, InfM_count = Survive(Population, sigma, 'R', feas_count, InfR_count, InfM_count)
    end
    if InfM_count >= mu + sigma
        feas_count, InfR_count, InfM_count = Survive(Population, sigma, 'M', feas_count, InfR_count, InfM_count)
    end
    if feas_count >= mu + sigma
        feas_count, InfR_count, InfM_count = Survive(Population, sigma, 'F', feas_count, InfR_count, InfM_count)
    end
    return feas_count, InfR_count, InfM_count
end

function Adjust_penalties(fractionFeasibleLoad::Float64, fractionInFeasibleRLoad::Float64, penaltyM::Float64,
     penaltyR::Float64, targetFeasible::Float64, flying_range::Float64)
    if fractionFeasibleLoad < targetFeasible - 0.05
        if fractionInFeasibleRLoad < 1 - fractionFeasibleLoad - fractionInFeasibleRLoad
            penaltyM = min(penaltyM * 1.1, 8)
        else
            penaltyR = min(penaltyR * 1.1, 5)
        end
    elseif fractionFeasibleLoad > targetFeasible + 0.05
        if flying_range == Inf
            penaltyM = max(penaltyM * 0.90909, 3)
            return penaltyM, penaltyR
        end
        if fractionInFeasibleRLoad < 1 - fractionFeasibleLoad - fractionInFeasibleRLoad
            penaltyR = max(penaltyR * 0.90909, 1.5)
        else
            penaltyM = max(penaltyM * 0.90909, 3)
        end
    end
    return penaltyM, penaltyR
end


function Generate_new_generation(TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, flying_range::Float64,
    popsize::Tuple{Int64,Int64}, k_tournament::Int64, targetFeasible::Float64, sR::Int, sL::Int, ClosenessT::Matrix{Int64},
    ClosenessD::Matrix{Int64}, initial_chrm::Vector{Int64}, Gen_num::Int64, feas_count::Int64, InfR_count::Int64, InfM_count::Int64,
    penaltyR::Float64, penaltyM::Float64, old_best::Float64, fractionFeasibleLoad::Float64, fractionInFeasibleRLoad::Float64,
    Population::Vector{Chromosome}, improve_count::Int64, problem_type::String)
    t1 = time()

    mu, sigma = popsize
    n_nodes = length(Population[1].genes)

    penaltyM, penaltyR = Adjust_penalties(fractionFeasibleLoad, fractionInFeasibleRLoad, penaltyM, penaltyR, targetFeasible, flying_range)


    if improve_count % 100 == 0
        Diversify(Population, mu, TT, DD, dEligible, n_nodes, "soft", flying_range, sR, sL, penaltyR, penaltyM, initial_chrm, problem_type)
    end
    Sort_based_on_power(Population)
    psize = length(Population)
    parent1, parent2 = Select_parents(Population, k_tournament, psize)

    child = Reproduce(parent1.genes, parent2.genes, n_nodes)
    child = mutate(child)

    fractionFeasibleLoad, fractionInFeasibleRLoad, feas_count, InfR_count, InfM_count = 
    process_child(Population, child, TT, DD, dEligible, flying_range, penaltyR, penaltyM, 
    fractionFeasibleLoad, fractionInFeasibleRLoad, feas_count, InfR_count, InfM_count, sR, sL, ClosenessT, ClosenessD, problem_type)

    sort!(Population, by=x -> x.fitness)
    if improve_count % 999 == 0
        Scape_local_optima(Population, TT, DD, dEligible, ClosenessT, ClosenessD, flying_range, sR, sL, penaltyR, Int(round(improve_count / 999)))
    end
    feas_count, InfR_count, InfM_count = Perform_Survival_plan(Population, mu, sigma, feas_count, InfR_count, InfM_count)

    new_best = best_objective(Population)
    if (old_best - new_best) / new_best > 0.0005
        old_best = new_best
        improve_count = 0
    else
        improve_count += 1
    end
    t2 = time()
    Gen_num += 1

    if Gen_num % 1000 == 0
        println("Generation ", Gen_num, " the best objective is: ", old_best, " Feasible fraction = ", round(fractionFeasibleLoad, digits=3), " penaltyR=", round(penaltyR, digits=2), " penaltyM=", round(penaltyM, digits=2))
    end
    return Gen_num, feas_count, InfR_count, InfM_count, penaltyR, penaltyM, old_best, fractionFeasibleLoad, fractionInFeasibleRLoad, Population, improve_count
end



function best_objective(Population::Vector{Chromosome})
    @inbounds for i in 1:length(Population)
        if Population[i].feasible == 'F'
            return Population[i].fitness
        end
    end
    return Inf # Is this correct?
end

function best_route(Population::Vector{Chromosome})
    @inbounds for i in 1:length(Population)
        if Population[i].feasible == 'F'
            @inbounds for j in Population[i].genes
                print(j, " ")
            end
            break
        end
    end
    # TODO: 
    # return ???
end

function Perform_Genetic_Algorithm(TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, h::Float64, popsize::Tuple{Int64,Int64},
    k_tournament::Int64, targetFeasible::Float64, sR::Int64, sL::Int64, num_iter::Int64, flying_range::Float64, initial_chrm::Vector{Int64}, problem_type::String)
    n_nodes = size(TT)[1] - 2
    t1 = time()
    ClosenessT, ClosenessD = Find_Closeness(TT, DD, h)
    mu, sigma = popsize
    improve_count = 0
    Gen_num = 1
    penaltyR = 2.0
    penaltyM = 2.0
    old_best = 0.0
    Population, old_best = Generate_initial_population(mu, TT, DD, dEligible, n_nodes, flying_range, penaltyR, penaltyM, sR, sL, initial_chrm, problem_type)
    count = 0
    feas_count = mu
    InfR_count = mu
    InfM_count = mu
    fractionFeasibleLoad = 0.2
    fractionInFeasibleRLoad = 0.4


    # TODO: 
    # static_info and dynamic_info should better be struct
    # Please create 
    # struct StaticInfo  and struct DynamicInfo

    @inbounds while improve_count < num_iter
        Gen_num, feas_count, InfR_count, InfM_count, penaltyR, penaltyM, old_best, fractionFeasibleLoad, fractionInFeasibleRLoad,
        Population, improve_count = Generate_new_generation(TT, DD, dEligible, flying_range, popsize, k_tournament, targetFeasible, sR, sL,
            ClosenessT, ClosenessD, initial_chrm, Gen_num, feas_count, InfR_count, InfM_count, penaltyR, penaltyM, old_best, fractionFeasibleLoad, fractionInFeasibleRLoad,
            Population, improve_count, problem_type)
        count += 1
    end
    t2 = time()

    println("The best objective achieved in ", Gen_num, " generations is: ", best_objective(Population), " and it took ", t2 - t1, " seconds.")
    println("And the best route is: ")
    best_route(Population)
    return Population
end



function run_GA(depot::Tuple{Float64}, Customers::Matrix{Float64}, tspeed::Int, dspeed::Int, h::Float64, flying_range::Float64, sR::Int, sL::Int,
     num_runs::Int64, num_generations::Int64, problem_type::String)
     n_nodes = size(T)[1] - 2
     T, D = Calculate_duration_matrices(tspeed, dspeed, depot, Customers)
    popsize = (15, 25)  #(mu,sigma)
    k_tournament = 5
    targetFeasible = 0.2

    objs = Float64[]
    times = Float64[]

    best_obj_all_time = Inf
    worst_obj = 0.0

    @inbounds for i in 1:num_runs
        initial_chrm = Build_Initial_chromosome(T, D, n_nodes, flying_range, sR, sL)
        t1 = time()
        println("Run ", i, ":")
        P = Perform_Genetic_Algorithm(T, D, dEligible, h, popsize, k_tournament, targetFeasible, sR, sL, 
        num_generations, flying_range, initial_chrm, problem_type)
        t2 = time()
        current_best = best_objective(P)
        push!(objs, current_best)
        push!(times, t2 - t1)
        println()
        if current_best < best_obj_all_time
            best_obj_all_time = current_best
        end
        if current_best > worst_obj
            worst_obj = current_best
        end
    end

    println()

    println("Genetic Algorithm results for solving instance ", string(sample_name), " for ", num_runs, " times:")
    println()
    println("The average objective found for this instance in ", num_runs, " runs is: ", round(sum(objs) / num_runs, digits=2),
        " ,the best found was: ", round(best_obj_all_time, digits=2), " and the worst found was: ", round(worst_obj, digits=2))
    println("The average time for each runs is: ", round(sum(times) / num_runs, digits=2))
    return worst_obj, best_obj_all_time, mean(objs), mean(times)
end


function Find_Closeness(TT::Matrix{Float64}, DD::Matrix{Float64}, h::Float64)
    n_nodes = size(TT)[1] - 2
    num = Int(ceil(h * n_nodes))
    ClosenessT = zeros(Int, n_nodes, num)
    ClosenessD = zeros(Int, n_nodes, num)
    @inbounds for i = 2:n_nodes+1
        a = copy(TT[i, 2:n_nodes+1])
        b = sortperm(a)
        ClosenessT[i-1, :] = b[2:num+1]
    end
    @inbounds for i = 2:n_nodes+1
        a = copy(DD[i, 2:n_nodes+1])
        b = sortperm(a)
        ClosenessD[i-1, :] = b[2:num+1]
    end
    return ClosenessT, ClosenessD
end

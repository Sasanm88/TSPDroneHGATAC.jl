using StatsBase

function N1(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)   #choose one free truck node and puts it randomly in truck tour     #OK
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    if length(free_tnodes) > 0
        remove_from = free_tnodes[rand(1:length(free_tnodes))]
        # TODO: specify Type
        tnodes_subs = Int[]
        @inbounds for i in tnodes_loc
            if c[i] in ClosenessT[c[remove_from], :]
                push!(tnodes_subs, i)
            end
        end
        if length(tnodes_subs) == 0
            return chrm
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs))]
        temp = c[remove_from]
        deleteat!(c, remove_from)
        insert!(c, insert_to, temp)

        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N2(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem) #choose two consecutive free truck nodes and puts them randomly in truck tour    #OK
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    consecutive_tnodes = Int[]
    f1 = 0.0
    f2 = 0.0
    @inbounds for i = 1:length(free_tnodes)-1
        if free_tnodes[i+1] == free_tnodes[i] + 1
            push!(consecutive_tnodes, free_tnodes[i])
        end
    end
    if length(consecutive_tnodes) > 0
        remove_from = consecutive_tnodes[rand(1:length(consecutive_tnodes))]
        tnodes_subs = Int[]
        @inbounds for i in tnodes_loc
            if i != remove_from && i != remove_from + 1
                if c[i] in ClosenessT[c[remove_from], :] || c[i] in ClosenessT[c[remove_from+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end
        if length(tnodes_subs) <= 1
            return chrm
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs)-1)]
        temp1 = c[remove_from]
        temp2 = c[remove_from+1]
        c1 = copy(c)
        deleteat!(c1, [remove_from, remove_from + 1])
        insert!(c1, insert_to, temp1)
        insert!(c1, insert_to, temp2)
        violating_drones = Is_feasibleR(c1, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f1, llc1, rllc1 = find_fitness(c1, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f1 < chrm.fitness
                    chrm.genes = c1
                    chrm.fitness = f1
                    chrm.LLnodes = llc1
                    chrm.Real_LLnodes = rllc1
                end
            end
        else
            f1, llc1, rllc1 = find_fitness(c1, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f1 < chrm.fitness
                chrm.genes = c1
                chrm.fitness = f1
                chrm.LLnodes = llc1
                chrm.Real_LLnodes = rllc1
            end
        end
        c2 = copy(c)
        deleteat!(c2, [remove_from, remove_from + 1])
        insert!(c2, insert_to, temp2)
        insert!(c2, insert_to, temp1)
        violating_drones = Is_feasibleR(c2, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f2, llc2, rllc2 = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f2 < chrm.fitness
                    chrm.genes = c2
                    chrm.fitness = f2
                    chrm.LLnodes = llc2
                    chrm.Real_LLnodes = rllc2
                end
            end
        else
            f2, llc2, rllc2 = find_fitness(c2, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f2 < chrm.fitness
                chrm.genes = c2
                chrm.fitness = f2
                chrm.LLnodes = llc2
                chrm.Real_LLnodes = rllc2
            end
        end

    end
    return chrm
end

function N3(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)   #swap two truck nodes 
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    idx1 = tnodes_loc[rand(1:length(tnodes_loc))]
    tnodes_subs = Int[]
    @inbounds for i in tnodes_loc
        if c[i] in ClosenessT[c[idx1], :]
            push!(tnodes_subs, i)
        end
    end
    if length(tnodes_subs) == 0
        return chrm
    end
    idx2 = tnodes_subs[rand(1:length(tnodes_subs))]
    c[idx1], c[idx2] = c[idx2], c[idx1]
    violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
    if chrm.feasible == 'F'
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N4(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)   #swap two consecutive free truck nodes with another truck node 
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    indices = Int[]
    @inbounds for node in free_tnodes
        if node > 1
            if c[node-1] > 0
                push!(indices, node)
            end
        end
    end
    if length(indices) > 0
        remove_from = indices[rand(1:length(indices))] - 1
        tnodes_subs = Int[]
        @inbounds for i in tnodes_loc
            if i != remove_from && i != remove_from + 1
                if c[i] in ClosenessT[c[remove_from], :] || c[i] in ClosenessT[c[remove_from+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end

        if length(tnodes_subs) == 0
            return chrm
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs))]
        c[remove_from], c[insert_to] = c[insert_to], c[remove_from]
        temp = c[remove_from+1]
        deleteat!(c, remove_from + 1)
        insert!(c, insert_to, temp)
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N5(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)   #swap two consecutive truck nodes with another truck node 
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    indices = Int[]
    @inbounds for i = 1:length(c)-1
        if c[i] > 0
            if c[i+1] > 0
                push!(indices, i)
            end
        end
    end
    if length(indices) > 0
        remove_from = indices[rand(1:length(indices))]
        tnodes_subs = Int[]
        @inbounds for i in tnodes_loc
            if i != remove_from && i != remove_from + 1
                if c[i] in ClosenessT[c[remove_from], :] || c[i] in ClosenessT[c[remove_from+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end

        if length(tnodes_subs) == 0
            return chrm
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs))]
        c[remove_from], c[insert_to] = c[insert_to], c[remove_from]
        temp = c[remove_from+1]
        deleteat!(c, remove_from + 1)
        insert!(c, insert_to, temp)
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N6(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)        #Truck swap 2-2
    c = copy(chrm.genes)
    indices = Int[]
    i = 1
    @inbounds while i < n_nodes
        if c[i] > 0
            if c[i+1] > 0
                push!(indices, i)
            end
            i = i + 2
        else
            i = i + 1
        end
    end
    if length(indices) > 1
        i1 = indices[rand(1:length(indices))]
        tnodes_subs = Int[]
        @inbounds for i in indices
            if i != i1
                if c[i] in ClosenessT[c[i1], :] || c[i] in ClosenessT[c[i1+1], :] || c[i+1] in ClosenessT[c[i1], :] || c[i+1] in ClosenessT[c[i1+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end

        if length(tnodes_subs) == 0
            return chrm
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        temp1 = c[i1]
        temp2 = c[i1+1]
        c[i1], c[i1+1] = c[i2], c[i2+1]
        c[i2] = temp1
        c[i2+1] = temp2
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N7(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)      #Truck swap 2-2
    c = copy(chrm.genes)
    indices = Int[]
    i = 1
    @inbounds while i < n_nodes
        if c[i] > 0
            if c[i+1] > 0
                push!(indices, i)
            end
            i = i + 2
        else
            i = i + 1
        end
    end
    if length(indices) > 1
        i1 = indices[rand(1:length(indices))]
        tnodes_subs = Int[]
        @inbounds for i in indices
            if i != i1
                if c[i] in ClosenessT[c[i1], :] || c[i] in ClosenessT[c[i1+1], :] || c[i+1] in ClosenessT[c[i1], :] || c[i+1] in ClosenessT[c[i1+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end
        if length(tnodes_subs) == 0
            return chrm
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        c[i2], c[i2+1] = c[i2+1], c[i2]
        c[i1+1], c[i2] = c[i2], c[i1+1]
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end


function N8(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)       #Truck swap 2-2
    c = copy(chrm.genes)
    indices = Int[]
    i = 1
    @inbounds while i < n_nodes
        if c[i] > 0
            if c[i+1] > 0
                push!(indices, i)
            end
            i = i + 2
        else
            i = i + 1
        end
    end
    if length(indices) > 1
        i1 = indices[rand(1:length(indices))]
        tnodes_subs = Int[]
        @inbounds for i in indices
            if i != i1
                if c[i] in ClosenessT[c[i1], :] || c[i] in ClosenessT[c[i1+1], :] || c[i+1] in ClosenessT[c[i1], :] || c[i+1] in ClosenessT[c[i1+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end
        if length(tnodes_subs) == 0
            return chrm
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        c[i1+1], c[i2+1] = c[i2+1], c[i1+1]
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N9(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)        #Truck swap 2-2
    c = copy(chrm.genes)
    indices = Int[]
    i = 1
    @inbounds while i < n_nodes
        if c[i] > 0
            if c[i+1] > 0
                push!(indices, i)
            end
            i = i + 2
        else
            i = i + 1
        end
    end
    if length(indices) > 1
        i1 = indices[rand(1:length(indices))]
        tnodes_subs = Int[]
        @inbounds for i in indices
            if i != i1
                if c[i] in ClosenessT[c[i1], :] || c[i] in ClosenessT[c[i1+1], :] || c[i+1] in ClosenessT[c[i1], :] || c[i+1] in ClosenessT[c[i1+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end
        if length(tnodes_subs) == 0
            return chrm
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        temp1 = c[i1]
        temp2 = c[i1+1]
        c[i1], c[i1+1] = c[i2+1], c[i2]
        c[i2] = temp1
        c[i2+1] = temp2
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end

    return chrm
end


function N10(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #swap a drone node with either one of randezvous nodes or any truck node between
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return chrm
    end
    dnode_idx = rand(1:length(dnodes_loc))
    tnode_indices = Int[]
    start = 0
    finish = 0
    if dnode_idx > 1
        @inbounds for i = 1:length(chrm.LLnodes)
            if chrm.LLnodes[i] > dnodes_loc[dnode_idx]
                if i > 1
                    start = chrm.LLnodes[i-1]
                else
                    start = 0
                end
                finish = chrm.LLnodes[i]
                break
            end
        end
    else
        if chrm.LLnodes[1] < dnodes_loc[dnode_idx]
            start = chrm.LLnodes[1]
            if length(chrm.LLnodes) == 1
                finish = n_nodes
            else
                finish = chrm.LLnodes[2]
            end
        else
            start = 0
            finish = chrm.LLnodes[1]
        end
    end

    if start == 0
        if finish == 0
            tnode_indices = tnodes_loc[findfirst(x -> x == chrm.LLnodes[length(chrm.LLnodes)], tnodes_loc):length(tnodes_loc)]
        else
            tnode_indices = tnodes_loc[1:findfirst(x -> x == finish, tnodes_loc)]
        end
    elseif finish == 0
        tnode_indices = tnodes_loc[findfirst(x -> x == start, tnodes_loc):length(tnodes_loc)]
    else
        tnode_indices = tnodes_loc[findfirst(x -> x == start, tnodes_loc):findfirst(x -> x == finish, tnodes_loc)]
    end
    tnode = tnode_indices[rand(1:length(tnode_indices))]
    c = copy(chrm.genes)
    c[tnode], c[dnodes_loc[dnode_idx]] = -c[dnodes_loc[dnode_idx]], -c[tnode]
    violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
    if chrm.feasible == 'F'
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N11(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #in drone route <i,j,k> swap i and j
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return chrm
    end
    dnode_idx = rand(1:length(dnodes_loc))
    tnode_indices = Int[]
    start = 0
    finish = 0
    if dnode_idx > 1
        @inbounds for i = 1:length(chrm.LLnodes)
            if chrm.LLnodes[i] > dnodes_loc[dnode_idx]
                if i > 1
                    start = chrm.LLnodes[i-1]
                else
                    start = 0
                end
                finish = chrm.LLnodes[i]
                break
            end
        end
    else
        if chrm.LLnodes[1] < dnodes_loc[dnode_idx] && length(chrm.LLnodes) > 1
            start = chrm.LLnodes[1]
            finish = chrm.LLnodes[2]
        else
            start = 0
            finish = chrm.LLnodes[1]
        end
    end
    if start == 0 && finish == 0
        start = chrm.LLnodes[length(chrm.LLnodes)]
    end
    if start > 0
        c[start], c[dnodes_loc[dnode_idx]] = -c[dnodes_loc[dnode_idx]], -c[start]
    end
    violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
    if chrm.feasible == 'F'
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N12(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #in drone route <i,j,k> swap k and j
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return chrm
    end
    dnode_idx = rand(1:length(dnodes_loc))
    tnode_indices = Int[]
    start = 0
    finish = 0
    if dnode_idx > 1
        @inbounds for i = 1:length(chrm.LLnodes)
            if chrm.LLnodes[i] > dnodes_loc[dnode_idx]
                if i > 1
                    start = chrm.LLnodes[i-1]
                else
                    start = 0
                end
                finish = chrm.LLnodes[i]
                break
            end
        end
    else
        if chrm.LLnodes[1] < dnodes_loc[dnode_idx] && length(chrm.LLnodes) > 1
            start = chrm.LLnodes[1]
            finish = chrm.LLnodes[2]
        else
            start = 0
            finish = chrm.LLnodes[1]
        end
    end

    if finish > 0
        c[finish], c[dnodes_loc[dnode_idx]] = -c[dnodes_loc[dnode_idx]], -c[finish]
    end
    violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
    if chrm.feasible == 'F'
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N13(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #in drone route <i,j,k> swap k and i
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return chrm
    end
    dnode_idx = rand(1:length(dnodes_loc))
    tnode_indices = Int[]
    start = 0
    finish = 0
    if dnode_idx > 1
        @inbounds for i = 1:length(chrm.LLnodes)
            if chrm.LLnodes[i] > dnodes_loc[dnode_idx]
                if i > 1
                    start = chrm.LLnodes[i-1]
                else
                    start = 0
                end
                finish = chrm.LLnodes[i]
                break
            end
        end
    else
        if chrm.LLnodes[1] < dnodes_loc[dnode_idx] && length(chrm.LLnodes) > 1
            start = chrm.LLnodes[1]
            finish = chrm.LLnodes[2]
        else
            start = 0
            finish = chrm.LLnodes[1]
        end
    end

    if finish > 0 && start > 0
        c[finish], c[start] = c[start], c[finish]
    end
    violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
    if chrm.feasible == 'F'
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N14(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #remove one free truck node and make it a drone node between two consecutive truck nodes
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    indices = Int[]

    if length(free_tnodes) > 0
        remove_from = free_tnodes[rand(1:length(free_tnodes))]
        setdiff!(tnodes_loc, remove_from)
        i = 1
        @inbounds while i < length(tnodes_loc) - 1
            if tnodes_loc[i+1] == tnodes_loc[i] + 1
                if c[tnodes_loc[i]] in ClosenessD[c[remove_from], :] || c[tnodes_loc[i+1]] in ClosenessD[c[remove_from], :]
                    push!(indices, tnodes_loc[i])
                end
            end
            i += 1
        end
        if length(indices) > 0
            insert_to = indices[rand(1:length(indices))]
            temp = c[remove_from]

            if remove_from < insert_to
                deleteat!(c, remove_from)
                insert!(c, insert_to, -temp)
            elseif remove_from > insert_to
                deleteat!(c, remove_from)
                insert!(c, insert_to + 1, -temp)
            end
            violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
            if chrm.feasible == 'F'
                if length(violating_drones) == 0
                    f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                    if f < chrm.fitness
                        chrm.genes = c
                        chrm.fitness = f
                        chrm.LLnodes = llc
                        chrm.Real_LLnodes = rllc
                    end
                end
            else
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        end
    end
    return chrm
end


function N14p(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #remove one middle truck node and make it a drone node
    c = copy(chrm.genes)
    middle_tnodes = Int[]
    one_positive = false
    two_positive = false
    i = 1
    @inbounds while i < n_nodes - 1
        if two_positive
            if c[i] > 0
                push!(middle_tnodes, i - 1)
            else
                one_positive = false
                two_positive = false
            end
        elseif one_positive
            if c[i] > 0
                two_positive = true
            else
                one_positive = false
            end
        else
            if c[i] > 0
                one_positive = true
            end
        end
        i = i + 1
    end

    if length(middle_tnodes) > 0
        tnode = middle_tnodes[rand(1:length(middle_tnodes))]
        c[tnode] = -c[tnode]
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end


function N15(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #remove a drone node and insert it randomly as a truck node between two consequtive truck nodes
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    dnode = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    i = 1
    @inbounds while i < n_nodes - 1
        if c[i] > 0
            if c[i+1] > 0
                if c[i] in ClosenessT[-c[dnode], :] || c[i+1] in ClosenessT[-c[dnode], :]
                    push!(indices, i)
                end
            end
            i = i + 2
        else
            i = i + 1
        end
    end

    if length(indices) > 0
        insert_to = indices[rand(1:length(indices))]
        temp = c[dnode]

        if dnode < insert_to
            deleteat!(c, dnode)
            insert!(c, insert_to, -temp)
        elseif dnode > insert_to
            deleteat!(c, dnode)
            insert!(c, insert_to + 1, -temp)
        end
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N19(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #remove a drone node and insert it randomly as a drone node between two consequtive truck nodes
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    dnode = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    i = 1
    @inbounds while i < n_nodes - 1
        if c[i] > 0
            if c[i+1] > 0
                if c[i] in ClosenessD[-c[dnode], :] || c[i+1] in ClosenessD[-c[dnode], :]
                    push!(indices, i)
                end
            end
            i = i + 2
        else
            i = i + 1
        end
    end

    if length(indices) > 0
        insert_to = indices[rand(1:length(indices))]
        temp = c[dnode]

        if dnode < insert_to
            deleteat!(c, dnode)
            insert!(c, insert_to, temp)
        elseif dnode > insert_to
            deleteat!(c, dnode)
            insert!(c, insert_to + 1, temp)
        end
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if chrm.feasible == 'F'
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N16(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #swap two drone nodes
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) > 1
        dnode1 = dnodes_loc[rand(1:length(dnodes_loc))]
        d_subs = Int[]
        @inbounds for i in dnodes_loc
            if -c[i] in ClosenessD[-c[dnode1], :]
                push!(d_subs, i)
            end
        end
        if length(d_subs) > 0
            dnode2 = d_subs[rand(1:length(d_subs))]
            c[dnode1], c[dnode2] = c[dnode2], c[dnode1]
            if chrm.feasible == 'F'
                violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
                if length(violating_drones) == 0
                    f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                    if f < chrm.fitness
                        chrm.genes = c
                        chrm.fitness = f
                        chrm.LLnodes = llc
                        chrm.Real_LLnodes = rllc
                    end
                end
            else
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        end
    end
    return chrm
end

function N17(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)   #Swap a truck node with a drone node
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return chrm
    end
    idx2 = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    @inbounds for i = 1:n_nodes
        if c[i] > 0 && c[i] in ClosenessT[-c[idx2], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return chrm
    end
    idx1 = indices[rand(1:length(indices))]
    c[idx1], c[idx2] = -c[idx2], -c[idx1]
    if chrm.feasible == 'F'
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N18(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #swap two truck arcs
    c = copy(chrm.genes)
    tnodes_loc = findall(x -> x > 0, c)
    idx1 = tnodes_loc[rand(1:length(tnodes_loc))]
    indices = Int[]
    @inbounds for i in tnodes_loc
        if c[i] in ClosenessT[c[idx1], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return chrm
    end
    idx2 = indices[rand(1:length(indices))]
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    c[idx1:idx2] = c[idx2:-1:idx1]
    if chrm.feasible == 'F'
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function N20(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #swap two drone nodes and make them truck nodes
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) > 1
        dnode1 = dnodes_loc[rand(1:length(dnodes_loc))]
        indices = Int[]
        @inbounds for i in dnodes_loc
            if -c[i] in ClosenessT[-c[dnode1], :]
                push!(indices, i)
            end
        end
        if length(indices) == 0
            return chrm
        end
        dnode2 = indices[rand(1:length(indices))]
        c[dnode1], c[dnode2] = -c[dnode2], -c[dnode1]
        if chrm.feasible == 'F'
            violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < chrm.fitness
                    chrm.genes = c
                    chrm.fitness = f
                    chrm.LLnodes = llc
                    chrm.Real_LLnodes = rllc
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    end
    return chrm
end

function N21(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  #
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return chrm
    end
    d1Loc = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    @inbounds for i in dnodes_loc
        if -c[i] in ClosenessD[-c[d1Loc], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return chrm
    end
    d2Loc = indices[rand(1:length(indices))]
    d1 = c[d1Loc]
    c[d2Loc] = -c[d2Loc]
    if d1Loc < d2Loc
        deleteat!(c, d1Loc)
        insert!(c, d2Loc, d1)
    elseif d1Loc > d2Loc
        deleteat!(c, d1Loc)
        insert!(c, d2Loc + 1, d1)
    end
    if chrm.feasible == 'F'
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end


function N22(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)  ##remove a drone node d, find a random drone pair <i,j,k> change the j to truck and put d in <i,j,d,k> or <i,d,j,k>
    c = copy(chrm.genes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return chrm
    end
    d1Loc = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    @inbounds for i in dnodes_loc
        if -c[i] in ClosenessD[-c[d1Loc], :]
            push!(indices, i)
        elseif i > 1 && c[i-1] in ClosenessD[-c[d1Loc], :]
            push!(indices, i)
        elseif i < n_nodes && c[i+1] in ClosenessD[-c[d1Loc], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return chrm
    end
    d2Loc = indices[rand(1:length(indices))]
    d1 = c[d1Loc]
    c[d2Loc] = -c[d2Loc]
    if d1Loc < d2Loc
        deleteat!(c, d1Loc)
        insert!(c, d2Loc - 1, d1)
    elseif d1Loc > d2Loc
        deleteat!(c, d1Loc)
        insert!(c, d2Loc, d1)
    end
    if chrm.feasible == 'F'
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < chrm.fitness
                chrm.genes = c
                chrm.fitness = f
                chrm.LLnodes = llc
                chrm.Real_LLnodes = rllc
            end
        end
    else
        f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
        if f < chrm.fitness
            chrm.genes = c
            chrm.fitness = f
            chrm.LLnodes = llc
            chrm.Real_LLnodes = rllc
        end
    end
    return chrm
end

function Improve_chromosome(chrm::Chromosome, n_nodes::Int64, TT::Matrix{Float64}, DD::Matrix{Float64}, 
    dEligible::Vector{Int64}, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64,
     sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, problem_type::problem)

    Search_methods = [N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N14p, N15, N16, N17, N18, N19, N20, N21, N22]
    # @show typeof(Search_methods)
    shuffle!(Search_methods)
    @inbounds for search in Search_methods
        chrm = search(chrm, n_nodes, TT, DD, dEligible, ClosenessT, ClosenessD, flying_range, sR, sL, penaltyR, penaltyM, problem_type)
    end
    return chrm
end


#remove one nonfree truck node (randezvous) and make it a drone node between two consecutive truck, change the corresponding drone to turck



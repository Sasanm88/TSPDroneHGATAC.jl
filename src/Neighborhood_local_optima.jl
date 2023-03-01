using StatsBase

function find_two_consecutive_truck_nodes_NF(c::Vector{Int64})  #finds the start points, the ones that can be removed without making infeasible
    # TODO: change the type
    locations = Int[]
    i = 1
    @inbounds while i <= length(c) - 2
        if c[i] < 0
            i += 1
        else
            if c[i+1] < 0
                i += 2
            else
                if c[i+2] < 0
                    i += 3
                else
                    push!(locations, i)
                    push!(locations, i + 1)
                    i += 1
                end
            end
        end
    end
    return unique(locations)
end


function NLO1(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #choose one free truck node and puts it randomly in truck tour
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    if length(free_tnodes) > 0
        remove_from = free_tnodes[rand(1:length(free_tnodes))]
        tnodes_subs = Int[]
        @inbounds for i in tnodes_loc
            if c[i] in ClosenessT[c[remove_from], :]
                push!(tnodes_subs, i)
            end
        end
        # tnodes_subs = findall(x->x in ClosenessT[c[remove_from],:], c)
        if length(tnodes_subs) == 0
            return c
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs))]
        temp = c[remove_from]
        deleteat!(c, remove_from)
        insert!(c, insert_to, temp)
    end
    return c
end

function NLO2(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #choose two consecutive free truck nodes and puts them randomly in truck tour
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    consecutive_tnodes = Int[]

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
            return c
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs)-1)]
        temp1 = c[remove_from]
        temp2 = c[remove_from+1]
        deleteat!(c, [remove_from, remove_from + 1])
        insert!(c, insert_to, temp1)
        insert!(c, insert_to, temp2)
    end
    return c
end

function NLO3(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #choose two consecutive truck nodes and puts them randomly in truck tour
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    consecutive_tnodes = find_two_consecutive_truck_nodes_NF(c)

    if length(consecutive_tnodes) > 0
        remove_from = consecutive_tnodes[rand(1:length(consecutive_tnodes))]
        # TODO: change the type
        tnodes_subs = Int[]
        @inbounds for i in tnodes_loc
            if i != remove_from && i != remove_from + 1
                if c[i] in ClosenessT[c[remove_from], :] || c[i] in ClosenessT[c[remove_from+1], :]
                    push!(tnodes_subs, i)
                end
            end
        end
        if length(tnodes_subs) <= 1
            return c
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs)-1)]
        temp1 = c[remove_from]
        temp2 = c[remove_from+1]

        deleteat!(c, [remove_from, remove_from + 1])
        insert!(c, insert_to, temp2)
        insert!(c, insert_to, temp1)
    end
    return c
end

function NLO4(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})   #swap two truck nodes 
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    idx1 = tnodes_loc[rand(1:length(tnodes_loc))]
    tnodes_subs = Int[]
    @inbounds for i in tnodes_loc
        if i != idx1 && c[i] in ClosenessT[c[idx1], :]
            push!(tnodes_subs, i)
        end
    end
    if length(tnodes_subs) == 0
        return c
    end
    idx2 = tnodes_subs[rand(1:length(tnodes_subs))]
    c[idx1], c[idx2] = c[idx2], c[idx1]
    return c
end

function NLO5(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})   #swap two consecutive truck nodes with another truck node 
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
            return c
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs))]
        c[remove_from], c[insert_to] = c[insert_to], c[remove_from]
        temp = c[remove_from+1]
        deleteat!(c, remove_from + 1)
        insert!(c, insert_to, temp)
    end
    return c
end

function NLO6(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})   #swap two consecutive truck nodes with another truck node 
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    free_tnodes = setdiff(tnodes_loc, chrm.LLnodes)
    indices = Int[]
    @inbounds for node in free_tnodes
        if node < n_nodes
            if c[node+1] > 0
                push!(indices, node)
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
            return c
        end
        insert_to = tnodes_subs[rand(1:length(tnodes_subs))]
        c[remove_from], c[insert_to] = c[insert_to], c[remove_from]
        temp = c[remove_from+1]
        deleteat!(c, remove_from + 1)
        insert!(c, insert_to, temp)
    end
    return c
end

function NLO7(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})        #Truck swap 2-2
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
            return c
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        temp1 = c[i1]
        temp2 = c[i1+1]
        c[i1], c[i1+1] = c[i2], c[i2+1]
        c[i2] = temp1
        c[i2+1] = temp2
    end
    return c
end

function NLO8(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})      #Truck swap 2-2
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
            return c
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        c[i2], c[i2+1] = c[i2+1], c[i2]
        c[i1+1], c[i2] = c[i2], c[i1+1]
    end
    return c
end


function NLO9(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})       #Truck swap 2-2
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
            return c
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        c[i1+1], c[i2+1] = c[i2+1], c[i1+1]
    end
    return c
end

function NLO10(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})        #Truck swap 2-2
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
            return c
        end
        i2 = tnodes_subs[rand(1:length(tnodes_subs))]
        temp1 = c[i1]
        temp2 = c[i1+1]
        c[i1], c[i1+1] = c[i2+1], c[i2]
        c[i2] = temp1
        c[i2+1] = temp2
    end

    return c
end


function NLO11(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #swap a drone node with either one of randezvous nodes or any truck node between
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) < 2
        return c
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
            finish = chrm.LLnodes[2]
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
    return c
end

function NLO12(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #in drone route <i,j,k> swap i and j
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    # tnodes_loc = findall(x->x>0,c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) < 2
        return c
    end
    dnode_idx = rand(1:length(dnodes_loc))
    # tnode_indices = Int[]
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
    return c
end

function NLO13(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #in drone route <i,j,k> swap k and j
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    # tnodes_loc = findall(x->x>0,c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) < 2
        return c
    end
    dnode_idx = rand(1:length(dnodes_loc))
    # tnode_indices = Int[]
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
            finish = chrm.LLnodes[2]
        else
            start = 0
            finish = chrm.LLnodes[1]
        end
    end

    if finish > 0
        c[finish], c[dnodes_loc[dnode_idx]] = -c[dnodes_loc[dnode_idx]], -c[finish]
    end
    return c
end

function NLO14(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #in drone route <i,j,k> swap k and i
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    # tnodes_loc = findall(x->x>0,c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) < 2
        return c
    end
    dnode_idx = rand(1:length(dnodes_loc))
    # tnode_indices = Int[]
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
            finish = chrm.LLnodes[2]
        else
            start = 0
            finish = chrm.LLnodes[1]
        end
    end

    if finish > 0 && start > 0
        c[finish], c[start] = c[start], c[finish]
    end
    return c
end

function NLO15(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #remove one free truck node and make it a drone node between two consecutive truck nodes
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
        end
    end
    return c
end


function NLO16(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #remove one middle truck node and make it a drone node
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
    end
    return c
end


function NLO17(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #remove a drone node and insert it randomly as a truck node between two consequtive truck nodes
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return c
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
    end
    return c
end

function NLO18(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #remove a drone node and insert it randomly as a drone node between two consequtive truck nodes
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return c
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
    end
    return c
end

function NLO19(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #swap two drone nodes
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
        end
    end
    return c
end

function NLO20(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #Swap a truck node with a drone node
    c = copy(chrm.genes)
    n_nodes = length(c)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0
        return c
    end
    idx2 = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    @inbounds for i = 1:n_nodes
        if c[i] > 0 && c[i] in ClosenessT[-c[idx2], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return c
    end
    idx1 = indices[rand(1:length(indices))]
    c[idx1], c[idx2] = -c[idx2], -c[idx1]
    return c
end

function NLO21(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #swap two truck arcs
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    tnodes_loc = findall(x -> x > 0, c)
    idx1 = tnodes_loc[rand(1:length(tnodes_loc))]
    indices = Int[]
    @inbounds for i in tnodes_loc
        if c[i] in ClosenessT[c[idx1], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return c
    end
    idx2 = indices[rand(1:length(indices))]
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    c[idx1:idx2] = c[idx2:-1:idx1]
    return c
end

function NLO22(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #swap two drone nodes and make them truck nodes
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
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
            return c
        end
        dnode2 = indices[rand(1:length(indices))]
        c[dnode1], c[dnode2] = -c[dnode2], -c[dnode1]
    end
    return c
end

function NLO23(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  #
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    # tnodes_loc = findall(x->x>0,c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return c
    end
    d1Loc = dnodes_loc[rand(1:length(dnodes_loc))]
    indices = Int[]
    @inbounds for i in dnodes_loc
        if -c[i] in ClosenessD[-c[d1Loc], :]
            push!(indices, i)
        end
    end
    if length(indices) == 0
        return c
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
    return c
end


function NLO24(chrm::Chromosome, ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64})  ##remove a drone node d, find a random drone pair <i,j,k> change the j to truck and put d in <i,j,d,k> or <i,d,j,k>
    c = copy(chrm.genes)
    # f= chrm.fitness
    # llc = copy(chrm.LLnodes)
    # tnodes_loc = findall(x->x>0,c)
    dnodes_loc = findall(x -> x < 0, c)
    if length(dnodes_loc) == 0 || length(chrm.LLnodes) == 0
        return c
    end
    d1Loc, d2Loc = sample(dnodes_loc, 2)
    d1 = c[d1Loc]
    c[d2Loc] = -c[d2Loc]
    if d1Loc < d2Loc
        deleteat!(c, d1Loc)
        insert!(c, d2Loc - 1, d1)
    elseif d1Loc > d2Loc
        deleteat!(c, d1Loc)
        insert!(c, d2Loc, d1)
    end
    return c
end


function Scape_local_optima(P::Vector{Chromosome}, TT::Matrix{Float64}, DD::Matrix{Float64}, dEligible::Vector{Int64},
    ClosenessT::Matrix{Int64}, ClosenessD::Matrix{Int64}, flying_range::Float64, sR::Float64, sL::Float64, penaltyR::Float64,
    penaltyM::Float64, problem_type::ProblemType, turn::Int64)

    bestF = 0
    bestR = 0
    substituesF = Vector{Chromosome}()
    substituesR = Vector{Chromosome}()
    fsF = Vector{Float64}()
    fsR = Vector{Float64}()
    @inbounds for i = 1:length(P)
        if P[i].feasible == 'F'
            bestF = i
            break
        end
    end
    chrmF = deepcopy(P[bestF])
    push!(fsF, chrmF.fitness)
    push!(substituesF, chrmF)
    best_fF = chrmF.fitness
    best_fR = 0.0

    if flying_range < Inf
        @inbounds for i = 1:length(P)
            if P[i].feasible == 'R'
                bestR = i
                break
            end
        end

        chrmR = deepcopy(P[bestR])


        push!(substituesR, chrmR)
        push!(fsR, chrmR.fitness)

        best_fR = chrmR.fitness
    end

    max_size = 20
    allowed_diff = 0.05

    methods = [NLO1, NLO2, NLO3, NLO4, NLO5, NLO6, NLO7, NLO8, NLO9, NLO10, NLO11, NLO12, NLO13, NLO14, NLO15, NLO16, NLO17, NLO18, NLO19, NLO20, NLO21, NLO22, NLO23, NLO24]




    @inbounds for i = 1:30000*turn
        chrmF = substituesF[rand(1:length(substituesF))]
        r = 0
        r = rand(1:length(methods))
        # c = @code_warntype NLO24(chrmF, ClosenessT, ClosenessD)
        c = methods[r](chrmF, ClosenessT, ClosenessD)
        violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
        if length(violating_drones) == 0
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
            if f < best_fF
                if length(substituesF) < max_size
                    push!(substituesF, Chromosome(c, llc, rllc, f, 'F', 0.0))
                    push!(fsF, f)
                else
                    sort!(substituesF, by=x -> x.fitness)
                    substituesF[max_size] = Chromosome(c, llc, rllc, f, 'F', 0.0)
                    sort!(fsF)
                    fsF[max_size] = f
                end
                best_fF = f
            elseif f > best_fF
                temp = Chromosome(c, llc, rllc, f, 'F', 0.0)
                if !(f in fsF)
                    if (f - best_fF) / best_fF < allowed_diff
                        if length(substituesF) < max_size
                            push!(substituesF, temp)
                            push!(fsF, f)
                        else
                            sort!(substituesF, by=x -> x.fitness)
                            substituesF[max_size] = temp
                            sort!(fsF)
                            fsF[max_size] = f
                        end
                    end
                end
            end
        else
            f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
            if f < best_fR
                if length(substituesR) < max_size
                    push!(substituesR, Chromosome(c, llc, rllc, f, 'R', 0.0))
                    push!(fsR, f)
                else
                    sort!(substituesR, by=x -> x.fitness)
                    substituesR[max_size] = Chromosome(c, llc, rllc, f, 'R', 0.0)
                    sort!(fsR)
                    fsR[max_size] = f
                end
                best_fR = f
            elseif f > best_fR
                temp = Chromosome(c, llc, rllc, f, 'R', 0.0)
                if !(f in fsR)
                    if (f - best_fR) / best_fR < allowed_diff
                        if length(substituesR) < max_size
                            push!(substituesR, temp)
                            push!(fsR, f)
                        else
                            sort!(substituesR, by=x -> x.fitness)
                            substituesR[max_size] = temp
                            sort!(fsR)
                            fsR[max_size] = f
                        end
                    end
                end
            end
        end
    end
    if flying_range < Inf
        @inbounds for i = 1:10000*turn
            chrmR = substituesR[rand(1:length(substituesR))]
            r = 0
            r = rand(1:length(methods))
            c = methods[r](chrmR, ClosenessT, ClosenessD)
            violating_drones = Is_feasibleR(c, DD, TT, dEligible, flying_range, sR, sL, problem_type)
            if length(violating_drones) == 0
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'F', problem_type)
                if f < best_fF
                    if length(substituesF) < max_size
                        push!(substituesF, Chromosome(c, llc, rllc, f, 'F', 0.0))
                        push!(fsF, f)
                    else
                        sort!(substituesF, by=x -> x.fitness)
                        substituesF[max_size] = Chromosome(c, llc, rllc, f, 'F', 0.0)
                        sort!(fsF)
                        fsF[max_size] = f
                    end
                    best_fF = f
                end
            else
                f, llc, rllc = find_fitness(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM, 'R', problem_type)
                if f < best_fR

                    if length(substituesR) < max_size
                        push!(substituesR, Chromosome(c, llc, rllc, f, 'R', 0.0))
                        push!(fsR, f)
                    else
                        sort!(substituesR, by=x -> x.fitness)
                        substituesR[max_size] = Chromosome(c, llc, rllc, f, 'R', 0.0)
                        sort!(fsR)
                        fsR[max_size] = f
                    end
                    best_fR = f
                elseif f > best_fR
                    temp = Chromosome(c, llc, rllc, f, 'R', 0.0)
                    if !(f in fsR)
                        if (f - best_fR) / best_fR < allowed_diff
                            if length(substituesR) < max_size
                                push!(substituesR, temp)
                                push!(fsR, f)
                            else
                                sort!(substituesR, by=x -> x.fitness)
                                substituesR[max_size] = temp
                                sort!(fsR)
                                fsR[max_size] = f
                            end
                        end
                    end
                end
            end
        end
    end

    sort!(substituesF, by=x -> x.fitness)
    @inbounds for i = 1:length(substituesF)
        chrm = substituesF[i]
        if chrm.fitness < P[bestF].fitness
            push!(P, chrm)
        else
            break
        end
    end
    if flying_range < Inf
        sort!(substituesR, by=x -> x.fitness)
        @inbounds for i = 1:length(substituesR)
            chrm = substituesR[i]
            if chrm.fitness < P[bestR].fitness
                push!(P, chrm)
            else
                break
            end
        end
    end
end




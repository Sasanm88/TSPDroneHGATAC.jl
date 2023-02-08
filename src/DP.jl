mutable struct drone_LL
    node::Int64
    luanch::Int64
    Land::Int64
end

mutable struct Tr_node
    tnode::Int64
    dnode::Int64
    before_drone::Vector{Int64}
    after_drone::Vector{Int64}
    f_value::Float64
    prev_action::Int64   #negative means drone flew from there to here visiting dnode, positive means truck moved
end

function find_tnodes(c::Vector{Int64}, n_nodes::Int64)
    tnodes = Vector{Int64}()
    num = 0
    @inbounds for i = 1:length(c)
        if c[i] > 0
            push!(tnodes, c[i])
            num += 1
        end
    end

    pushfirst!(tnodes, 0)
    push!(tnodes, n_nodes + 1)
    return tnodes, num + 2
end


function truck_time(TT::Matrix{Float64}, tnodes::Vector{Int64}, inv_tnodes::Dict{Int64,Int64}, start::Int64, finish::Int64)
    summ = 0.0
    @inbounds for i = inv_tnodes[start]:inv_tnodes[finish]-1
        summ += TT[tnodes[i]+1, tnodes[i+1]+1]
    end
    return summ
end


function find_fitness(c::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, flying_range::Float64,
     sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64, feasibility::Char, problem_type::problem)
    if problem_type == "TSPD"
        if feasibility == 'F'
            return find_fitness_F_TSPD(c, TT, DD, flying_range)
        elseif feasibility == 'R'
            return find_fitness_infR_TSPD(c, TT, DD, flying_range, penaltyR)
        else
            return find_fitness_infM_TSPD(c, TT, DD, flying_range, penaltyR, penaltyM)
        end
    else 
        if feasibility == 'F'
            return find_fitness_F_FSTSP(c, TT, DD, flying_range, sR, sL)
        elseif feasibility == 'R'
            return find_fitness_infR_FSTSP(c, TT, DD, flying_range, sR, sL, penaltyR)
        else
            return find_fitness_infM_FSTSP(c, TT, DD, flying_range, sR, sL, penaltyR, penaltyM)
        end
    end
end



function find_fitness_F_FSTSP(c::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, flying_range::Float64,
     sR::Int64, sL::Int64)
    n_nodes = length(c)
    tnodes, num_t_nodes = find_tnodes(c, n_nodes)
    inv_tnodes = Dict{Int64,Int64}()
    @inbounds for i = 1:num_t_nodes
        inv_tnodes[tnodes[i]] = i
    end

    M = Vector{Tr_node}()
    @inbounds for i = 1:num_t_nodes
        push!(M, Tr_node(tnodes[i], 0, Vector{Int64}(), Vector{Int64}(), Inf, 0))
    end
    M[num_t_nodes].f_value = 0
    tnode_loc = n_nodes + 1
    @inbounds for ii = 1:num_t_nodes-1    #Start from here
        i = num_t_nodes - ii
        if i == 1
            tnode_loc = 0
        else
            while true
                tnode_loc -= 1
                if c[tnode_loc] == tnodes[i]
                    break
                end
            end
        end
        #         tnode_loc = findfirst(x->x==tnodes[i], c)
        if tnode_loc == n_nodes
            M[i].before_drone = [n_nodes + 1]
        else
            if c[tnode_loc+1] < 0
                M[i].after_drone = copy(M[i+1].before_drone)
                M[i].dnode = -c[tnode_loc+1]
                push!(M[i].after_drone, tnodes[i+1])
            else
                M[i].before_drone = copy(M[i+1].before_drone)
                M[i].after_drone = copy(M[i+1].after_drone)
                M[i].dnode = M[i+1].dnode
                push!(M[i].before_drone, tnodes[i+1])
            end
        end


        @inbounds for j in M[i].after_drone
            dtime = DD[j+1, M[i].dnode+1] + DD[M[i].dnode+1, tnodes[i]+1] + sR
            ttime = truck_time(TT, tnodes, inv_tnodes, tnodes[i], j) + sR
            if M[inv_tnodes[j]].prev_action < 0
                ttime += sL
            end
            if dtime < flying_range && ttime < flying_range
                temp = M[inv_tnodes[j]].f_value + max(ttime, dtime)
                if temp < M[i].f_value
                    M[i].f_value = temp
                    M[i].prev_action = -j

                end
            end
        end
        @inbounds for j in M[i].before_drone
            temp = M[inv_tnodes[j]].f_value + truck_time(TT, tnodes, inv_tnodes, tnodes[i], j)
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = j
            end
        end
    end

    # Landnodes = []
    # Launchnodes = []
    LLnodesLoc = Vector{Int64}()
    LLnodes = Vector{Int64}()
    k = 0
    last_landnode = 0
    @inbounds for i = 1:num_t_nodes
        if M[i].tnode == k
            if M[i].prev_action < 0

                k = -M[i].prev_action

                #             push!(Landnodes, k)
                #             push!(Launchnodes,M[i].tnode)

                if M[i].tnode != last_landnode && last_landnode != 0
                    push!(LLnodes, last_landnode)
                end
                if M[i].tnode != 0
                    push!(LLnodes, M[i].tnode)
                end
                x = findfirst(x -> x == M[i].tnode, c)
                if !isnothing(x)
                    push!(LLnodesLoc, x)
                end
                if k == n_nodes + 1
                    if M[i].tnode != last_landnode
                        push!(LLnodes, last_landnode)
                    end
                    break
                end
                last_landnode = k
            else
                k = M[i].prev_action
                if k == n_nodes + 1
                    if M[i].tnode != last_landnode
                        push!(LLnodes, last_landnode)
                    end
                    break
                end
            end
        end
    end

    # LLnodesLoc = Vector{Int64}()
    nLL = length(LLnodes)
    l = 1
    if length(LLnodes) > 0
        @inbounds for i = 1:n_nodes
            if c[i] == LLnodes[l]
                push!(LLnodesLoc, i)
                l += 1
            end
            if l > nLL
                break
            end
        end
    end
    return M[1].f_value, LLnodesLoc
end


function find_fitness_infR_FSTSP(c::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, 
    flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64)
    n_nodes = length(c)
    tnodes, num_t_nodes = find_tnodes(c, n_nodes)
    inv_tnodes = Dict{Int64,Int64}()
    @inbounds for i = 1:num_t_nodes
        inv_tnodes[tnodes[i]] = i
    end

    M = Vector{Tr_node}()
    @inbounds for i = 1:num_t_nodes
        push!(M, Tr_node(tnodes[i], 0, Vector{Int64}(), Vector{Int64}(), Inf, 0))
    end
    M[num_t_nodes].f_value = 0
    tnode_loc = n_nodes + 1
    @inbounds for ii = 1:num_t_nodes-1    #Start from here
        i = num_t_nodes - ii
        if i == 1
            tnode_loc = 0
        else
            while true
                tnode_loc -= 1
                if c[tnode_loc] == tnodes[i]
                    break
                end
            end
        end
        #         tnode_loc = findfirst(x->x==tnodes[i], c)
        if tnode_loc == n_nodes
            M[i].before_drone = [n_nodes + 1]
        else
            if c[tnode_loc+1] < 0
                M[i].after_drone = copy(M[i+1].before_drone)
                M[i].dnode = -c[tnode_loc+1]
                push!(M[i].after_drone, tnodes[i+1])
            else
                M[i].before_drone = copy(M[i+1].before_drone)
                M[i].after_drone = copy(M[i+1].after_drone)
                M[i].dnode = M[i+1].dnode
                push!(M[i].before_drone, tnodes[i+1])
            end
        end

        for j in M[i].after_drone
            dtime = DD[j+1, M[i].dnode+1] + DD[M[i].dnode+1, tnodes[i]+1] + sL
            ttime = truck_time(TT, tnodes, inv_tnodes, tnodes[i], j) + sR
            if M[inv_tnodes[j]].prev_action < 0
                ttime += sL
            end

            temp = M[inv_tnodes[j]].f_value + max(ttime, dtime) + max(0, dtime - flying_range) * penaltyR + max(0, ttime - flying_range) * penaltyR
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = -j

            end
        end
        for j in M[i].before_drone
            temp = M[inv_tnodes[j]].f_value + truck_time(TT, tnodes, inv_tnodes, tnodes[i], j)
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = j
            end
        end
    end

    # Landnodes = []
    # Launchnodes = []
    LLnodesLoc = Vector{Int64}()
    LLnodes = Vector{Int64}()
    k = 0
    last_landnode = 0
    @inbounds for i = 1:num_t_nodes
        if M[i].tnode == k
            if M[i].prev_action < 0

                k = -M[i].prev_action

                #             push!(Landnodes, k)
                #             push!(Launchnodes,M[i].tnode)

                if M[i].tnode != last_landnode && last_landnode != 0
                    push!(LLnodes, last_landnode)
                end
                if M[i].tnode != 0
                    push!(LLnodes, M[i].tnode)
                end
                x = findfirst(x -> x == M[i].tnode, c)
                if !isnothing(x)
                    push!(LLnodesLoc, x)
                end
                if k == n_nodes + 1
                    if M[i].tnode != last_landnode
                        push!(LLnodes, last_landnode)
                    end
                    break
                end
                last_landnode = k
            else
                k = M[i].prev_action
                if k == n_nodes + 1
                    if M[i].tnode != last_landnode
                        push!(LLnodes, last_landnode)
                    end
                    break
                end
            end
        end
    end

    # LLnodesLoc = []
    nLL = length(LLnodes)
    l = 1
    if nLL > 0
        @inbounds for i = 1:n_nodes
            if c[i] == LLnodes[l]
                push!(LLnodesLoc, i)
                l += 1
            end
            if l > nLL
                break
            end
        end
    end
    return M[1].f_value, LLnodesLoc
end

mutable struct Tr_node_inf
    tnode::Int64
    dnode::Vector{Int64}
    before_drone::Vector{Int64}
    after_drone::Vector{Int64}
    f_value::Float64
    prev_action::Int64   #negative means drone flew from there to here visiting dnode, positive means truck moved
end

function find_tdnodes(c::Vector{Int64}, n_nodes::Int64)
    tnodes = Vector{Int64}()
    dnodes = Vector{Vector{Int64}}()
    num = 0
    dnode = Vector{Int64}()
    @inbounds for i = 1:length(c)
        if c[i] > 0
            push!(tnodes, c[i])
            num += 1
            if length(dnode) > 0
                push!(dnodes, dnode)
                dnode = Vector{Int64}()
            end
        else
            push!(dnode, -c[i])
        end
    end
    if length(dnode) > 0
        push!(dnodes, dnode)
    end

    pushfirst!(tnodes, 0)
    push!(tnodes, n_nodes + 1)
    return tnodes, dnodes, num + 2
end

function drone_time(DD::Matrix{Float64}, dnodes::Vector{Int64}, start::Int64, finish::Int64, penalty1::Float64)
    summ = DD[start+1, dnodes[1]+1]
    real_sum = summ
    if length(dnodes) > 1
        @inbounds for i = 2:length(dnodes)
            summ += DD[dnodes[i-1]+1, dnodes[i]+1] * penalty1^(i - 1)
            real_sum += DD[dnodes[i-1]+1, dnodes[i]+1]
        end
    end
    summ += DD[dnodes[length(dnodes)]+1, finish+1] * penalty1^(length(dnodes) - 1)
    real_sum += DD[dnodes[length(dnodes)]+1, finish+1]
    return summ, real_sum
end

function find_fitness_infM_FSTSP(c::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, 
    flying_range::Float64, sR::Int64, sL::Int64, penaltyR::Float64, penaltyM::Float64)
    n_nodes = length(c)
    tnodes, dnodes, num_t_nodes = find_tdnodes(c, n_nodes)
    inv_tnodes = Dict{Int64,Int64}()
    @inbounds for i = 1:num_t_nodes
        inv_tnodes[tnodes[i]] = i
    end

    M = Vector{Tr_node_inf}()
    @inbounds for i = 1:num_t_nodes
        push!(M, Tr_node_inf(tnodes[i], Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), Inf, 0))
    end
    M[num_t_nodes].f_value = 0
    tnode_loc = n_nodes + 1
    d_counter = length(dnodes)

    @inbounds for ii = 1:num_t_nodes-1    #Start from here
        i = num_t_nodes - ii
        if i == 1
            tnode_loc = 0
        else
            while true
                tnode_loc -= 1
                if c[tnode_loc] == tnodes[i]
                    break
                end
            end
        end
        #         tnode_loc = findfirst(x->x==tnodes[i], c)
        if tnode_loc == n_nodes
            M[i].before_drone = [n_nodes + 1]
        else
            if c[tnode_loc+1] < 0
                M[i].after_drone = copy(M[i+1].before_drone)
                M[i].dnode = dnodes[d_counter]
                d_counter -= 1
                push!(M[i].after_drone, tnodes[i+1])
            else
                M[i].before_drone = copy(M[i+1].before_drone)
                M[i].after_drone = copy(M[i+1].after_drone)
                M[i].dnode = M[i+1].dnode
                push!(M[i].before_drone, tnodes[i+1])
            end
        end
        @inbounds for j in M[i].after_drone
            dtime, real_dtime = drone_time(DD, M[i].dnode, j, tnodes[i], penaltyM)
            ttime = truck_time(TT, tnodes, inv_tnodes, tnodes[i], j) + sR
            if M[inv_tnodes[j]].prev_action < 0
                ttime += sL
            end

            temp = M[inv_tnodes[j]].f_value + max(ttime, dtime) + max(0, real_dtime - flying_range) * penaltyR + max(0, ttime - flying_range) * penaltyR
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = -j
            end
        end
        @inbounds for j in M[i].before_drone
            temp = M[inv_tnodes[j]].f_value + truck_time(TT, tnodes, inv_tnodes, tnodes[i], j)
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = j
            end
        end
    end

    # Landnodes = []
    # Launchnodes = []
    LLnodesLoc = Vector{Int64}()
    LLnodes = Vector{Int64}()
    k = 0
    last_landnode = 0
    @inbounds for i = 1:num_t_nodes
        if M[i].tnode == k
            if M[i].prev_action < 0

                k = -M[i].prev_action

                #             push!(Landnodes, k)
                #             push!(Launchnodes,M[i].tnode)

                if M[i].tnode != last_landnode && last_landnode != 0
                    push!(LLnodes, last_landnode)
                end
                if M[i].tnode != 0
                    push!(LLnodes, M[i].tnode)
                end
                x = findfirst(x -> x == M[i].tnode, c)
                if !isnothing(x)
                    push!(LLnodesLoc, x)
                end
                if k == n_nodes + 1
                    if M[i].tnode != last_landnode
                        push!(LLnodes, last_landnode)
                    end
                    break
                end
                last_landnode = k
            else
                k = M[i].prev_action
                if k == n_nodes + 1
                    if M[i].tnode != last_landnode
                        push!(LLnodes, last_landnode)
                    end
                    break
                end
            end
        end
    end

    # LLnodesLoc = []
    nLL = length(LLnodes)
    l = 1
    @inbounds for i = 1:n_nodes
        if c[i] == LLnodes[l]
            push!(LLnodesLoc, i)
            l += 1
        end
        if l > nLL
            break
        end
    end
    return M[1].f_value, LLnodesLoc
end




function find_fitness_F_TSPD(cc::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, flying_range::Float64)
    c = copy(cc)
    n_nodes = length(c)
    tnodes, num_t_nodes = find_tnodes(c, n_nodes)
    inv_tnodes = Dict{Int64, Int64}()
    for i=1:num_t_nodes
        inv_tnodes[tnodes[i]] = i
    end

    M = Vector{Tr_node}()
    for i=1:num_t_nodes
        push!(M,Tr_node(tnodes[i],0,Int[],Int[],Inf,0))
    end
    M[1].f_value = 0

    for i=2:num_t_nodes
        if i==num_t_nodes
            tnode_loc = length(c)+1
        else
            tnode_loc = findfirst(x->x==tnodes[i], c)
        end
        if tnode_loc == 1 
            M[i].after_drone = [0]
        else
            if c[tnode_loc-1]<0
                M[i].before_drone = copy(M[i-1].after_drone)
                M[i].dnode = -c[tnode_loc-1]
                push!(M[i].before_drone, tnodes[i-1])
            else
                M[i].before_drone = copy(M[i-1].before_drone)
                M[i].after_drone = copy(M[i-1].after_drone)
                M[i].dnode = M[i-1].dnode
                push!(M[i].after_drone, tnodes[i-1])
            end
        end
        
        for j in M[i].before_drone
            dtime = DD[j+1,M[i].dnode+1]+DD[M[i].dnode+1,tnodes[i]+1]
            if dtime < flying_range
                temp = M[inv_tnodes[j]].f_value + max(truck_time(TT, tnodes, inv_tnodes, j, tnodes[i]), dtime)
                if temp < M[i].f_value
                    M[i].f_value = temp
                    M[i].prev_action = -(j + (n_nodes+2)*Int(j==0))
                end
            end
        end
        for j in M[i].after_drone
            temp = M[inv_tnodes[j]].f_value + truck_time(TT, tnodes, inv_tnodes, j, tnodes[i])
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = j
            end
        end

    end

    LLnodes = Int[]
    LLnodesLoc = Int[]
    k = n_nodes+1
    for i=1:num_t_nodes
        j = num_t_nodes-i+1
        if M[j].tnode == k
            if M[j].prev_action<0
                k = -M[j].prev_action
                if M[j].tnode != n_nodes+1
                    push!(LLnodes,M[j].tnode)
                end
                x = findfirst(x->x==M[j].tnode,c)
                if !isnothing(x)
                    push!(LLnodesLoc, x)
                end
                if k==n_nodes+2 || k==0
                    break
                end
            else
                k = M[j].prev_action
            end
        end
    end
    reverse!(LLnodes)
    reverse!(LLnodesLoc)
    if n_nodes+1 in LLnodes
        pop!(LLnodes)
        pop!(LLnodesLoc)
    end
    return M[num_t_nodes].f_value, LLnodesLoc
end
            

function find_fitness_infR_TSPD(cc::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, flying_range::Float64, penaltyR::Float64)
    c = copy(cc)
    n_nodes = length(c)
    tnodes, num_t_nodes = find_tnodes(c, n_nodes)
    inv_tnodes = Dict{Int64, Int64}()
    for i=1:num_t_nodes
        inv_tnodes[tnodes[i]] = i
    end

    M = Vector{Tr_node}()
    for i=1:num_t_nodes
        push!(M,Tr_node(tnodes[i],0,[],[],Inf,0))
    end
    M[1].f_value = 0

    for i=2:num_t_nodes
        if i==num_t_nodes
            tnode_loc = length(c)+1
        else
            tnode_loc = findfirst(x->x==tnodes[i], c)
        end
        if tnode_loc == 1 
            M[i].after_drone = [0]
        else
            if c[tnode_loc-1]<0
                M[i].before_drone = copy(M[i-1].after_drone)
                M[i].dnode = -c[tnode_loc-1]
                push!(M[i].before_drone, tnodes[i-1])
            else
                M[i].before_drone = copy(M[i-1].before_drone)
                M[i].after_drone = copy(M[i-1].after_drone)
                M[i].dnode = M[i-1].dnode
                push!(M[i].after_drone, tnodes[i-1])
            end
        end

        for j in M[i].before_drone
            dtime = DD[j+1,M[i].dnode+1]+DD[M[i].dnode+1,tnodes[i]+1]
            temp = M[inv_tnodes[j]].f_value + max(truck_time(TT, tnodes, inv_tnodes, j, tnodes[i]), dtime) + penaltyR*max(0,dtime-flying_range)
            if temp  < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = -(j + (n_nodes+2)*Int(j==0))
            end
        end

        for j in M[i].after_drone
            temp = M[inv_tnodes[j]].f_value + truck_time(TT, tnodes, inv_tnodes, j, tnodes[i])
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = j
            end
        end

    end

    LLnodes = Int[]
    LLnodesLoc = Int[]
    k = n_nodes+1
    for i=1:num_t_nodes
        j = num_t_nodes-i+1
        if M[j].tnode == k
            if M[j].prev_action<0
                k = -M[j].prev_action
                if M[j].tnode != n_nodes+1
                    push!(LLnodes,M[j].tnode)
                end
                x = findfirst(x->x==M[j].tnode,c)
                if !isnothing(x)

                    push!(LLnodesLoc,x)
                end
                if k==n_nodes+2 || k==0
                    break
                end
            else
                k = M[j].prev_action
            end
        end
    end
    reverse!(LLnodes)
    reverse!(LLnodesLoc)
    if n_nodes+1 in LLnodes
        pop!(LLnodes)
        pop!(LLnodesLoc)
    end
    return M[num_t_nodes].f_value, LLnodesLoc
end

function find_fitness_infM_TSPD(c::Vector{Int64}, TT::Matrix{Float64}, DD::Matrix{Float64}, flying_range::Float64,
     penaltyR::Float64, penaltyM::Float64)
    n_nodes = length(c)
    tnodes,dnodes, num_t_nodes = find_tdnodes(c, n_nodes)
    inv_tnodes = Dict{Int64, Int64}()
    for i=1:num_t_nodes
        inv_tnodes[tnodes[i]] = i
    end

    M = Vector{Tr_node_inf}()
    for i=1:num_t_nodes
        push!(M,Tr_node_inf(tnodes[i],[],[],[],Inf,0))
    end
    M[1].f_value = 0
    d_counter = 1
    for i=2:num_t_nodes
        if i==num_t_nodes
            tnode_loc = length(c)+1
        else
            tnode_loc = findfirst(x->x==tnodes[i], c)
        end
        if tnode_loc == 1 
            M[i].after_drone = [0]
        else
            if c[tnode_loc-1]<0
                M[i].before_drone = copy(M[i-1].after_drone)
                M[i].dnode = dnodes[d_counter]
                d_counter += 1
                push!(M[i].before_drone, tnodes[i-1])
            else
                M[i].before_drone = copy(M[i-1].before_drone)
                M[i].after_drone = copy(M[i-1].after_drone)
                M[i].dnode = M[i-1].dnode
                push!(M[i].after_drone, tnodes[i-1])
            end
        end

        for j in M[i].before_drone
            dtime, real_dtime = drone_time(DD,M[i].dnode,j,tnodes[i],penaltyM)
            temp = M[inv_tnodes[j]].f_value + max(truck_time(TT, tnodes, inv_tnodes, j, tnodes[i]),dtime) +  penaltyR*max(0,real_dtime-flying_range)
            if temp  < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = -(j + (n_nodes+2)*Int(j==0))
            end
        end

        for j in M[i].after_drone
            temp = M[inv_tnodes[j]].f_value + truck_time(TT, tnodes, inv_tnodes, j, tnodes[i])
            if temp < M[i].f_value
                M[i].f_value = temp
                M[i].prev_action = j
            end
        end

    end
    LLnodes = Int[]
    LLnodesLoc = Int[]
    k = n_nodes+1
    for i=1:num_t_nodes
        j = num_t_nodes-i+1
        if M[j].tnode == k
            if M[j].prev_action<0
                k = -M[j].prev_action
                if M[j].tnode != n_nodes+1
                    push!(LLnodes,M[j].tnode)
                end
                x = findfirst(x->x==M[j].tnode,c)
                if !isnothing(x)
                    push!(LLnodesLoc, x)
                end
                if k==n_nodes+2 || k==0
                    break
                end
            else
                k = M[j].prev_action
            end
        end
    end
    reverse!(LLnodes)
    reverse!(LLnodesLoc)
    if n_nodes+1 in LLnodes
        pop!(LLnodes)
        pop!(LLnodesLoc)
    end
    return M[num_t_nodes].f_value, LLnodesLoc
end
@enum problem begin
    TSPD
    FSTSP
end

mutable struct Chromosome
    genes::Vector{Int64}      #Represents the sequence and the types of nodes, positive number means the node is visited by truck, while negative is for drone 
    LLnodes::Vector{Int64}    #Represents the location of luanch and land nodes so along with genes it gives us the tour
    fitness::Float64          #This is the Makespan of the given sequence found by JOIN algorithm (smaller, better)
    feasible::Char            #'F': feasible, 'R':infeasible type range, 'M':infeasible type multiple drone nodes   (In TSPD with unlimited flight range, there would be no 'R')
    power::Float64            # This is the fitness after considering the diversity contribution (smaller, better)
end

function Is_feasibleM(c::Vector{Int64})
    prev_negative = false
    @inbounds for i = 1:length(c)
        if c[i] < 0
            if prev_negative
                return false
            else
                prev_negative = true
            end
        else
            prev_negative = false
        end
    end
    return true
end

function fix_negatives(chromosome::Vector{Int64}, start::Int, count::Int)
    if count > 1
        if count % 2 == 1
            @inbounds for i = 1:Int((count - 1) / 2)
                chromosome[start+2*(i-1)+1] *= -1
            end
        else
            if rand() < 0.5
                @inbounds for i = 1:Int(count / 2)
                    chromosome[start+2*(i-1)] *= -1
                end
            else
                @inbounds for i = 1:Int(count / 2)
                    chromosome[start+2*(i-1)+1] *= -1
                end
            end
        end
    end
    return chromosome
end


function make_feasibleM(chromosome::Vector{Int64})
    gene = 1
    start = 0
    count = 0
    @inbounds while gene <= length(chromosome)
        if chromosome[gene] < 0
            if count == 0
                start = gene
            end
            count += 1
            if gene == length(chromosome) && count > 1
                chromosome = fix_negatives(chromosome, start, count)
            end
        else
            if count > 1
                chromosome = fix_negatives(chromosome, start, count)
            end
            count = 0
        end
        gene += 1
    end
    return chromosome
end

mutable struct Dr_node    #Not required for unlimited TSPD
    node::Int64
    before::Vector{Int64}
    after::Vector{Int64}
    Ttimes::Matrix{Float64}
end

function find_beforeANDafter_nodes(c::Vector{Int64}, TT::Matrix{Float64}, problem_type::problem)   #Not required for unlimited TSPD
    n_nodes = length(c)
    DrNodes = Vector{Dr_node}()
    dnodes_loc = findall(x -> x < 0, c)
    @inbounds for i in dnodes_loc
        push!(DrNodes, Dr_node(i, Int[], Int[], zeros(1, 1)))
    end
    if length(dnodes_loc) == 1
        if dnodes_loc[1] == 1
            DrNodes[1].before = [0]
        else
            DrNodes[1].before = [0; c[1:dnodes_loc[1]-1]]
        end
        if dnodes_loc[1] == n_nodes
            DrNodes[1].after = [n_nodes + 1]
        else
            DrNodes[1].after = [c[dnodes_loc[1]+1:n_nodes]; n_nodes + 1]
        end
    else
        @inbounds for i = 1:length(dnodes_loc)
            if i == 1
                if dnodes_loc[i] == 1
                    DrNodes[i].before = [0]
                else
                    DrNodes[i].before = [0; c[1:dnodes_loc[i]-1]]
                end
                DrNodes[i].after = c[dnodes_loc[i]+1:dnodes_loc[i+1]-1]
            elseif i == length(dnodes_loc)
                DrNodes[i].before = copy(DrNodes[i-1].after)
                if dnodes_loc[i] == n_nodes
                    DrNodes[i].after = [n_nodes + 1]
                else
                    DrNodes[i].after = [c[dnodes_loc[i]+1:n_nodes]; n_nodes + 1]
                end
            else
                DrNodes[i].before = copy(DrNodes[i-1].after)
                DrNodes[i].after = c[dnodes_loc[i]+1:dnodes_loc[i+1]-1]
            end
        end
    end
    if problem_type == "TSPD"
        return DrNodes
    end
    @inbounds for drnode in DrNodes
        m = length(drnode.before)
        n = length(drnode.after)
        drnode.Ttimes = zeros(m, n)
        drnode.Ttimes[m, 1] = TT[drnode.before[m]+1, drnode.after[1]+1]
        @inbounds for i = 1:m-1
            drnode.Ttimes[m-i, 1] = drnode.Ttimes[m-i+1, 1] + TT[drnode.before[m-i]+1, drnode.before[m-i+1]+1]
        end
        @inbounds for j = 1:n-1
            @inbounds for i = 1:m
                drnode.Ttimes[i, j+1] = drnode.Ttimes[i, j] + TT[drnode.after[j]+1, drnode.after[j+1]+1]
            end
        end
    end
    return DrNodes
end

function Is_feasibleR(c::Vector{Int64}, DD::Matrix{Float64}, TT::Matrix{Float64}, dEligible::Vector{Int64},
     flying_range::Float64, sR::Int, sL::Int, problem_type::problem)   #Not required for unlimited TSPD
    
    violating_nodes = Vector{Int64}()
    if flying_range == Inf        #For some problems, flying_range is not Inf but technically is. We should manually take care of this
        return violating_nodes
    end

    Dr_nodes = find_beforeANDafter_nodes(c, TT, problem_type)

    if problem_type == "TSPD"
        for drnode in Dr_nodes
            violates = true
            for i in drnode.before
                for j in drnode.after
                    if DD[i+1,-c[drnode.node]+1] + DD[-c[drnode.node]+1, j+1] < flying_range
                        violates = false
                    end
                end
            end
            if violates
                push!(violating_nodes, drnode.node)
            end
        end
        return violating_nodes
    end

    first_feasible = 1
    next_first_feasible = 1
    for drnode in Dr_nodes
        next_first_feasible = length(drnode.after) + 1
        @inbounds for i = first_feasible:length(drnode.before)
            @inbounds for j = 1:length(drnode.after)
                if DD[drnode.before[i]+1, -c[drnode.node]+1] + DD[-c[drnode.node]+1, drnode.after[j]+1] + sR < flying_range
                    if drnode.Ttimes[i, j] + sR < flying_range && !(-c[drnode.node] in dEligible)
                        next_first_feasible = min(j + 1, next_first_feasible)
                    end
                    if drnode.Ttimes[i, j] + sR + sL < flying_range && !(-c[drnode.node] in dEligible)
                        next_first_feasible = min(j, next_first_feasible)
                    end
                end
            end
        end
        if next_first_feasible <= length(drnode.after)
            first_feasible = next_first_feasible
        else
            first_feasible = 1
            push!(violating_nodes, drnode.node)
        end
    end
    return violating_nodes
end


function make_feasibleR(c::Vector{Int64}, violating_nodes::Vector{Int64})   #Not required for unlimited TSPD
    @inbounds for i in violating_nodes
        c[i] = -c[i]
    end
    return c
end




function best_objective(Population::Vector{Chromosome})
    @inbounds for i in 1:length(Population)
        if Population[i].feasible == 'F'
            return Population[i].fitness
        end
    end
    return Inf # Is this correct?
end

function Print_best_route(Population::Vector{Chromosome})
    @inbounds for i in 1:length(Population)
        if Population[i].feasible == 'F'
            @inbounds for j in Population[i].genes
                print(j, " ")
            end
            break
        end
    end
end

function Return_best_route(Population::Vector{Chromosome})
    @inbounds for i in 1:length(Population)
        if Population[i].feasible == 'F'
            return Population[i].genes, Population[i].LLnodes
        end
    end
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
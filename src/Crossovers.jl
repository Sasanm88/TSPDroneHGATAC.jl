
function Crossover_DX1(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)
    child = zeros(Int64, n_nodes)
    idx1 = rand(1:n_nodes)
    idx2 = rand(1:n_nodes)
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    p2 = copy(parent2)
    @inbounds for i = idx1:idx2
        if parent1[i] in p2
            child[i] = parent1[i]
            deleteat!(p2, findfirst(x -> x == parent1[i], p2))
        end
    end
    child[findall(x -> x == 0, child)] = p2
    return child
end

function Crossover_DX2(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #drone crossover
    child = zeros(Int64, n_nodes)
    p1_idx = Vector{Int64}()
    if rand() < 0.5
        p1_idx = findall(x -> x > 0, parent1)
    else
        p1_idx = findall(x -> x < 0, parent1)
    end
    if length(p1_idx) < 2
        return parent1
    end
    p1_idx2 = sample(1:length(p1_idx), Weights(ones(length(p1_idx))), 2, replace=false)
    idx1 = p1_idx2[1]
    idx2 = p1_idx2[2]
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    @inbounds for idx in p1_idx[idx1:idx2]
        child[idx] = parent1[idx]
    end
    i = 1
    j = 1
    @inbounds while i <= n_nodes
        if child[i] != 0
            i += 1
        else
            if abs(parent2[j]) in abs.(child)
                j += 1
            else
                child[i] = parent2[j]
                j += 1
                i += 1
            end
        end
    end
    return child
end

function Crossover_DX3(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)
    signs1 = zeros(Int, n_nodes)
    signs2 = zeros(Int, n_nodes)
    @inbounds for node in parent1
        if node > 0
            signs1[node] = 1
        else
            signs1[-node] = -1
        end
    end
    @inbounds for node in parent2
        if node > 0
            signs2[node] = 1
        else
            signs2[-node] = -1
        end
    end
    idx1, idx2 = sample(1:n_nodes, 2, replace=false)
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end

    child = zeros(Int64, n_nodes)
    p1_copy = abs.(parent1[idx1:idx2])
    p2 = abs.(parent2)

    child[idx1:idx2] = p1_copy
    index = 1
    @inbounds for gene in p2
        if index == idx1
            if idx2 == n_nodes
                break
            else
                index = idx2 + 1
            end
        end
        if !(gene in p1_copy)
            child[index] = gene
            index += 1
        end
    end
    @inbounds for i = 1:n_nodes
        if i < idx1 || i > idx2
            child[i] *= signs1[child[i]]
        else
            child[i] *= signs2[child[i]]
        end
    end
    return child
end


function Crossover_OX1(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #order crossover
    child = zeros(Int64, n_nodes)
    idx1 = rand(2:n_nodes-1)
    idx2 = rand(2:n_nodes-1)
    if idx1 > idx2
        temp = idx1
        idx1 = idx2
        idx2 = temp
    end
    copy_from_P1 = parent1[idx1:idx2]
    child[idx1:idx2] = copy_from_P1
    index = 1
    abs_p1 = abs.(copy_from_P1)
    @inbounds for gene in parent2
        if index == idx1
            if idx2 == n_nodes
                break
            else
                index = idx2 + 1
            end
        end
        if !(abs(gene) in abs_p1)
            child[index] = gene
            index += 1
        end
    end
    return child
end

function Crossover_OX2(parent1::Vector{Int64}, parent2::Vector{Int64}, n_nodes::Int64)   #Order-based crossover
    child = zeros(Int64, n_nodes)
    num_pos = rand(1:n_nodes)
    selected_pos2 = sample(1:n_nodes, Weights(ones(n_nodes)), num_pos, replace=false)
    abs_p2 = abs.(parent2[selected_pos2])

    selected_pos1 = findall(x -> abs(x) in abs_p2, parent1)
    child[setdiff(1:n_nodes, selected_pos1)] = parent1[setdiff(1:n_nodes, selected_pos1)]
    child[selected_pos1] = parent2[sort(selected_pos2)]
    return child
end

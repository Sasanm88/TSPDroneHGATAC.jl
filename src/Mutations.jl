
function sign_mutation(cr::Vector{Int64})
    c = copy(cr)
    @inbounds for i = 1:length(c)
        if rand() < 0.1
            c[i] = -c[i]
        end
    end
    return c
end

function tour_mutation(cr::Vector{Int64})
    c = copy(cr)
    idx = sample(1:length(c), Int(round(0.2 * length(c))), replace=false)
    c[idx] = cr[shuffle(idx)]
    return c
end
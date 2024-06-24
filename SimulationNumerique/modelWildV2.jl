abstract type mathematicalModel end

struct modelWildV2 <: mathematicalModel
    rF::Float64
    f::Float64
    LS::Float64

    rV::Float64
    KV::Float64
    alpha::Float64
        
    function modelWildV2(modelParam::Dict{String, Float64})
        rF, f, LS = modelParam["rF"], modelParam["f"], modelParam["LS"]
        rV, KV, alpha = modelParam["rV"], modelParam["KV"], modelParam["alpha"]
        Veq = KV / (1 + KV * alpha * f / rV)
        # println(f * Veq, " and ", Veq)
        new(rF, f, LS, rV, KV, alpha)
    end
end

function equationModel(model::modelWildV2, variables::Vector{Float64})
    F, V = variables[1], variables[2]
    dF = model.rF * V / (V + model.LS) * (1 - F / (model.f * V)) * F
    dV = model.rV * (1 - V / model.KV) * V - model.alpha * V * F
    return [dF, dV]
end
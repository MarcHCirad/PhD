abstract type mathematicalModel end

struct modelWildV1 <: mathematicalModel
    rF::Float64
    KF::Float64
    LS::Float64

    rV::Float64
    KV::Float64
    alpha::Float64
        
    function modelWildV1(modelParam::Dict{String, Float64})
        rF, KF, LS = modelParam["rF"], modelParam["KF"], modelParam["LS"]
        rV, KV, alpha = modelParam["rV"], modelParam["KV"], modelParam["alpha"]
        Veq = KV / (1 + KV * alpha * KF / rV)
        new(rF, KF, LS, rV, KV, alpha)
    end
end

function equationModel(model::modelWildV1, variables::Vector{Float64})
    F, V = variables[1], variables[2]
    dF = model.rF * V / (V + model.LS) * (1 - F / (model.KF)) * F
    dV = model.rV * (1 - V / model.KV) * V - model.alpha * V * F
    return [dF, dV]
end
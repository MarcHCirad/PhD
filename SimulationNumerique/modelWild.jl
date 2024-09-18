struct modelWild <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    LW::Float64
    f::Float64

    rV::Float64
    KV::Float64
    alpha::Float64
        
    function modelWild(modelParam::Dict{String, Float64})
        rF, f, LW = modelParam["rF"], modelParam["f"], modelParam["LW"]
        rV, KV, alpha = modelParam["rV"], modelParam["KV"], modelParam["alpha"]

        variablesNames = ["F_W", "V_W"]

        new(variablesNames, rF, LW, f, rV, KV, alpha)
    end
end

function equationModel(model::modelWild, variables::Vector{Float64})
    F, V = variables[1], variables[2]
    dF = model.rF * V / (V + model.LW) * (1 - F / (model.f * V)) * F
    dV = model.rV * (1 - V / model.KV) * V - model.alpha * V * F
    return [dF, dV]
end

function equilibriumFV(model::modelWild)
    V = model.KV / (1 + model.KV * model.alpha * model.f / model.rV)
    return [model.f * V, V]
end
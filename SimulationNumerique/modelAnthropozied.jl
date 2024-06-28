abstract type mathematicalModel end

struct modelAlleeEffect <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    KF::Float64
    LF::Float64
    muF::Float64
    lambdaFH::Float64

    rV::Float64
    KV::Float64
    LV::Float64
    muV::Float64
    lambdaVH::Float64

    rH::Float64
    a::Float64
    b::Float64
    c::Float64
    beta::Float64

    thresholdsName::Vector{String}
    
    function modelAlleeEffect(modelParam::Dict{String, Float64})
        rF, KF, LF, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["LF"],
                        modelParam["muF"], modelParam["lambdaFH"]
        rV, KV, LV, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["LV"],
                        modelParam["muV"], modelParam["lambdaVH"]
        rH, a, b, c = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"]
        beta = modelParam["beta"]

        variablesNames = ["F", "V", "H"]
       
        new(variablesNames,
            rF, KF, LF, muF, lambdaFH, rV, KV, LV, muV, 
            lambdaVH, rH, a, b, c, beta,
            thresholdsName)
    end
end

function equationModel(model::modelAlleeEffect, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = (model.rF * H / (H + model.LF) * (1 - F / model.KF) * F - 
                model.muF * F - model.lambdaFH * F * H)

    dV = (model.rV * H / (H + model.LV) * (1 - V / model.KV) * V - 
                model.muV * V - model.lambdaVH * V * H)

    dH = (model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * H * (H/model.beta - 1))

    return [dF, dV, dH]
end

function thresholdH(model::modelAlleeEffect ; F = 0, V =0)
    return (model.a * F + model.b * V + model.c) / model.beta
end

function thresholdF(model::modelAlleeEffect ; H = 0)
    fact1 = H / (model.LF + H)
    fact2 = model.rF / (model.muF + H * model.lambdaFH)
    return fact1 * fact2
end

function thresholdV(model::modelAlleeEffect, H::Float64 ; F = 0)
    fact1 = H / (model.LV + H)
    fact2 = model.rV / (model.lambdaVH * H + model.muV)
    return fact1 * fact2
end


abstract type Model end

struct modelEcoService <: Model
    rF::Float64
    KF::Float64
    omega::Float64
    f::Float64
    muF::Float64
    lambdaFH::Float64

    rV::Float64
    KV::Float64
    alpha::Float64
    muV::Float64
    lambdaVH::Float64
    gamma::Float64
    betaF::Float64
    betaH::Float64

    rH::Float64
    a::Float64
    b::Float64
    c::Float64
    g::Float64

    
    function modelEcoService(modelParam::Dict{String, Float64})
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"], modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        rV, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["alpha"], modelParam["muV"], modelParam["lambdaVH"]
        gamma, betaF, betaH = modelParam["gamma"], modelParam["betaF"], modelParam["betaH"]
        rH, a, b, c, g = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"], modelParam["g"]
       
        new(rF, KF, omega, f, muF, lambdaFH, rV, KV, alpha, muV, lambdaVH, gamma, betaF, betaH, rH, a, b, c, g)
    end
end

function equationModel(model::modelEcoService, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = model.rF * (1 - F / model.KF) * F - model.omega * model.f * H * F 
                        - model.muF * F - model.lambdaFH * model.g * F * H
    dV = model.rV * (1 - model.gamma*exp(-model.betaF * F - model.betaH * H)) * (1 - V / model.KV) * V
             - model.alpha * V * F - model.muV * V - model.lambdaVH * g * V * H
    dH = model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * H

    return [dF, dV, dH]
end   
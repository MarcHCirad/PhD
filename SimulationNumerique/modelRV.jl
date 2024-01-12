
struct modelRV <: Model
    rV::Float64
    H0::Float64
    KV::Float64
    alpha::Float64
    muV::Float64
    lambdaVH::Float64
    rF::Float64
    KF::Float64
    omega::Float64
    f::Float64
    muF::Float64
    lambdaFH::Float64
    e::Float64
    muH::Float64

    
    function modelRV(modelParam::Dict{String, Float64})
        rV, H0, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["H0"], modelParam["KV"], modelParam["alpha"], modelParam["muV"], modelParam["lambdaVH"]
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"], modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        e, muH = modelParam["e"], modelParam["muH"]
       
        new(rV, H0, KV, alpha, muV, lambdaVH, rF, KF, omega, f, muF, lambdaFH, e, muH)
    end
end

function equationModel(model::modelRV, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    
    dF = model.rF * (1 - F / model.KF) * F - model.omega * model.f * F - model.muF * F - model.lambdaFH * F * H
    dV = model.rV * H / (H + model.H0) * (1 - V / model.KV) * V - model.alpha * V * F - model.muV * V - model.lambdaVH * V * H
    dH = model.e * (model.lambdaFH * F + model.lambdaVH * V) * H - model.muH * H^2
    
    return [dF, dV, dH]
end

function bio(model::modelRV)
    println("Test")
    println(model.rV)
end

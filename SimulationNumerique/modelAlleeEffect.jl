abstract type mathematicalModel end

struct modelAlleeEffect <: mathematicalModel
    rF::Float64
    KF::Float64
    omega::Float64
    f::Float64
    muF::Float64
    lambdaFH::Float64

    rV::Float64
    KV::Float64
    LV::Float64
    alpha::Float64
    muV::Float64
    lambdaVH::Float64

    rH::Float64
    a::Float64
    b::Float64
    c::Float64
    beta::Float64
    
    function modelAlleeEffect(modelParam::Dict{String, Float64})
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"],
                        modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        rV, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["alpha"], 
                        modelParam["muV"], modelParam["lambdaVH"]
        LV = modelParam["LV"]
        rH, a, b, c = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"]
        beta = modelParam["beta"]
       
        new(rF, KF, omega, f, muF, lambdaFH, rV, KV, LV, alpha, muV, 
                    lambdaVH, rH, a, b, c, beta)
    end
end

function equationModel(model::modelAlleeEffect, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = (model.rF * (1 - F / model.KF) * F - model.omega * model.f * H * F - model.muF * F - model.lambdaFH * F * H)

    dV = (model.rV * H / (H + model.LV) * 
            (1 - V / model.KV) * V - model.alpha * V * F - model.muV * V -
            model.lambdaVH * V * H)

    dH = (model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * H * (H/model.beta - 1))

    return [dF, dV, dH]
end

function NF(model::modelAlleeEffect ; H = 0)
    return model.rF / (model.muF + H * (model.omega * model.f + model.lambdaFH))
end

function thresholdV(model::modelAlleeEffect, H::Float64 ; F = 0)
    fact1 = H / (model.LV + H)
    fact2 = model.rV / (model.lambdaVH * H + model.alpha * F + model.muV)
    return fact1 * fact2
end

function thresholdVH1(model::modelAlleeEffect)
    fact1 = (model.b * model.KV + model.c) / (model.b * model.KV)
    fact2 = model.rV / (model.lambdaVH * model.LV + model.muV)
    return fact1 * fact2
end

function thresholdVH2(model::modelAlleeEffect)
    A = model.lambdaVH + model.rV / model.KV * model.b
    B = model.rV * (1 + model.c / (model.KV * model.b)) - model.lambdaVH * model.LV - model.muV
    C = model.muV * model.LV

    return B^2 / (4 * C * A)
end

function eqTrivial(modelAlleeEffect)
    return (0,0,0)
end

function eqHBeta(model::modelAlleeEffect)
    return (0,0, model.beta)
end

function eqF(model::modelAlleeEffect)
    Feq = model.KF * ( 1 - model.muF / model.rF)
    return (Feq, 0, 0)
end

function eqFHBeta(model::modelAlleeEffect)
    Feq = model.KF * (1 - 
        (model.muF + model.beta*(model.omega*model.f + model.lambdaFH) / model.rF))

    return (Feq, 0, model.beta)
end

function eqVHBeta(model::modelAlleeEffect)
    Veq = model.KV * (1 - 
        (model.beta + model.LV)/model.beta * (model.muV + model.beta*model.LV)/model.rV)

    return (0, Veq, model.beta)
end

function eqFVHBeta(model::modelAlleeEffect)
    Feq = model.KF * (1 - 
    (model.muF + model.beta*(model.omega*model.f + model.lambdaFH) / model.rF))
    Veq = model.KV * (1 - (model.beta + model.LV)/model.beta *
     (model.muV + model.beta*model.LV + model.alpha * Feq)/model.rV)

    return (Feq, Veq, model.beta)
end

function eqH(model::modelAlleeEffect)
    return (0,0,model.c)
end

function eqFH(model::modelAlleeEffect)
    numFeq = 1 - model.muF/model.rF - model.c * (model.omega*model.f + model.lambdaFH)/model.rF
    denFeq = 1 + model.a * model.KF * (model.omega*model.f + model.lambdaFH)/model.rF
    Feq = model.KF * numFeq / denFeq
    return (Feq, 0, model.a * Feq + model.c)
end

function eqVH1(model::modelAlleeEffect)
    A = model.lambdaVH + model.rV / (model.KV * model.b)
    B = model.rV * (1 + model.c / (model.KV * model.b)) - model.lambdaVH * model.LV - model.muV
    C = model.muV * model.LV

    B = B/A
    C = C/A

    Heq = 1/2 * B - 1/2 * sqrt(B^2 - 4*C)
    Veq = (Heq - model.c) / model.b

    return (0, Veq, Heq)
end

function eqVH2(model::modelAlleeEffect)
    A = model.lambdaVH + model.rV / (model.KV * model.b)
    B = model.rV * (1 + model.c / (model.KV * model.b)) - model.lambdaVH * model.LV - model.muV
    C = model.muV * model.LV

    B = B/A
    C = C/A

    Heq = 1/2 * B + 1/2 * sqrt(B^2 - 4*C)
    Veq = (Heq - model.c) / model.b

    return (0, Veq, Heq)
end
struct modelAnthropized <: mathematicalModel
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
    
    function modelAnthropized(modelParam::Dict{String, Float64})
        rF, KF, LF, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["LF"],
                        modelParam["muF"], modelParam["lambdaFH"]
        rV, KV, LV, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["LV"],
                        modelParam["muV"], modelParam["lambdaVH"]
        rH, a, b, c = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"]
        beta = modelParam["beta"]

        variablesNames = ["F", "V", "H"]
       
        new(variablesNames,
            rF, KF, LF, muF, lambdaFH, rV, KV, LV, muV, 
            lambdaVH, rH, a, b, c, beta)
    end
end

function equationModel(model::modelAnthropized, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = (model.rF * H / (H + model.LF) * (1 - F / model.KF) * F - 
                model.muF * F - model.lambdaFH * F * H)

    dV = (model.rV * H / (H + model.LV) * (1 - V / model.KV) * V - 
                model.muV * V - model.lambdaVH * V * H)

    dH = (model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * H * (H/model.beta - 1))

    return [dF, dV, dH]
end

function thresholdH(model::modelAnthropized ; F = 0, V =0)
    return (model.a * F + model.b * V + model.c) / model.beta
end

function thresholdF(model::modelAnthropized, H::Float64)
    fact1 = H / (model.LF + H)
    fact2 = model.rF / (model.muF + H * model.lambdaFH)
    return fact1 * fact2
end

function thresholdV(model::modelAnthropized, H::Float64)
    fact1 = H / (model.LV + H)
    fact2 = model.rV / (model.lambdaVH * H + model.muV)
    return fact1 * fact2
end

function equilibriumTE(model::modelAnthropized)
    return [0,0,0]
end

function equilibriumF(model::modelAnthropized)
    if model.LF != 0 ##error
        return [0,0,0]
    else
        return [model.KF, 0, 0]
    end
end

function equilibriumHBeta(model::modelAnthropized)
    return [0, 0, model.beta]
end

function equilibriumFHBeta(model::modelAnthropized)
    TFBeta = thresholdF(model, model.beta)
    if TFBeta > 1
        return [model.KF * (1 - 1/TFBeta), 0, model.beta]
    else
        return [0,0,0] ##error
    end
end

function equilibriumVHBeta(model::modelAnthropized)
    TVBeta = thresholdV(model, model.beta)
    if TVBeta > 1
        return [0, model.KV * (1 - 1/TVBeta), model.beta]
    else
        return [0,0,0] ##error
    end
end

function equilibriumFVHBeta(model::modelAnthropized)
    TFBeta = thresholdF(model, model.beta)
    TVBeta = thresholdV(model, model.beta)
    if (TVBeta > 1) && (TFBeta > 1)
        return [model.KF * (1 - 1/TFBeta), model.KV * (1 - 1/TVBeta), model.beta]
    else
        return [0,0,0] ##error
    end
end

function equilibriumH(mode::modelAnthropized)
    return [0,0, model.c]
end

function equilibriumFH(model::modelAnthropized)
    A = 1 + model.a * model.KF * model.lambdaFH / model.rF
    B = model.a * model.KF * (1 - (model.muF + model.LF * model.lambdaFH) / model.rF) +
        model.c
    C = model.a * model.KF * model.LF * model.muF / model.rF

    Delta = B^2 - 4 * A * C
    if (Delta >= 0) && (B > 0)
        H = B / (2 * A) + sqrt(Delta) / (2 * A)
        TF = thresholdF(model, H)
        if TF > 1
            F = model.KF * (1 - TF)
            return (F, 0, H)
        end
    end
    ##error
    return (0,0,-5)
end

function equilibriumFVH(model::modelAnthropized)
    A = 1 + model.a * model.KF * model.lambdaFH / model.rF +
             model.b * model.KV * model.lambdaVH / model.rV
    B = model.a * model.KF * (1 - (model.muF + model.LF * model.lambdaFH) / model.rF) +
        model.b * model.KV * (1 - (model.muV + model.LV * model.lambdaVH) / model.rV) +
        model.c
    C = model.a * model.KF * model.LF * model.muF / model.rF +
        model.b * model.KV * model.LV * model.muV / model.rV

    Delta = B^2 - 4 * A * C
    if (Delta >= 0) && (B > 0)
        H = B / (2 * A) + sqrt(Delta) / (2 * A)
        TF = thresholdF(model, H)
        TV = thresholdV(model, H)
        if (TF > 1) && (TV > 1)
            F = model.KF * (1 - TF)
            V = model.KV * (1 - TV)
            return (F, V, H)
        end
    end
    ##error
    return (0,0,-5)       
end


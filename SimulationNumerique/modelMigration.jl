struct modelMigration <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    LW::Float64
    f::Float64
    lambdaFHW::Float64

    rV::Float64
    KV::Float64
    lambdaW::Float64
    lambdaVHW::Float64

    rH::Float64
    aW::Float64
    bW::Float64
    c::Float64
    beta::Float64
        
    function modelMigration(modelParam::Dict{String, Float64})
        rF, f, LW, lambdaFHW = modelParam["rF"], modelParam["f"], modelParam["LW"], modelParam["lambdaFHW"]
        rV, KV, lambdaW, lambdaVHW = modelParam["rV"], modelParam["KV"], modelParam["lambdaW"], modelParam["lambdaVHW"]
        rH, aW, bW, c = modelParam["rH"], modelParam["aW"], modelParam["bW"], modelParam["c"]
        beta = modelParam["beta"]

        variablesNames = ["HA", "FW", "VW"]

        new(variablesNames, rF, LW, f, lambdaFHW, 
                            rV, KV, lambdaW, lambdaVHW,
                            rH, aW, bW, c, beta)
    end
end

function equationModel(model::modelMigration, variables::Vector{Float64})
    HA, FW, VW = variables[1], variables[2], variables[3]
    dHA = model.rH * (1 - HA / (model.aW * FW + model.bW * VW + model.c)) * (HA/model.beta - 1) * HA
    dFW = model.rF * VW / (VW + model.LW) * (1 - FW / (model.f * VW)) * FW - model.lambdaFHW * FW * HA 
    dVW = model.rV * (1 - VW / model.KV) * VW - model.lambdaW * VW * FW - model.lambdaVHW * HA * VW
    return [dHA, dFW, dVW]
end

function thresholdF(model::modelMigration, H::Float64, V::Float64)
    return model.rF / (model.lambdaFHW * H) * V / (model.LW + V)
end

function equilibriumFV(model::modelMigration)
    V = model.KV / (1 + model.KV * model.lambdaW * model.f / model.rV)
    return [0, model.f * V, V]
end

function equilibriumHV(model::modelMigration)
    if (model.b == 0) && (model.lambdaVHW == 0)
        return (model.c, 0, model.KV)
    else
        error("try to compute equilibriumHFV while not implemented")
    end
end


function equilibriumHFV(model::modelMigration)
    if (model.bW == 0) && (model.lambdaVHW == 0)
        A = model.aW * model.lambdaFHW * model.lambdaW / (model.rF * model.rV)
        B = -(1/(model.f * model.KV) + model.lambdaW / model.rV 
                + model.aW * model.lambdaFHW / model.rF * (model.KV + model.LW) / model.KV
                - model.c * model.lambdaFHW * model.lambdaW / (model.rF * model.rV))
        C = 1 - model.c * model.lambdaFHW / model.rF * (model.KV + model.LW) / model.KV
        Delta = B^2 - 4*A*C
        F1 = -B / (2*A) - sqrt(Delta) / (2 * A)
        return (model.aW * F1 + model.c, F1, model.KV * (1 - model.lambdaW * F1 / model.rV))
    else
        error("try to compute equilibriumHFV while not implemented")
    end
end
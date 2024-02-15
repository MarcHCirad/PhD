abstract type MathematicalModel end

struct modelEcoService <: MathematicalModel
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
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"],
                        modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        rV, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["alpha"], 
                        modelParam["muV"], modelParam["lambdaVH"]
        gamma, betaF, betaH = modelParam["gamma"], modelParam["betaF"], modelParam["betaH"]
        rH, a, b, c, g = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"], 
                        modelParam["g"]
       
        new(rF, KF, omega, f, muF, lambdaFH, rV, KV, alpha, muV, 
                    lambdaVH, gamma, betaF, betaH, rH, a, b, c, g)
    end
end

function equationModel(model::modelEcoService, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = (model.rF * (1 - F / model.KF) * F - model.omega * model.f * H * F )
                - model.muF * F - model.lambdaFH * model.g * F * H

    dV = (model.rV * (1 - model.gamma*exp(-model.betaF * F - model.betaH * H)) * 
            (1 - V / model.KV) * V - model.alpha * V * F - model.muV * V 
            - model.lambdaVH * model.g * V * H)

    dH = model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * H

    return [dF, dV, dH]
end

function existenceTresholds(model::modelEcoService)
    println("WARNING : the thresolds for FVH existence are not implemented")
    N0F = model.rF / model.muF
    N0V = model.rV * (1-model.gamma) / model.muV

    F_FV = model.KF * (1-1/N0F)
    N1V = model.rV*(1-model.gamma*exp(-model.betaF*F_FV)) / (model.alpha*F_FV + model.muV)

    TFH = N0F/(1 + model.c * (model.omega*model.f + model.g*model.lambdaFH)/model.muF)

    T2VH = model.rV / (model.muV + model.lambdaVH*model.g)
    T1VH = (1-model.gamma*exp(-model.betaH*model.c))*T2VH

    return Dict([("N0F", N0F), ("N0V", N0V), ("N1V", N1V), ("TFH", TFH), 
                ("T1VH", T1VH), ("T2VH", T2VH)])

end

function interpretExistenceTresholds(model::modelEcoService)
    thresolds = existenceTresholds(model)
    N0F, N0V, N1V  = thresolds["N0F"], thresolds["N0V"], thresolds["N1V"]
    TFH, T1VH, T2VH = thresolds["TFH"], thresolds["T1VH"], thresolds["T2VH"]
    existence = Dict([("TE", 1), ("F", undef), ("V", undef), 
            ("FV", undef), ("FH", undef), ("VH", undef), ("FVH", undef) ])

    if N0F < 1
        existence["F"] = 0
        existence["FV"] = 0
        existence["FH"] = 0
        existence["FVH"] = 0
    else
        existence["F"] = 1
        if N1V > 1
            existence["FV"] = 1
        else
            existence["FV"] = 0
        end
        if TFH > 1
            existence["FH"] = 1
        else
            existence["FH"] = 0
        end
    end

    if N0V < 1
        existence["V"] = 0
    else
        existence["V"] = 1
    end

    if T1VH < 1
        if T2VH < 1
            existence["VH"] = 0
        else
            existence["VH"] = 1
        end
    else
        existence["VH"] = 1
    end

    return existence
end

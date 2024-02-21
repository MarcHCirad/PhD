abstract type mathematicalModel end

struct modelEcoService <: mathematicalModel
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
    dF = (model.rF * (1 - F / model.KF) * F - model.omega * model.f * H * F
                - model.muF * F - model.lambdaFH * model.g * F * H)

    dV = (model.rV * (1 - model.gamma*exp(-model.betaF * F - model.betaH * H)) * 
            (1 - V / model.KV) * V - model.alpha * V * F - model.muV * V 
            - model.lambdaVH * model.g * V * H)

    dH = model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * H

    return [dF, dV, dH]
end

function tresholdTF(model::modelEcoService, V::Float64)
    den = model.muF + (model.b * V + model.c) * (model.omega*model.f + model.g*model.lambdaFH)
    return model.rF/den
end

function tresholdT1V(model::modelEcoService, F::Float64)
    rslt1 = model.rV / (model.muV + model.alpha*F + model.lambdaVH*model.g * (model.a*F + model.c))
    return (1 - model.gamma * exp(-model.betaH * (model.c + model.a * F) - model.betaF * F)) * rslt1
end

function tresholdT2V(model::modelEcoService, F::Float64)
    rslt1 = model.rV / (model.muV + model.alpha*F + model.lambdaVH*model.g * (model.a*F + model.c))
    return rslt1
end

function tresholdSVH(model::modelEcoService, V::Float64)
    p = model.rV*model.gamma*exp(-model.betaH*model.c)
    q = model.betaH * model.b
    fV = p * exp(-q * V) * (1 - V / model.KV)
    return fV - model.lambdaVH * model.g / model.betaH - (model.rV - p*exp(-q * V)) / (q * model.KV)
end

function zerosEqTranscendental(p::Float64, q::Float64, w::Float64, z::Float64, KV::Float64)
    eqTranscendental = V -> p * exp(-q * V) * (1 - V / KV) + w * V - z
    noSolution = Float64[]
    if p <= 0
        println("Value of p must be positive; stop solving")
        return noSOlution
    else
        if (0 < q) && (0 < w)
            if (z < 0) || (max(p, KV * w) < z)
                println("Case1")
                return noSolution
            elseif (w * KV < z < p)
                println("Case2")
                return [Roots.find_zero(eqTranscendental, (0,KV))]
            elseif (p < z < w*KV)
                println("Case3")
                # println(eqTranscendental(0), eqTranscendental(z/w))
                return [Roots.find_zero(eqTranscendental, (0, z/w+0.4))]
            elseif (0 < z < min(p, w*KV))
                println("Case4")
                return Roots.find_zeros(eqTranscendental, 0, KV)
            end
        elseif (0 < q) && (w < 0)
            if (p < z) || (z < w * KV)
                println("Case5")
                return noSolution
            else
                println("Case6")
                return [Roots.find_zero(eqTranscendental, (0, KV))]
            end
        elseif (q < 0) && (0 < w)
            if (z < 0) || (0 < z < min(p, KV*w))
                println("Case7")
                return noSolution
            elseif (w*KV < z < p)
                Vm = KV + 1/q
                if eqTranscendental(Vm) < 0
                    Vm = 0
                end
                println("Case8")
                println(eqTranscendental(Vm), eqTranscendental(KV))
                return [Roots.find_zero(eqTranscendental, (Vm, KV))]
            elseif (p < z < w*KV)
                println("Case9")
                return [Roots.find_zero(eqTranscendental, (0, w*KV))]
            else
                Vm = KV + 1/q
                if eqTranscendental(Vm) > 0
                    println("Case10")
                    firstZero = Roots.find_zero(eqTranscendental, (0, Vm))
                    secondZero = Roots.find_zero(eqTranscendental, (Vm, Kv))
                    return append!(firstZero, secondZero)
                else
                    println("Case11")
                    return Roots.find_zeros(eqTranscendental, 0, KV)
                end
            end
        else
            println("Not yet done")
            return noSolution
        end
    end
end

function zerosEqTranscendentalVH(model::modelEcoService)
    p, q = model.rV * model.gamma * exp(-model.betaH * model.c), model.b * model.betaH
    w = model.rV/model.KV + model.lambdaVH * model.b * model.g
    z = model.rV - model.muV - model.lambdaVH * model.c * model.g
    return zerosEqTranscendental(p, q, w, z, model.KV)
end

function existenceTresholds(model::modelEcoService)
    println("WARNING : the thresolds for FVH existence are not implemented")
    N0F = model.rF / model.muF
    N0V = model.rV * (1-model.gamma) / model.muV

    TF0 = tresholdTF(model, 0.)

    FFV = model.KF * (1-1/N0F)
    TFV = model.rV*(1-model.gamma*exp(-model.betaF*FFV))/(model.alpha*FFV + model.muV)

    T1V0 = tresholdT1V(model, 0.)
    T2V0 = tresholdT2V(model, 0.)

    return Dict([("N0F", N0F), ("N0V", N0V), ("T1V0", T1V0), ("T2V0", T2V0), 
                ("TF0", TF0), ("TFV", TFV)])
end

function interpretExistenceTresholds(model::modelEcoService)
    thresholds = existenceTresholds(model)
    N0F, N0V, TFV = thresholds["N0F"], thresholds["N0V"], thresholds["TFV"]
    TF0, T1V0, T2V0 = thresholds["TF0"], thresholds["T1V0"], thresholds["T2V0"]
    existence::Dict{String, Vector{Vector{Float64}}} = Dict([("TE", [zero(Vector{Float64}(undef,3))]), ("F", []), ("V", []), 
            ("FV", []), ("FH", []), ("VH", []), ("FVH", []) ])

    if N0F < 1
        existence["F"] = []
        existence["FV"] = []
        existence["FH"] = []
        existence["FVH"] = []
    else
        F = model.KF * (1 - 1/N0F)
        existence["F"] = [[F, 0. ,0.]]
        if TFV > 1
            V = model.KV * (1 - (model.alpha * F + model.muV) / 
                        (model.rV * (1 - model.gamma*exp(-model.betaF*F))))
            existence["FV"] = [[F, V, 0.]]
        else
            existence["FV"] = []
        end
        if TF0 > 1
            num = (N0F - model.c * (model.omega*model.f + model.lambdaFH*model.g) / model.muF - 1.0)
            den = (N0F + model.KF * model.a * (model.omega*model.f + model.lambdaFH*model.g) / model.muF)
            F = (model.KF * num) / den
            existence["FH"] = [[F, 0., model.a * F + model.c]]
        else
            existence["FH"] = []
        end
    end

    if N0V < 1
        existence["V"] = []
    else
        existence["V"] = [[0., model.KV * (1 - 1 / N0V), 0.]]
    end

    if T1V0 < 1
        if T2V0 < 1
            existence["VH"] = []
        else
            ### There may be 0,1 or 2 equilibrium VH
            valuesV = zerosEqTranscendentalVH(model)
            equilibrium = []
            for V in valuesV
                push!(equilibrium, [0., V, model.b * V + model.c])
            end
            existence["VH"] = equilibrium
        end
    else
        ### There is exactly one equilibria VH
        valuesV = zerosEqTranscendentalVH(model)
        equilibrium = []
        for V in valuesV
            push!(equilibrium, [0., V, model.b * V + model.c])
        end
        existence["VH"] = equilibrium
    end

    return existence
end

function interpretStability(model::modelEcoService)   
    equilibrium = interpretExistenceTresholds(model)
    stableEquilibrium = Dict{String, Vector{Vector{Float64}}}()

    if !(isempty(equilibrium["FH"]))
        localEq = equilibrium["FH"][1]
        TV1F = tresholdT1V(model, localEq[1])
        if TV1F < 1
            stableEquilibrium["FH"] = [localEq]
        end
    end

    if !(isempty(equilibrium["VH"]))
        stableVH = Vector{Vector{Float64}}()
        for localEq in equilibrium["VH"]
            localV = localEq[2]
            sV = tresholdSVH(model, localV)
            TF = tresholdTF(model, localV)
            if (sV < 0) && (TF < 1)
                push!(stableVH, localEq)
            end
        end
        stableEquilibrium["VH"] = stableVH
    end

    return stableEquilibrium
end

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

    thresholdsName::Vector{String}
    
    function modelAlleeEffect(modelParam::Dict{String, Float64})
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"],
                        modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        rV, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["alpha"], 
                        modelParam["muV"], modelParam["lambdaVH"]
        LV = modelParam["LV"]
        rH, a, b, c = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"]
        beta = modelParam["beta"]

        thresholdsName = ["TH0", #
            "THFH", #
            "THVVH1",
            "THVVH2",
            "TF0", #
            "TFBeta", #
            "TFHVH1",
            "TFHVH2",
            "TVBeta0", #
            "TVBetaFBeta", #
            "TFc", #
            "TVc0", #
            "TVcFFH", #
            "DeltaVH", #
            "TVi",
            "THVH1"
            ] #
       
        new(rF, KF, omega, f, muF, lambdaFH, rV, KV, LV, alpha, muV, 
                    lambdaVH, rH, a, b, c, beta, thresholdsName)
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

function thresholdH(model::modelAlleeEffect ; F = 0, V =0)
    return (model.a * F + model.b * V + model.c) / model.beta
end

function computeHmin(model::modelAlleeEffect)
    A = model.rV / model.rH * model.beta / (model.KV * model.b) + 1
    B = model.LV - model.beta - model.c * model.rV / model.rH * model.beta / (model.KV * model.b)
    C = - model.beta * model.LV
    Hmin = (-B + sqrt(B^2 - 4*A*C)) / (2 * A)
    return Hmin
end

function computeBetaMax(model::modelAlleeEffect)
    return (model.rF - model.muF) / (model.lambdaFH + model.omega * model.f)
end

function thresholdHMin(model::modelAlleeEffect, H::Float64)
    A = model.rV / model.rH * model.beta / (model.KV * model.b) + 1
    B = model.LV - model.beta - model.c * model.rV / model.rH * model.beta / (model.KV * model.b)
    C = - model.beta * model.LV
    Hmin = (-B + sqrt(B^2 - 4*A*C)) / (2 * A)
    return H / Hmin
end

function thresholdF(model::modelAlleeEffect ; H = 0)
    return model.rF / (model.muF + H * (model.omega * model.f + model.lambdaFH))
end

function thresholdV(model::modelAlleeEffect, H::Float64 ; F = 0)
    fact1 = H / (model.LV + H)
    fact2 = model.rV / (model.lambdaVH * H + model.alpha * F + model.muV)
    return fact1 * fact2
end

function thresholdVi(model::modelAlleeEffect)
    num = 1 - (model.muV + model.lambdaVH * model.LV) / model.rV
    den = 1 + 2 * model.lambdaVH * model.b * model.KV / model.rV
    fraction1 = num / den
    fraction2 = model.b * model.KV / model.c
    return fraction1 * fraction2
end

function thresholdDeltaVH(model::modelAlleeEffect)
    frac1 = model.b / (4 * model.KV )
    frac2 = (1 + model.b * model.KV * model.lambdaVH / model.rV)
    frac3 = (model.c + model.LV) * (model.muV + model.lambdaVH * model.c) / model.rV - model.c
    frac4 = model.KV * ((model.lambdaVH * (2 * model.c + model.LV) + model.muV) / model.rV -1) + model.c / model.b
    return frac1 * frac2^(-1) * frac3^(-1) * frac4^2
end

function thresholdFVH1(model::modelAlleeEffect)
    FHterm = model.KF * (model.omega * model.f + model.lambdaFH) / model.rF
    frac1 = (1 - model.alpha / model.lambdaFH * FHterm) / (1 + model.a * FHterm)
    frac2 = model.b * model.KV * model.lambdaVH / model.rV

    return -frac2 * frac1
end

function thresholdDeltaFVH(model::modelAlleeEffect)
    FFH, _, HFH = eqFH(model)
    F = model.KF * (1 - model.muF / model.rF)
    H = model.a * F + model.c
    FHterm = model.KF * (model.omega * model.f + model.lambdaFH) / model.rF

    term1 = 4 * model.KV * H / model.b
    term2 = 1 + model.b * model.KV * model.lambdaVH / model.rV * 
            (1 - model.alpha / model.lambdaFH * FHterm) / (1 + model.a * FHterm)
    term3 = (HFH + model.LV) / HFH * (model.muV + model.alpha * FFH + model.lambdaVH * HFH) / rV - 1

    term4 = H / model.b
    term5 = model.KV / model.rV * ((model.lambdaVH - model.alpha * FHterm) * (2 * HFH + model.LV) + 
                model.muV + model.alpha * F) - model.KV
    
    term5 += term4
    return term5^2 / term1 / term2 / term3 
end
function thresholdFVHi(model::modelAlleeEffect)
    _, _, HFH = eqFH(model)
    F = model.KF * (1 - model.muF / model.rF)
    H = model.a * F + model.c
    FHterm = model.KF * (model.omega * model.f + model.lambdaFH) / model.rF

    num = model.b * model.KV - model.b * model.KV / model.rV * ((model.lambdaVH - model.alpha * FHterm) * (2 * HFH + model.LV) + 
                model.muV + model.alpha * F)
    
    return num / H
end

function computeThresholds(model::modelAlleeEffect)
    dictThresholds = Dict{String, Float64}()
    
    dictThresholds["TH0"] = thresholdH(model)
    dictThresholds["TF0"] = thresholdF(model)
    dictThresholds["TVBeta0"] = thresholdV(model, model.beta)
    dictThresholds["TVi"] = thresholdVi(model)
    
    
    TFc = thresholdF(model ; H = model.c)
    dictThresholds["TFc"] = TFc
    if TFc > 1
        FFH = eqFH(model)[1]
        dictThresholds["THFH"] = thresholdH(model, F = FFH)
        dictThresholds["TVcFFH"] = thresholdV(model, model.c, F = FFH)
    else
        dictThresholds["THFH"] = 0
        dictThresholds["TVcFFH"] = 0
    end
    
    TFBeta = thresholdF(model ; H = model.beta)
    dictThresholds["TFBeta"] = TFBeta
    if TFBeta > 1
        FBeta = eqFHBeta(model)[1]
        dictThresholds["TVBetaFBeta"] = thresholdV(model, model.beta ; F = FBeta)
    else
        dictThresholds["TVBetaFBeta"] = 0
    end

    TVc0 = thresholdV(model, model.c)
    dictThresholds["TVc0"] = TVc0
    TDeltaVH = thresholdDeltaVH(model)
    dictThresholds["DeltaVH"] = TDeltaVH
    if TVc0 > 1
        _, VVH2, HVH2 = eqVH2(model)
        dictThresholds["TFHVH2"] = thresholdF(model; H = HVH2)
        dictThresholds["THVVH2"] = thresholdH(model; V = VVH2)
    else
        if TDeltaVH < 1
        dictThresholds["TFHVH2"] = 0
        dictThresholds["THVVH2"] = 0
        dictThresholds["TFHVH1"] = 0
        dictThresholds["THVVH1"] = 0
        dictThresholds["THHVH1"] = 0
        else
            _, VVH2, HVH2 = eqVH2(model)
            dictThresholds["TFHVH2"] = thresholdF(model; H = HVH2)
            dictThresholds["THVVH2"] = thresholdH(model; V = VVH2)
            _, VVH1, HVH1 = eqVH1(model)
            dictThresholds["TFHVH1"] = thresholdF(model; H = HVH1)
            dictThresholds["THVVH1"] = thresholdH(model; V = VVH1)
            dictThresholds["THVH1"] = thresholdHMin(model, HVH1)
        end
    end

    TFVH1 = thresholdFVH1(model)
    dictThresholds["TFVH1"] = TFVH1
    FFH, _, HFH = eqFH(model)
    TVFFHHFH = thresholdV(model, HFH ; F = FFH)
    dictThresholds["TVFFHHFH"] = TVFFHHFH
    DeltaFVH = thresholdDeltaFVH(model)
    dictThresholds["DeltaFVH"] = DeltaFVH
    if TFVH1 < 1  ### Cas 0 < 1 + ..(1-alpha t)/(1+at)
        if 1 < TVFFHHFH
            VFVH2 = eqFVH2(model)[2]
            dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
        else
            if DeltaFVH < 1
                dictThresholds["TFFVH2"] = 0
                dictThresholds["TFFVH1"] = 0
            else
                VFVH2 = eqFVH2(model)[2]
                dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
                VFVH1 = eqFVH1(model)[2]
                dictThresholds["TFFVH1"] = thresholdF(model; H = model.b * VFVH1 + model.c)
            end
        end
    else ### Cas 1 + ..(1-alpha t)/(1+at) < 0
        if TVFFHHFH < 1
            VFVH2 = eqFVH2(model)[2]
            dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
        else
            if DeltaFVH < 1
                dictThresholds["TFFVH2"] = 0
                dictThresholds["TFFVH1"] = 0
            else
                VFVH2 = eqFVH2(model)[2]
                dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
                VFVH1 = eqFVH1(model)[2]
                dictThresholds["TFFVH1"] = thresholdF(model; H = model.b * VFVH1 + model.c)
            end
        end
    end

    return dictThresholds
end

function interpretThresholds(model::modelAlleeEffect, thresholdsValues)
    listEq = ""
    if thresholdsValues["TF0"] < 1
        listEq = listEq * "(TE)"
    else
        listEq = listEq * "(F)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] < 1 && thresholdsValues["TVBeta0"] < 1
        listEq = listEq * "(H_β)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] > 1 && thresholdsValues["TVBetaFBeta"] < 1
        listEq = listEq * "(FH_β)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] < 1 && thresholdsValues["TVBeta0"] > 1
        listEq = listEq * "(VH_β)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] > 1 && thresholdsValues["TVBetaFBeta"] > 1
        listEq = listEq * "(FVH_β)"
    end
    if thresholdsValues["TH0"] > 1 && thresholdsValues["TFc"] < 1 && thresholdsValues["TVc0"] < 1
        listEq = listEq * "(H)"
    end
    if thresholdsValues["THFH"] > 1 && thresholdsValues["TFc"] > 1 && thresholdsValues["TVcFFH"] < 1
        listEq = listEq * "(FH)"
    end
    if thresholdsValues["TVc0"] > 1 && thresholdsValues["THVVH2"] > 1 && thresholdsValues["TFHVH2"] < 1 
        listEq = listEq * "(VH_2)"
    end
    if (thresholdsValues["TVc0"] < 1 && thresholdsValues["DeltaVH"] >= 1 && thresholdsValues["TVi"] > 1 &&
        thresholdsValues["THVVH2"] > 1 && thresholdsValues["TFHVH2"] < 1 )
        listEq = listEq * "(VH_2)"
    end
    if (thresholdsValues["TVc0"] < 1 && thresholdsValues["DeltaVH"] >= 1 && thresholdsValues["TVi"] > 1 &&
        thresholdsValues["THVVH1"] < 1 && thresholdsValues["TFHVH1"] < 1 && thresholdsValues["THVH1"] > 1)
    listEq = listEq * "(VH_1)"
end
    return listEq
end

function eqTrivial(model::modelAlleeEffect)
    return (0, 0, 0)
end

function eqHBeta(model::modelAlleeEffect)
    return (0, 0, model.beta)
end

function eqH(model::modelAlleeEffect)
    return (0, 0, model.c)
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

function eqFH(model::modelAlleeEffect)
    numFeq = 1 - model.muF/model.rF - model.c * (model.omega*model.f + model.lambdaFH)/model.rF
    denFeq = 1 + model.a * model.KF * (model.omega*model.f + model.lambdaFH)/model.rF
    Feq = model.KF * numFeq / denFeq
    return (Feq, 0, model.a * Feq + model.c)
end

function eqVHBeta(model::modelAlleeEffect)
    Veq = model.KV * (1 - 
        (model.beta + model.LV)/model.beta * (model.muV + model.beta * model.lambdaVH) / model.rV)
    return (0, Veq, model.beta)
end

function eqVH1(model::modelAlleeEffect)
    A = 1 + model.b * model.KV * model.lambdaVH / model.rV
    B = model.KV * (model.lambdaVH * (2 * model.c + model.LV) + model.muV)/model.rV - model.KV + 
            model.c / model.b
    C = model.KV * model.c / model.b *
        (model.c + model.LV) / model.c * (model.lambdaVH * model.c + model.muV) / model.rV -
        model.KV * model.c / model.b
    
    B = B / A
    C = C / A

    Veq = -Inf
    Heq = -Inf
    if B^2 - 4*C > 0
        Veq = -B / 2 - 1/2 * sqrt(B^2 - 4*C)
        Heq = model.b * Veq + model.c
    end 

    return (0, Veq, Heq)
end

function eqVH2(model::modelAlleeEffect)
    A = 1 + model.b * model.KV * model.lambdaVH / model.rV
    B = model.KV * (model.lambdaVH * (2 * model.c + model.LV) + model.muV)/model.rV - model.KV + 
            model.c / model.b
    C = model.KV * model.c / model.b *
        (model.c + model.LV) / model.c * (model.lambdaVH * model.c + model.muV) / model.rV -
        model.KV * model.c / model.b
    
    B = B / A
    C = C / A

    Veq = -Inf
    Heq = -Inf
    if B^2 - 4*C > 0
        Veq = -B / 2 + 1/2 * sqrt(B^2 - 4*C)
        Heq = model.b * Veq + model.c
    end 

    return (0, Veq, Heq)
end

function eqFVHBeta(model::modelAlleeEffect)
    Feq = model.KF * (1 - 
        (model.muF + model.beta * (model.omega * model.f + model.lambdaFH) / model.rF))
    Veq = model.KV * (1 - (model.beta + model.LV)/model.beta *
        (model.muV + model.beta * model.LV + model.alpha * Feq) / model.rV)
    return (Feq, Veq, model.beta)
end

function eqFVH1(model::modelAlleeEffect)
    FFH, _, HFH = eqFH(model) 

    FHterm = model.KF * (model.omega * model.f + model.lambdaFH) / model.rF
    A = 1 + model.KV * model.lambdaVH * model.b / model.rV * 
            (1 - model.alpha / model.lambdaVH * FHterm) /
            (1 + model.a * FHterm)

    B = model.KV / model.rV * ((model.lambdaVH - model.alpha * FHterm) * (model.LV + 2 * HFH) + 
                        model.muV + model.alpha * model.KF * (1 - model.muF / model.rF)) - model.KV
    B = B + (model.a * model.KF * (1 - model.muF / model.rF) + model.c) / model.b

    C = model.KV * (model.a * model.KF * (1 - model.muF / model.rF) + model.c) / model.b * 
                ((HFH + model.LV) / HFH * (model.muV + model.lambdaVH * HFH + model.alpha * FFH) / model.rV - 1)

    if B^2 - 4 * A * C >= 0
        Veq = -B / (2 * A) - sqrt(B^2 - 4 * A * C) / (2 * A)
    else
        Veq = -Inf
    end    
    Feq = FFH - model.b * FHterm / (1 + model.a * FHterm) * Veq
    Heq = HFH + model.b / (1 + model.a * FHterm) * Veq
    return (Feq, Veq, Heq)
end

function eqFVH2(model::modelAlleeEffect)
    FFH, _, HFH = eqFH(model) 

    FHterm = model.KF * (model.omega * model.f + model.lambdaFH) / model.rF
    A = 1 + model.KV * model.lambdaVH * model.b / model.rV * 
            (1 - model.alpha / model.lambdaVH * FHterm) /
            (1 + model.a * FHterm)

    B = model.KV * ((model.lambdaVH - model.alpha * FHterm) * (model.LV + 2 * HFH) / model.rV + 
                        model.muV/model.rV + model.alpha * model.KF * (1 - model.muF / model.rF)/model.rV - 1)
    B = B + (model.a * model.KF * (1 - model.muF / model.rF) + model.c) / model.b

    C = model.KV * (model.a * model.KF * (1 - model.muF / model.rF) + model.c) / model.b * 
                ((HFH + model.LV) / HFH * (model.muV + model.lambdaVH * HFH + model.alpha * FFH) / model.rV - 1)

    if B^2 - 4 * A * C >= 0
        Veq = -B / (2 * A) + sqrt(B^2 - 4 * A * C) / (2 * A)
    else
        Veq = -Inf
    end    
    Feq = FFH - model.b * FHterm / (1 + model.a * FHterm) * Veq
    Heq = HFH + model.b / (1 + model.a * FHterm) * Veq
    return (Feq, Veq, Heq)
end

function createTable(tresholdsName::Vector{String})
    table = DataFrame([[] for _ = tresholdsName], tresholdsName)
    for line in (2^length(tresholdsName)-1):-1:0
        newRow = digits(line, base =2)
        completeNewRow = [0 for _ in 1:(length(tresholdsName)-length(newRow))]
        append!(newRow, completeNewRow)
        push!(table, 2 .*newRow)
    end
    return table
end

function modelTable(model::modelAlleeEffect; biologicalModel = false)

    table = createTable(model.thresholdsName)

    function implication(TH0, THFH, THVVH1, THVVH2, TF0, TFBeta, 
        TFHVH1, TFHVH2, TVBeta0, TVBetaFBeta, TFc, TVc0, TVcFFH)::Bool

        cond1 = (TF0 >= TFc) && (TFHVH1 >= TFHVH2) && (TF0 >= TFBeta)  #TF is decreasing
        cond2 = ((TH0 < 1) && TFBeta <= TFc) || ((1 < TH0) && TFc <= TFBeta)
        cond3 = ((THVVH1 < 1) && TFBeta <= TFHVH1) || ((1 < THVVH1) && TFHVH1 <= TFBeta)
        cond4 = ((THVVH2 < 1) && TFBeta <= TFHVH2) || ((1 < THVVH2) && TFHVH2 <= TFBeta)
        cond5 = (TH0 <= THFH) && (TH0 <= THVVH1) && (THVVH1 <= THVVH2)      #TH is increasing
        cond6 = (TVBetaFBeta <= TVBeta0) && (TVcFFH <= TVc0) #TV(H, F) is decreasing wrt F
        return cond1 && cond2 && cond3 && cond4 && cond5 && cond6
    end

    function implicationBiological(TH0, THFH, THVVH1, THVVH2, TF0, TFBeta,
                        TFHVH1, TFHVH2, TVBeta0, TVBetaFBeta, TFc, TVc0, TVcFFH)::Bool    
        condBio1 = TFBeta >= 1                      # Beta < Bmax
        cond1 = (TF0 >= TFc) && (TFHVH1 >= TFHVH2) && (TF0 >= TFBeta)  #TF is decreasing
        cond2 = ((TH0 < 1) && TFBeta <= TFc) || ((1 < TH0) && TFc <= TFBeta)
        cond3 = ((THVVH1 < 1) && TFBeta <= TFHVH1) || ((1 < THVVH1) && TFHVH1 <= TFBeta)
        cond4 = ((THVVH2 < 1) && TFBeta <= TFHVH2) || ((1 < THVVH2) && TFHVH2 <= TFBeta)
        cond5 = (TH0 <= THFH) && (TH0 <= THVVH1) && (THVVH1 <= THVVH2)      #TH is increasing
        cond6 = (TVBetaFBeta <= TVBeta0) && (TVcFFH <= TVc0) #TV(H, F) is decreasing wrt F

        return condBio1 && cond1 && cond2 && cond3 && cond4 && cond5 && cond6
    end

    if biologicalModel
        modelTable = filter([:TH0, :THFH, :THVVH1, :THVVH2, :TF0, :TFBeta, :TFHVH1, :TFHVH2, 
            :TVBeta0, :TVBetaFBeta, :TFc, :TVc0, :TVcFFH] => implicationBiological, table)
    else
        modelTable = filter([:TH0, :THFH, :THVVH1, :THVVH2, :TF0, :TFBeta, :TFHVH1, :TFHVH2,
            :TVBeta0, :TVBetaFBeta, :TFc, :TVc0, :TVcFFH] => implication, table)
    end

    modelTable[!, :Equilibrium] .= ""
    for row in eachrow(modelTable)
        row[:Equilibrium] = interpretThresholds(model, row)
    end
    names = propertynames(modelTable)
    for ind in 1:length(names)-1
        columnSelected = names[1:end .!= ind]
        unique!(modelTable, columnSelected)
    end

    return modelTable

end

function computeBifurcationDiagram(model::modelAlleeEffect, 
        nameParam1::String, listParam1::Vector{Float64}, 
        nameParam2::String, listParam2::Vector{Float64})
    
    bifurcationMatrix = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)
    
    localParam = Dict([
                ("rF",model.rF), ("KF",model.KF), ("omega",model.omega), ("f",model.f), ("muF",model.muF), ("lambdaFH", model.lambdaFH),
                ("rV",model.rV), ("KV",model.KV), ("LV", model.LV), ("alpha",model.alpha), ("muV",model.muV), ("lambdaVH", model.lambdaVH),
                ("rH",model.rH), ("a", model.a), ("b", model.b), ("c", model.c), ("beta", model.beta)
                ])

    bifurcationMatrix[2:end,1] = listParam1

    for ind_col in 2:length(listParam2)+1
        bifurcationMatrix[1, ind_col] = listParam2[ind_col-1]
        localParam[nameParam2] = listParam2[ind_col-1]

        for ind_row in 2:length(listParam1)+1
            localParam[nameParam1] = listParam1[ind_row-1]
            localModel = modelAlleeEffect(localParam)
            thresholdsValues = computeThresholds(localModel)
            bifurcationMatrix[ind_row, ind_col] = interpretThresholds(localModel, thresholdsValues)
        end
    end
    bifurcationMatrix[1,1] = nameParam1 * "\\" * nameParam2
    return bifurcationMatrix
end
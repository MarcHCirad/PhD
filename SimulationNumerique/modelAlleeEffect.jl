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

function findEqNames(myString::String)
    openBracket = findall("(", myString)
    eqNames = []
    for ind in 1:length(openBracket) - 1
        beg = openBracket[ind][1] + 1
        en = openBracket[ind+1][1] - 2
        push!(eqNames, myString[beg:en])
    end
    en = openBracket[end][1] + 1
    push!(eqNames, myString[en:end-1])
    return eqNames
end

function computeEq(model::modelAlleeEffect, eqName::String)
    F, V, H = 0, 0, 0
    if eqName == "F"
        F, V, H = eqF(model)
    elseif eqName == "H_\beta"
        F, V, H = eqHBeta(model)
    elseif eqName == "FH_\beta"
        F, V, H = eqFHBeta(model)
    elseif eqName == "VH_\beta"
        F, V, H = eqVHBeta(model)
    elseif eqName == "FVH_\beta"
        F, V, H = eqFVHBeta(model)
    elseif eqName == "H"
        F, V, H = eqH(model)
    elseif eqName == "FH"
        F, V, H = eqFH(model)
    elseif eqName == "VH_2"
        F, V, H = eqVH2(model)
    elseif eqName == "VH1"
        F, V, H = eqVH1(model)
    elseif eqName == "FVH2"
        F, V, H = eqFVH2(model)
    elseif eqName == "FVH1"
        F, V, H = eqFVH1(model)
    end

    return F, V, H
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

function thresholdFVHAlpha(model::modelAlleeEffect)
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
    term3 = (HFH + model.LV) / HFH * (model.muV + model.alpha * FFH + model.lambdaVH * HFH) / model.rV - 1

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

function thresholdFVHq2(model::modelAlleeEffect, numberEq::Int)
    F,V,H = 0,0,0
    if numberEq == 2
        F, V, H = eqFVH2(model)
    elseif numberEq == 1
        F, V, H = eqFVH1(model)
    end
    rslt = model.rF * F / model.KF + model.rV * H / (H + model.LV) * V / model.KV + model.rH * (H / model.beta - 1)
    return rslt
end

function thresholdFVH2ndCompound(model::modelAlleeEffect, numberEq::Int)
    F,V,H = 0,0,0
    if numberEq == 2
        F, V, H = eqFVH2(model)
    elseif numberEq == 1
        F, V, H = eqFVH1(model)
    end

    FHterm = model.omega * model.f + model.lambdaFH
    HBeta = (H / model.beta - 1)

    term1 = (model.rF * F / model.KF + model.rV * H / (H + model.LV) * V / model.KV) * 
            (model.rF * F / model.KF + model.rH * HBeta) * 
            ( model.rV * H / (H + model.LV) * V / model.KV + model.rH * HBeta)
    term2 = model.alpha * FHterm * model.b * model.rH * HBeta * F * V
    term3 = FHterm * model.a * model.rH * F * HBeta * (model.rF * F / model.KF + model.rH * HBeta)
    term4 = model.b * model.rH * HBeta * (model.rH * HBeta + model.rV * H / (H + model.LV) * V / model.KV) * 
            ((model.alpha * F + model.muV)*model.LV - model.lambdaVH * H^2) * V / (H * (H + model.LV))

    return -term1 - term2 - term3 + term4

end

function computeThresholds(model::modelAlleeEffect)
    dictThresholds = Dict{String, Float64}()
    
    dictThresholds["TH0"] = thresholdH(model)
    dictThresholds["TF0"] = thresholdF(model)
    dictThresholds["TVBeta0"] = thresholdV(model, model.beta)
    dictThresholds["TVi"] = thresholdVi(model)
    dictThresholds["TFVHi"] = thresholdFVHi(model)
    
    
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
        dictThresholds["TFHVH1"] = 0
        dictThresholds["THVVH1"] = 0
        dictThresholds["THHVH1"] = 0
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

    TFVHAlpha = thresholdFVHAlpha(model)
    dictThresholds["TFVHAlpha"] = TFVHAlpha
    FFH, _, HFH = eqFH(model)
    TVFFHHFH = thresholdV(model, HFH ; F = FFH)
    dictThresholds["TVFFHHFH"] = TVFFHHFH
    DeltaFVH = thresholdDeltaFVH(model)
    dictThresholds["DeltaFVH"] = DeltaFVH
    if TFVHAlpha < 1  ### Cas 0 < 1 + ..(1-alpha t)/(1+at)
        if 1 < TVFFHHFH ##Seulement FVH2 existe
            FFVH2, VFVH2, _ = eqFVH2(model)
            dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
            dictThresholds["q22"] = thresholdFVHq2(model, 2)
            dictThresholds["2ndCompound2"] = thresholdFVH2ndCompound(model, 2)
            dictThresholds["THFVH2"] = thresholdH(model; F = FFVH2, V = VFVH2)
            dictThresholds["TFFVH1"] = 0
            dictThresholds["q21"] = -1
            dictThresholds["2ndCompound1"] = 1
            dictThresholds["THFVH1"] = 0
        else
            if DeltaFVH < 1 ##Aucun des deux n'existe
                dictThresholds["TFFVH2"] = 0
                dictThresholds["q22"] = -1
                dictThresholds["2ndCompound2"] = 1
                dictThresholds["THFVH2"] = 0
                dictThresholds["TFFVH1"] = 0
                dictThresholds["q21"] = -1
                dictThresholds["2ndCompound1"] = 1
                dictThresholds["THFVH1"] = 0
            else  ## FVH1,2 existent
                FFVH2, VFVH2, _ = eqFVH2(model)                
                dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
                dictThresholds["q22"] = thresholdFVHq2(model, 2)
                dictThresholds["2ndCompound2"] = thresholdFVH2ndCompound(model, 2)
                dictThresholds["THFVH2"] = thresholdH(model; F = FFVH2, V = VFVH2)
                FFVH1, VFVH1, _ = eqFVH1(model)               
                dictThresholds["TFFVH1"] = thresholdF(model; H = model.b * VFVH1 + model.c)
                dictThresholds["q21"] = thresholdFVHq2(model, 1)
                dictThresholds["2ndCompound1"] = thresholdFVH2ndCompound(model, 1)
                dictThresholds["THFVH1"] = thresholdH(model; F = FFVH1, V = VFVH1)
            end
        end
    else ### Cas 1 + ..(1-alpha t)/(1+at) < 0
        if TVFFHHFH < 1 ##Seulement FVH2 existe
            FFVH2, VFVH2, _ = eqFVH2(model)
            dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
            dictThresholds["q22"] = thresholdFVHq2(model, 2)
            dictThresholds["2ndCompound2"] = thresholdFVH2ndCompound(model, 2)
            dictThresholds["THFVH2"] = thresholdH(model; F = FFVH2, V = VFVH2)
            dictThresholds["TFFVH1"] = 0
            dictThresholds["q21"] = -1
            dictThresholds["2ndCompound1"] = 1
            dictThresholds["THFVH1"] = 0
        else
            if DeltaFVH < 1 ##Aucun des deux n'existe
                dictThresholds["TFFVH2"] = 0
                dictThresholds["q22"] = -1
                dictThresholds["2ndCompound2"] = 1
                dictThresholds["THFVH2"] = 0
                dictThresholds["TFFVH1"] = 0
                dictThresholds["q21"] = -1
                dictThresholds["2ndCompound1"] = 1
                dictThresholds["THFVH1"] = 0
            else ## FVH1,2 existent
                FFVH2, VFVH2, _ = eqFVH2(model)
                dictThresholds["TFFVH2"] = thresholdF(model; H = model.b * VFVH2 + model.c)
                dictThresholds["q22"] = thresholdFVHq2(model, 2)
                dictThresholds["2ndCompound2"] = thresholdFVH2ndCompound(model, 2)
                dictThresholds["THFVH2"] = thresholdH(model; F = FFVH2, V = VFVH2)
                FFVH1, VFVH1, _ = eqFVH1(model)
                dictThresholds["TFFVH1"] = thresholdF(model; H = model.b * VFVH1 + model.c)
                dictThresholds["q21"] = thresholdFVHq2(model, 1)
                dictThresholds["2ndCompound1"] = thresholdFVH2ndCompound(model, 1)
                dictThresholds["THFVH1"] = thresholdH(model; F = FFVH1, V = VFVH1)
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
        listEq = listEq * "(H_\beta)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] > 1 && thresholdsValues["TVBetaFBeta"] < 1
        listEq = listEq * "(FH_\beta)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] < 1 && thresholdsValues["TVBeta0"] > 1
        listEq = listEq * "(VH_\beta)"
    end
    if thresholdsValues["TH0"] < 1 && thresholdsValues["TFBeta"] > 1 && thresholdsValues["TVBetaFBeta"] > 1
        listEq = listEq * "(FVH_\beta)"
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
    if (thresholdsValues["TFVHAlpha"] < 1 && thresholdsValues["TVFFHHFH"] > 1 && thresholdsValues["TFFVH2"] > 1
        && thresholdsValues["2ndCompound2"] < 0 && thresholdsValues["q22"] > 0 && thresholdsValues["THFVH2"] > 1)
        listEq = listEq * "(FVH_2)"
    end
    if (thresholdsValues["TFVHAlpha"] < 1 && thresholdsValues["TVFFHHFH"] < 1 && thresholdsValues["DeltaFVH"] >= 1 && 
        thresholdsValues["TFVHi"] > 1 && thresholdsValues["TFFVH2"] > 1 && 
        thresholdsValues["2ndCompound2"] < 0 && thresholdsValues["q22"] > 0 && thresholdsValues["THFVH2"] > 1)
        listEq = listEq * "(FVH_2)"
    end
    if (thresholdsValues["TFVHAlpha"] < 1 && thresholdsValues["TVFFHHFH"] < 1 && thresholdsValues["DeltaFVH"] >= 1 && 
        thresholdsValues["TFVHi"] > 1 && thresholdsValues["TFFVH1"] > 1 && 
        thresholdsValues["2ndCompound1"] < 0 && thresholdsValues["q21"] > 0 && thresholdsValues["THFVH1"] < 1)
        listEq = listEq * "(FVH_1)"
    end
    if (thresholdsValues["TFVHAlpha"] > 1 && thresholdsValues["TVFFHHFH"] < 1 && thresholdsValues["TFFVH2"] > 1
        && thresholdsValues["2ndCompound2"] > 0 && thresholdsValues["q22"] > 0 && thresholdsValues["THFVH2"] < 1)
        listEq = listEq * "(FVH_2)"
    end
    if (thresholdsValues["TFVHAlpha"] > 1 && thresholdsValues["TVFFHHFH"] > 1 && thresholdsValues["DeltaFVH"] >= 1 && 
        thresholdsValues["TFVHi"] < 1 && thresholdsValues["TFFVH2"] > 1 && 
        thresholdsValues["2ndCompound2"] < 0 && thresholdsValues["q22"] > 0 && thresholdsValues["THFVH2"] < 1)
        listEq = listEq * "(FVH_2)"
    end
    if (thresholdsValues["TFVHAlpha"] > 1 && thresholdsValues["TVFFHHFH"] > 1 && thresholdsValues["DeltaFVH"] >= 1 && 
        thresholdsValues["TFVHi"] < 1 && thresholdsValues["TFFVH1"] > 1 && 
        thresholdsValues["2ndCompound1"] < 0 && thresholdsValues["q21"] > 0 && thresholdsValues["THFVH1"] > 1)
        listEq = listEq * "(FVH_1)"
    end

    listEq = replace(listEq, ")("=>"), (")

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
        (model.muF + model.beta*(model.omega*model.f + model.lambdaFH)) / model.rF)
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
        (model.muF + model.beta * (model.omega * model.f + model.lambdaFH)) / model.rF)
    Veq = model.KV * (1 - (model.beta + model.LV)/model.beta *
        (model.muV + model.beta * model.lambdaVH + model.alpha * Feq) / model.rV)
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

    Veq = -Inf
    if B^2 - 4 * A * C >= 0
        Veq = -B / (2 * A) + sqrt(B^2 - 4 * A * C) / (2 * A)
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
        nameParam2::String, listParam2::Vector{Float64} ;
        eqValues = false,
        Values = "max")
    
    bifurcationMatrix = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)
    if eqValues
        Fval = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)
        Vval = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)
        Hval = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)
    end
    
    localParam = Dict([
                ("rF",model.rF), ("KF",model.KF), ("omega",model.omega), ("f",model.f), ("muF",model.muF), ("lambdaFH", model.lambdaFH),
                ("rV",model.rV), ("KV",model.KV), ("LV", model.LV), ("alpha",model.alpha), ("muV",model.muV), ("lambdaVH", model.lambdaVH),
                ("rH",model.rH), ("a", model.a), ("b", model.b), ("c", model.c), ("beta", model.beta)
                ])

    bifurcationMatrix[2:end,1] = listParam1
    if eqValues
        Fval[2:end,1] = listParam1
        Vval[2:end,1] = listParam1
        Hval[2:end,1] = listParam1
    end

    for ind_col in 2:length(listParam2)+1
        bifurcationMatrix[1, ind_col] = listParam2[ind_col-1]
        if eqValues
            Fval[1, ind_col] = listParam2[ind_col-1]
            Vval[1, ind_col] = listParam2[ind_col-1]
            Hval[1, ind_col] = listParam2[ind_col-1]
        end

        localParam[nameParam2] = listParam2[ind_col-1]

        for ind_row in 2:length(listParam1)+1
            localParam[nameParam1] = listParam1[ind_row-1]
            localModel = modelAlleeEffect(localParam)
            thresholdsValues = computeThresholds(localModel)
            eqNames = interpretThresholds(localModel, thresholdsValues)
            bifurcationMatrix[ind_row, ind_col] = eqNames

            if eqValues
                Feq, _, _ = eqF(localModel)
                listEq = findEqNames(eqNames)
                allEqs = [computeEq(localModel, eq) for eq in listEq]
                # Fval[ind_row, ind_col] = (sum([eq[1] for eq in allEqs]) - Feq) / (length(allEqs)-1)
                Fval[ind_row, ind_col] = sum([eq[1] for eq in allEqs]) / length(allEqs)
                Vval[ind_row, ind_col] = sum([eq[2] for eq in allEqs]) / length(allEqs)
                Hval[ind_row, ind_col] = sum([eq[3] for eq in allEqs]) / length(allEqs)
            end

        end
    end
    bifurcationMatrix[1,1] = nameParam1 * "\\" * nameParam2
    if eqValues
        Fval[1,1] = nameParam1 * "\\" * nameParam2
        Vval[1,1] = nameParam1 * "\\" * nameParam2
        Hval[1,1] = nameParam1 * "\\" * nameParam2
        return bifurcationMatrix, Fval, Vval, Hval
    end
    
    return bifurcationMatrix

end
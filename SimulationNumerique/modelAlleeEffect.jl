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
            "THVVH2",
            "TF0", #
            "TFBeta", #
            "TFHVH2",
            "TVBeta0", #
            "TVBetaFBeta", #
            "TFc", #
            "TVc0", #
            "TVcFFH", #
            "DeltaVH", #
            "TVi"
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

function thresholdF(model::modelAlleeEffect ; H = 0)
    return model.rF / (model.muF + H * (model.omega * model.f + model.lambdaFH))
end

function thresholdV(model::modelAlleeEffect, H::Float64 ; F = 0)
    fact1 = H / (model.LV + H)
    fact2 = model.rV / (model.lambdaVH * H + model.alpha * F + model.muV)
    return fact1 * fact2
end

function thresholdVi(model::modelAlleeEffect)
    fraction1 = (model.b * model.KV - model.c) / (model.b * model.KV)
    fraction2 = model.rV / (model.lambdaVH * (2 * model.c + model.LV) + model.muV)
    return fraction1 * fraction2
end

function thresholdDeltaVH(model::modelAlleeEffect)
    TVc = thresholdV(model, model.c)
    TVi = thresholdVi(model)
    fraction1 = 1 / (1 + model.b * model.KV * model.lambdaVH / model.rV)
    fraction2 = (model.b * model.KV - model.c)^2 / (4 * model.b * model.c * model.KV)
    return fraction1 * fraction2 / (1 / TVc - 1) * (1 / TVi - 1)^2
end

function computeThresholds(model::modelAlleeEffect)
    dictThresholds = Dict{String, Float64}()
    
    dictThresholds["TH0"] = thresholdH(model)
    dictThresholds["TF0"] = thresholdF(model)
    dictThresholds["TVc0"] = thresholdV(model, model.c)
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
    
    TDeltaVH = thresholdDeltaVH(model)
    TVc0 = thresholdV(model, model.c)
    dictThresholds["DeltaVH"] = TDeltaVH
    if TDeltaVH > 1 || TVc0 > 1
        _, VVH2, HVH2 = eqVH2(model)
        dictThresholds["TFHVH2"] = thresholdF(model; H = HVH2)
        dictThresholds["THVVH2"] = thresholdH(model; V = VVH2)
    else
        dictThresholds["TFHVH2"] = 0
        dictThresholds["THVVH2"] = 0
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
    if thresholdsValues["THVVH2"] > 1 && thresholdsValues["TFHVH2"] < 1 && thresholdsValues["TVc0"] > 1
        listEq = listEq * "(VH_2)"
    end
    if (thresholdsValues["THVVH2"] > 1 && thresholdsValues["TFHVH2"] < 1 && thresholdsValues["TVc0"] < 1 &&  
            thresholdsValues["DeltaVH"] > 1 && thresholdsValues["TVi"] > 1)
        listEq = listEq * "(VH_2)"
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
    B = (model.b * model.KV - model.c) / model.b * 
            (model.b * model.KV / (model.b * model.KV - model.c) * 
            (model.lambdaVH * (2*model.c + model.LV) + model.muV) / model.rV - 1)
    B = B / A 
    C = model.KV * model.c / model.b *
        ((model.c + model.LV) / model.c * (model.lambdaVH * model.c + model.muV) / model.rV - 1)
    C = C / A

    Veq = -1/2 * B - 1/2 * sqrt(B^2 - 4*C)
    Heq = model.b * Veq + model.c

    return (0, Veq, Heq)
end

function eqVH2(model::modelAlleeEffect)
    A = 1 + model.b * model.KV * model.lambdaVH / model.rV
    B = (model.b * model.KV - model.c) / model.b * 
            (model.b * model.KV / (model.b * model.KV - model.c) * 
            (model.lambdaVH * (2*model.c + model.LV) + model.muV) / model.rV - 1)
    B = B / A 
    C = model.KV * model.c / model.b *
        ((model.c + model.LV) / model.c * (model.lambdaVH * model.c + model.muV) / model.rV - 1)
    C = C / A

    Veq = -1/2 * B + 1/2 * sqrt(B^2 - 4*C)
    Heq = model.b * Veq + model.c

    return (0, Veq, Heq)
end

function eqFVHBeta(model::modelAlleeEffect)
    Feq = model.KF * (1 - 
        (model.muF + model.beta * (model.omega * model.f + model.lambdaFH) / model.rF))
    Veq = model.KV * (1 - (model.beta + model.LV)/model.beta *
        (model.muV + model.beta * model.LV + model.alpha * Feq) / model.rV)
    return (Feq, Veq, model.beta)
end

function createTable(tresholdsName::Vector{String})
    table = DataFrame([[] for _ = tresholdsName], tresholdsName)
    for line in (2^length(tresholdsName)-1):-1:0
    # for line in (2^4):-1:0
        newRow = digits(line, base =2)
        completeNewRow = [0 for _ in 1:(length(tresholdsName)-length(newRow))]
        append!(newRow, completeNewRow)
        push!(table, 2 .*newRow)
    end
    return table
end

function modelTable(model::modelAlleeEffect)

    table = createTable(model.thresholdsName)
    function implication(TH0, THFH, THVVH2, TF0, TFBeta, TFHVH2, TVBeta0, TVBetaFBeta, TFc, TVc0, TVcFFH)::Bool
        
        condOpt = (TFBeta > 1)
        cond1 = (TF0 >= TFBeta) && (TF0 >= TFc) && (TF0 >= TFHVH2) && (TFc >= TFHVH2)#TF is decreasing
        cond2 = !((TH0 < 1) && !(TFBeta <= TFc))   #If TH(0) < 1, TF(Beta) < TF(c)
        cond3 = !((1 < TH0) && !(TFc <= TFBeta))   #If 1 < TH(0), TF(c) < TF(Beta)
        cond4 = (TH0 <= THFH) && (TH0 <= THFH) && (TH0 <= THVVH2)      #TH is increasing
        cond5 = (TVBetaFBeta <= TVBeta0) && (TVcFFH <= TVc0) #TV(H, F) is decreasing wrt F
        cond6 = !((TH0 < 1) && !(TVBeta0 <= TVc0))    #If TH(0) < 1, TV(Beta, 0) < TV(c, 0)
        cond7 = !((1 < TH0) && !(TVc0 <= TVBeta0))    #If 1 < TH(0), TV(c) < TV(Beta)

        return cond1 && cond2 && cond3 && cond4 && cond5 && cond6 && cond7
    end
    modelTable = filter([:TH0, :THFH, :THVVH2, :TF0, :TFBeta, :TFHVH2, 
            :TVBeta0, :TVBetaFBeta, :TFc, :TVc0, :TVcFFH] => implication, table)

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
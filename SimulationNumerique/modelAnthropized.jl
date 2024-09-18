struct modelAnthropized <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    delta0F::Float64
    delta1::Float64
    LF::Float64
    muF::Float64
    lambdaFH::Float64

    rV::Float64
    delta0V::Float64
    LV::Float64
    muV::Float64
    lambdaVH::Float64

    rH::Float64
    a::Float64
    b::Float64
    c::Float64
    beta::Float64
    
    function modelAnthropized(modelParam::Dict{String, Float64})
        rF, delta0F, delta1, LF, muF, lambdaFH = modelParam["rF"], modelParam["delta0F"], modelParam["delta1"], modelParam["LF"],
                        modelParam["muF"], modelParam["lambdaFH"]
        rV, delta0V, LV, muV, lambdaVH = modelParam["rV"], modelParam["delta0V"], modelParam["LV"],
                        modelParam["muV"], modelParam["lambdaVH"]
        rH, a, b, c = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"]
        beta = modelParam["beta"]
        if c <= beta
            error("parameter c is lower than parameter β")
        end

        variablesNames = ["F_A", "V_A", "H_A"]
       
        new(variablesNames,
            rF, delta0F, delta1, LF, muF, lambdaFH,
            rV, delta0V, LV, muV, lambdaVH, 
            rH, a, b, c, beta)
    end
end

function equationModel(model::modelAnthropized, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = (model.rF * H / (H + model.LF) * F - model.delta0F / (1 + model.delta1*H) * F^2 - 
                model.muF * F - model.lambdaFH * F * H)

    dV = (model.rV * H / (H + model.LV) * V - model.delta0V * V*V - 
                model.muV * V - model.lambdaVH * V * H)

    dH = (model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * (H/model.beta - 1) * H)

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
    if (model.LV == 0)
        fact1 = 1
    else
        fact1 = H / (model.LV + H)
    end
    fact2 = model.rV / (model.lambdaVH * H + model.muV)
    return fact1 * fact2
end

function thresholdDeltaFH(model::modelAnthropized)
    TFc = thresholdF(model, model.c)
    a0 = model.c * (1 + model.c * model.delta1) * (1/TFc - 1)
    a1 = -model.a * (1 + 2*model.c * model.delta1 
            - model.muF / model.rF * (1 + 2 * model.delta1 * model.c + model.LF * model.delta1)
            - model.lambdaFH / model.rF * (model.LF * (1 + 2*model.delta1 * model.c) + model.c * (2 + 3 * model.delta1 * model.c))
            - model.delta0F * (model.c + model.LV) / (model.rF * model.a))
    a2 = -model.a^2 * (model.delta1 - 
            (model.muF * model.delta1 + model.lambdaFH * (1 + 3 * model.c * model.delta1 + model.LF * model.delta1)) / model.rF
            - model.delta0F / (model.a * model.rF))
    a3 = model.a^3 * model.delta1 * model.lambdaFH / model.rF
    
    Delta = 18 * a0 * a1 * a2 * a3 - 4 * a2^3 * a0 + a2^2 * a1^2 - 4 * a1^3 * a3 - 27 * a0^2 * a3^2

    return Delta   
end

function thresholdA1F(model::modelAnthropized)
    den = model.muF * (1 + model.delta1 * model.LF / (1 + 2 * model.c * model.delta1))
    den += model.lambdaFH * (model.LF + model.c * (2 + 3 * model.delta1 * model.c) / (1 + 2 * model.c * model.delta1))
    den += model.delta0F * (model.c + model.LF) / (model.a * (1 + 2 * model.c * model.delta1))
    
    return model.rF / den
end

function thresholdA2F(model::modelAnthropized)
    den = (model.muF + model.lambdaFH * (3 * model.c + model.LF + 1/model.delta1) + 
                model.delta0F / (model.a * model.delta1))
    
    return model.rF / den
end

function thresholdDeltaVH(model::modelAnthropized)
    A = 1 + model.b * model.lambdaVH / model.delta0V
    B = (model.rV / model.delta0V * (
        (model.lambdaVH * (2*model.c + model.LV) + model.muV) / model.rV -1) +
        (model.c + model.LV) / model.b)
    C = (model.rV * model.c / (model.delta0V * model.b) * 
        ((model.c + model.LV)/model.c * (model.muV + model.c * model.lambdaVH) / model.rV - 1))

    Delta = B^2 - 4 * A * C
    if (B < 0) && (Delta >= 0)
        return 2
    else
        return 0
    end
end

function equilibriumTE(model::modelAnthropized)
    return [0,0,0]
end

function equilibriumV(model::modelAnthropized)
    if model.LV == 0
        V = model.rV / model.delta0V * (1 - model.muV/model.rV)
        return [0, V, 0]
    else
        error("try to compute equilibriumV while does not exist")
    end
end

function equilibriumHBeta(model::modelAnthropized)
    return [0, 0, model.beta]
end

function equilibriumFHBeta(model::modelAnthropized)
    TFBeta = thresholdF(model, model.beta)
    if TFBeta > 1
        factor1 = model.rF * (1 + model.delta1 * model.beta)/model.delta0F
        factor2 = model.beta / (model.LF + model.beta)
        return [factor1 * factor2 * (1 - 1/TFBeta), 0, model.beta]
    else
        error("try to compute equilibriumFHBeta while does not exist")
    end
end

function equilibriumVHBeta(model::modelAnthropized)
    TVBeta = thresholdV(model, model.beta)
    if TVBeta > 1
        factor1 = model.rV / model.delta0V
        factor2 = model.beta / (model.LV + model.beta)
        return [0, factor1 * factor2 * (1 - 1/TVBeta), model.beta]
    else
        error("try to compute equilibriumVHBeta while does not exist")
    end
end

function equilibriumFVHBeta(model::modelAnthropized)
    TFBeta = thresholdF(model, model.beta)
    TVBeta = thresholdV(model, model.beta)
    if (TVBeta > 1) && (TFBeta > 1)
        factor1 = model.rF * (1 + model.delta1 * model.beta)/model.delta0F
        factor2 = model.beta / (model.LF + model.beta)
        factor3 = model.rF * (1 + model.delta1 * model.beta)/model.delta0F
        factor4 = model.beta / (model.LF + model.beta)
        return [factor1 * factor2 * (1 - 1/TFBeta), factor3 * factor4 * (1 - 1/TVBeta), model.beta]
    else
        error("try to compute equilibriumFVHBeta while does not exist")
    end
end

function equilibriumH(mode::modelAnthropized)
    return [0,0, model.c]
end

function equilibriumFH(model::modelAnthropized)
    TFc = thresholdF(model, model.c)
    a0 = model.c * (1 + model.c * model.delta1) * (1/TFc - 1)
    a1 = -model.a * (1 + 2*model.c * model.delta1 
            - model.muF / model.rF * (1 + 2 * model.delta1 * model.c + model.LF * model.delta1)
            - model.lambdaFH / model.rF * (model.LF * (1 + 2*model.delta1 * model.c) + model.c * (2 + 3 * model.delta1 * model.c))
            - model.delta0F * (model.c + model.LV) / (model.rF * model.a))
    a2 = -model.a^2 * (model.delta1 - 
            (model.muF * model.delta1 + model.lambdaFH * (1 + 3 * model.c * model.delta1 + model.LF * model.delta1)) / model.rF
            - model.delta0F / (model.a * model.rF))
    a3 = model.a^3 * model.delta1 * model.lambdaFH / model.rF
    
    Delta = 18 * a0 * a1 * a2 * a3 - 4 * a2^3 * a0 + a2^2 * a1^2 - 4 * a1^3 * a3 - 27 * a0^2 * a3^2
    polynomialF = Polynomial([a0, a1, a2, a3])
    rootF = PolynomialRoots.roots(coeffs(polynomialF))

    eqFH = []
    for F in rootF
        if isreal(F) && (0 < real(F))
            H = model.a * real(F) + model.c
            push!(eqFH, (real(F), 0, H))
        end
    end

    return eqFH
end

function equilibriumVH(model::modelAnthropized)
    A = 1 + model.b * model.lambdaVH / model.delta0V
    B = (model.rV / model.delta0V * (
        (model.lambdaVH * (2*model.c + model.LV) + model.muV) / model.rV -1) +
        (model.c + model.LV) / model.b)
    C = (model.rV * model.c / (model.delta0V * model.b) * 
        ((model.c + model.LV)/model.c * (model.muV + model.c * model.lambdaVH) / model.rV - 1))

    Delta = B^2 - 4 * A * C
    if (Delta >= 0)
        V = -B / (2 * A) + sqrt(Delta) / (2 * A)
        if V > 0
            H = model.b * V + model.c
            return (0, V, H)
        end
    end
    error("try to compute equilibriumVH while it does not exist")
end

# function equilibriumFVH(model::modelAnthropized)
#     A = 1 + model.a * model.KF * model.lambdaFH / model.rF +
#              model.b * model.KV * model.lambdaVH / model.rV
#     B = model.a * model.KF * (1 - (model.muF + model.LF * model.lambdaFH) / model.rF) +
#         model.b * model.KV * (1 - (model.muV + model.LV * model.lambdaVH) / model.rV) +
#         model.c
#     C = model.a * model.KF * model.LF * model.muF / model.rF +
#         model.b * model.KV * model.LV * model.muV / model.rV

#     Delta = B^2 - 4 * A * C
#     if (Delta >= 0) && (B > 0)
#         H = B / (2 * A) + sqrt(Delta) / (2 * A)
#         TF = thresholdF(model, H)
#         TV = thresholdV(model, H)
#         if (TF > 1) && (TV > 1)
#             F = model.KF * (1 - 1 / TF)
#             V = model.KV * (1 - 1 / TV)
#             return (F, V, H)
#         end
#     end
#     error("try to compute equilibriumFVH while it does not exist")       
# end

function computeThresholdsVH(model::modelAnthropized)
    dictThresholds = Dict{String, Float64}()
    dictThresholds["TV0"] = thresholdV(model, 0.)
    dictThresholds["TVc"] = thresholdV(model, model.c)
    dictThresholds["DeltaVH"] = thresholdDeltaVH(model)

    return dictThresholds
end

function computeThresholdsFH(model::modelAnthropized)
    dictThresholds = Dict{String, Float64}()
    dictThresholds["TFc"] = thresholdF(model, model.c)
    dictThresholds["DeltaFH"] = thresholdDeltaFH(model)
    dictThresholds["TFA1"] = thresholdA1F(model)
    dictThresholds["TFA2"] = thresholdA2F(model)

    return dictThresholds
end   

function interpretThresholdsVH(model::modelAnthropized, thresholdsValues::Dict{String, Float64})
    
    if model.LV > 0
        listEq = "(TE)"
        if thresholdsValues["TVc"] > 1
            listEq = listEq * "(V_AH_A)"
        else
            listEq = listEq * "(H_A)"
            if thresholdsValues["DeltaVH"] > 1
                listEq = listEq * "(V_AH_A)"
            end
        end
    else
        listEq = ""
        if thresholdsValues["TV0"] < 1
            listEq = listEq * "(TE)" * "(H_A)"
        else
            listEq = listEq * "(V_A)"
            if thresholdsValues["TVc"] < 1
                listEq = listEq * "(H_A)"
            else
                listEq = listEq * "(V_AH_A)"
            end
        end
    end
    return listEq
end

function interpretThresholdsFH(model::modelAnthropized, thresholdsValues::Dict{String, Float64})
    listEq = "(TE)"
    if thresholdsValues["TFc"] < 1
        if (thresholdsValues["DeltaFH"] > 0) 
            if (thresholdsValues["TFA1"] < 1) && (thresholdsValues["TFA2"] < 1)
                listEq = listEq * "(H_A)"
            else
                listEq = listEq * "(H_A)" * "(V_AH_A)"
            end
        else
            listEq = listEq * "(H_A)"
        end
    else
        if (thresholdsValues["DeltaFH"] < 0)
            listEq = listEq * "(V_AH_A)"
        else
            if (thresholdsValues["TFA1"] < 1) && (1 < thresholdsValues["TFA2"])
                listEq = listEq * "(V_AH_A_1)" * "(V_AH_A)"
            else
                listEq = listEq * "(V_AH_A)"
            end
        end
    end
    return listEq
end

function computeThresholds(model::modelAnthropized)
    dictThresholds = Dict{String, Float64}()

    dictThresholds["LF"] = model.LF
    dictThresholds["TF0"] = thresholdF(model, 0.)
    TFc = thresholdF(model, model.c)
    dictThresholds["TFc"] = TFc
    TVc = thresholdV(model, model.c)
    dictThresholds["TVc"] = TVc

    deltaFH = 0
    if TFc < 1
        deltaFH = thresholdDeltaFH(model)
    end
    dictThresholds["DeltaFH"] = deltaFH

    TVH_FH = 2
    if (1 < TFc) || (TFc < 1 && 1 <= deltaFH)
        H_FH = equilibriumFH(model)[3]
        TVH_FH = thresholdV(model, H_FH)
    end
    dictThresholds["TVH_FH"] = TVH_FH

    deltaVH = 0
    if TVc < 1
        deltaVH = thresholdDeltaVH(model)
    end
    dictThresholds["DeltaVH"] = deltaVH

    TFH_VH = 2
    if (1 < TVc) || (TVc < 1 && 1 <= deltaVH)
        H_VH = equilibriumVH(model)[3]
        TFH_VH = thresholdF(model, H_VH)
    end
    dictThresholds["TFH_VH"] = TFH_VH
    return dictThresholds
end


function interpretThresholds(model::modelAnthropized, thresholdsValues::Dict{String, Float64})
    listEq = ""

    if thresholdsValues["LF"] == 0
        if thresholdsValues["TF0"] < 1
            listEq = listEq * "(TE)"
        elseif 1 < thresholdsValues["TF0"]
            listEq = listEq * "(F_A)"
        end
    else
        listEq = listEq * "(TE)"
    end

    if thresholdsValues["TFc"] < 1 && thresholdsValues["TVc"] < 1
        listEq = listEq * "(H_A)"
    end

    if 1 < thresholdsValues["TFc"] && thresholdsValues["TVH_FH"] < 1
        listEq = listEq * "(F_AH_A)"
    end

    if thresholdsValues["TFc"] < 1 && 1 <= thresholdsValues["DeltaFH"] &&
            thresholdsValues["TVH_FH"] < 1
        listEq = listEq * "(F_AH_A)"
    end

    if 1 < thresholdsValues["TVc"] && thresholdsValues["TFH_VH"] < 1
        listEq = listEq * "(V_AH_A)"
    end

    if thresholdsValues["TVc"] < 1 && 1 <= thresholdsValues["DeltaVH"] &&
            thresholdsValues["TFH_VH"] < 1
        listEq = listEq * "(V_AH_A)"
    end

    return listEq
end

function computeBifurcationDiagramVH(model::modelAnthropized, 
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
                ("rF",model.rF), ("LF", model.LF), ("delta0F",model.delta0F), ("delta1",model.delta1), ("muF",model.muF), ("lambdaFH", model.lambdaFH),
                ("rV",model.rV), ("LV", model.LV), ("delta0V",model.delta0V), ("muV",model.muV), ("lambdaVH", model.lambdaVH),
                ("rH",model.rH), ("a", model.a), ("b", model.b), ("c", model.c), ("beta", model.beta)
                ])

    bifurcationMatrix[2:end,1] = listParam1


    for ind_col in 2:length(listParam2)+1
        bifurcationMatrix[1, ind_col] = listParam2[ind_col-1]
        
        localParam[nameParam2] = listParam2[ind_col-1]

        for ind_row in 2:length(listParam1)+1
            localParam[nameParam1] = listParam1[ind_row-1]
            localModel = modelAnthropized(localParam)
            thresholdsValues = computeThresholdsVH(localModel)
            eqNames = interpretThresholdsVH(localModel, thresholdsValues)
            bifurcationMatrix[ind_row, ind_col] = eqNames
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

function computeBifurcationDiagramFH(model::modelAnthropized, 
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
                ("rF",model.rF), ("LF", model.LF), ("delta0F",model.delta0F), ("delta1",model.delta1), ("muF",model.muF), ("lambdaFH", model.lambdaFH),
                ("rV",model.rV), ("LV", model.LV), ("delta0V",model.delta0V), ("muV",model.muV), ("lambdaVH", model.lambdaVH),
                ("rH",model.rH), ("a", model.a), ("b", model.b), ("c", model.c), ("beta", model.beta)
                ])

    bifurcationMatrix[2:end,1] = listParam1


    for ind_col in 2:length(listParam2)+1
        bifurcationMatrix[1, ind_col] = listParam2[ind_col-1]
        
        localParam[nameParam2] = listParam2[ind_col-1]

        for ind_row in 2:length(listParam1)+1
            localParam[nameParam1] = listParam1[ind_row-1]
            localModel = modelAnthropized(localParam)
            thresholdsValues = computeThresholdsFH(localModel)
            eqNames = interpretThresholdsFH(localModel, thresholdsValues)
            bifurcationMatrix[ind_row, ind_col] = eqNames
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
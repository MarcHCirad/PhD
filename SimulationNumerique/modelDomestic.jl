struct modelDomestic <: mathematicalModel
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
    
    function modelDomestic(modelParam::Dict{String, Float64})
        rF, delta0F, delta1, LF, muF, lambdaFH = modelParam["rF"], modelParam["delta0F"], modelParam["delta1"], modelParam["LF"],
                        modelParam["muF"], modelParam["lambdaFH"]
        rV, delta0V, LV, muV, lambdaVH = modelParam["rV"], modelParam["delta0V"], modelParam["LV"],
                        modelParam["muV"], modelParam["lambdaVH"]
        rH, a, b, c = modelParam["rH"], modelParam["a"], modelParam["b"], modelParam["c"]
        beta = modelParam["beta"]
        if c <= beta
            error("parameter c is lower than parameter β")
        end

        variablesNames = ["F_D", "V_D", "H_D"]
       
        new(variablesNames,
            rF, delta0F, delta1, LF, muF, lambdaFH,
            rV, delta0V, LV, muV, lambdaVH, 
            rH, a, b, c, beta)
    end
end

function equationModel(model::modelDomestic, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = ((model.rF - model.lambdaFH) * H / (H + model.LF) * F - model.delta0F / (1 + model.delta1*H) * F^2 - 
                model.muF * F)

    dV = ((model.rV - model.lambdaVH) * H / (H + model.LV) * V - model.delta0V * V^2 - 
                model.muV * V)

    dH = (model.rH * (1 - H / (model.a * F + model.b * V + model.c)) * (H/model.beta - 1) * H)

    return [dF, dV, dH]
end

function thresholdF(model::modelDomestic, H::Float64)
    fact1 = (model.LF + H) / H
    threshold = model.rF / (model.muF * fact1 + model.lambdaFH)
    return threshold
end

function thresholdV(model::modelDomestic, H::Float64)
    fact1 = (model.LV + H) / H
    threshold = model.rV / (model.muV * fact1 + model.lambdaVH)
    return threshold
end

function thresholdDeltaFH(model::modelDomestic)
    TFc = thresholdF(model, model.c)
    A = 1 + model.rF * ((model.lambdaFH + model.muF) / model.rF - 1) * model.delta1 * model.a / model.delta0F
    C = (model.c / model.a * model.rF * (1 + model.delta1 * model.c) / model.delta0F * 
            (1/TFc - 1))
    B = (model.rF * ((model.lambdaFH + model.muF) / model.rF - 1) * (1 + 2 * model.delta1 * model.c) / model.delta0F +
            (model.muF * model.LF * model.delta1) / model.delta0F + (model.LF + model.c) / model.a)
    Delta = B^2 - 4 * A * C

    if A < 0
        if C > 0
            return 2
        else
            if (Delta >= 0) && (B > 0)
                return 2
            end
        end
    else
        if C < 0
            return 2
        else
            if (Delta >= 0) && (B < 0)
                return 2
            end
        end
    end
    return 0   
end

function thresholdDeltaVH(model::modelDomestic)
    TVc = thresholdV(model, model.c)
    A = 1

    if TVc > 1
        return 2
    else
        C = (model.c / model.b * model.rV / model.delta0V * (1/TVc - 1))
        B = (model.rV * ((model.lambdaVH + model.muV) / model.rV - 1) / model.delta0V +
                (model.LV + model.c) / model.b)
        Delta = B^2 - 4 * A * C

        if (Delta >= 0) && (B < 0)
            return 2
        end
    end
    return 0   
end

function thresholdEqFVH(model::modelDomestic)
    q3 = (model.delta0F * model.delta0V * (1 - model.a * model.rF * model.delta1 / model.delta0F * 
            (1 - (model.lambdaFH + model.muF)/ model.rF)) )
    q0 = model.LV * model.LF * model.delta0F * model.delta0V * (model.a * model.muF / model.delta0F +
            model.b * model.muV / model.delta0V - model.c)
    
    q2 = (model.delta0F * model.delta0V * (model.LF + model.LV - model.c )
            - model.delta0V * model.a * model.rF * (1 - (model.lambdaFH + model.muF) / model.rF) * (1 + model.LV * model.delta1)
            + model.delta0V * model.a * model.muF * model.LF * model.delta1
            - model.delta0F * model.b * model.rV * (1 - (model.lambdaVH + model.muV) / model.rV)
    )

    q1 = (model.delta0F * model.delta0V * model.LF * model.LV * (1 - model.c / model.LF - model.c / model.LV)
        - model.delta0V * model.a * model.rF * (1 - (model.lambdaFH + model.muF) / model.rF) * model.LV
        + model.delta0V * model.a * model.muF * (1 + model.delta1 * model.LV) * model.LF
        - model.delta0F * model.b * model.rV * (1 - (model.lambdaVH + model.muV) / model.rV) * model.LF
        + model.delta0F * model.b * model.muV * model.LV)
    
    polynomeH = Polynomial([q0, q1, q2, q3])

    eq = 0

    rootH = PolynomialRoots.roots(coeffs(polynomeH))
    for H in rootH
        if isreal(H) && (real(H) > 0)
            if (thresholdF(model, real(H)) > 1) && (thresholdV(model, real(H)) > 1)
                eq = eq + 1
            end
        end
    end

    if eq == 3
        return 3
    elseif eq > 0
        return 2
    else
        return 0
    end
end

function equilibriumTE(model::modelDomestic)
    return [0,0,0]
end

function equilibriumH(model::modelDomestic)
    return [0,0, model.c]
end

function equilibriumFH(model::modelDomestic)
    TFc = thresholdF(model, model.c)
    A = 1 + model.rF * ((model.lambdaFH + model.muF) / model.rF - 1) * model.delta1 * model.a / model.delta0F
    C = (model.c / model.a * model.rF * (1 + model.delta1 * model.c) / model.delta0F * 
            (1/TFc - 1))
    B = (model.rF * ((model.lambdaFH + model.muF) / model.rF - 1) * (1 + 2 * model.delta1 * model.c) / model.delta0F +
            (model.muF * model.LF * model.delta1) / model.delta0F + (model.LF + model.c) / model.a)
    Delta = B^2 - 4 * A * C

    exist = thresholdDeltaFH(model)
    if exist > 1
        Feq = (-B + sqrt(Delta)) / (2*A)
        Heq = model.a * Feq + model.c
        return [Feq, 0, Heq]
    else
        return [-10,-10,-10]
    end
end

function equilibriumVH(model::modelDomestic)
    TVc = thresholdV(model, model.c)
    A = 1

    C = (model.c / model.b * model.rV / model.delta0V * (1/TVc - 1))
    B = (model.rV * ((model.lambdaVH + model.muV) / model.rV - 1) / model.delta0V +
            (model.LV + model.c) / model.b)
    Delta = B^2 - 4 * A * C

    exist = thresholdDeltaVH(model)
    if exist > 1
        Veq = (-B + sqrt(Delta)) / (2*A)
        Heq = model.b * Veq + model.c
        return [0, Veq, Heq]
    else
        return [-10,-10,-10]
    end
end

function equilibriumFVH(model::modelDomestic)
    q3 = (model.delta0F * model.delta0V * (1 - model.a * model.rF * model.delta1 / model.delta0F * 
            (1 - (model.lambdaFH + model.muF)/ model.rF)) )
    q0 = model.LV * model.LF * model.delta0F * model.delta0V * (model.a * model.muF / model.delta0F +
            model.b * model.muV / model.delta0V - model.c)
    
    q2 = (model.delta0F * model.delta0V * (model.LF + model.LV - model.c )
            - model.delta0V * model.a * model.rF * (1 - (model.lambdaFH + model.muF) / model.rF) * (1 + model.LV * model.delta1)
            + model.delta0V * model.a * model.muF * model.LF * model.delta1
            - model.delta0F * model.b * model.rV * (1 - (model.lambdaVH + model.muV) / model.rV)
    )

    q1 = (model.delta0F * model.delta0V * model.LF * model.LV * (1 - model.c / model.LF - model.c / model.LV)
        - model.delta0V * model.a * model.rF * (1 - (model.lambdaFH + model.muF) / model.rF) * model.LV
        + model.delta0V * model.a * model.muF * (1 + model.delta1 * model.LV) * model.LF
        - model.delta0F * model.b * model.rV * (1 - (model.lambdaVH + model.muV) / model.rV) * model.LF
        + model.delta0F * model.b * model.muV * model.LV)
    
    polynomeH = Polynomial([q0, q1, q2, q3])

    nbrEq = 0
    eq = []

    rootH = PolynomialRoots.roots(coeffs(polynomeH))
    for H in rootH
        if isreal(H) && (real(H) > 0)
            TFH = thresholdF(model, real(H))
            TVH = thresholdV(model, real(H))
            if (TFH > 1) && (TVH > 1)
                Heq = real(H)
                Feq = model.rF * (1 + model.delta1*Heq) / model.delta0F * Heq / (Heq + model.LF) * (1 - 1 / TFH)
                Veq = model.rV / model.delta0V * Heq / (Heq + model.LV) * (1 - 1 / TVH)
                eq = push!(eq, [Feq, Veq, Heq])
            end
        end
    end
    sort!(eq, by = x -> x[3])
    if length(eq) == 3
        return [eq[1], eq[3]]
    elseif length(eq) > 0
        return [eq[end]]
    end
    return [[-10,-10,-10]]
end


function computeThresholds(model::modelDomestic)
    dictThresholds = Dict{String, Float64}()
    dictThresholds["TFc"] = thresholdF(model, model.c)
    dictThresholds["TVc"] = thresholdV(model, model.c)

    dictThresholds["DeltaFH"] = thresholdDeltaFH(model)
    if dictThresholds["DeltaFH"] > 1
        H_FH = equilibriumFH(model)[3]
        dictThresholds["TVH_FH"] = thresholdV(model, H_FH)
    else
        dictThresholds["TVH_FH"] = 2
    end

    dictThresholds["DeltaVH"] =  thresholdDeltaVH(model)
    if dictThresholds["DeltaVH"] > 1
        H_VH = equilibriumVH(model)[3]
        dictThresholds["TFH_VH"] = thresholdF(model, H_VH)
    else
        dictThresholds["TFH_VH"] = 2
    end

    dictThresholds["FVH"] = thresholdEqFVH(model)

    return dictThresholds
end


function interpretThresholds(model::modelDomestic, thresholdsValues::Dict{String, Float64})
    listEq = "(TE)"
    if (thresholdsValues["TFc"] < 1) && (thresholdsValues["TVc"] < 1)
        listEq = listEq * "(H_D)"
    end
    if (thresholdsValues["DeltaFH"] > 1) && (thresholdsValues["TVH_FH"] < 1)
        listEq = listEq * "(F_DH_D)"
    end
    if (thresholdsValues["DeltaVH"] > 1) && (thresholdsValues["TFH_VH"] < 1)
        listEq = listEq * "(V_DH_D)"
    end
    if (thresholdsValues["FVH"] > 1)
        if thresholdsValues["FVH"] == 3
            listEq = listEq * "(F_DV_DH_D,1)" * "(F_DV_DH_D,3)"
        else
            listEq = listEq * "(F_DV_DH_D,3)"
        end
    end
    return listEq
end

function computeBifurcationDiagram(model::modelDomestic, 
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
            localModel = modelDomestic(localParam)
            thresholdsValues = computeThresholds(localModel)
            eqNames = interpretThresholds(localModel, thresholdsValues)
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
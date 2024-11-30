struct modelHunter <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    KF::Float64
    lambdaFWH::Float64

    rH::Float64
    aW::Float64
    c::Float64
    beta::Float64
    mA::Float64
    mW::Float64

    m::Float64
        
    function modelHunter(modelParam::Dict{String, Float64})
        rF, lambdaFWH, KF = modelParam["rF"], modelParam["lambdaFWH"], modelParam["KF"]
        rH, aW, c = modelParam["rH"], modelParam["aW"], modelParam["c"]
        beta, mA, mW = modelParam["beta"], modelParam["mA"], modelParam["mW"]

        variablesNames = ["H_A", "F_W", "H_W"]
        m  = mA / mW
        new(variablesNames, rF, KF, lambdaFWH, rH, aW, c, beta, mA, mW, m)
    end
end

function equationModel(model::modelHunter, variables::Vector{Float64})
    HA, FW, HW = variables[1], variables[2], variables[3]
    dHA = model.rH * (1 - HA / (model.aW * FW + model.c)) * (HA/model.beta - 1) * HA - model.mA * HA + model.mW * HW
    dFW = model.rF * (1 - FW / model.KF) * FW - model.lambdaFWH * FW * HW
    dHW = model.mA * HA - model.mW * HW
    return [dHA, dFW, dHW]
end

function thresholdFW(model::modelHunter, H::Float64)
    return model.rF / (model.m * model.lambdaFWH * H) 
end

function thresholdStability(model::modelHunter)
    HA, FW, _ = equilibriumHAFWHW(model)
    RH = model.rH * (HA / model.beta - 1)
    RF = model.rF * FW / model.KF

    Tr = -(RH + model.mA + RF + model.mW)
    det = - model.mW * RH * RF * (1 + model.aW * model.lambdaFWH * model.KF * model.m / model.rF)
    q1 = (RH + model.mA) * RF + (RH + RF) * model.mW

    return -Tr * q1 + det
    
end

function equilibriumFW(model::modelHunter)
    return [0, model.KF, 0]
end

function equilibriumHAHW(model::modelHunter)
    return [model.c, 0, model.m * model.c]
end

function equilibriumHAFWHW(model::modelHunter)
    numFW = 1 - model.m * model.lambdaFWH * model.c / model.rF
    denFW = 1 + model.aW * model.KF * model.m * model.lambdaFWH / model.rF
    FW = model.KF * numFW / denFW
    return [model.aW * FW + model.c, FW, model.m * (model.aW * FW + model.c)]
end

function computeThresholds(model::modelHunter)
    dictThresholds = Dict{String, Float64}()
    dictThresholds["TFc"] = thresholdFW(model, model.c)
    dictThresholds["Stab"] = thresholdStability(model)
    dictThresholds["TFaKV"] = thresholdFW(model, model.aW * model.KF)

    return dictThresholds
end

function interpretThresholds(model::modelHunter, thresholdsValues::Dict{String, Float64})
    
    listEq = "(F_W)"
    if thresholdsValues["TFc"] < 1
        listEq = listEq * "(H_AH_W)"
    else
        if thresholdsValues["Stab"] > 0
            listEq = listEq * "(H_AF_WH_W)"
        end
    end
    return listEq
end

function computeBifurcationDiagram(model::modelHunter, 
    nameParam1::String, listParam1::Vector{Float64}, 
    nameParam2::String, listParam2::Vector{Float64})
    println("model m : ", model.m)
    bifurcationMatrix = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)

    localParam = Dict([
                ("rF",model.rF), ("KF", model.KF), ("lambdaFWH",model.lambdaFWH),
                ("rH",model.rH), ("aW", model.aW), ("c", model.c), 
                ("beta", model.beta), ("mA", model.mA), ("mW", model.mW)
                ])

    bifurcationMatrix[2:end,1] = listParam1


    for ind_col in 2:length(listParam2)+1
        bifurcationMatrix[1, ind_col] = listParam2[ind_col-1]
        
        localParam[nameParam2] = listParam2[ind_col-1]

        for ind_row in 2:length(listParam1)+1
            localParam[nameParam1] = listParam1[ind_row-1]
            localModel = modelHunter(localParam)
            thresholdsValues = computeThresholds(localModel)
            eqNames = interpretThresholds(localModel, thresholdsValues)
            bifurcationMatrix[ind_row, ind_col] = eqNames
        end
    end
    bifurcationMatrix[1,1] = nameParam1 * "\\" * nameParam2

    return bifurcationMatrix
end
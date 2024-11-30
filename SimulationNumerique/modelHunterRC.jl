struct modelHunter <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    KF::Float64
    alpha::Float64
    lambdaFWH::Float64

    e::Float64
    c::Float64
    muD::Float64
    rD::Float64
    mD::Float64
    mW::Float64

    m::Float64
        
    function modelHunter(modelParam::Dict{String, Float64})
        rF, lambdaFWH, KF, alpha = modelParam["rF"], modelParam["lambdaFWH"], modelParam["KF"], modelParam["alpha"]
        e, c = modelParam["e"], modelParam["c"]
        muD, rD, mD, mW = modelParam["muD"], modelParam["rD"], modelParam["mD"], modelParam["mW"]
        alpha = modelParam["alpha"]
        variablesNames = ["H_D", "F_W", "H_W"]
        m  = mD / mW
        new(variablesNames, rF, KF, alpha, lambdaFWH, e, c, muD, rD, mD, mW, m)
    end
end

function equationModel(model::modelHunter, variables::Vector{Float64})
    HD, FW, HW = variables[1], variables[2], variables[3]
    dHD = model.c + model.e * model.lambdaFWH * HW * FW + (model.rD- model.muD) * HD - model.mD * HD + model.mW * HW
    dFW = model.rF * (1 - FW / (model.KF * (1 - model.alpha))) * FW - model.lambdaFWH * FW * HW
    dHW = model.mD * HD - model.mW * HW
    return [dHD, dFW, dHW]
end

function thresholdHD(model::modelHunter)
    return (model.muD - model.rD) * model.rF / (model.lambdaFWH * model.m * model.c)
end


function equilibriumHDFW(model::modelHunter)
    A = model.e * model.rF / (model.KF * (1 - model.alpha))
    B = -(model.e * model.rF + (model.muD - model.rD) * model.rF / (model.lambdaFWH * model.m * model.KF * (1-model.alpha)))
    C = (model.muD - model.rD) * model.rF / (model.lambdaFWH * model.m) - model.c

    Delta = B^2 - 4 *A *C
    if Delta >= 0
        F1 = -B/(2*A) - sqrt(Delta) / (2*A)
        H1 = model.rF / (model.lambdaFWH * model.m) * (1 - F1 / (model.KF * (1 - model.alpha)))
        return [(H1, F1, model.m * H1)]
    else
        return [(0,0,0)]
    end

end

function computeDeltaStab(model::modelHunter ; lambdaFWH = -1., alpha = -1.0)
    if lambdaFWH < 0
        lambdaFWH = model.lambdaFWH
    end
    if alpha < 0
        alpha = model.alpha
    end

    KFalpha = model.KF * (1-alpha)
    A = model.e * model.rF / KFalpha
    B = -(model.e * model.rF + (model.muD - model.rD) * model.rF / (lambdaFWH * model.m * KFalpha))
    C = (model.muD - model.rD) * model.rF / (lambdaFWH * model.m) - model.c
    
    DeltaF = B^2 - 4 *A *C


    @assert(DeltaF > 0 && C > 0)
    Feq = -B / (2*A) - sqrt(DeltaF) / (2*A)

    Tr = model.rD - model.muD - model.mD - model.mW - model.rF * Feq / KFalpha
    det = -model.mD * lambdaFWH * Feq * sqrt(DeltaF)
    a1 = ((model.muD - model.rD + model.mD + model.mW)*model.rF * Feq / KFalpha + 
            model.mD * model.e * lambdaFWH * ((model.muD - model.rD) / (model.m * lambdaFWH * model.e) - Feq) )

    @assert (a1 > 0 && det < 0 && Tr < 0)

    return -Tr * a1 + det, Feq
end


function thresholdStabHDFW(model::modelHunter)
    eq = equilibriumHDFW(model)
    if eq == [(0,0,0)]
        println("Pas d'équilibre de coexistence")
        return 0
    else
        eq1 = eq[1]
        Feq = eq1[2]
        if Feq < 0
            println("lambda trop grand, le petit equilibre est négatif")
            return 0
        else
            DeltaStab = computeDeltaStab(model)
        end
    end
end

function computeDeltaStabMatrix(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})

    deltaStabMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    deltaStabMatrix[2:end,1] = listAlpha

    FeqMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    FeqMatrix[2:end,1] = listAlpha


    for ind_col in 2:length(listLambdaFWH)+1
        deltaStabMatrix[1, ind_col] = listLambdaFWH[ind_col-1]
        FeqMatrix[1, ind_col] = listLambdaFWH[ind_col-1]

        lambdaFWH = listLambdaFWH[ind_col-1]
        for ind_row in 2:length(listAlpha)+1
            alpha = listAlpha[ind_row-1]
            Delta, Feq = computeDeltaStab(model; lambdaFWH=lambdaFWH, alpha=alpha)
            rsltDelta = Delta
            # if Delta < 0
            #     rsltDelta = 0
            #     Feq = 0.001
            # elseif Delta == 0
            #     rsltDelta = 1
            # else
            #     rsltDelta = 2
            # end

            deltaStabMatrix[ind_row, ind_col] = rsltDelta
            FeqMatrix[ind_row, ind_col] = Feq
        end
    end

    deltaStabMatrix[1,1] = "alpha" * "\\" * "lambdaFWH"
    FeqMatrix[1,1] = "alpha" * "\\" * "lambdaFWH"

    return deltaStabMatrix, FeqMatrix
end


# function thresholdFW(model::modelHunter, H::Float64)
#     return model.rF / (model.m * model.lambdaFWH * H) 
# end

# function thresholdStability(model::modelHunter)
#     HD, FW, _ = equilibriumHDFWHW(model)
#     RH = model.rH * (HD / model.muD - 1)
#     RF = model.rF * FW / model.KF

#     Tr = -(RH + model.mD + RF + model.mW)
#     det = - model.mW * RH * RF * (1 + model.e * model.lambdaFWH * model.KF * model.m / model.rF)
#     q1 = (RH + model.mD) * RF + (RH + RF) * model.mW

#     return -Tr * q1 + det
    
# end

# function equilibriumFW(model::modelHunter)
#     return [0, model.KF, 0]
# end

# function equilibriumHDHW(model::modelHunter)
#     return [model.c, 0, model.m * model.c]
# end

# function equilibriumHDFWHW(model::modelHunter)
#     """
#     When alpha = 0
#     """
#     numFW = 1 - model.m * model.lambdaFWH * model.c / model.rF
#     denFW = 1 + model.e * model.KF * model.m * model.lambdaFWH / model.rF
#     FW = model.KF * numFW / denFW
#     return [model.e * FW + model.c, FW, model.m * (model.e * FW + model.c)]
# end

# function equilibriumHDFWHW1(model::modelHunter)
#     """
#     When alpha > 0
#     """
    
#     A = model.e^2 * model.KF * model.alpha * model.lambdaFWH * model.m / model.rF
#     B = (-1 -model.e * model.KF * (model.alpha + model.lambdaFWH * model.m / model.rF) +
#         2 * model.e * model.KF * model.alpha * model.lambdaFWH * model.m * model.c / model.rF)
#     C = model.KF * (1 - model.c * model.alpha) * (1 - model.c * model.lambdaFWH * model.m / model.rF)

#     Delta = B^2 - 4*A*C
#     if Delta >= 0
#         FW1 = -B / (2*A) - sqrt(Delta) / (2*A)
#         FW2 = -B / (2*A) + sqrt(Delta) / (2*A)

#         return [[model.e * FW1 + model.c, FW1, model.m * (model.e * FW1 + model.c)],
#                     [model.e * FW2 + model.c, FW2, model.m * (model.e * FW2 + model.c)]]
#     else
#         return [-10,-10,-10]
#     end
# end

# function computeThresholds(model::modelHunter)
#     dictThresholds = Dict{String, Float64}()
#     dictThresholds["TFc"] = thresholdFW(model, model.c)
#     dictThresholds["Stab"] = thresholdStability(model)
#     dictThresholds["TFaKV"] = thresholdFW(model, model.e * model.KF)

#     return dictThresholds
# end

# function interpretThresholds(model::modelHunter, thresholdsValues::Dict{String, Float64})
    
#     listEq = "(F_W)"
#     if thresholdsValues["TFc"] < 1
#         listEq = listEq * "(H_DH_W)"
#     else
#         if thresholdsValues["Stab"] > 0
#             listEq = listEq * "(H_DF_WH_W)"
#         end
#     end
#     return listEq
# end

# function computeBifurcationDiagram(model::modelHunter, 
#     nameParam1::String, listParam1::Vector{Float64}, 
#     nameParam2::String, listParam2::Vector{Float64})
#     bifurcationMatrix = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)

#     localParam = Dict([
#                 ("rF",model.rF), ("KF", model.KF), ("lambdaFWH",model.lambdaFWH),
#                 ("alpha", model.alpha),
#                 ("rH",model.rH), ("e", model.e), ("c", model.c), 
#                 ("muD", model.muD), ("mD", model.mD), ("mW", model.mW)
#                 ])

#     bifurcationMatrix[2:end,1] = listParam1


#     for ind_col in 2:length(listParam2)+1
#         bifurcationMatrix[1, ind_col] = listParam2[ind_col-1]
        
#         localParam[nameParam2] = listParam2[ind_col-1]

#         for ind_row in 2:length(listParam1)+1
#             localParam[nameParam1] = listParam1[ind_row-1]
#             localModel = modelHunter(localParam)
#             thresholdsValues = computeThresholds(localModel)
#             eqNames = interpretThresholds(localModel, thresholdsValues)
#             bifurcationMatrix[ind_row, ind_col] = eqNames
#         end
#     end
#     bifurcationMatrix[1,1] = nameParam1 * "\\" * nameParam2

#     return bifurcationMatrix
# end
abstract type Model end

struct modelFVH <: Model
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

    e::Float64
    muH::Float64

    
    function modelFVH(modelParam::Dict{String, Float64})
        rV, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["KV"], modelParam["alpha"], modelParam["muV"], modelParam["lambdaVH"]
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"], modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        e, muH = modelParam["e"], modelParam["muH"]
       
        new(rF, KF, omega, f, muF, lambdaFH, rV, KV, alpha, muV, lambdaVH, e, muH)
    end
end

function equationModel(model::modelFVH, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    dF = model.rF * (1 - F / model.KF) * F - model.omega * model.f * F - model.muF * F - model.lambdaFH * F * H
    dV = model.rV * (1 - V / model.KV) * V - model.alpha * V * F - model.muV * V - model.lambdaVH * V * H
    dH = model.e * (model.lambdaFH * F + model.lambdaVH * V) * H - model.muH * H^2
    
    return [dF, dV, dH]
end

function stabilityTresholds(model::modelFVH)
    """
    Compute the stabilities tresholds for equationModel
    """
    R0V = model.rV / model.muV
    R0F = model.rF / (model.muF + model.omega*model.f)
    TF = (model.lambdaVH / model.muV) * (model.omega*model.f + model.muF) / model.lambdaFH * (R0F - 1) / (R0V - 1) * (1 + model.muH / (model.e*model.lambdaVH^2)*model.rV/model.KV)
    TV = (model.muV / model.lambdaVH) * model.lambdaFH / (model.omega*model.f + model.muF) * (R0V - 1) / (R0F - 1) * (1 + model.muH / (model.e*model.lambdaFH^2)*model.rF/model.KF) / (1 + model.alpha*model.muH/(model.e * model.lambdaFH * model.lambdaVH))

    return Dict([("R_0^V",R0V), ("R_0^F",R0F), ("T^F",TF), ("T^V",TV)])
end

function interpretStabilityTresholds(model::modelFVH)
    """
    Return a string of stable equilibrium
    """
    tresholds = stabilityTresholds(model)
    R0V, R0F, TV, TF = tresholds["R_0^V"], tresholds["R_0^F"], tresholds["T^V"], tresholds["T^F"]
    if (R0V < 1) & (R0F < 1)
        return "TE"
    elseif (R0V > 1) & (R0F > 1)
        if (TF < 1) & (TV < 1)
            return "VH & FH"
        elseif (TF > 1) & (TV > 1)
            return "FVH"
        elseif (TF < 1) & (TV > 1)
            return "VH"
        elseif (TF > 1) & (TV < 1)
            return "FH"
        end
    elseif (R0V > 1) & (R0F < 1)
        if (TF < 1)
            return "VH"
        end
    elseif (R0V < 1) & (R0F > 1)
        if (TV < 1)
            return "FH"
        end
    end
    return "_"
end

function computeEquilibria(model::modelFVH)
    """
    Compute all possible stable equilibria of equationModel
    """
    equilibriaName = interpretStabilityTresholds(model)
    rslt = Dict{String, Vector{Float64}}()

    deltaF = model.rF - model.muF - model.omega*model.f
    deltaV = model.rV - model.muV

    if (equilibriaName == "VH") || (equilibriaName == "VH & FH")
            V_VH = 1 / (model.rV / model.KV + model.e*model.lambdaVH^2/model.muH)*deltaV
            H_VH = model.e * model.lambdaVH / model.muH * V_VH
            rslt["VH"] = [0., V_VH, H_VH]

    elseif (equilibriaName == "FH") || (equilibriaName == "VH & FH")
            F_FH = 1 / (model.rF / model.KF + model.e * model.lambdaFH^2 / model.muH)*deltaF
            H_FH = model.e * model.lambdaFH / model.muH * F_FH
            rslt["FH"] = [F_FH, 0., H_FH]

    elseif equilibriaName == "FVH"
            V_0 = model.KV / model.rV * (deltaV - model.lambdaVH / model.lambdaFH * deltaF)
            V_1 = model.lambdaVH / model.lambdaFH * model.rF / model.rV * model.KV / model.KF - model.alpha * model.KV / model.rV
            H_0 = deltaF / model.lambdaFH
            H_1 = model.rF / (model.lambdaFH * model.KF)
            F_FVH = (model.muH / model.e * H_0 - model.lambdaVH * V_0) / (model.lambdaFH + model.lambdaVH * V_1 + model.muH / model.e *H_1)
            V_FVH = V_0 + V_1 * F_FVH
            H_FVH = H_0 - H_1 * F_FVH
            rslt["FVH"] = [F_FVH, V_FVH, H_FVH]
    end
    return rslt
end

function computeBifurcationDiagram(model::modelFVH, listLambdaFH::Vector{Float64}, listLambdaVH::Vector{Float64}, nameFile::String)
    
    stability = Matrix{Any}(undef, length(listLambdaFH)+1, length(listLambdaVH)+1)
    localParam = Dict([("rV",model.rV), ("KV",model.KV), ("alpha",model.alpha), ("muV",model.muV),
                    ("rF",model.rF), ("KF",model.KF), ("omega",model.omega), ("f",model.f), ("muF",model.muF), 
                    ("e",model.e), ("muH", model.muH)])

    stability[2:end,1] = listLambdaFH

    for ind_col in 2:length(listLambdaVH)+1
        stability[1, ind_col] = listLambdaVH[ind_col-1]
        localParam["lambdaVH"] = listLambdaVH[ind_col-1]

        for ind_row in 2:length(listLambdaFH)+1
            localParam["lambdaFH"] = listLambdaFH[ind_row-1]
            localModel = modelFVH(localParam)
            stability[ind_row, ind_col] = interpretStabilityTresholds(localModel)
        end
    end
    stability[1,1] = "FH \\ VH"
    CSV.write(nameFile, DataFrame(stability, :auto))
end

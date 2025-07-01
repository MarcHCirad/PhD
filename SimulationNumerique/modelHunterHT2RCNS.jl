include("modelHunterHT2RC.jl")

struct hunterHT2NS <: numericalModel
    variablesNames::Vector{String}
    mathModel::modelHunterHT2

    t0::Float64
    tf::Float64
    dt::Float64
    n::Int64

    phiDt::Float64

    result::Matrix{Float64}

    function hunterHT2NS(mathModel::modelHunterHT2, numericalParam::Dict{String, Float64}, 
        initialValues::Dict{String, Float64})
        t0, tf, dt = numericalParam["t0"], numericalParam["tf"], numericalParam["dt"]
        n = Int((tf-t0)/dt)

        Q = max(mathModel.muD - mathModel.fD + mathModel.mD, mathModel.mW)
        phiDt = (1- exp(-Q * dt)) / Q

        result = Matrix{Float64}(undef, 7, n+1)
        result[:,1] = [t0, initialValues["HD0"], initialValues["FW0"], mathModel.m * initialValues["HD0"], 
                initialValues["HD0"], initialValues["FW0"], mathModel.m * initialValues["HD0"]]

        
        variablesNames = mathModel.variablesNames
        pushfirst!(variablesNames, "time")

        new(variablesNames, mathModel, t0, tf, dt, n, phiDt, result)
        end
end

function numericalScheme(model::hunterHT2NS, variables::Vector{Float64})
    mModel = model.mathModel

    HD2D, FW2D, HW2D = variables[1], variables[2], variables[3]

    FW2Dn = 1 + model.phiDt * mModel.rF * (1 - mModel.alpha) * (1 + mModel.beta * mModel.m * HD2D)
    denFW2Dn = (1 + model.phiDt * (mModel.rF / mModel.KF * (1 + mModel.beta * mModel.m * HD2D) * FW2D) + 
                model.phiDt * mModel.lambdaFWH * mModel.m * HD2D / (1+ mModel.theta * mModel.lambdaFWH * FW2D))
    FW2Dn = FW2Dn / denFW2Dn * FW2D

    HD2Dn = ((1- model.phiDt * (mModel.muD - mModel.fD)) * HD2D + 
        model.phiDt * (mModel.I + 
            (mModel.e * mModel.lambdaFWH * FW2Dn / (1+ mModel.theta * mModel.lambdaFWH * FW2D)) * mModel.m * HD2D))
    
    HW2Dn = mModel.m * HD2D


    HD, FW, HW = variables[4], variables[5], variables[6]

    FWn = 1 + model.phiDt * mModel.rF * (1 - mModel.alpha) * (1 + mModel.beta * HW)
    denFWn = (1 + model.phiDt * (mModel.rF / mModel.KF * (1 + mModel.beta * HW) * FW) + 
                model.phiDt * mModel.lambdaFWH * HW / (1+ mModel.theta * mModel.lambdaFWH * FW))
    FWn = FWn / denFWn * FW

    HDn = ((1- model.phiDt * (mModel.muD - mModel.fD + mModel.mD)) * HD + 
        model.phiDt * (mModel.I + (mModel.mW + 
            mModel.e * mModel.lambdaFWH * FWn / (1+ mModel.theta * mModel.lambdaFWH * FW)) * HW))
    
    HWn = (1 - model.phiDt * mModel.mW) * HW + model.phiDt * mModel.mD * HD

    return [HD2Dn, FW2Dn, HW2Dn, HDn, FWn, HWn]
end

function solveModel(model::hunterHT2NS)
    for ind in 1:model.n
        model.result[1, ind+1] = model.result[1, ind] + model.dt
        model.result[2:7, ind+1] = numericalScheme(model, model.result[2:7, ind])
    end
end
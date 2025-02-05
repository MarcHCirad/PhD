include("modelHunterRC.jl")

struct hunterNS <: numericalModel
    variablesNames::Vector{String}
    mathModel::modelHunter

    t0::Float64
    tf::Float64
    dt::Float64
    n::Int64

    phiDt::Float64

    result::Matrix{Float64}

    function hunterNS(mathModel::modelHunter, numericalParam::Dict{String, Float64}, 
        initialValues::Dict{String, Float64})
        t0, tf, dt = numericalParam["t0"], numericalParam["tf"], numericalParam["dt"]
        n = Int((tf-t0)/dt)

        Q = max(mathModel.muD + mathModel.mD - mathModel.fD, mathModel.mW)
        phiDt = (1- exp(-Q * dt)) / Q
        

        result = Matrix{Float64}(undef, 4, n+1)
        result[:,1] = [t0, initialValues["HD0"], initialValues["FW0"], initialValues["HW0"]]
        
        variablesNames = mathModel.variablesNames
        pushfirst!(variablesNames, "time")

        new(variablesNames, mathModel, t0, tf, dt, n, phiDt,result)
        end
end

function numericalScheme(model::hunterNS, variables::Vector{Float64})
    mModel = model.mathModel
    HD, FW, HW = variables[1], variables[2], variables[3]

    dFW = (1 + mModel.rF * model.phiDt) * FW / 
            (1 + model.phiDt * (mModel.rF / (mModel.KF * (1-mModel.alpha)) * FW + mModel.lambdaFWH * HW))
    dHW = (1-model.phiDt * mModel.mW) * HW + model.phiDt * mModel.mD * HD
    dHD = (1 - model.phiDt * (mModel.muD - mModel.fD + mModel.mD))*HD + 
            model.phiDt * (mModel.c  + (mModel.e * mModel.lambdaFWH * dFW + mModel.mW)* HW)
    
    return [dHD, dFW, dHW]
end

function solveModel(model::hunterNS)
    for ind in 1:model.n
        model.result[1, ind+1] = model.result[1, ind] + model.dt
        model.result[2:4, ind+1] = numericalScheme(model, model.result[2:4, ind])
    end
end

function computeCharacteristicLC(model::hunterNS)
    mModel = model.mathModel
    @assert mModel.c == 0
    eq = equilibriumHFW(mModel)
    model.result[:,1] = [0, eq[1]+10, eq[2], eq[3]]
    solveModel(model)
    (maxHD, lineMaxHD) = findmax(model.result[2,:])
    (minHD, lineMinHD) = findmin(model.result[2,:])

    timeMaxHD = model.result[1, lineMaxHD]
    timeMinHD = model.result[1, lineMinHD]

    return (maxHD, timeMaxHD, minHD, timeMinHD)
end
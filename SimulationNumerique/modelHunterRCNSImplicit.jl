include("modelHunterRC.jl")

struct hunterNSImplicit <: numericalModel
    variablesNames::Vector{String}
    mathModel::modelHunter

    t0::Float64
    tf::Float64
    dt::Float64
    n::Int64

    phiDt::Float64

    result::Matrix{Float64}

    function hunterNSImplicit(mathModel::modelHunter, numericalParam::Dict{String, Float64}, 
        initialValues::Dict{String, Float64})
        t0, tf, dt = numericalParam["t0"], numericalParam["tf"], numericalParam["dt"]
        n = Int((tf-t0)/dt)

        Q = mathModel.rF
        phiDt = (1- exp(-Q * dt)) / Q
        

        result = Matrix{Float64}(undef, 4, n+1)
        result[:,1] = [t0, initialValues["HD0"], initialValues["FW0"], initialValues["HW0"]]
        
        variablesNames = mathModel.variablesNames
        pushfirst!(variablesNames, "time")

        new(variablesNames, mathModel, t0, tf, dt, n, phiDt,result)
        end
end

function numericalScheme(model::hunterNSImplicit, variables::Vector{Float64})
    mModel = model.mathModel
    HD, FW, HW = variables[1], variables[2], variables[3]
    schemeMatrix = Matrix{Float64}(undef, 3,3)
    schemeMatrix[1,:] = [1 - model.phiDt * (mModel.fD - mModel.muD - mModel.mD),
                    -model.phiDt * mModel.e * mModel.lambdaFWH * HW, - model.phiDt * mModel.mW]
    schemeMatrix[2,:] = [0, 1 - model.phiDt * (mModel.rF * (1 - FW / ((1-mModel.alpha)*mModel.KF)) - mModel.lambdaFWH * HW), 0]
    schemeMatrix[3,:] = [- model.phiDt * mModel.mD, 0, 1 + model.phiDt * mModel.mW]
    
    problem = LinearProblem(schemeMatrix, variables + [mModel.c*model.phiDt, 0.,0.])
    variablesN1 = solve(problem)
    
    return variablesN1
end

function solveModel(model::hunterNSImplicit)
    for ind in 1:model.n
        model.result[1, ind+1] = model.result[1, ind] + model.dt
        model.result[2:4, ind+1] = numericalScheme(model, model.result[2:4, ind])
    end
end
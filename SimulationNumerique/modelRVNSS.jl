include("modelRV.jl")

struct RVNSS <: Model

    mathModel::modelRV

    t0::Float64
    tf::Float64
    dt::Float64
    n::Int64

    phiF::Float64
    phiV::Float64
    phiH::Float64
    phi::Float64

    result::Matrix{Float64}

    function RVNSS(modelParam::Dict{String, Float64}, numericalParam::Dict{String, Float64}, initialValues::Dict{String, Float64})
        mathModel = modelRV(modelParam)

        t0, tf, dt = numericalParam["t0"], numericalParam["tf"], numericalParam["dt"]
        n = Int((tf-t0)/dt)
        
        result = Matrix{Float64}(undef, 4, n+1)
        result[:,1] = [t0, initialValues["F0"], initialValues["V0"], initialValues["H0"]]
        
        deltaF = mathModel.rF - mathModel.muF - mathModel.omega*mathModel.f
        deltaV = mathModel.rV - mathModel.muV

        phiF = (exp(deltaF * dt) - 1)/deltaF
        phiV = (exp(deltaV * dt) - 1)/deltaV
        phiH = dt
        phi = max(phiF, phiV)

        new(mathModel, t0, tf, dt, n, phiF, phiV, phiH, phi, result)
    end
end

function numericalScheme(model::RVNSS, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    rVH = model.mathModel.rV * H/(H+model.mathModel.H0)

    deltaF = model.mathModel.rF - model.mathModel.muF - model.mathModel.omega*model.mathModel.f
    deltaVH = rVH - model.mathModel.muV

    Fdt = F*(1 + model.phi*deltaF) / (1 + model.mathModel.rF/model.mathModel.KF * model.phi * F + model.mathModel.lambdaFH * model.phi * H)
    Vdt = V*(1 + model.phi*deltaVH) / (1 + rVH/model.mathModel.KV * model.phi * V + model.mathModel.alpha * model.phi * F + model.mathModel.lambdaVH * model.phi * H)
    Hdt = H*(1 + model.phi * model.mathModel.e * (model.mathModel.lambdaFH * Fdt + model.mathModel.lambdaVH * Vdt)) / (1 + model.mathModel.muH * model.phi * H)

    return [Fdt, Vdt, Hdt]
end

function solveModel(model::RVNSS)
    for ind in 1:model.n
        model.result[1, ind+1] = model.result[1, ind] + model.dt
        model.result[2:4, ind+1] = numericalScheme(model, model.result[2:4, ind])
    end
end
include("modelEcoServiceFV.jl")

struct ecoServiceFVRK4 <: numericalModel

    mathModel::modelEcoServiceFV

    t0::Float64
    tf::Float64
    dt::Float64
    n::Int64

    result::Matrix{Float64}

    function ecoServiceFVRK4(modelParam::Dict{String, Float64}, numericalParam::Dict{String, Float64}, 
                            initialValues::Dict{String, Float64})
        mathModel = modelEcoServiceFV(modelParam)

        t0, tf, dt = numericalParam["t0"], numericalParam["tf"], numericalParam["dt"]
        n = Int((tf-t0)/dt)
        
        result = Matrix{Float64}(undef, 4, n+1)
        if initialValues["H0"] != (mathModel.a * initialValues["F0"] + mathModel.b * initialValues["V0"] + mathModel.c)
            println("In this model, H must be equal at \$H = a F + bV + c\$, even for initial values")
            println("Initial values of H has been changed to ensure this requirement")
            initialValues["H0"] = mathModel.a * initialValues["F0"] + mathModel.b * initialValues["V0"] + mathModel.c
        end

        result[:,1] = [t0, initialValues["F0"], initialValues["V0"], initialValues["H0"]]

        new(mathModel, t0, tf, dt, n, result)
    end

    function ecoServiceFVRK4(mathModel::modelEcoServiceFV, numericalParam::Dict{String, Float64}, 
        initialValues::Dict{String, Float64})

        t0, tf, dt = numericalParam["t0"], numericalParam["tf"], numericalParam["dt"]
        n = Int((tf-t0)/dt)

        result = Matrix{Float64}(undef, 4, n+1)
        if initialValues["H0"] != (mathModel.a * initialValues["F0"] + mathModel.b * initialValues["V0"] + mathModel.c)
            println("In this model, H must be equal at \$H = a F + bV + c\$, even for initial values")
            println("Initial values of H has been changed to ensure this requirement")
            initialValues["H0"] = mathModel.a * initialValues["F0"] + mathModel.b * initialValues["V0"] + mathModel.c
        end
        result[:,1] = [t0, initialValues["F0"], initialValues["V0"], initialValues["H0"]]

        new(mathModel, t0, tf, dt, n, result)
        end
end

function numericalScheme(model::ecoServiceFVRK4, variables::Vector{Float64})
    k1 = equationModel(model.mathModel, variables)
    k2 = equationModel(model.mathModel, variables + model.dt/2 * k1)
    k3 = equationModel(model.mathModel, variables + model.dt/2 * k2)
    k4 = equationModel(model.mathModel, variables + model.dt * k3)

    return variables + model.dt/6*(k1 + 2*k2 + 2*k3 + k4)
end

function solveModel(model::ecoServiceFVRK4)
    for ind in 1:model.n
        model.result[1, ind+1] = model.result[1, ind] + model.dt
        model.result[2:3, ind+1] = numericalScheme(model, model.result[2:3, ind])
        model.result[4, ind+1] = (model.mathModel.a * model.result[2,ind+1] + 
                                    model.mathModel.b * model.result[3,ind+1] + model.mathModel.c)
    end
end
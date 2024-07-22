struct modelMigration <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    LW::Float64
    f::Float64
    lambdaFHW::Float64

    rV::Float64
    KV::Float64
    alpha::Float64
    lambdaVHW::Float64

    rH::Float64
    aW::Float64
    bW::Float64
    c::Float64
    beta::Float64
        
    function modelMigration(modelParam::Dict{String, Float64})
        rF, f, LW = modelParam["rF"], modelParam["f"], modelParam["LW"]
        rV, KV, alpha = modelParam["rV"], modelParam["KV"], modelParam["alpha"]
        
        variablesNames = ["F", "V"]

        new(variablesNames, rF, LW, f, rV, KV, alpha)
    end
end
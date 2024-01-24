
struct modelRV <: Model
    rF::Float64
    KF::Float64
    omega::Float64
    f::Float64
    muF::Float64
    lambdaFH::Float64

    rV::Float64
    H0::Float64
    KV::Float64
    alpha::Float64
    muV::Float64
    lambdaVH::Float64

    e::Float64
    muH::Float64

    
    function modelRV(modelParam::Dict{String, Float64})
        rV, H0, KV, alpha, muV, lambdaVH = modelParam["rV"], modelParam["H0"], modelParam["KV"], modelParam["alpha"], modelParam["muV"], modelParam["lambdaVH"]
        rF, KF, omega, f, muF, lambdaFH = modelParam["rF"], modelParam["KF"], modelParam["omega"], modelParam["f"], modelParam["muF"], modelParam["lambdaFH"]
        e, muH = modelParam["e"], modelParam["muH"]
       
        new(rF, KF, omega, f, muF, lambdaFH, rV, H0, KV, alpha, muV, lambdaVH, e, muH)
    end
end

function equationModel(model::modelRV, variables::Vector{Float64})
    F, V, H = variables[1], variables[2], variables[3]
    
    dF = model.rF * (1 - F / model.KF) * F - model.omega * model.f * F - model.muF * F - model.lambdaFH * F * H
    dV = model.rV * H / (H + model.H0) * (1 - V / model.KV) * V - model.alpha * V * F - model.muV * V - model.lambdaVH * V * H
    dH = model.e * (model.lambdaFH * F + model.lambdaVH * V) * H - model.muH * H^2
    
    return [dF, dV, dH]
end

# function stabilityTresholds(model::modelRV)
#     R0F = model.rF / (model.muF + model.omega*model.f)
    
#     palha = model.alpha*model.muH/(model.e * model.lambdaFH * model.lambdaVH)
#     pf = model.muH * model.rF/ (model.e*model.lambdaFH^2 *model.KF)
#     pv = model.muH * model.rV/ (model.e*model.lambdaVH^2 *model.KV)
    
#     TF = (model.lambdaVH / model.muV) * (model.omega*model.f + model.muF) / model.lambdaFH *
#         (R0F - 1) / (model.rV/model.muV - 1) * (1 + palha) / (1 + pf)
#         + model.H0 * (1+palha) * model.lambdaVH / model.muV / (model.rV - model.muV)
#         + model.H0 * model.lambdaFH / (model.omega*model.f + model.muF) * (1 + pf) / (R0F - 1)

#     TV11 = model.muV * model.lambdaFH / (model.lambdaVH * (model.omega * model.f + model.muF)) / (1 + pv)
#             * ((model.rV / model.muV - 1 - model.lambdaVH/model.muV * model.H0) - 
#                 sqrt((model.rV / model.muV - 1 - model.lambdaVH/model.muV * model.H0)^2 - 4 *(1 + pv)*model.lambdaVH/model.rV * model.H0)
#             )

#     TV12 = model.muV * model.lambdaFH / (model.lambdaVH * (model.omega * model.f + model.muF)) / (1 + pv)
#             * ((model.rV / model.muV - 1 - model.lambdaVH/model.muV * model.H0) + 
#                 sqrt((model.rV / model.muV - 1 - model.lambdaVH/model.muV * model.H0)^2 - 4 *(1 + pv)*model.lambdaVH/model.rV * model.H0)
#             )
    
# end

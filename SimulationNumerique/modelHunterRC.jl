struct modelHunter <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    KF::Float64
    alpha::Float64
    beta::Float64
    lambdaFWH::Float64

    e::Float64
    I::Float64
    muD::Float64
    fD::Float64
    mD::Float64
    mW::Float64

    m::Float64
        
    function modelHunter(modelParam::Dict{String, Any})

        rF, lambdaFWH, KF = modelParam["rF"], modelParam["lambdaFWH"], modelParam["KF"]
        alpha, beta = modelParam["alpha"], modelParam["beta"]
        e, I = modelParam["e"], modelParam["I"]
        muD, fD, mD, mW = modelParam["muD"], modelParam["fD"], modelParam["mD"], modelParam["mW"]
        variablesNames = ["H_D", "F_W", "H_W"]
        m  = mD / mW

        @assert(alpha < 1)
        bound = 4 * (muD - fD) / (m * e * rF * KF * (1-alpha)^2)
        @assert (bound > beta) "bound for beta is $bound"

        new(variablesNames, rF, KF, alpha, beta, lambdaFWH, e, I, muD, fD, mD, mW, m)
    end
end

function equationModel(model::modelHunter, variables::Vector{Float64})
    """
    Return the right hand side of the model
    """
    HD, FW, HW = variables[1], variables[2], variables[3]
    dHD = model.I + model.e * model.lambdaFWH * HW * FW + (model.fD- model.muD) * HD - model.mD * HD + model.mW * HW
    dFW = ((1 - model.alpha) * (1 + model.beta * HW) * model.rF * (1 - FW / (model.KF * (1 - model.alpha))) * FW -
            model.lambdaFWH * FW * HW)
    dHW = model.mD * HD - model.mW * HW
    return [dHD, dFW, dHW]
end

function computeBetaMax(model::modelHunter ; alpha = -1.0)
    if alpha < 0
        alpha = model.alpha
    end
    betaMax = 4 * (model.muD - model.fD) / (model.m * model.e * model.rF * model.KF * (1-alpha)^2)
    return betaMax
end

function thresholdHD(model::modelHunter)
    """
    When I > 0 ;
    If < 1, EE^H is GAS ; if > 1 EE^HFW exists
    """
    rslt = (((model.muD - model.fD) / (model.m * model.I) + model.beta) * 
                model.rF * (1 - model.alpha) / model.lambdaFWH)
    return rslt
end

function thresholdFW(model::modelHunter)
    """
    When I = 0 ;
    If < 1, EE^F is GAS ; if > 1 EE^HFW exists
    """
    return (model.lambdaFWH * model.m * model.e * model.KF * (1 - model.alpha)) / (model.muD - model.fD)
end

function printAllCharacteristicValues(model::modelHunter)
    betaMinMax = computeBetaMax(model ; alpha = 0)
    betaMax = computeBetaMax(model)

    println("Beta max pour alpha = 0 : ", betaMinMax)
    println("Beta max pour valeur de alpha donnée : ", betaMax)

    if model.I == 0
        lambdaMin = computeLambdaMinI0(model)
        println("Cas I = 0 ; lambdaMin = ", lambdaMin)
        if model.beta == 0
            lambdaMax = computeLambdaMaxI0(model)
            println("Cas I = beta = 0 ; lambdaMax = ", lambdaMax)
        end
    end

    if model.I > 0
        lambdaMaxMax = computeLambdaMax(model; alpha = 0)
        lambdaMax = computeLambdaMax(model)
        println("Cas I > 0 ; lambdaMax pour alpha = 0 : ", lambdaMaxMax)
        println("Cas I > 0 ; lambdaMax pour alpha donné : ", lambdaMax)
    end
end


function equilibriumFW(model::modelHunter)
    """
    When I == 0 ;
    return eq EE^F
    """
    return [0, model.KF * (1 - model.alpha), 0]
end

function equilibriumH(model::modelHunter)
    """
    When I > 0 ;
    return eq EE^H
    """
    HD = model.I / (model.muD - model.fD)
    return [HD, 0, model.m * HD]
end

function equilibriumHFW(model::modelHunter)
    """
    For I >= 0
    Compute the values of EE^HFW
    Rise a warning if negatif (threshold not respected)
    """
    
    if model.I == 0
        Feq = (model.muD - model.fD) / (model.m * model.e * model.lambdaFWH)
        thresoldExistence = thresholdFW(model)
    else
        A = model.e * model.rF / model.KF
        B = -(model.e * (1 - model.alpha) * model.rF + (model.muD - model.fD) * model.rF / (model.lambdaFWH * model.m * model.KF)
                + model.I * model.beta * model.rF / (model.lambdaFWH * model.KF))
        C = ((model.muD - model.fD) * (1 - model.alpha) * model.rF / (model.lambdaFWH * model.m) - 
                model.I * (1 - (1 - model.alpha) * model.beta * model.rF / model.lambdaFWH))

        Delta = B^2 - 4 *A *C # Delta is always positive

        thresoldExistence = thresholdHD(model)

        Feq = -B/(2*A) - sqrt(Delta) / (2*A)        
    end

    if !(thresoldExistence > 1)
        @warn "Equilibrium EE^HFW computed while it does not exists ; negative values return "
    end
    Heq = (1 - model.alpha) * model.rF * (1 - Feq / (model.KF * (1 - model.alpha)))
    Heq = Heq / (model.m * (model.lambdaFWH - model.beta * (1 - model.alpha) * model.rF + 
                    model.beta * model.rF * Feq / model.KF))

    return [Heq, Feq, model.m * Heq]
    
end

function computeMaxVal(model::modelHunter)
    """
    Compute the bound for the invariant region Omega 
        {H_D + H_W + e F_W < Smax, F_W < Fmax, H_W < m_D/m_D+m_W Smax}
    """
    Fmax = model.KF * (1 - model.alpha)
    Smax = (1 + model.m) * (model.I + 
        (model.muD - model.fD + model.rF * (1-model.alpha) / 4) * model.e * model.KF * (1-model.alpha))
    Smax = Smax / ((model.muD - model.fD) / model.m - model.e * model.rF * (1-model.alpha)^2 * model.KF * model.beta / 4)
    Hwmax = model.mD / (model.mD + model.mW) * Smax
    
    return (Smax, Fmax, Hwmax)
end

function computeLambdaMax(model::modelHunter ; alpha = -1.0)
    """
    When I > 0, compute lambda^max such that lambda < lambda^max iff EE^HFW exists
    """
    @assert model.I > 0
    if alpha < 0
        alpha = model.alpha
    end
    lambdaMax = (1 - alpha) * model.rF * ( model.beta + 
        (model.muD - model.fD) / (model.m * model.I))
    return lambdaMax
end

function computeLambdaMinI0(model::modelHunter ; alpha = -1.0, e = -1.0, mW = -1.0)
    """
    When I = 0, compute lambda^min such that lambda > lambda^min iff EE^HFW exists
    """
    @assert model.I == 0
    if alpha < 0
        alpha = model.alpha
    end
    if e < 0
        e = model.e
    end
    if mW < 0
        mW = model.mW
    end

    m = model.mD / mW
    lambdaMin = (model.muD - model.fD) / (m * e * model.KF * (1-alpha))
    return lambdaMin
end

function computeLambdaMaxI0(model::modelHunter ; alpha = -1.0, e = -1.0, mW = -1.0)
    """
    When I = \beta = 0, compute lambda^max such that lambda < lambda^max iff EE^HFW is AS
    """    
    @assert (model.I == 0) && (model.beta == 0)  
    if alpha < 0
        alpha = model.alpha
    end

    if e < 0
        e = model.e
    end

    if mW < 0
        mW = model.mW
    end

    m = model.mD / mW

    a0 = (-(model.muD - model.fD + model.mD + mW) * 
            model.rF * (model.muD - model.fD) / (e * m * model.KF)
            )
    a1 = -(mW * (model.muD - model.fD) + 
            (model.muD - model.fD + model.mD + mW)^2)
    a2 = (1 - alpha) * model.KF * e * model.mD

    PDeltaStab = Polynomial([a0, a1, a2])
    rootsPDeltaStab = real(PolynomialRoots.roots(coeffs(PDeltaStab)))
    lambdaMaxStab = maximum(rootsPDeltaStab)

    return lambdaMaxStab
end

function computeDeltaStab(model::modelHunter ; lambdaFWH = -1., alpha = -1.0, rsltFeq = true)
    """
    For I >= 0
    Compute the value of DeltaStab = a2a1 - a0
    If positive, eq EE^HFW is GAS
    Return DeltaStab, and eventually Feq
    Error when EE^HFW does not exist
    """
    
    if lambdaFWH < 0
        lambdaFWH = model.lambdaFWH
    end
    if alpha < 0
        alpha = model.alpha
    end

    ## Compute the values at equilibrium with eventually
    ## new parameters values 

    A = model.e * model.rF / model.KF
    B = -(model.e * (1 - alpha) * model.rF + (model.muD - model.fD) * model.rF / (lambdaFWH * model.m * model.KF)
            + model.I * model.beta * model.rF / (lambdaFWH * model.KF))
    C = ((model.muD - model.fD) * (1 - alpha) * model.rF / (lambdaFWH * model.m) - 
            model.I * (1 - (1 - alpha) * model.beta * model.rF / lambdaFWH))

    DeltaF = B^2 - 4 *A *C # Delta is always positive
    ## Check that equilibrium is well defined
    # @assert(C > 0)
    
    Feq = -B / (2*A) - sqrt(DeltaF) / (2*A)
    HWeq = (1 - alpha) * model.rF * (1 - Feq / ((1 - alpha) * model.KF))
    HWeq = HWeq / (lambdaFWH - model.beta * (1 - alpha) * model.rF + 
                    model.beta * model.rF * Feq / model.KF)
    if C <= 0
        println("C : ", C)
        println("lambda : ", lambdaFWH)
        println("Alpha : ", alpha)
        println("Feq : ", Feq)
    end

    ## Compute a2, a1 and a0             
    a2 = -(model.fD - model.muD - model.mD - model.mW - 
        (1 + model.beta * HWeq) * model.rF * Feq / model.KF)

    a0 = model.e * lambdaFWH * model.mD * model.rF * (
        (model.muD - model.fD) / (model.e * model.m * model.KF * lambdaFWH) - 2 * Feq / model.KF + 
        1 - alpha + model.beta * HWeq * 
                ((model.muD - model.fD) / (model.e * model.m * model.KF * lambdaFWH) - Feq / model.KF)
        ) * Feq
    
    a1 = ((model.muD - model.fD + model.mD + model.mW) * model.rF * (1 + model.beta * HWeq) *
            Feq / model.KF + model.mD * model.e * lambdaFWH * 
                ((model.muD - model.fD) / (model.m * lambdaFWH * model.e) - Feq) )

    @assert (a1 > 0 && a2 > 0 && a0 > 0)

    if rsltFeq
        return a2 * a1 - a0, Feq
    else
        return a2 * a1 - a0
    end
end

function LambdaLC(model::modelHunter)

    alpha = model.alpha
    function Delta(lambdaFWH)
        A = model.e * model.rF / model.KF
        B = -(model.e * (1 - alpha) * model.rF + (model.muD - model.fD) * model.rF / (lambdaFWH * model.m * model.KF)
                + model.I * model.beta * model.rF / (lambdaFWH * model.KF))
        C = ((model.muD - model.fD) * (1 - alpha) * model.rF / (lambdaFWH * model.m) - 
                model.I * (1 - (1 - alpha) * model.beta * model.rF / lambdaFWH))
    
        DeltaF = B^2 - 4 *A *C # Delta is always positive
        ## Check that equilibrium is well defined
                
        Feq = -B / (2*A) - sqrt(DeltaF) / (2*A)
        HWeq = (1 - alpha) * model.rF * (1 - Feq / ((1 - alpha) * model.KF))
        HWeq = HWeq / (lambdaFWH - model.beta * (1 - alpha) * model.rF + 
                        model.beta * model.rF * Feq / model.KF)
    
        ## Compute a2, a1 and a0             
        a2 = -(model.fD - model.muD - model.mD - model.mW - 
            (1 + model.beta * HWeq) * model.rF * Feq / model.KF)
    
        a0 = model.e * lambdaFWH * model.mD * model.rF * (
            (model.muD - model.fD) / (model.e * model.m * model.KF * lambdaFWH) - 2 * Feq / model.KF + 
            1 - alpha + model.beta * HWeq * 
                    ((model.muD - model.fD) / (model.e * model.m * model.KF * lambdaFWH) - Feq / model.KF)
            ) * Feq
        
        a1 = ((model.muD - model.fD + model.mD + model.mW) * model.rF * (1 + model.beta * HWeq) *
                Feq / model.KF + model.mD * model.e * lambdaFWH * 
                    ((model.muD - model.fD) / (model.m * lambdaFWH * model.e) - Feq) )
        return a2 * a1 - a0
    end

    if model.I == 0
        lambdaMin = computeLambdaMinI0(model)
        zerosDelta = find_zeros(Delta, lambdaMin+0.01, 100)
    end

    if model.I > 0
        lambdaMax = computeLambdaMax(model)
        zerosDelta = find_zeros(Delta, 0.0001, lambdaMax)
    end

    for lambda in zerosDelta
        println("Lambda = ", lambda, "; Delta = ", Delta(lambda))
    end
    return zerosDelta
end

function printLongTermDynamic(model::modelHunter)
    if model.I == 0
        lambdaMin = computeLambdaMinI0(model)
        if model.lambdaFWH < lambdaMin
            eq = equilibriumFW(model)
            return "I = 0 ; EE^FW = " * string(eq) * " is GAS"
        elseif model.lambdaFWH > lambdaMin
            eq = equilibriumHFW(model)
            DeltaStab = computeDeltaStab(model, rsltFeq = false)
            if DeltaStab > 0
                return "I = 0 ; EE^HFW = " * string(eq) * " is GAS"
            elseif DeltaStab < 0
                return "I = 0 ; Limit Cycle around EE^HFW = " * string(eq)
            else
                return "I = 0 ; DeltaStab = 0"
            end
        else
            return "I = 0 ; lambda = lambdaMin"
        end
    else
        lambdaMax = computeLambdaMax(model)
        if model.lambdaFWH > lambdaMax
            eq = equilibriumH(model)
            return "I > 0 ; EE^H = " * string(eq) * " is GAS"
        elseif model.lambdaFWH < lambdaMax
            eq = equilibriumHFW(model)
            DeltaStab = computeDeltaStab(model, rsltFeq = false)
            println(DeltaStab)
            if DeltaStab > 0
                return "I > 0 ; EE^HFW = " * string(eq) * " is GAS"
            elseif DeltaStab < 0
                return "I > 0 ; Limit Cycle around EE^HFW = " * string(eq)
            else
                return "I > 0 ; DeltaStab = 0"
            end
        else
            return "I > 0 ; lambda = lambdaMax"
        end
    end
end

function longTermDynamic(model::modelHunter)
    if model.I == 0
        lambdaMin = computeLambdaMinI0(model)
        if model.lambdaFWH < lambdaMin
            eq = equilibriumFW(model)
            return push!(eq, 0.)
        else
            eq = equilibriumHFW(model)
            DeltaStab = computeDeltaStab(model, rsltFeq = false)
            if DeltaStab > 0
                return push!(eq, 0.)
            else
                return push!(eq, 1.)
            end
        end
    else
        lambdaMax = computeLambdaMax(model)
        if model.lambdaFWH >= lambdaMax
            eq = equilibriumH(model)
            return push!(eq, 0)
        else
            eq = equilibriumHFW(model)
            DeltaStab = computeDeltaStab(model, rsltFeq = false)
            if DeltaStab > 0
                return push!(eq, 0.)
            else DeltaStab < 0
                return push!(eq, 1.)
            end
        end
    end
end


function bifurcationEquilibriumValues(model::modelHunter, paramName::String,
                listParam::Vector{Float64})

            bifurcationTable = Matrix{Any}(undef, length(listParam)+1, 5)
            bifurcationTable[1,:] = [paramName, "H_D", "F_W", "H_W", "LC"]
            bifurcationTable[2:end,1] = listParam

            localParam = Dict(string(name) => getfield(model, name) for name in fieldnames(modelHunter))
            delete!(localParam, "variablesNames")
            for (ind, paramVal) in enumerate(listParam)
                localParam[paramName] = paramVal
                localModel = modelHunter(localParam)
                bifurcationTable[ind+1, 2:end] .= longTermDynamic(localModel)
            end
            return bifurcationTable
    
end

function computeDeltaStabMatrix(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})
    """
    For I >= 0
    Compute the value of DeltaStab = a2a1 - a0 for values of alpha and lambda
    Return 
        - matrix of DeltaStab values, with alpha in first line, lambda in first column
        - matrix of Feq values, with alpha in first line, lambda in first column
    Values of alpha and lambda must be such EE^HFW exists
    """

    deltaStabMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    deltaStabMatrix[2:end,1] = model.KF * (1 .- listAlpha)

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
            deltaStabMatrix[ind_row, ind_col] = rsltDelta
            FeqMatrix[ind_row, ind_col] = Feq
        end
    end

    deltaStabMatrix[1,1] = "alpha" * "\\" * "lambdaFWH"
    FeqMatrix[1,1] = "alpha" * "\\" * "lambdaFWH"

    return deltaStabMatrix, FeqMatrix


end

function computeBifurcationDiagram(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})
    """
    For I >= 0
    Return the bifurcation matrix in function of alpha and lambda
        First line : alpha ; first column : lambda
    """
    if model.I == 0
        bifurcationMatrix = computeBifurcationDiagramI0(model, listAlpha, listLambdaFWH)
    else
        bifurcationMatrix =computeBifurcationDiagramI(model, listAlpha, listLambdaFWH)
    end
    return bifurcationMatrix
end

function computeBifurcationDiagramI(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})
    
    """
    For I > 0
    Return the bifurcation matrix in function of alpha and lambda
        First line : alpha ; first column : lambda
    """

    @assert model.I > 0

    bifurcationMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    bifurcationMatrix[2:end,1] = listAlpha

    for ind_row in 2:length(listAlpha)+1
        alpha = listAlpha[ind_row-1]
        lambdaMax = computeLambdaMax(model; alpha = alpha)
        
        for ind_col in 2:length(listLambdaFWH)+1
            bifurcationMatrix[1, ind_col] = listLambdaFWH[ind_col-1]
            lambdaFWH = listLambdaFWH[ind_col-1]

            if lambdaFWH >= lambdaMax
                bifurcationMatrix[ind_row, ind_col] = "H" ## Specific for I>0
            else
                Delta, _ = computeDeltaStab(model; lambdaFWH=lambdaFWH, alpha=alpha)

                if Delta > 0
                    bifurcationMatrix[ind_row, ind_col] = "HFW"
                else
                    bifurcationMatrix[ind_row, ind_col] = "LC"
                end
            end
        end
    end

    bifurcationMatrix[1,1] = "alpha" * "\\" * "lambdaFWH"

    return bifurcationMatrix
end

function computeBifurcationDiagramI0(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})

    """
    For I == 0
    Return the bifurcation matrix in function of alpha and lambda
        First line : alpha ; first column : lambda
    """

    @assert model.I == 0

    bifurcationDiagram = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    bifurcationDiagram[2:end,1] = listAlpha

    for ind_row in 2:length(listAlpha)+1
        alpha = listAlpha[ind_row-1]
        lambdaMin = computeLambdaMinI0(model; alpha)
        lambdaMaxStab = computeLambdaMaxI0(model; alpha)

        for ind_col in 2:length(listLambdaFWH)+1
            bifurcationDiagram[1, ind_col] = listLambdaFWH[ind_col-1]
            lambdaFWH = listLambdaFWH[ind_col-1]

            if lambdaFWH <= lambdaMin ## Specific for I = 0
                rslt = "F"
            else
                deltaStab,_ = computeDeltaStab(model; alpha=alpha, lambdaFWH=lambdaFWH)
                if deltaStab <= 0
                    rslt = "LC"
                else
                    rslt = "HFW"
                end
            end
            bifurcationDiagram[ind_row, ind_col] = rslt
        end
    end

    bifurcationDiagram[1,1] = "alpha" * "\\" * "lambdaFWH"

    return bifurcationDiagram
end

function computeLambdaMatrixAlphaE(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listE::Vector{Float64})
    """
    For I = 0
    Compute the values of lambda^max and lambda^min as function of parameters alpha and e
    Return two matrixs, which contains the values
    """

    lambdaMaxMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listE)+1)
    lambdaMaxMatrix[2:end,1] = model.KF * (1 .- listAlpha)

    lambdaMinMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listE)+1)
    lambdaMinMatrix[2:end,1] = model.KF * (1 .- listAlpha)

    for ind_col in 2:length(listE)+1
        lambdaMaxMatrix[1, ind_col] = listE[ind_col-1]
        lambdaMinMatrix[1, ind_col] = listE[ind_col-1]
        
        e = listE[ind_col-1]
        for ind_row in 2:length(listAlpha)+1
            alpha = listAlpha[ind_row-1]
            lambdaStab = computeLambdaMaxI0(model; e=e, alpha=alpha)
            lambdaMaxMatrix[ind_row, ind_col] = lambdaStab

            lambdaExistence = computeLambdaMinI0(model; e=e, alpha=alpha)
            lambdaMinMatrix[ind_row, ind_col] = lambdaExistence
        end
    end

    lambdaMaxMatrix[1,1] = "alpha" * "\\" * "e"
    lambdaMinMatrix[1,1] = "alpha" * "\\" * "e"

    return lambdaMinMatrix, lambdaMaxMatrix
end

function computeLambdaMatrixAlphaMW(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listMW::Vector{Float64})
    """
    For I = 0
    Compute the values of lambda^max and lambda^min as function of parameters alpha and mW
    Return two matrixs, which contains the values
    """

    lambdaMaxMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listMW)+1)
    lambdaMaxMatrix[2:end,1] = model.KF * (1 .- listAlpha)

    lambdaMinMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listMW)+1)
    lambdaMinMatrix[2:end,1] = model.KF * (1 .- listAlpha)

    for ind_col in 2:length(listMW)+1
        lambdaMaxMatrix[1, ind_col] = model.mD / listMW[ind_col-1]
        lambdaMinMatrix[1, ind_col] = model.mD / listMW[ind_col-1]
        
        mW = listMW[ind_col-1]
        for ind_row in 2:length(listAlpha)+1
            alpha = listAlpha[ind_row-1]
            lambdaStab = computeLambdaMaxI0(model; mW=mW, alpha=alpha)
            lambdaMaxMatrix[ind_row, ind_col] = lambdaStab

            lambdaExistence = computeLambdaMinI0(model; mW=mW, alpha=alpha)
            lambdaMinMatrix[ind_row, ind_col] = lambdaExistence
        end
    end

    lambdaMaxMatrix[1,1] = "alpha" * "\\" * "m"
    lambdaMinMatrix[1,1] = "alpha" * "\\" * "m"

    return lambdaMinMatrix, lambdaMaxMatrix
end


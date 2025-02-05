struct modelHunter <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    KF::Float64
    alpha::Float64
    lambdaFWH::Float64

    e::Float64
    c::Float64
    muD::Float64
    fD::Float64
    mD::Float64
    mW::Float64

    m::Float64
        
    function modelHunter(modelParam::Dict{String, Float64})
        rF, lambdaFWH, KF, alpha = modelParam["rF"], modelParam["lambdaFWH"], modelParam["KF"], modelParam["alpha"]
        e, c = modelParam["e"], modelParam["c"]
        muD, fD, mD, mW = modelParam["muD"], modelParam["fD"], modelParam["mD"], modelParam["mW"]
        alpha = modelParam["alpha"]
        variablesNames = ["H_D", "F_W", "H_W"]
        m  = mD / mW
        new(variablesNames, rF, KF, alpha, lambdaFWH, e, c, muD, fD, mD, mW, m)
    end
end

function equationModel(model::modelHunter, variables::Vector{Float64})
    """
    Return the right hand side of the model
    """
    HD, FW, HW = variables[1], variables[2], variables[3]
    dHD = model.c + model.e * model.lambdaFWH * HW * FW + (model.fD- model.muD) * HD - model.mD * HD + model.mW * HW
    dFW = model.rF * (1 - FW / (model.KF * (1 - model.alpha))) * FW - model.lambdaFWH * FW * HW
    dHW = model.mD * HD - model.mW * HW
    return [dHD, dFW, dHW]
end

function thresholdHD(model::modelHunter)
    """
    When cI > 0 ;
    If < 1, EE^H is GAS ; if > 1 EE^HFW exists
    """
    @assert model.c > 0
    return (model.muD - model.fD) * model.rF / (model.lambdaFWH * model.m * model.c)
end

function thresholdFW(model::modelHunter)
    """
    When cI = 0 ;
    If > 1, EE^F is GAS ; if < 1 EE^HFW exists
    """
    @assert model.c == 0
    return (model.muD - model.fD) / (model.lambdaFWH * model.m * model.e * model.KF * (1 - model.alpha))
end

function equilibriumFW(model::modelHunter)
    """
    When cI == 0 ;
    return eq EE^F
    """
    return (0, model.KF * (1 - model.alpha), 0)
end

function equilibriumH(model::modelHunter)
    """
    When cI > 0 ;
    return eq EE^H
    """
    HD = model.c / (model.muD - model.fD)
    return (HD, 0, model.m * HD)
end

function equilibriumHFW(model::modelHunter)
    """
    For cI >= 0
    Compute the values of EE^HFW
    Return (0,0,0) if does not exist (+warning)
    """
    
    if model.c == 0
        Feq = (model.muD - model.fD) / (model.m * model.e * model.lambdaFWH)
        if Feq < model.KF * (1 - model.alpha)
            Heq = model.rF / (model.m * model.lambdaFWH) * (1 - Feq / (model.KF * (1 - model.alpha)))
            return (Heq, Feq, model.m * Heq)
        else
            @warn "Equilibrium EE^HFW computed while it does not exists; return 0"
            return (0,0,0)
        end

    else
        A = model.e * model.rF / (model.KF * (1 - model.alpha))
        B = -(model.e * model.rF + (model.muD - model.fD) * model.rF / (model.lambdaFWH * model.m * model.KF * (1-model.alpha)))
        C = (model.muD - model.fD) * model.rF / (model.lambdaFWH * model.m) - model.c

        Delta = B^2 - 4 *A *C
        if Delta >= 0
            Feq = -B/(2*A) - sqrt(Delta) / (2*A)
            Heq = model.rF / (model.lambdaFWH * model.m) * (1 - Feq / (model.KF * (1 - model.alpha)))
            println(-B/(2*A) + sqrt(Delta) / (2*A))
            return (Heq, Feq, model.m * Heq)
        else
            @warn "Equilibrium EE^HFW computed while it does not exists; return 0"
            return (0,0,0)
        end
    end
end

function computeValMax(model::modelHunter)
    """
    Compute the bound for the invariant region Omega 
        {H_D + H_W + e F_W < Smax, F_W < Fmax, H_W < m_D/m_D+m_W Smax}
    """
    Fmax = model.KF * (1 - model.alpha)
    Smax = (1 + model.m) / (model.muD - model.fD) *
         (model.c + model.e * model.KF * (1 - model.alpha)*(model.rF + model.muD - model.fD))
    Hwmax = model.mD / (model.mD + model.mW) * Smax
    
    return (Smax, Fmax, Hwmax)
end

function computeLambdaMax(model::modelHunter)
    """
    When cI > 0, compute lambda^max such that lambda < lambda^max iff EE^HFW exists
    """
    @assert model.c > 0
    lambdaMax = (model.muD - model.fD) * model.rF / (model.m * model.c)
    return lambdaMax
end

function computeLambdaMincI0(model::modelHunter ; alpha = -1.0, e = -1.0, mW = -1.0)
    """
    When cI = 0, compute lambda^min such that lambda > lambda^min iff EE^HFW exists
    """
    @assert model.c == 0
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

function computeLambdaMaxcI0(model::modelHunter ; alpha = -1.0, e = -1.0, mW = -1.0)
    """
    When cI = 0, compute lambda^max such that lambda < lambda^max iff EE^HFW is AS
    """    
    @assert model.c == 0
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
            model.rF * (model.muD - model.fD) / (e * m)
            )
    a1 = -(mW * (model.muD - model.fD) + 
            (model.muD- model.fD + model.mD + mW)^2)
    a2 = e * model.mD

    PDeltaStab = Polynomial([a0, a1, a2])
    rootsPDeltaStab = real(PolynomialRoots.roots(coeffs(PDeltaStab)))
    root = maximum(rootsPDeltaStab)

    lambdaMaxStab = root / (model.KF * (1 - alpha))
    return lambdaMaxStab
end

function computeDeltaStab(model::modelHunter ; lambdaFWH = -1., alpha = -1.0, rsltFeq = true)
    """
    For cI >= 0
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

    KFalpha = model.KF * (1-alpha)
    A = model.e * model.rF / KFalpha
    B = -(model.e * model.rF + (model.muD - model.fD) * model.rF / (lambdaFWH * model.m * KFalpha))
    C = (model.muD - model.fD) * model.rF / (lambdaFWH * model.m) - model.c
    
    DeltaF = B^2 - 4 *A *C


    @assert(DeltaF >= 0 && C > 0)
    Feq = -B / (2*A) - sqrt(DeltaF) / (2*A)

    Tr = model.fD - model.muD - model.mD - model.mW - model.rF * Feq / KFalpha
    det = -model.mD * lambdaFWH * Feq * sqrt(DeltaF)
    a1 = ((model.muD - model.fD + model.mD + model.mW)*model.rF * Feq / KFalpha + 
            model.mD * model.e * lambdaFWH * ((model.muD - model.fD) / (model.m * lambdaFWH * model.e) - Feq) )

    if !(DeltaF > 0 && C > 0)
        println(alpha)
        println("C :", C, "Delta : ", DeltaF)
        println("det", det, "a1, ", a1, "Tr", Tr)
    end

    @assert (a1 > 0 && det < 0 && Tr < 0)

    if rsltFeq
        return -Tr * a1 + det, Feq
    else
        return -Tr * a1 + det
    end
end

function longTermDynamic(model::modelHunter)
    if model.c == 0
        lambdaMin = computeLambdaMincI0(model)
        lambdaMax = computeLambdaMaxcI0(model)
        if model.lambdaFWH < lambdaMin
            eq = equilibriumFW(model)
            return "cI = 0 ; EE^FW = " * string(eq) * " is GAS"
        elseif model.lambdaFWH > lambdaMax
            eq = equilibriumHFW(model)
            return "cI = 0 ; Limit Cycle around EE^HFW = " * string(eq)
        elseif lambdaMin < model.lambdaFWH < lambdaMax
            eq = equilibriumHFW(model)
            return "cI = 0 ; EE^HFW = " * string(eq) * " is GAS"
        else
            return "cI = 0 ; lambda equals to lambdaMin or lambdaMax"
        end
    else
        lambdaMax = computeLambdaMax(model)
        if model.lambdaFWH > lambdaMax
            eq = equilibriumH(model)
            return "cI > 0 ; EE^H = " * string(eq) * " is GAS"
        elseif model.lambdaFWH < lambdaMax
            eq = equilibriumHFW(model)
            DeltaStab = computeDeltaStab(model, rsltFeq = false)
            if DeltaStab > 0
                return "cI > 0 ; EE^HFW = " * string(eq) * " is GAS"
            elseif DeltaStab < 0
                return "cI > 0 ; Limit Cycle around EE^HFW = " * string(eq)
            else
                return "cI > 0 ; DeltaStab = 0"
            end
        else
            return "cI > 0 ; lambda = lambdaMax"
        end
    end
end

function computeDeltaStabMatrix(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})
    """
    For cI >= 0
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
    For cI > 0
    Return the bifurcation matrix in function of alpha and lambda
        First line : K_F(1-alpha) ; first column : lambda
    """

    @assert model.c > 0

    bifurcationMatrix = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    bifurcationMatrix[2:end,1] = model.KF * (1 .- listAlpha)

    lambdaMax = computeLambdaMax(model)

    for ind_col in 2:length(listLambdaFWH)+1
        bifurcationMatrix[1, ind_col] = listLambdaFWH[ind_col-1]
        lambdaFWH = listLambdaFWH[ind_col-1]

        if lambdaFWH >= lambdaMax
            bifurcationMatrix[2:end, ind_col] .= "H"
        else
            for ind_row in 2:length(listAlpha)+1
                alpha = listAlpha[ind_row-1]
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

function computeBifurcationDiagramcI0(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listLambdaFWH::Vector{Float64})

    """
    For cI == 0
    Return the bifurcation matrix in function of alpha and lambda
        First line : K_F(1-alpha) ; first column : lambda
    """

    @assert model.c == 0

    bifurcationDiagram = Matrix{Any}(undef, length(listAlpha)+1, length(listLambdaFWH)+1)
    bifurcationDiagram[2:end,1] = model.KF * (1 .- listAlpha)

    lambdaMinTab = Matrix{Any}(undef, length(listAlpha)+1, 2)
    lambdaMinTab[2:end,1] = model.KF * (1 .- listAlpha)

    lambdaMaxTab = Matrix{Any}(undef, length(listAlpha)+1, 2)
    lambdaMaxTab[2:end,1] = model.KF * (1 .- listAlpha)

    for ind_row in 2:length(listAlpha)+1
        alpha = listAlpha[ind_row-1]
        lambdaMin = computeLambdaMincI0(model; alpha)
        lambdaMinTab[ind_row, 2] = lambdaMin
        lambdaMaxStab = computeLambdaMaxcI0(model; alpha)
        lambdaMaxTab[ind_row, 2] = lambdaMaxStab

        for ind_col in 2:length(listLambdaFWH)+1
            bifurcationDiagram[1, ind_col] = listLambdaFWH[ind_col-1]
            lambdaFWH = listLambdaFWH[ind_col-1]

            if lambdaFWH <= lambdaMin
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
    lambdaMaxTab[1,1] = "alpha"
    lambdaMaxTab[1,2] =  "lambdaMax"
    lambdaMinTab[1,1] = "alpha"
    lambdaMinTab[1,2] =  "lambdaMin"

    return bifurcationDiagram, lambdaMinTab, lambdaMaxTab
end

function computeLambdaMatrixAlphaE(model::modelHunter, 
    listAlpha::Vector{Float64}, 
    listE::Vector{Float64})
    """
    For cI = 0
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
            lambdaStab = computeLambdaMaxcI0(model; e=e, alpha=alpha)
            lambdaMaxMatrix[ind_row, ind_col] = lambdaStab

            lambdaExistence = computeLambdaMincI0(model; e=e, alpha=alpha)
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
    For cI = 0
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
            lambdaStab = computeLambdaMaxcI0(model; mW=mW, alpha=alpha)
            lambdaMaxMatrix[ind_row, ind_col] = lambdaStab

            lambdaExistence = computeLambdaMincI0(model; mW=mW, alpha=alpha)
            lambdaMinMatrix[ind_row, ind_col] = lambdaExistence
        end
    end

    lambdaMaxMatrix[1,1] = "alpha" * "\\" * "m"
    lambdaMinMatrix[1,1] = "alpha" * "\\" * "m"

    return lambdaMinMatrix, lambdaMaxMatrix
end
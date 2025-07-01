struct modelHunterHT2 <: mathematicalModel
    variablesNames::Vector{String}

    rF::Float64
    KF::Float64
    alpha::Float64
    beta::Float64
    lambdaFWH::Float64
    theta::Float64
    e::Float64
    I::Float64
    muD::Float64
    fD::Float64
    mD::Float64
    mW::Float64

    m::Float64
        
    function modelHunterHT2(modelParam::Dict{String, Any})

        rF, lambdaFWH, KF = modelParam["rF"], modelParam["lambdaFWH"], modelParam["KF"]
        alpha, beta = modelParam["alpha"], modelParam["beta"]
        theta, e, I = modelParam["theta"], modelParam["e"], modelParam["I"]
        muD, fD, mD, mW = modelParam["muD"], modelParam["fD"], modelParam["mD"], modelParam["mW"]
        variablesNames = ["H_D2d", "F_W2d", "H_W2d", "H_D", "F_W", "H_W"]
        m  = mD / mW

        @assert(alpha < 1)
        bound = 4 * (muD - fD) / (m * e * rF * KF * (1-alpha)^2)
        @assert (bound > beta) "bound for beta is $bound"

        new(variablesNames, rF, KF, alpha, beta, lambdaFWH, theta, e, I, muD, fD, mD, mW, m)
    end
end

function equationModel(model::modelHunterHT2, variables::Vector{Float64})
    """
    Return the right hand side of the model
    """
    HD2d, FW2d, HW2d = variables[1], variables[2], variables[3]
    dHD2d = (model.I + model.e * model.lambdaFWH * FW2d / (1 + model.theta * model.lambdaFWH * FW2d) * model.m * HD2d + 
        (model.fD - model.muD) * HD2d)
    dFW2d = ((1 - model.alpha) * (1 + model.beta * model.m * HD2d) * model.rF * (1 - FW2d / (model.KF * (1 - model.alpha))) * FW2d -
            model.lambdaFWH * FW2d / (1 + model.theta * model.lambdaFWH * FW2d) * model.m * HD2d)

    HD, FW, HW = variables[4], variables[5], variables[6]
    dHD = (model.I + model.e * model.lambdaFWH * FW / (1 + model.theta * model.lambdaFWH * FW) * HW + 
        (model.fD- model.muD) * HD - model.mD * HD + model.mW * HW)
    dFW = ((1 - model.alpha) * (1 + model.beta * HW) * model.rF * (1 - FW / (model.KF * (1 - model.alpha))) * FW -
            model.lambdaFWH * FW / (1 + model.theta * model.lambdaFWH * FW) * HW)
    dHW = model.mD * HD - model.mW * HW

    return [dHD2d, dFW2d, model.m * dHD2d, dHD, dFW, dHW]
end

function computeBetaMax(model::modelHunterHT2 ; alpha = -1.0)
    if alpha < 0
        alpha = model.alpha
    end
    betaMax = 4 * (model.muD - model.fD) / (model.m * model.e * model.rF * model.KF * (1-alpha)^2)
    return betaMax
end

function computeMaxVal(model::modelHunterHT2)
    """
    Compute the bound for the invariant region Omega 
        {H_D + H_W + e F_W < Smax, F_W < Fmax, H_W < m_D/(m_D+m_W) Smax}
    """
    Fmax = model.KF * (1 - model.alpha)
    Smax = (1 + model.m) * (model.I + 
        (model.muD - model.fD + model.rF * (1-model.alpha) / 4) * model.e * model.KF * (1-model.alpha))
    Smax = Smax / ((model.muD - model.fD) / model.m - model.e * model.rF * (1-model.alpha)^2 * model.KF * model.beta / 4)
    
    Hwmax = model.mD / (model.mD + model.mW) * Smax
    
    return (Smax, Fmax, Hwmax)
end

function equilibriumFW(model::modelHunterHT2)
    """
    When I == 0 ;
    return eq EE^F
    """
    return [0, model.KF * (1 - model.alpha), 0]
end

function thresholdNI0(model::modelHunterHT2)
    """
    When I == 0 ;
    return N_{I = 0}
    """
    @assert (model.I ==0) 
    NI = (model.lambdaFWH * model.KF * (1 - model.alpha)) * 
        (model.m * model.e - model.theta * (model.muD - model.fD)) / (model.muD - model.fD)
    return NI
end

function equilibriumH(model::modelHunterHT2)
    """
    When I > 0 ;
    return eq EE^H
    """
    HD = model.I / (model.muD - model.fD)
    return [HD, 0, model.m * HD]
end

function thresholdN(model::modelHunterHT2)
    """
    When I > 0
    return N_{I>0}
    """
    @assert (model.I > 0)
    NI = model.rF * (1-model.alpha) * 
            ((model.muD - model.fD) / (model.m * model.I) + model.beta) / model.lambdaFWH
    return NI
end

function equilibriumHFW(model::modelHunterHT2; nbrExist = false)
    if model.I == 0
        return equilibriumHFWI0(model; nbrExist)
    else
        return equilibriumHFWI(model; nbrExist)
    end
    
end

function equilibriumHFWI0(model::modelHunterHT2; nbrExist = false)
    """
    For I = 0
    Compute the values of EE^HFW
    Rise a warning if negatif (threshold not respected)
    """
    @assert model.I == 0
    if thresholdNI0(model) > 1

        Feq = (model.muD - model.fD) / ((model.m * model.e - model.theta * (model.muD - model.fD ))* model.lambdaFWH)
        
        Heq = (1 - model.alpha) * model.rF * (1 - Feq / (model.KF * (1 - model.alpha)))
        Heq = Heq / (model.m * (model.lambdaFWH - model.beta * (1 - model.alpha) * model.rF + 
                    model.beta * model.rF * Feq / model.KF))
        eq = [[Heq, Feq, model.m * Heq]]
        exist = 1
    else
        @warn "Equilibrium EE^HFW computed while it does not exists ; empty list return "
        eq = []
        exist = 0
    end
    if nbrExist
        return eq, exist
    else
        return eq
    end
end

function equilibriumHFWI(model::modelHunterHT2; nbrExist = false)
    """
    For I > 0
    Compute the values of EE^HFW
    Rise a warning if negatif (threshold not respected)
    """
    @assert model.I > 0
    numA = (model.e * model.m - model.theta * (model.muD - model.fD + model.beta * model.m * model.I))
    A = model.rF / model.KF * numA
    B = -((1 - model.alpha) * model.rF * numA  + 
        (model.muD - model.fD + model.beta * model.m * model.I) * model.rF / (model.lambdaFWH * model.KF))
    C = ((model.muD - model.fD + model.m * model.beta * model.I) * (1 - model.alpha) * model.rF / model.lambdaFWH - 
            model.m * model.I)

    Delta = B^2 - 4 * A * C
    if C > 0 ### equivalent to N_{I > 0} > 1
        if A > 0
            Feq = -B / (2*A) - sqrt(Delta) / (2*A)
        elseif A == 0 ### impossible when theta = 0
            Feq = -C / B
        else ### impossible when theta = 0
            Feq = -B / (2*A) + sqrt(Delta) / (2*A)
        end
        Heq = (1 - model.alpha) * model.rF * (1 - Feq / (model.KF * (1 - model.alpha)))
        Heq = Heq / (model.m * (model.lambdaFWH - model.beta * (1 - model.alpha) * model.rF + 
                model.beta * model.rF * Feq / model.KF))
        eq = [[Heq, Feq, model.m * Heq]]
        exist = 1

    elseif C == 0 ### equivalent to N_{I > 0} = 1
        if (A < 0) && (B > 0) ### A < 0 impossible when theta = 0
            Feq = -B / A
            exist = 1
            Heq = (1 - model.alpha) * model.rF * (1 - Feq / (model.KF * (1 - model.alpha)))
            Heq = Heq / (model.m * (model.lambdaFWH - model.beta * (1 - model.alpha) * model.rF + 
                model.beta * model.rF * Feq / model.KF))
            eq = [[Heq, Feq, model.m * Heq]]
        else
            # @warn "Equilibrium EE^HFW computed while it does not exists ; empty list return "
            eq = []
            exist = 0
        end
        

    else ### equivalent to N_{I > 0} < 1
        if (A < 0) && (B > 0) && (Delta >= 0)  ### A < 0 impossible when theta = 0
            condition5 = model.e * model.m - model.theta * (model.muD - model.fD)
            Flim = (model.muD - model.fD) / (model.lambdaFWH * (model.e * model.m - model.theta * (model.muD - model.fD)))
            
            existence2Eq = ((condition5 <= 0) | (2*A*Flim + B < 0))
            if existence2Eq
                Feq1 = -B / (2*A) - sqrt(Delta) / (2*A)
                Feq2 = -B / (2*A) + sqrt(Delta) / (2*A)

                Heq1 = (model.I) / (model.muD - model.fD - (model.e * model.lambdaFWH * model.m) / (1+model.theta * model.lambdaFWH * Feq1) * Feq1)
                Heq2 = (model.I) / (model.muD - model.fD - (model.e * model.lambdaFWH * model.m) / (1+model.theta * model.lambdaFWH * Feq2) * Feq2)

                eq = [[Heq2, Feq2, model.m * Heq2], [Heq1, Feq1, model.m * Heq1]]
                exist = 2
            end
        else
            # @warn "Equilibrium EE^HFW computed while it does not exists ; empty list return "
            eq = []
            exist = 0
        end
    end

    if nbrExist
        return eq, exist
    else
        return eq
    end
end

function computeLambdaMinI0(model::modelHunterHT2 ; alpha = -1.0, e = -1.0, mW = -1.0)
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
    lambdaMin = (model.muD - model.fD) / ((1 - alpha) * model.KF * (m * e - model.theta * (model.muD - model.fD)))
    
    return lambdaMin
end

function computelambdaMaxI0(model::modelHunterHT2)
    """
    When I = 0 and theta > 0
    or I = theta = beta = 0, 
    compute lambda^max such that lambda < lambda^max iff EE^HFW is AS
    """
    if model.theta > 0
        num = (model.m * model.e + model.theta * (model.muD - model.fD))
        den = (model.m * model.e - model.theta * (model.muD - model.fD))
        return num / den / (model.KF * (1-model.alpha) * model.theta)
    else
        if model.beta == 0
            a0 = (-(model.muD - model.fD + model.mD + mW) * 
                    model.rF * (model.muD - model.fD) / (e * m * model.KF))
            a1 = -(mW * (model.muD - model.fD) + 
            (model.muD - model.fD + model.mD + mW)^2)
            a2 = (1 - alpha) * model.KF * e * model.mD

            PDeltaStab = Polynomial([a0, a1, a2])
            rootsPDeltaStab = real(PolynomialRoots.roots(coeffs(PDeltaStab)))
            lambdaMaxStab = maximum(rootsPDeltaStab)

            return lambdaMaxStab
        else
            @warn "lambda max I = 0 computed while not defined; infinite value return"
            return Inf64
        end
    end
end

function computeLambdaMax(model::modelHunterHT2 ; alpha = -1.0, mD = -1.0, e = -1.0)
    """
    When theta = 0, I > 0, compute lambda^max such that lambda < lambda^max iff EE^HFW exists
    """
    if model.I > 0
        if alpha < 0
            alpha = model.alpha
        end
        if mD < 0
            mD = model.mD
        end
        if e < 0
            e = model.e
        end
        m = mD / model.mW
        lambdaMax = (1 - alpha) * model.rF * ( model.beta + 
            (model.muD - model.fD) / (m * model.I))
    else
        @warn "lambda max I > 0 computed while not defined; infinite value return"
        lambdaMax = Inf64
    end
    return lambdaMax
end

function computeLambdaMaxI0(model::modelHunterHT2 ; alpha = -1.0, e = -1.0, mW = -1.0)
    """
    When I = \beta = 0, compute lambda^max such that lambda < lambda^max iff EE^HFW is AS
    """    
    @assert (model.I == 0) && (model.beta == 0) && (model.theta == 0) "beta : $(model.beta) and I : $(model.I)" 
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

function computeEigenVal2D(model::modelHunterHT2, val::Vector{Float64}; eigenvector = false,)
    HD, FW = val[1], val[2]
    Jac = Matrix{Float64}(undef, 2, 2)
    Jac[1,1] = -(model.muD - model.fD) + model.e * model.m * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) * FW
    Jac[1,2] = model.e * model.m * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HD
    Jac[2,1] = (-model.m * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) + 
                model.m * model.beta * (1 - model.alpha) * model.rF * (1 - FW / ((1 - model.alpha) * model.KF))) * FW
    Jac[2,2] = ((1 - model.alpha) * model.rF * (1 + model.beta * model.m * HD) * (1 - 2 * FW / ((1 - model.alpha) * model.KF)) -
                model.m * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HD)

    if eigenvector
        EV = eigvecs(Jac)
        EVA = eigvals(Jac)
        println("eigenvec : ", EV[:, 2])
        println("test : ", norm(EV[:, 2] - Jac * EV[:, 2] / EVA[2]))
        return eigvals(Jac), eigvecs(Jac)
    else
        return eigvals(Jac)
    end
end

function computeTrace2D(model::modelHunterHT2, val::Vector{Float64})
    HD, FW = val[1], val[2]
    J11 = -(model.muD - model.fD) + model.e * model.m * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) * FW
    J22 = ((1 - model.alpha) * model.rF * (1 + model.beta * model.m * HD) * (1 - 2 * FW / ((1 - model.alpha) * model.KF)) -
                model.m * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HD)
    return J11 + J22
end

function computeEigenVal3D(model::modelHunterHT2, val::Vector{Float64})
    HD, FW, HW = val[1], val[2], val[3]
    Jac = Matrix{Float64}(undef, 3, 3)
    Jac[1,1] = -(model.muD - model.fD + model.mD)
    Jac[1,2] = model.e * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HW
    Jac[1,3] = model.e * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) * FW + model.mW
    
    Jac[2,1] = 0
    Jac[2,2] = ((1 - model.alpha) * model.rF * (1 + model.beta * HW) * (1 - 2 * FW / ((1 - model.alpha) * model.KF)) -
                model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HW)
    Jac[2,3] = (-model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) + 
                model.beta * (1 - model.alpha) * model.rF * (1 - FW / ((1 - model.alpha) * model.KF))) * FW
    
    Jac[3,1] = model.mD
    Jac[3,2] = 0
    Jac[3,3] = model.mW

    return eigvals(Jac)
end

function computeDeltaStab(model::modelHunterHT2 ; lambdaFWH = -1., alpha = -1.0, mD = -1.0, e = -1.0, rsltFeq = true)
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
    if mD < 0
        mD = model.mD
    end
    if e < 0
        e = model.e
    end

    m = mD / model.mW

    ## Compute the values at equilibrium with eventually
    ## new parameters values 

    A = e * model.rF / model.KF
    B = -(e * (1 - alpha) * model.rF + (model.muD - model.fD) * model.rF / (lambdaFWH * m * model.KF)
            + model.I * model.beta * model.rF / (lambdaFWH * model.KF))
    C = ((model.muD - model.fD) * (1 - alpha) * model.rF / (lambdaFWH * m) - 
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
    a2 = -(model.fD - model.muD - mD - model.mW - 
        (1 + model.beta * HWeq) * model.rF * Feq / model.KF)

    a0 = e * lambdaFWH * mD * model.rF * (
        (model.muD - model.fD) / (e * m * model.KF * lambdaFWH) - 2 * Feq / model.KF + 
        1 - alpha + model.beta * HWeq * 
                ((model.muD - model.fD) / (e * m * model.KF * lambdaFWH) - Feq / model.KF)
        ) * Feq
    
    a1 = ((model.muD - model.fD + mD + model.mW) * model.rF * (1 + model.beta * HWeq) *
            Feq / model.KF + mD * e * lambdaFWH * 
                ((model.muD - model.fD) / (m * lambdaFWH * e) - Feq) )

    @assert (a1 > 0 && a2 > 0 && a0 > 0) "alpha : $alpha, mD : $mD, lambda : $lambdaFWH,  a0 : $a0, a1 : $a1, a2 : $a2"
    println("DeltaSatb : a2 :" , a2, " a1 : ", a1, " a0 : ", a0)
    println(Feq)

    if rsltFeq
        return a2 * a1 - a0, Feq
    else
        return a2 * a1 - a0
    end
end

function computeDeltaStab3D(model::modelHunterHT2, val::Vector{Float64})
    HD, FW, HW = val[1], val[2], val[3]
    Jac = Matrix{Float64}(undef, 3, 3)
    Jac[1,1] = -(model.muD - model.fD + model.mD)
    Jac[1,2] = model.e * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HW
    Jac[1,3] = model.e * model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) * FW + model.mW
    
    Jac[2,1] = 0
    Jac[2,2] = ((1 - model.alpha) * model.rF * (1 + model.beta * HW) * (1 - 2 * FW / ((1 - model.alpha) * model.KF)) -
                model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW)^2 * HW)
    Jac[2,3] = (-model.lambdaFWH / (1 + model.lambdaFWH * model.theta * FW) + 
                model.beta * (1 - model.alpha) * model.rF * (1 - FW / ((1 - model.alpha) * model.KF))) * FW
    
    Jac[3,1] = model.mD
    Jac[3,2] = 0
    Jac[3,3] = -model.mW

    a0 = -det(Jac)
    a2 = -(Jac[1,1] + Jac[2,2] + Jac[3,3])
    a1 = Jac[1,1] * Jac[2,2] + Jac[1,1] * Jac[3,3] + Jac[2,2] * Jac[3,3] - Jac[1,3] * Jac[3,1]
    # println("DeltaStab3D : a2 :" , a2, "a1 : ", a1, "a0 : ", a0)
    return a2 * a1 - a0
end

function longTermDynamicStr(model::modelHunterHT2)
    if model.I == 0
        return longTermDynamicStrI0(model)
    else
        return longTermDynamicStrI(model)

    end
end

function longTermDynamicStrI0(model::modelHunterHT2)
    @assert model.I == 0
    lambdaMin = computeLambdaMinI0(model)
    if model.lambdaFWH <= lambdaMin
        return "F"
    else ## coexistence exist, F is unstable
        if model.theta > 0 ## ie 2D case
            lambdaMax = computelambdaMaxI0(model)
            if model.lambdaFWH < lambdaMax
                return "HFW"
            else
                return "LC"
            end
        else ## 3D case
            listEq = equilibriumHFW(model)
            DeltaStab = computeDeltaStab3D(model, listEq[1])
            if DeltaStab > 0
                return "HFW"
            else
                return "LC"
            end
        end
    end
end

function longTermDynamicStrI(model::modelHunterHT2)
    @assert model.I > 0
    NI = thresholdN(model)
    if NI > 1 ## H is unstable ; one HFW equilibrium
        listEq = equilibriumHFW(model)
        if model.theta > 0
            Tr = computeTrace2D(model, listEq[1])
            if Tr < 0
                return "HFW"
            else
                return "LC"
            end
        else
            DeltaStab = computeDeltaStab3D(model, listEq[1])
            if DeltaStab > 0
                return "HFW"
            else
                return "LC"
            end
        end
    else ## H is Astable ; zero or two HFW equilibrium, one being unstable
        rslt = "H"
        if model.theta > 0
            listEq = equilibriumHFW(model)
            if length(listEq) > 0
                Tr = computeTrace2D(model, listEq[2])
                if Tr < 0
                    rslt = rslt * "HFW2" * "?" #dont know if there is a LC
                else
                    rslt = rslt * "??"
                end
            end
        end  # no need to look theta = 0 case since there is 0 HFW equilibrium in this case
        return rslt
    end
end

function computeBifurcationDiagram(model::modelHunterHT2, 
    paramName1::String, listParam1::Vector{Float64},
    paramName2::String, listParam2::Vector{Float64})
    
    """
    For I >= 0
    Return the bifurcation matrix in function of param 1 & param 2
        First line : param 1 ; first column : param 2
    """

    bifurcationMatrix = Matrix{Any}(undef, length(listParam1)+1, length(listParam2)+1)
    bifurcationMatrix[2:end,1] = listParam1

    localParam = Dict(string(name) => getfield(model, name) for name in fieldnames(modelHunterHT2))
    delete!(localParam, "variablesNames")

    for (indRow, param1Val) in enumerate(listParam1)
        localParam[paramName1] = param1Val
       
        for (indCol, param2Val) in enumerate(listParam2)
            
            bifurcationMatrix[1, indCol + 1] = param2Val
            
            localParam[paramName2] = param2Val
            localModel = modelHunterHT2(localParam)
            
            bifurcationMatrix[indRow + 1, indCol + 1] = longTermDynamicStr(localModel)
        end
    end

    bifurcationMatrix[1,1] = paramName1 * "\\" * paramName2

    return bifurcationMatrix
end

function computeLambdaMatrixAlpha(model::modelHunterHT2, 
    listAlpha::Vector{Float64})
    """
    For I >= 0
    Return the bifurcation matrix in function of alpha
        First column : alpha
    """
    if model.I == 0
        bifurcationMatrixs = computeLambdaMatrixAlphaI0(model, listAlpha)
    else
        bifurcationMatrixs = computeLambdaMatrixAlphaI(model, listAlpha)
    end
    return bifurcationMatrixs
end

function computeLambdaMatrixAlphaI(model::modelHunterHT2, 
    listAlpha::Vector{Float64})
    """
    For I > 0
    Compute the values of lambda^max and lambda^min as function of parameters alpha
    Return two matrixs, which contains the values
    """
    lambdaMinMatrix = Matrix{Any}(undef, 1, 2) #lambda min is not defined for I > 0 -> empty matrix

    lambdaMaxMatrix = Matrix{Any}(undef, length(listAlpha)+1, 2)
    lambdaMaxMatrix[2:end,1] = listAlpha

    for ind_row in 2:length(listAlpha)+1
        alpha = listAlpha[ind_row-1]
        lambdaExistence = computeLambdaMax(model; alpha = alpha)
        lambdaMaxMatrix[ind_row, 2] = lambdaExistence
    end

    lambdaMaxMatrix[1,1] = "alpha"
    lambdaMaxMatrix[1,2] = "lambdaMAx"
    lambdaMinMatrix[1,1] = "alpha"
    lambdaMinMatrix[1,2] = "lambdaMin"

    return lambdaMinMatrix, lambdaMaxMatrix
end

function computeLambdaMatrixAlphaI0(model::modelHunterHT2, 
    listAlpha::Vector{Float64})
    """
    For I = beta = 0
    Compute the values of lambda^max and lambda^min as function of parameters alpha
    Return two matrixs, which contains the values
    """

    if (model.beta) == 0 && (model.theta == 0)
        lambdaMaxMatrix = Matrix{Any}(undef, length(listAlpha)+1, 2)
        lambdaMaxMatrix[2:end,1] = listAlpha
    else     #else lambda max is not defined, empty matrix return
        lambdaMaxMatrix = Matrix{Any}(undef, 1, 2)
    end

    lambdaMinMatrix = Matrix{Any}(undef, length(listAlpha)+1, 2)
    lambdaMinMatrix[2:end,1] = listAlpha

    for ind_row in 2:length(listAlpha)+1
        alpha = listAlpha[ind_row-1]
        if (model.beta) == 0 && (model.theta == 0)
            lambdaStab = computeLambdaMaxI0(model; alpha = alpha)
            lambdaMaxMatrix[ind_row, 2] = lambdaStab
        end

        lambdaExistence = computeLambdaMinI0(model; alpha = alpha)
        lambdaMinMatrix[ind_row, 2] = lambdaExistence
    end

    lambdaMaxMatrix[1,1] = "alpha"
    lambdaMaxMatrix[1,2] = "lambdaMAx"
    lambdaMinMatrix[1,1] = "alpha"
    lambdaMinMatrix[1,2] = "lambdaMin"

    return lambdaMinMatrix, lambdaMaxMatrix
end
import numpy as np
import matplotlib.pyplot as plt
import csv

def equationTest(variables, param):
    """
    variables is a 3 length np.array, param is a dict
    Return the right hand side of the test model
    """
    dF = 1
    dV = 2
    dH = 3
    return np.array([dF, dV, dH])

def equationModel(variables, param):
    """
    variables is a 3 length np.array, param is a dict
    Return the right hand side of the model
    """
    F, V, H = variables[0], variables[1], variables[2]
    dF = param["rF"] * (1-F/param["KF"])*F - param["omega"]*param["f"]*F - param["muF"]*F - param["lambdaFH"]*F*H
    dV = param["rV"] * (1-V/param["KV"])*V - param["alpha"]*V*F - param["muV"]*V - param["lambdaVH"]*V*H
    dH = param["e"] * (param["lambdaFH"]*F + param["lambdaVH"]*V)*H - param["muH"]*H**2
    return np.array([dF, dV, dH])

def rungeKutta4(equation, variables, dt, param):
    """
    Apply RK4 method over dt using equationModel function
    """
    k1 = equation(variables, param)
    k2 = equation(variables + dt/2*k1, param)
    k3 = equation(variables + dt/2*k2, param)
    k4 = equation(variables + dt*k3, param)

    return variables + dt/6*(k1 + 2*k2 + 2*k3 + k4)

def solveModel(equation, init, n, dt, param):
    result = np.zeros((n+1, 4))
    result[0] = init
    for ind in range(n):
        result[ind+1][0] = result[ind][0] + dt
        result[ind+1][1:] = rungeKutta4(equation, result[ind][1:], dt, param)
    return result

def computeEquilibria(param, equilibria=["VH", "FH", "FVH"]):
    """
    Compute all possible stable equilibria of equationModel
    """
    rV, KV, alpha, muV, lambdaVH = param["rV"], param["KV"], param["alpha"], param["muV"], param["lambdaVH"]
    rF, KF, omega, f, muF, lambdaFH = param["rF"], param["KF"], param["omega"], param["f"], param["muF"], param["lambdaFH"]
    e, muH = param["e"], param["muH"]

    dV = (rV - muV)
    dF = (rF - muF - omega*f)

    rslt = {}
    for eq in equilibria:
        if eq=="VH":
            V_VH = 1/(rV/KV + e*lambdaVH**2/muH)*dV
            H_VH = e*lambdaVH/muH*V_VH
            rslt[eq] = [0, V_VH, H_VH]
        elif eq=="FH":
            F_FH = 1/(rF/KF + e*lambdaFH**2/muF)*dF
            H_FH = e*lambdaFH/muH*F_FH
            rslt[eq] = [F_FH, 0, H_FH]
        elif eq=="FVH":
            V_0 = KV/rV*(dV - lambdaVH/lambdaFH*dF)
            V_1 = lambdaVH/lambdaFH*rF/rV*KV/KF - alpha*KV/rV
            H_0 = dF/lambdaFH
            H_1 = rF/(lambdaFH*KF)
            F_FVH = (muH/e*H_0 - lambdaVH*V_0) / (lambdaFH + lambdaVH*V_1 + muH/e*H_1)
            V_FVH = V_0 + V_1*F_FVH
            H_FVH = H_0 - H_1*F_FVH
            rslt[eq] = [F_FVH, V_FVH, H_FVH]

    return rslt

def stabilityCondition(param):
    """
    Compute the stabilities conditions for equationModel
    """
    rV, KV, alpha, muV, lambdaVH = param["rV"], param["KV"], param["alpha"], param["muV"], param["lambdaVH"]
    rF, KF, omega, f, muF, lambdaFH = param["rF"], param["KF"], param["omega"], param["f"], param["muF"], param["lambdaFH"]
    e, muH = param["e"], param["muH"]
    R0V = rV/muV
    R0F = rF/(muF + omega*f)
    TF = (lambdaVH/muV) * (omega*f + muF)/lambdaFH * (R0F - 1) / (R0V - 1) * (1 + muH/(e*lambdaVH**2)*rV/KV)
    TV = (muV/lambdaVH) * lambdaFH/(omega*f + muF) * (R0V - 1) / (R0F - 1) * (1 + muH/(e*lambdaFH**2)*rF/KF) / (1 + alpha*muH/(e*lambdaFH*lambdaVH))

    return {"R_0^V":R0V, "R_0^F":R0F, "T^F":TF, "T^V":TV}

def printStabilityCondition(stabilityCondition):
    """
    Print which equilibrium is stable
    """
    R0V, R0F, TV, TF = stabilityCondition["R_0^V"], stabilityCondition["R_0^F"], stabilityCondition["T^V"], stabilityCondition["T^F"]
    if (R0V < 1) and (R0F < 1) :
        print("TE is stable")
    elif (R0V > 1) and (R0F > 1) :
        if (TF < 1) and (TV < 1):
            print("VH and FH are asymptotically stable")
        elif (TF > 1) and (TV > 1):
            print("FVH is asymptotically stable")
        elif (TF < 1) and (TV > 1):
            print("VH is asymptotically stable")
        elif (TF > 1) and (TV < 1):
            print("FH is asymptotically stable")
    elif (R0V > 1) and (R0F < 1):
        if (TF < 1):
            print("VH is asymptotically stable")
    elif (R0V < 1) and (R0F > 1):
        if (TV < 1):
            print("FH is asymptotically stable")
    else:
        print("No equilibrium is stable")

def writeResult(filename, resultSimu, param, paramSimu, writeStabilityCondition=False):
    """
    Write the results on a csv file
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        writer.writerow(['####_Global_Param_####'])
        writer.writerow(list(param.keys()))
        writer.writerow(list(param.values()))
        writer.writerow(["####_Simulation_Param_####"])
        writer.writerow(list(paramSimu.keys()))
        writer.writerow(list(paramSimu.values()))

        if writeStabilityCondition:
            writer.writerow(["####_Stability_Condition_####"])
            stabilityCond = stabilityCondition(param)
            writer.writerow(list(stabilityCond.keys()))
            writer.writerow(list(stabilityCond.values()))

        writer.writerow(["############"])
        writer.writerow(["Time", "F", "V", "H"])
        for result in resultSimu:
            writer.writerow([result[0], result[1], result[2], result[3]])

def plotVectorField(equilibrium, nameEquilibrium, param, boxSize, step):
    [F, V, H] = equilibrium
    coord_F = np.arange(max(F-boxSize, 0), F+boxSize, step)
    coord_V = np.arange(max(V-boxSize, 0), V+boxSize, step)
    coord_H = np.arange(max(H-boxSize, 0), H+boxSize, step)
    x, y, z = np.meshgrid(coord_F,
                            coord_V,
                            coord_H, indexing='ij')
    dF = np.zeros((len(coord_F), len(coord_V), len(coord_H)))
    dV = np.zeros((len(coord_F), len(coord_V), len(coord_H)))
    dH = np.zeros((len(coord_F), len(coord_V), len(coord_H)))

    for ind_F in range(len(coord_F)):
        for ind_V in range(len(coord_V)):
            for ind_H in range(len(coord_H)):
                val = equationModel([coord_F[ind_F], coord_V[ind_V], coord_H[ind_H]], param)
                dF[ind_F, ind_V, ind_H] = val[0]
                dV[ind_F, ind_V, ind_H] = val[1]
                dH[ind_F, ind_V, ind_H] = val[2]
    
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(x, y, z, dF, dV, dH, normalize=True, length=1)
    ax.plot(F, V, H, label=r'$%s$'%nameEquilibrium, marker='x', color = 'red', linestyle='')
    ax.set(xlabel='F', ylabel='V', zlabel='H')
    ax.legend()
    plt.show()

def plotResult(equilibrium, nameEquilibrium, listResult):
    [F, V, H] = equilibrium
    ax = plt.figure().add_subplot(projection="3d")
    ax.plot(F, V, H, label=r'$%s$'%nameEquilibrium, marker='x', color = 'red', linestyle='')

    for result in listResult:
        ax.plot(result[:,1], result[:,2], result[:,3])
    
    
    ax.set(xlabel='F', ylabel='V', zlabel='H')
    ax.legend()
    plt.show()
        
def main():

    ## Param given by Yatat, 2021
    # param = {"rV":1.8, "KV":19.9, "alpha":0.01, "muV":0.1, "rF":0.71, "KF":429.2, "omega":0.1, "f":1, "muF":0.1, "e":0.8}

    ## Param for FVH stable
    paramFVH = {"rV":1.8, "KV":19.9, "alpha":0.01, "muV":0.1, "rF":0.71, "KF":429.2, "omega":0.1, "f":1, "muF":0.1, "e":0.8}
    paramFVH["muH"] = 0.01
    paramFVH["lambdaVH"] = 0.001
    paramFVH["lambdaFH"] = 0.1

    ## Param for VH stable
    # paramVH = {"rV":1.8, "KV":19.9, "alpha":0.01, "muV":0.1, "rF":0.71, "KF":429.2, "omega":0.1, "f":1, "muF":0.1, "e":0.8}
    # paramVH["muH"] = 0.01
    # paramVH["lambdaVH"] = 0.1
    # paramVH["lambdaFH"] = 0.1

    # print(stabilityCondition(paramVH))
    printStabilityCondition(stabilityCondition(paramFVH))

    t0, tf = 0., 100.
    n = 1000
    dt = (tf-t0)/n
    [F, V, H] = computeEquilibria(paramFVH, ["FVH"])["FVH"]
    F0, V0, H0 = F +3, V + 5, H+10
    F1, V1, H1 = 0.1, V -5, H
    
    
    init = np.array([t0, F0, V0, H0])
    init1 = np.array([t0, F1, V1, H1])
    resultSimu = solveModel(equationModel, init, n, dt, paramFVH)
    resultSimu1 = solveModel(equationModel, init1, n, dt, paramFVH)

    paramSimu = {"t0":t0, "tf":tf,"dt":dt,"F0":F0, "V0":V0, "H0":H0}
    plotResult([F,V,H], "EE^{FVH}", [resultSimu, resultSimu1])
    writeResult("ModelFVH.csv", resultSimu, paramFVH, paramSimu, True)

    return

if __name__ == "__main__":
    main()

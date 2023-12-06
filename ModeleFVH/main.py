import numpy as np
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
            stabilityCondition = stabilityCondition(param)
            writer.writerow(list(stabilityCondition.keys()))
            writer.writerow(list(stabilityCondition.values()))
            
        writer.writerow(["############"])
        writer.writerow(["Time", "F", "V", "H"])
        for result in resultSimu:
            writer.writerow([result[0], result[1], result[2], result[3]])
        
def main():
    t0, tf = 0., 10.
    n = 10
    dt = (tf-t0)/n

    F0, V0, H0 = 1, 0, 0
    init = np.array([t0, F0, V0, H0])
    param = {}
    resultSimu = solveModel(equationTest, init, n, dt, param)
    writeResult("ModelTest.csv", resultSimu, param, {"t0":t0, "tf":tf,"dt":dt,"F0":F0, "V0":V0, "H0":H0})


if __name__ == "__main__":
    main()

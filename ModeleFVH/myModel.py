import numpy as np

def equationTest(variables, param):
    """
    variables is a 3 length np.array, param is a dict
    Return the right hand side of the test model
    """
    dF = 1
    dV = 2
    dH = 3
    return np.array([dF, dV, dH])

def equationFVH(variables, param):
    """
    variables is a 3 length np.array, param is a dict
    Return the right hand side of the model F,V,H
    """
    F, V, H = variables[0], variables[1], variables[2]
    dF = param["rF"] * (1-F/param["KF"])*F - param["omega"]*param["f"]*F - param["muF"]*F - param["lambdaFH"]*F*H
    dV = param["rV"] * (1-V/param["KV"])*V - param["alpha"]*V*F - param["muV"]*V - param["lambdaVH"]*V*H
    dH = param["e"] * (param["lambdaFH"]*F + param["lambdaVH"]*V)*H - param["muH"]*H**2
    return np.array([dF, dV, dH])

def equationRV(variables, param):
    """
    variables is a 3 length np.array, param is a dict
    Return the right hand side of the model F,V,H with rV = fct(H)
    """
    F, V, H = variables[0], variables[1], variables[2]
    dF = param["rF"] * (1-F/param["KF"])*F - param["omega"]*param["f"]*F - param["muF"]*F - param["lambdaFH"]*F*H
    dV = param["rV"] * (H / (H+param["H0"])) * (1-V/param["KV"])*V - param["alpha"]*V*F - param["muV"]*V - param["lambdaVH"]*V*H
    dH = param["e"] * (param["lambdaFH"]*F + param["lambdaVH"]*V)*H - param["muH"]*H**2
    return np.array([dF, dV, dH])
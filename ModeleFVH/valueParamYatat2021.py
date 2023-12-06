import numpy as np

"""
This file compute the parameters value using the ones given in Yatat, 2021, table 4.
When a value was not available (because it depends on the variables of the model), it is taken equal at the average of its possible range.
"""

# Variables in table 4
cg, ct, bg, bt, ag, at = 20, 430, 500, 1100, 0.0029, 0.004
dg, dt, gg, gt, dg, dt = 14.73, 107, 2.7, 1.5, 0.1, 0.1
lambdamin, lambdamax, a, b, c, d = 0.005, 0.4, 0.01, 600, 120, 0.0045

# Average of other variables
W = 1000
vT = (lambdamax+lambdamin)/2
omegaG = 1/2

# Variables for model FVH
rV, rF = gg*W/(bg + W), gt*W/(bt + W)
KV, KF = cg/(1+dg*np.exp(-ag*W)), ct/(1+dt*np.exp(-at*W))
muV, muF = dg, dt
alpha = a*np.tanh((W-b)/c) + d
omega, f = vT*omegaG, 1

param = {"rV":rV, "KV":KV, "alpha":alpha, "muV":muV, "rF":rF, "KF":KF, "omega":omega, "f":f, "muF":muF}
print(param)
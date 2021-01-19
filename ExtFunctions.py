from scipy import interpolate
import numpy as np
def interpCores(valMin, valMax, d_in, x_in, DeltaInput = False, DeltaIn = 0):


    d = d_in[(d_in >= valMin) & (d_in <= valMax)]
    x = x_in[(d_in >= valMin) & (d_in <= valMax)]

    if DeltaInput:
        Delta = DeltaIn
    else:
        diff = np.diff(d)
        Delta = round(min(diff), 3)

    d_min = Delta * np.ceil(d.values[0]/Delta)
    d_max = Delta * np.floor(d.values[-1]/Delta)

    n = int(1 + (d_max - d_min)/Delta)

    j_arr = np.linspace(0,n,n)
    dhat = d_min + (j_arr - 1)*Delta

    f = interpolate.CubicSpline(d,x)

    xhat = f(dhat)

    return dhat, xhat, Delta

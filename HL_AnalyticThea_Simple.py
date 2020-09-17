import numpy as np
import matplotlib.pyplot as plt
import time
'''
    Function to model density and age profiles for a given depth, using Herron-Langway dynamics.
    Models bubble close-off at density 0.550. Uses m in water equivalent (w.e.) per yea (a^-1)

    Input parameters:
        :param z_vec    : [array of floats] Depth vector, consisting of depths between top and bottom depths. In [m]
        :param rhoSurf  : [float] Surface density (e.g. average of first 15 m of firn). In [g/m^3]
        :param bdot0    : [float] Accumulation rate, surface. In [kg m^-2 a^-1]
        :param Temp0    : [float] Surface temperature. In [K]

    Output:
        :rho_h          : [array of floats] Modeled density profile, corresponding to depth. [g/m^3]
        :t              : [array of floats] Modeled age profile, corresponding to depth. [y]

    Function based on M. Herron & C. Langway "Firn Densification: An Empirical Model" (1980)
'''


def HL_analyticThea(z_vec, rhoSurf, bdot0, Temp0):
    # Defining relevant constants
    R           = 8.314             # Gas constant for computing Arrhenius rate constants
    rho_ice     = 0.917             # Density of ice in [g m^-3]
    rho055      = 0.550             # Density at close-off
#    sPerYear    = 31557600.0        # Seconds per year
    Acc         = bdot0 * rho_ice   # Accumulation from accumulation rate
    zSize       = np.size(z_vec)    # Length of depth array

    # Calculating rate constants (eq. 6a, 6b)
    k0 = 11.0 * np.exp(-10160 / (R * Temp0))
    k1 = 575.0 * np.exp(-21400 / (R * Temp0))

    # Computing depth and age of critical density, bubble close-off, 0.55 (eq. 8, 9)
    h055 = (1 / (rho_ice * k0)) * (np.log(rho055 / (rho_ice - rho055)) - np.log(rhoSurf / (rho_ice - rhoSurf)))
    t055 = (1 / (k0 * Acc)) * np.log((rho_ice - rhoSurf) / (rho_ice - rho055))

    # From critical depth h055, rate constant k1, constants and depth vector, compute Z0 and Z1 densification stage constants
    Z0 = np.exp(rho_ice * k0 * z_vec + np.log(rhoSurf / (rho_ice - rhoSurf)))
    Z1 = np.exp(((rho_ice * k1 * (z_vec - h055)) / (Acc**(0.5))) + np.log(rho055 / (rho_ice - rho055)))
    # concatenate Z0 and Z1 to mae vector to describe both stage 1 and 2 densification
    Z = np.concatenate((Z0[z_vec < h055], Z1[z_vec > h055]))

    # From Z dens. stage constants, compute the density (eq. 7, 10)
    rho_h = (rho_ice * Z) / (1 + Z)
    rho_h_age = np.copy(rho_h)

    # If any value in the density vector is larger than the density of ice,
    # then set all values hereafter to last density smaller than ice density
    if max(rho_h) >= rho_ice - 0.001:
        idx_Ice_rho = len(rho_h[rho_h <= rho_ice - 0.001])
        idx_Ice_age = len(rho_h[rho_h < rho_ice - 0.01])
        rho_hMax = rho_h[idx_Ice_rho]
        rho_h_ageMax = rho_h[idx_Ice_age]
        rho_h[idx_Ice_rho:] = rho_hMax * np.ones(len(rho_h[idx_Ice_rho:]))
        rho_h_age[idx_Ice_age:] = rho_h_ageMax * np.ones(len(rho_h_age[idx_Ice_age:]))

    # Compute the age profile through (eq. 9, 11)
    t0 = (1 / (k0 * Acc)) * np.log((rho_ice - rhoSurf) / (rho_ice - rho_h_age))
    t1 = (1 / (k1 * np.sqrt(Acc))) * np.log((rho_ice - rho055) / (rho_ice - rho_h_age)) + t055
    # Choose t0, if depth above crit. dens. boundary, choose t1 if above
    t = np.concatenate((Z0[z_vec < h055], Z1[z_vec > h055]))

    # Return density and age profiles
    return rho_h, t




# H_bot       = 3400
# H_base      = 3100
# bdot0       = 0.13100
# stpsPerYear = 1
# rhoSurf     = 350 / 1000.0
# Temp0       = 242
# rho_ice     = 0.917
#
# num_gridPoints = int((H_bot - H_base) / (bdot0 / stpsPerYear))
# H_vec = np.linspace(H_bot, H_base, num_gridPoints)
# z_vec = H_bot - H_vec
# zSize = np.size(z_vec)
# print(z_vec)
# start_time = time.time()
# dens, age = HL_analyticThea(z_vec, rhoSurf, bdot0, Temp0)
# print("--- %s seconds ---" % (time.time() - start_time))
#
#
# fig, ax = plt.subplots(1, 2, figsize=(8,8))
# ax[0].set(xlabel='density', ylabel='depth')
# ax[0].plot(dens*1000, z_vec)
# ax[0].invert_yaxis()
# ax[0].axvline(x=rho_ice*1000,color='k')
# ax[1].set(xlabel='age', ylabel='depth')
# ax[1].plot(age, z_vec)
# ax[1].invert_yaxis()
# #ax[1].axhline(y=z_vec[idx_Ice_age], color='k')
# fig.tight_layout()
# fig.savefig('densTest.png', dpi=600)

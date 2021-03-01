import sys
import copy
import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import optimize

import DiffusionProfiles_calculations



class sigma_Solver():
    '''
        Class to make temperature/accumulation rate estimates based on data and theory
        diffusion length estimates.

        Numerically finds the roots of Eq. (15) in [Holme, 2018].
    '''

    def __init__(self, rho_ice = 917., rho_Cr = 550., rho_surf = 330.):
        """
            Initialization of sigma_Solver class.
        """

        self.rho_ice = rho_ice
        self.rho_Cr = rho_Cr
        self.rho_surf = rho_surf

        return

    def norm_function(self, temp = 243.15, accum = 0.18, P_atm = 0.7, rho_CO = np.array([804.3]),\
                    f0 = 1, f1 = 1, fractFact_18O = "Majoube", fractFact_D = "Merlivat", fractFact_17O = "Majoube", \
                    sigma_data = 0.075, isotope = 2):

        """
            Function definition to find roots of, (rho_CO/rho)^2 * sigma(T(z), Accum)^2 - sigma_data^2 = 0.
            Uses an analytical Herron-Langway solution to the diffusion equation.
        """

        sigma_model = DiffusionProfiles_calculations.DiffusionLength(P = P_atm, rho_surf = self.rho_surf, f0 = f0, \
                                  f1 = f1, fractFact_18O = fractFact_18O, fractFact_D = fractFact_D, \
                                  fractFact_17O = fractFact_17O).analytical_HL(rho_arr = rho_CO, T = temp, \
                                  accum = accum)[isotope - 1]

        sigma_model = sigma_model * (rho_CO/self.rho_ice)

        return sigma_model - sigma_data

    def solveTemp(self, sigma_data = 0.08, accum = 0.18, P_atm = 0.7, rho_CO = np.array([804.3]),\
                f0 = 1, f1 = 1, fractFact_18O = "Majoube", fractFact_D = "Merlivat", \
                fractFact_17O = "Majoube", isotope = 2):

        """
            Method to solve temperature estimate from theoretical and data determined
            diffusion length estimates.
            Accumulation is in water equivalent, diffusion length is in ice equivalent.
        """

        T_est = sp.optimize.newton(self.norm_function, 243.5, \
                                    args = (accum, P_atm, rho_CO, f0, f1, fractFact_18O, \
                                            fractFact_D, fractFact_17O, sigma_data, isotope), \
                                    maxiter = 100)

        return T_est

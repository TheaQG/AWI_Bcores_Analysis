import sys
import numpy as np
import scipy as sp
from scipy import interpolate
from scipy import integrate
from scipy import stats
import copy


import diffusivity
import HL_AnalyticThea_class


class DiffusionLength():
    '''
        Performs diffusion length calculations. Analytical, semi-analytical and numerical.
    '''
    def __init__(self, P = 1., rho_surf = 330., rho_ice = 917., rho_Cr = 550., rho_CO = 804.3,\
                f0 = 1, f1 = 1, fractFact_18O = "Majoube", fractFact_D = "Merlivat", \
                fractFact_17O = "Majoube"):

        """

        Sigma class instance

        Arguments:
        ------------------
            P = 1.: Pressure in atm
            rho_surf = 330: density at z = 0 [kgrm-3]
            rho_ice = 917., Ice density [kgrm-3]
            rho_Cr = 550., Densification transition density [kgrm-3]
            rho_CO = 804., Close off density [kgrm-3]
            f0 = 1: HL scaling factor for ko
            f1 = 1: HL scaling factor for k1
            fractFact_18O: Version of fract. factor for o18. "Majoube" or "Ellehoj"
            fractFact_17O: Version of fract. factor for o17. "Majoube" or "Ellehoj"
            fractFact_D: Version of fract. factor for deuterium. "Merlivat" or "Ellehoj"
        """
        self.P = P
        self.rho_surf = rho_surf/1000 #Convert to Mgrm-3
        self.rho_ice = rho_ice/1000 #Convert to Mgrm-3
        self.rho_Cr = rho_Cr/1000 #Convert to Mgrm-3
        self.rho_CO = rho_CO/1000 #Convert to Mgrm-3
        self.f0 = f0
        self.f1 = f1
        self.fractFact_18O = fractFact_18O
        self.fractFact_17O = fractFact_17O
        self.fractFact_D = fractFact_D
        self.R = 8.314472 #Gas constant JK-1mol-1

        return

    def analytical_HL(self, rho_arr = np.array((804.3,)), T = 218.15, accum = 0.025):
        '''
            Diffusion length calculations based on analytical solution to diffusion
            equation (Fick's 2nd law).
            Assumes simple strain rate (MEANING?) and uses Herron-Langway Densification
            model.

            Arguments:
            ----------
                rho_arr:        [1d array of floats] Default np.arr of 1 float, rho_CO.
                                Array of densities to be evaluated. [kg/m^3]
                T:              [float] Temperature [K]
                accum:          [float] Acummulation rate in W.E. [m/yr]

            Returns:
            --------
                Tuple with diffusion length estimates:
                tuple[0]:       [1d array of floats] Deuterium diffusion lengths [m]
                tuple[1]:       [1d array of floats] O18 diffusion lengths [m]
                tuple[2]:       [1d array of floats] O17 diffusion lengths [m]
        '''

            # First, divide rho array into sections of rho <= rho_Cr and rho > rho_Cr
        rhos = np.reshape(rho_arr, -1)/1000 # Convert to Mg/m^3
        rhos_clip = np.clip(rhos, self.rho_surf, self.rho_CO) # set all values below rho_surf to rho_surf, all values above rho_CO to rho_CO

            # Define sigma arrays of same size as rho array
        sigma2_D = np.zeros(np.size(rhos))
        sigma2_o18 = np.zeros(np.size(rhos))
        sigma2_o17 = np.zeros(np.size(rhos))

            # Find indices corresponding to upper (rho <= rho_Cr) and lower (rho > rho_Cr) densification zones.
        try:
            idx_up = np.where(rhos_clip <= self.rho_Cr)[0]
        except:
            idx_up = np.array(())
        try:
            idx_bot = np.where(rhos_clip > self.rho_Cr)[0]
        except:
            idx_bot = np.array(())

            # Set values of used parameters
        m = 18.e-6 # Molar weight
        P_sat = diffusivity.P_Ice().clausius_clapeyron_Lt(T) # Saturation vapour pressure over ice
        k0 = self.f0 * 11. * np.exp(-10160. / (self.R * T))
        k1 = self.f1 * 575. * np.exp(-21400. / (self.R * T))
        a = 1. # HL accumulation dependency, zone 1
        b = 0.5 # HL accumulation dependency, zone 2

            ###################
            #### Deuterium ####
            ###################

            # Calculate D fractionation factor and air diffusivity as sepcified in __init__
        alpha_D = diffusivity.FractionationFactor(T).deuterium()[self.fractFact_D]
        airDiff_D = diffusivity.AirDiffusivity(T, self.P).deuterium()*3600*24*365.25 # In [m^2/yr]
        Z_D = (m * P_sat * airDiff_D) / (self.R * T * alpha_D)

            # Compute sigma_sq_D for upper densification zone:
        if np.size(idx_up) != 0:
            sigma2_D[idx_up] = Z_D / (self.rho_ice * k0 * (accum**a) * (rhos_clip[idx_up]**2)) * \
            ((rhos_clip[idx_up]**2) - (self.rho_surf**2) - 1.3 / (2 * (self.rho_ice**2)) * ((rhos_clip[idx_up]**4) - (self.rho_surf**4)))
        else:
            pass

            # Compute sigma_sq_D for bottom densification zone:
        if np.size(idx_bot) != 0:
            sigma2_D[idx_bot] = Z_D / (self.rho_ice * k1 * (accum**b) * (rhos_clip[idx_bot]**2)) * \
            ((rhos_clip[idx_bot]**2) - (self.rho_Cr**2) - 1.3 / (2 * (self.rho_ice**2)) * ((rhos_clip[idx_bot]**4) - (self.rho_Cr**4))) + \
            Z_D / (self.rho_ice * k0 * (accum**a) * (rhos_clip[idx_bot]**2)) * \
            ((self.rho_Cr**2) - (self.rho_surf**2) - 1.3 / (2 * (self.rho_ice**2)) * ((self.rho_Cr**4) - (self.rho_surf**4)))
        else:
            pass



            ###################
            ####### o18 #######
            ###################

            # Calculate o18 fractionation factor and air diffusivity as sepcified in __init__
        alpha_o18 = diffusivity.FractionationFactor(T).o18()[self.fractFact_18O]
        airDiff_o18 = diffusivity.AirDiffusivity(T, self.P).o18()*3600*24*365.25 # In [m^2/yr]
        Z_o18 = (m * P_sat * airDiff_o18) / (self.R * T * alpha_o18)

            # Compute sigma_sq_o18 for upper densification zone:
        if np.size(idx_up) != 0:
            sigma2_o18[idx_up] = Z_o18 / (self.rho_ice * k0 * (accum**a) * (rhos_clip[idx_up]**2)) * \
            ((rhos_clip[idx_up]**2) - (self.rho_surf**2) - 1.3 / (2 * (self.rho_ice**2)) * ((rhos_clip[idx_up]**4) - (self.rho_surf**4)))
        else:
            pass

            # Compute sigma_sq_D for bottom densification zone:
        if np.size(idx_bot) != 0:
            sigma2_o18[idx_bot] = Z_o18 / (self.rho_ice * k1 * (accum**b) * (rhos_clip[idx_bot]**2)) * \
            ((rhos_clip[idx_bot]**2) - (self.rho_Cr**2) - 1.3 / (2 * (self.rho_ice**2)) * ((rhos_clip[idx_bot]**4) - (self.rho_Cr**4))) + \
            Z_o18 / (self.rho_ice * k0 * (accum**a) * (rhos_clip[idx_bot]**2)) * \
            ((self.rho_Cr**2) - (self.rho_surf**2) - 1.3 / (2 * (self.rho_ice**2)) * ((self.rho_Cr**4) - (self.rho_surf**4)))
        else:
            pass



            ###################
            ####### o17 #######
            ###################

            # Calculate o17 fractionation factor and air diffusivity as sepcified in __init__
        alpha_o17 = diffusivity.FractionationFactor(T).o17()[self.fractFact_17O]
        airDiff_o17 = diffusivity.AirDiffusivity(T, self.P).o17()*3600*24*365.25 # In [m^2/yr]
        Z_o17 = (m * P_sat * airDiff_o17) / (self.R * T * alpha_o17)

            # Compute sigma_sq_o18 for upper densification zone:
        if np.size(idx_up) != 0:
            sigma2_o17[idx_up] = Z_o17 / (self.rho_ice * k0 * (accum**a) * (rhos_clip[idx_up]**2)) * \
            ((rhos_clip[idx_up]**2) - (self.rho_surf**2) - 1.3 / (2 * (self.rho_ice**2)) * ((rhos_clip[idx_up]**4) - (self.rho_surf**4)))
        else:
            pass

            # Compute sigma_sq_D for bottom densification zone:
        if np.size(idx_bot) != 0:
            sigma2_o17[idx_bot] = Z_o17 / (self.rho_ice * k1 * (accum**b) * (rhos_clip[idx_bot]**2)) * \
            ((rhos_clip[idx_bot]**2) - (self.rho_Cr**2) - 1.3 / (2 * (self.rho_ice**2)) * ((rhos_clip[idx_bot]**4) - (self.rho_Cr**4))) + \
            Z_o17 / (self.rho_ice * k0 * (accum**a) * (rhos_clip[idx_bot]**2)) * \
            ((self.rho_Cr**2) - (self.rho_surf**2) - 1.3 / (2 * (self.rho_ice**2)) * ((self.rho_Cr**4) - (self.rho_surf**4)))
        else:
            pass


            # Compute final diffusion length estimate as sqrt of sigma2 times the densities
        sigma_D = np.sqrt(sigma2_D) * rhos_clip / rhos
        sigma_o18 = np.sqrt(sigma2_o18) * rhos_clip / rhos
        sigma_o17 = np.sqrt(sigma2_o17) * rhos_clip / rhos

        return sigma_D, sigma_o18, sigma_o17

    def semi_analytical_HL(self, rho = 804.3, T = 218.15, accum = 0.025):
        '''
            Diffusion length calculations based on semi-analytical solution to diffusion
            equation.
            Assuming a simple vertical strain rate and uses Herron-Langway densification
            model.
            Semi-analytical by integrating the diffusivity from t = 0 to t = t_rho.


            Arguments:
            ----------
                rho:        [float] Density to create diffusion length up to. [kg/m^3]
                T:          [float] Temperature [K]
                accum:      [float] Accumulation rate. [m/yr]

            Returns:
            --------
                tuple[0]:   [array of floats] Deuterium diffusion lengths [m]
                tuple[1]:   [array of floats] o18 diffusion lengths [m]
                tuple[2]:   [array of floats] o17 diffusion lengths [m]
        '''

            # Define a HL density model instance
        hl_inst = HL_AnalyticThea_class.HL_Thea(Temp_0 = T, Acc_0 = accum, f0_init = self.f0, f1_init = self.f1, rho_0 = self.rho_surf*1000)

            # Compute a HL density model for a depth array
        z_model = np.arange(0,200,0.05)
        hl_model = hl_inst.model(z_model)

            # Compute HL timescale for depth density model
        age_hl = hl_inst.time_scale_HL(z_model, hl_model['rhoHL'], [self.f0, self.f1])[1]
        rho_hl = hl_model['rhoHL']*1000

            # Create interpolation instance and age scale w. one year intervals
        age_1yr = np.arange(0, age_hl[-1], 1)
        interp_inst = interpolate.interp1d(age_hl, rho_hl)

            # Make age-density interpolation model w. 1 yr btw. points
        rho_hl_1yr = interp_inst(age_1yr)

            # Use 1yr timescale for rho <= rho_input as final timescale
        age_fin = age_1yr[rho_hl_1yr <= rho]

            # Use final timescale as input for diffusivity instance
        firnDiff_inst = diffusivity.FirnDiffusivity(rho_hl_1yr, rho_co = self.rho_CO*1000, \
            T = T, P = self.P)

            ###################
            #### Deuterium ####
            ###################

            # Compute firn diffusivity
        firnDiff_D = firnDiff_inst.deuterium(f_factor_version = self.fractFact_D) * 3600 * 24 * 365.25 # Diffusivity in [m^2/yr]

            # Calculate sigma_D through numerical integration using Simpson's rule.
        sigma2_D = (1 / rho**2) * sp.integrate.simps(2 * firnDiff_D[rho_hl_1yr <= rho] * \
                    (rho_hl_1yr[rho_hl_1yr <= rho])**2, age_fin)
        sigma_D = np.sqrt(sigma2_D)

            ###################
            ####### o18 #######
            ###################

            # Compute firn diffusivity
        firnDiff_o18 = firnDiff_inst.o18(f_factor_version = self.fractFact_18O) * 3600 * 24 * 365.25 # Diffusivity in [m^2/yr]

            # Calculate sigma_o18 through numerical integration using Simpson's rule.
        sigma2_o18 = (1 / rho**2) * sp.integrate.simps(2 * firnDiff_o18[rho_hl_1yr <= rho] * \
                    (rho_hl_1yr[rho_hl_1yr <= rho])**2, age_fin)
        sigma_o18 = np.sqrt(sigma2_o18)

            ###################
            ####### o17 #######
            ###################

            # Compute firn diffusivity
        firnDiff_o17 = firnDiff_inst.o17(f_factor_version = self.fractFact_17O) * 3600 * 24 * 365.25 # Diffusivity in [m^2/yr]

            # Calculate sigma_o18 through numerical integration using Simpson's rule.
        sigma2_o17 = (1 / rho**2) * sp.integrate.simps(2 * firnDiff_o17[rho_hl_1yr <= rho] * \
                    (rho_hl_1yr[rho_hl_1yr <= rho])**2, age_fin)
        sigma_o17 = np.sqrt(sigma2_o17)

        return sigma_D, sigma_o18, sigma_o17

    def numerical_HL(self, rho = 804.3, drho = 0.1, T = 218.15, accum = 0.025):
        '''
            Diffusion length calculations using a numerical ordinary differential equation
            solver to solve the diffusion equation.
            Uses a simple vertical strain rate and a Herron Langway densification model.

            Arguments:
            ---------
                rho:        [float] Density to evaluate up to [kg/m^3]
                drho:       [float] Step size for density model.
                T:          [float] Temperature [K]
                accum:      [float] Accumulation rate [m^2/yr]

            Returns:
            --------
                tuple[0]:   [array of floats] Deuterium diffusion lengths [m]
                tuple[1]:   [array of floats] o18 diffusion lengths [m]
                tuple[2]:   [array of floats] o17 diffusion lengths [m]

        '''

            ###################
            #### Deuterium ####
            ###################

            # Define diffusion length time evolution differential equation for upper densification zone
        def dsigma2_dtD_up(sig, rhos):
                # Define diffusivity instance
            firnDiff_inst = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_CO*1000,\
                            T = T, P = self.P)
                # Set firn diffusivity coefficients
            firnDiff_D = firnDiff_inst.deuterium(f_factor_version = self.fractFact_D)*3600*24*365.25
                # Compute density time evolution
            drho_dt = 1000 * self.f0 * 11 * np.exp(-10160/(self.R * T)) * accum * (self.rho_ice - rhos/1000)
                # Define ODE to solve
            dsigma2_dt = (2./drho_dt) * firnDiff_D - (2./rhos) * sig
            return dsigma2_dt

            # Define diffusion length time evolution differential equation for lower densification zone
        def dsigma2_dtD_bot(sig, rhos, accum):
                # Define diffusivity instance
            firnDiff_inst = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_CO*1000,\
                            T = T, P = self.P)
                # Set firn diffusivity coefficients
            firnDiff_D = firnDiff_inst.deuterium(f_factor_version = self.fractFact_D)*3600*24*365.25
                # Compute density time evolution
            drho_dt = 1000 * self.f1 * 575 * np.exp(-21400/(self.R * T)) * np.sqrt(accum) * (self.rho_ice - rhos/1000)
                # Define ODE to solve
            dsigma2_dt = (2./drho_dt) * firnDiff_D - (2./rhos) * sig
            return dsigma2_dt


            ###################
            ####### o18 #######
            ###################

            # Define diffusion length time evolution differential equation for upper densification zone
        def dsigma2_dtO18_up(sig, rhos):
                # Define diffusivity instance
            firnDiff_inst = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_CO*1000,\
                            T = T, P = self.P)
                # Set firn diffusivity coefficients
            firnDiff_o18 = firnDiff_inst.o18(f_factor_version = self.fractFact_18O)*3600*24*365.25
                # Compute density time evolution
            drho_dt = 1000 * self.f0 * 11 * np.exp(-10160/(self.R * T)) * accum * (self.rho_ice - rhos/1000)
                # Define ODE to solve
            dsigma2_dt = (2./drho_dt) * firnDiff_o18 - (2./rhos) * sig
            return dsigma2_dt

            # Define diffusion length time evolution differential equation for lower densification zone
        def dsigma2_dtO18_bot(sig, rhos, accum):
                # Define diffusivity instance
            firnDiff_inst = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_CO*1000,\
                            T = T, P = self.P)
                # Set firn diffusivity coefficients
            firnDiff_o18 = firnDiff_inst.o18(f_factor_version = self.fractFact_18O)*3600*24*365.25
                # Compute density time evolution
            drho_dt = 1000 * self.f1 * 575 * np.exp(-21400/(self.R * T)) * np.sqrt(accum) * (self.rho_ice - rhos/1000)
                # Define ODE to solve
            dsigma2_dt = (2./drho_dt) * firnDiff_o18 - (2./rhos) * sig
            return dsigma2_dt




            ###################
            ####### o17 #######
            ###################

            # Define diffusion length time evolution differential equation for upper densification zone
        def dsigma2_dtO17_up(sig, rhos):
                # Define diffusivity instance
            firnDiff_inst = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_CO*1000,\
                            T = T, P = self.P)
                # Set firn diffusivity coefficients
            firnDiff_o17 = firnDiff_inst.o17(f_factor_version = self.fractFact_17O)*3600*24*365.25
                # Compute density time evolution
            drho_dt = 1000 * self.f0 * 11 * np.exp(-10160/(self.R * T)) * accum * (self.rho_ice - rhos/1000)
                # Define ODE to solve
            dsigma2_dt = (2./drho_dt) * firnDiff_o17 - (2./rhos) * sig
            return dsigma2_dt

            # Define diffusion length time evolution differential equation for lower densification zone
        def dsigma2_dtO17_bot(sig, rhos, accum):
                # Define diffusivity instance
            firnDiff_inst = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_CO*1000,\
                            T = T, P = self.P)
                # Set firn diffusivity coefficients
            firnDiff_o17 = firnDiff_inst.o17(f_factor_version = self.fractFact_17O)*3600*24*365.25
                # Compute density time evolution
            drho_dt = 1000 * self.f1 * 575 * np.exp(-21400/(self.R * T)) * np.sqrt(accum) * (self.rho_ice - rhos/1000)
                # Define ODE to solve
            dsigma2_dt = (2./drho_dt) * firnDiff_o17 - (2./rhos) * sig
            return dsigma2_dt

            # Solve numerical diffusion equations for upper and lower densification zone
        if rho <= self.rho_Cr*1000:
            rhos = np.arange(self.rho_surf*1000, rho, drho)
            sigma_surf = 0
            sigma2_D = sp.integrate.odeint(dsigma2_dtD_up, sigma_surf, rhos)
            sigma2_o18 = sp.integrate.odeint(dsigma2_dtO18_up, sigma_surf, rhos)
            sigma2_o17 = sp.integrate.odeint(dsigma2_dtO17_up, sigma_surf, rhos)
        elif rho > self.rho_Cr*1000:
                # Calculate sigma2 at critical depth (densification zone boundary)
            rhos_Cr = np.array([self.rho_surf*1000, self.rho_Cr*1000])
            sigma_surf = 0
            sigma2_D_Cr = sp.integrate.odeint(dsigma2_dtD_up, sigma_surf, rhos_Cr)[1]
            sigma2_o18_Cr = sp.integrate.odeint(dsigma2_dtO18_up, sigma_surf, rhos_Cr)[1]
            sigma2_o17_Cr = sp.integrate.odeint(dsigma2_dtO17_up, sigma_surf, rhos_Cr)[1]
                # Calculate sigma2 for lower zone, given boundary values at rho_Cr
            rhos = np.arange(self.rho_Cr*1000, rho, drho)
            sigma2_D = sp.integrate.odeint(dsigma2_dtD_bot, sigma2_D_Cr, rhos, args = (accum,))
            sigma2_o18 = sp.integrate.odeint(dsigma2_dtO18_bot, sigma2_o18_Cr, rhos, args = (accum,))
            sigma2_o17 = sp.integrate.odeint(dsigma2_dtO17_bot  , sigma2_o17_Cr, rhos, args = (accum,))

        return rhos, np.sqrt(sigma2_D), np.sqrt(sigma2_o18), np.sqrt(sigma2_o17)



def sampling_sigma(dx = 0.05):
    """
    Returns the diffusion length due to discrete sampling at a resolution dx in [cm]
    """

    samp_sigma = np.sqrt(2*np.log(np.pi/(2))/(np.pi**2)*dx**2)

    return samp_sigma

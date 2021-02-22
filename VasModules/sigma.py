"""
#################################################################################
History for sigma.py
09092014: Added O17 in analytical HL. It is the 3rd tuple element returned.
10092014: Added O17 in semi_analytical_HL. It is the 3rd tuple element returned.
16102014: Added control on parametrization of fractionation factors using
          "Majoube", "Merlivat" and "Ellehoj" versions
22102014: Sampling sigma added
29042015: In the sigma_ice function in the conversion to m2yr-1
          I corrected 365 with 365.25.
          !!!Note that the ice diffusivity calculation has been changed and now
          there is a possibility for deifferent parametrizations (Ramseier,
          Delibaltas, Blicks and Sigfus). Here as default I use Ramseier.
12072016: Added O17 calculation in experiment 2 in . Also in experiment_2 of
          SigmaToolbox the first plot created is now O18_analy and O18 num as
          a function of depth instead of density.
20082016: In function sigma_ice I added functionality for choosing between
          different ice diffusivity parametrizations. Argument is ice_diffusivity
          and can be one of strings "ramseier", "blicks", "sigfus", "delibaltas"
          "itagaki100"
          Also added a numerical calculation of diffusion length integrating on
          densities. Still in progress..!!
15122016: Change in Sigma.analytical_HL. There was an error before that resulted in
          too low values of sigma because of negative diffusivities past the close-off
          density. The modification is such thatthe input density array is clipped at
          the self.rho_co value. The output diffusion lengths are scaled as
          sigma*rho_array_clipped/rho_array
#################################################################################
"""

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, FixedFormatter
from scipy import interpolate
from scipy import integrate
from scipy import stats
import copy
#sys.path.append("/Users/vasilis/vaspy")
#sys.path.append("/home/vasileios")
import diffusivity
import herron_lang
#plt.rc('text', usetex=True)
#plt.rc('font', family='sans-serif', size = 10)
#plt.rc('axes', labelsize = 12, titlesize = 12)
#plt.rc("xtick", labelsize = 12, direction = "in")
#plt.rc("ytick", labelsize = 12, direction = "in")
#plt.ion()
#plt.close("all")



class Sigma():
    """
  Class Sigma
  Performs calculations of diffusion lengths
    """
    def __init__(self, P = 1., rho_o = 330., rho_i = 917., rho_c = 550.,\
        rho_co = 804.3, fo = 1, f1 = 1, f_factor_o18 = "Majoube", f_factor_deuterium = "Merlivat", \
        f_factor_o17 = "Majoube"):
        """

        Sigma class instance

        Arguments:
        ------------------
            P = 1.: Pressure in atm

            rho_o = 330: density at z = 0 [kgrm-3]

            rho_i = 917., Ice density [kgrm-3]

            rho_c = 550., Densification transition density [kgrm-3]

            rho_co = 804., Close off density [kgrm-3]

            fo = 1: HL scaling factor for ko

            f1 = 1: HL scaling factor for k1

            f_factor_o18: Version of fract. factor for o18. "Majoube" or "Ellehoj"

            f_factor_o17: Version of fract. factor for o17. "Majoube" or "Ellehoj"

            f_factor_deuterium: Version of fract. factor for deuterium. "Merlivat" or "Ellehoj"


        Examples
        --------
        """
        self.P = P
        self.rho_o = rho_o/1000 #Convert to Mgrm-3
        self.rho_i = rho_i/1000 #Convert to Mgrm-3
        self.rho_c = rho_c/1000 #Convert to Mgrm-3
        self.rho_co = rho_co/1000 #Convert to Mgrm-3
        self.fo = fo
        self.f1 = f1
        self.f_factor_o18 = f_factor_o18
        self.f_factor_o17 = f_factor_o17
        self.f_factor_deuterium = f_factor_deuterium
        self.R = 8.314472 #Gas constant JK-1mol-1

        return



    def analytical_HL(self, rho_array = np.array((804.3,)), T = 218.15, accum = 0.025):
        """
        Calculation of  diffusion lengths for a certain density\n
        Uses an analytical solution of the diffusion length equation (Johnsen2000 eq.1)\n
        assuming simple strain rate and applying a Herron=Langway firn model\n

        Arguments:\n
        ----------

        rho_array: Array of densitites to be evaluated [kgrm-3]\n
        T: Temperature [K]\n
        accum: Accumulation **water equivalent** [myr-1]\n

        Returns:\n
        --------
        A tuple where\n
        tuple[0]: Deuterium diffusion lengths [m]\n
        tuple[1]: O18 diffusion lengths [m]\n
        tuple[2]: O17 diffusion lengths [m]\n
        """

        rho_array = np.reshape(rho_array, -1)/1000 #Convert to Mgrm-3
        rho_array_clipped = np.clip(rho_array, self.rho_o, self.rho_co)
        sigma_sq_D = np.zeros(np.size(rho_array))
        sigma_sq_o18 = np.zeros(np.size(rho_array))
        sigma_sq_o17 = np.zeros(np.size(rho_array))


        try:
            i_up_part = np.where(rho_array_clipped <= self.rho_c)[0]
        except:
            i_up_part = np.array(())
        try:
            i_bot_part = np.where(rho_array_clipped > self.rho_c)[0]
        except:
            i_bot_part = np.array(())

        m = 18.e-6 #molar weight in Mgr
        p_sat = diffusivity.P_Ice().clausius_clapeyron_Lt(T) # p_sat in Pa
        ko = self.fo*11.*np.exp(-10160./(self.R*T))
        k1 = self.f1*575.*np.exp(-21400./(self.R*T))
        alpha = 1.#HL accumulation depandancy
        beta = 0.5 ##HL accumulation depandancy


        #### Deuterium Block  #####
        alpha_D = diffusivity.FractionationFactor(T).deuterium()[self.f_factor_deuterium]
        air_diffusivity_D = diffusivity.AirDiffusivity(T, self.P).deuterium()*3600*24*365.25 # Convert to m2yr-1
        Z_D = m*p_sat*air_diffusivity_D/(self.R*T*alpha_D)

        if np.size(i_up_part) != 0:
            sigma_sq_D[i_up_part] = Z_D/(self.rho_i*ko*accum**alpha*rho_array_clipped[i_up_part]**2)*\
                (rho_array_clipped[i_up_part]**2- self.rho_o**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_up_part]**4 - self.rho_o**4))
        else:
            pass

        if np.size(i_bot_part) != 0:
            sigma_sq_D[i_bot_part] = Z_D/(k1*accum**beta*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (rho_array_clipped[i_bot_part]**2 - self.rho_c**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_bot_part]**4\
                - self.rho_c**4)) + Z_D/(ko*accum**alpha*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (self.rho_c**2 - self.rho_o**2 - 1.3/(2*self.rho_i**2)*(self.rho_c**4 - self.rho_o**4))
        else:
            pass


        #### Oxygen 18 Block  ####

        alpha_o18 = diffusivity.FractionationFactor(T).o18()[self.f_factor_o18]
        air_diffusivity_o18 = diffusivity.AirDiffusivity(T, self.P).o18()*3600*24*365.25
        Z_o18 = m*p_sat*air_diffusivity_o18/(self.R*T*alpha_o18)

        if np.size(i_up_part) != 0:
            sigma_sq_o18[i_up_part] = Z_o18/(self.rho_i*ko*accum**alpha*rho_array_clipped[i_up_part]**2)*\
                (rho_array_clipped[i_up_part]**2- self.rho_o**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_up_part]**4 - self.rho_o**4))
        else:
            pass

        if np.size(i_bot_part) != 0:
            sigma_sq_o18[i_bot_part] = Z_o18/(k1*accum**beta*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (rho_array_clipped[i_bot_part]**2 - self.rho_c**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_bot_part]**4\
                - self.rho_c**4)) + Z_o18/(ko*accum**alpha*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (self.rho_c**2 - self.rho_o**2 - 1.3/(2*self.rho_i**2)*(self.rho_c**4 - self.rho_o**4))
        else:
            pass

        #### Oxygen 17 Block  ####

        alpha_o17 = diffusivity.FractionationFactor(T).o17()[self.f_factor_o17]
        air_diffusivity_o17 = diffusivity.AirDiffusivity(T, self.P).o17()*3600*24*365.25
        Z_o17 = m*p_sat*air_diffusivity_o17/(self.R*T*alpha_o17)

        if np.size(i_up_part) != 0:
            sigma_sq_o17[i_up_part] = Z_o17/(self.rho_i*ko*accum**alpha*rho_array_clipped[i_up_part]**2)*\
                (rho_array_clipped[i_up_part]**2- self.rho_o**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_up_part]**4 - self.rho_o**4))
        else:
            pass

        if np.size(i_bot_part) != 0:
            sigma_sq_o17[i_bot_part] = Z_o17/(k1*accum**beta*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (rho_array_clipped[i_bot_part]**2 - self.rho_c**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_bot_part]**4\
                - self.rho_c**4)) + Z_o17/(ko*accum**alpha*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (self.rho_c**2 - self.rho_o**2 - 1.3/(2*self.rho_i**2)*(self.rho_c**4 - self.rho_o**4))
        else:
            pass

        sigma_D = np.sqrt(sigma_sq_D)*rho_array_clipped/rho_array
        sigma_o18 = np.sqrt(sigma_sq_o18)*rho_array_clipped/rho_array
        sigma_o17 = np.sqrt(sigma_sq_o17)*rho_array_clipped/rho_array
        return sigma_D, sigma_o18, sigma_o17




    def semi_analytical_HL(self, rho = 804.3, T = 218.15, accum = 0.025):
        """
        Calculation of  diffusion lengths for a certain density\n
        Uses a semi analytical solution to the diffusion length equation\n
        (Johnsen2000 eq.1) assuming a simple vertical strain rate\n
        It integrates the diffusivity from t = 0 until t = t_rho and \n
        uses a Herron-Langway densification model


        Arguments:\n
        ----------

        rho: Density to be evaluated **float..!!** [kgrm-3]\n
        T: Temperature [K]\n
        accum: Accumulation **!water equivalent!** [myr-1]\n

        Returns:\n
        --------
        A tuple where\n
        tuple[0]: Deuterium diffusion length [m]\n
        tuple[1]: O18 diffusion lengths [m]\n
        tuple[2]: O17 diffusion lengths [m]\n
        """
        hl_instance = herron_lang.HL(temp = T, accu = accum, fo = self.fo, f1 = self.f1, rho_o = self.rho_o*1000)
        hl_model = hl_instance.model(np.arange(0,200, 0.05))
        age_herron = hl_instance.time_scale_hl(np.arange(0,200,0.05), hl_model["rho_herron"],\
                                               [self.fo, self.f1])[1]
        rho_herron = 1000*hl_model["rho_herron"]

        age_one_yr = np.arange(0, age_herron[-1], 1)
        interp_instance = interpolate.interp1d(age_herron, rho_herron)
        rho_herron_one_yr = interp_instance(age_one_yr)
        t_final = age_one_yr[rho_herron_one_yr <= rho]


        firn_diffusivity_instance = diffusivity.FirnDiffusivity(rho_herron_one_yr, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
        #Deuterium Block
        firn_diffusivity_D = firn_diffusivity_instance.deuterium(f_factor_version = self.f_factor_deuterium)\
        *3600*24*365.25 ##diffusivity in m2yr-1
        sigma_sq_D = 1./(rho**2)*sp.integrate.simps(2*firn_diffusivity_D[rho_herron_one_yr <= rho]*\
            (rho_herron_one_yr[rho_herron_one_yr<=rho]**2), \
            t_final)

        sigma_D_at_rho = np.sqrt(sigma_sq_D)

        #O18 Block
        firn_diffusivity_o18 = firn_diffusivity_instance.o18(f_factor_version = self.f_factor_o18)\
        *3600*24*365.25 ##diffusivity in m2yr-1
        sigma_sq_o18 = 1./(rho**2)*sp.integrate.simps(2*firn_diffusivity_o18[rho_herron_one_yr <= rho]* \
            (rho_herron_one_yr[rho_herron_one_yr<=rho]**2), \
            t_final)
        sigma_o18_at_rho = np.sqrt(sigma_sq_o18)

        #O17 Block
        firn_diffusivity_o17 = firn_diffusivity_instance.o17(f_factor_version = self.f_factor_o17)\
        *3600*24*365.25 ##diffusivity in m2yr-1
        sigma_sq_o17 = 1./(rho**2)*sp.integrate.simps(2*firn_diffusivity_o17[rho_herron_one_yr <= rho]* \
            (rho_herron_one_yr[rho_herron_one_yr<=rho]**2), \
            t_final)
        sigma_o17_at_rho = np.sqrt(sigma_sq_o17)

        return sigma_D_at_rho, sigma_o18_at_rho, sigma_o17_at_rho


    def numerical_HL(self, rho = 804.3, drho = 0.1, T = 218.15, accum = 0.025):
        """
        Calculation of  diffusion lengths for a certain density\n
        Uses a semi analytical solution to the diffusion length equation\n
        (Johnsen2000 eq.1) assuming a simple vertical strain rate\n
        It integrates the diffusivity from t = 0 until t = t_rho and \n
        uses a Herron-Langway densification model


        Arguments:\n
        ----------

        rho: Density to be evaluated **float..!!** [kgrm-3]\n
        T: Temperature [K]\n
        accum: Accumulation **!water equivalent!** [myr-1]\n

        Returns:\n
        --------
        A tuple where\n
        tuple[0]: Deuterium diffusion length [m]\n
        tuple[1]: O18 diffusion lengths [m]\n
        tuple[2]: O17 diffusion lengths [m]\n
        """

        rho = rho
        rhos = np.arange(self.rho_o*1000, rho, drho)
        # rhos = np.array((self.rho_o*1000, rho))

        sigma_o = 0.

        def dsigma2_dtD_upper(sig, rhos):
            firn_diffusivity_instance = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
            firn_diffusivity_deuterium = firn_diffusivity_instance.deuterium(f_factor_version = self.f_factor_deuterium)\
        *3600*24*365.25 ##diffusivity in m2yr-1
            drho_dt = 1000*self.fo*11*np.exp(-10160/(self.R*T))*accum*(self.rho_i - rhos/1000)
            dsigma2dt = 2./drho_dt*firn_diffusivity_deuterium - 2./rhos*sig
            return dsigma2dt


        def dsigma2_dt18_upper(sig, rhos):
            firn_diffusivity_instance = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
            firn_diffusivity_o18 = firn_diffusivity_instance.o18(f_factor_version = self.f_factor_o18)\
        *3600*24*365.25 ##diffusivity in m2yr-1
            drho_dt = 1000*self.fo*11*np.exp(-10160/(self.R*T))*accum*(self.rho_i - rhos/1000)
            dsigma2dt = 2./drho_dt*firn_diffusivity_o18 - 2./rhos*sig
            return dsigma2dt


        def dsigma2_dt17_upper(sig, rhos):
            firn_diffusivity_instance = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
            firn_diffusivity_o17 = firn_diffusivity_instance.o17(f_factor_version = self.f_factor_o17)\
        *3600*24*365.25 ##diffusivity in m2yr-1
            drho_dt = 1000*self.fo*11*np.exp(-10160/(self.R*T))*accum*(self.rho_i - rhos/1000)
            dsigma2dt = 2./drho_dt*firn_diffusivity_o17 - 2./rhos*sig
            return dsigma2dt



        def dsigma2_dtD_lower(sig, rhos, accum):
            firn_diffusivity_instance = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
            firn_diffusivity_deuterium = firn_diffusivity_instance.deuterium(f_factor_version = self.f_factor_deuterium)\
        *3600*24*365.25 ##diffusivity in m2yr-1
            # drho_dt = 1000*self.fo*11*np.exp(-10160/(self.R*T))*accum*(self.rho_i - rhos/1000)
            drho_dt = 1000*self.f1*575*np.exp(-21400/(self.R*T))*np.sqrt(accum)*(self.rho_i - rhos/1000.)
            dsigma2dt = 2./drho_dt*firn_diffusivity_deuterium - 2./rhos*sig
            return dsigma2dt


        def dsigma2_dt18_lower(sig, rhos, accum):
            firn_diffusivity_instance = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
            firn_diffusivity_o18 = firn_diffusivity_instance.o18(f_factor_version = self.f_factor_o18)\
        *3600*24*365.25 ##diffusivity in m2yr-1
            # drho_dt = 1000*self.fo*11*np.exp(-10160/(self.R*T))*accum*(self.rho_i - rhos/1000)
            drho_dt = 1000*self.f1*575*np.exp(-21400/(self.R*T))*np.sqrt(accum)*(self.rho_i - rhos/1000.)
            dsigma2dt = 2./drho_dt*firn_diffusivity_o18 - 2./rhos*sig
            return dsigma2dt


        def dsigma2_dt17_lower(sig, rhos, accum):
            firn_diffusivity_instance = diffusivity.FirnDiffusivity(rhos, rho_co = self.rho_co*1000, \
            T = T, P = self.P)
            firn_diffusivity_o17 = firn_diffusivity_instance.o17(f_factor_version = self.f_factor_o17)\
        *3600*24*365.25 ##diffusivity in m2yr-1
            # drho_dt = 1000*self.fo*11*np.exp(-10160/(self.R*T))*accum*(self.rho_i - rhos/1000)
            drho_dt = 1000*self.f1*575*np.exp(-21400/(self.R*T))*np.sqrt(accum)*(self.rho_i - rhos/1000.)
            dsigma2dt = 2./drho_dt*firn_diffusivity_o17 - 2./rhos*sig
            return dsigma2dt


        if rho<=550.:
            rhos = np.arange(self.rho_o*1000, rho, drho)
            sigma2_deuterium = sp.integrate.odeint(dsigma2_dtD_upper, 0, rhos)
            sigma2_18 = sp.integrate.odeint(dsigma2_dt18_upper, 0, rhos)
            sigma2_17 = sp.integrate.odeint(dsigma2_dt17_upper, 0, rhos)

        elif rho>550.:
            rhos = np.array((self.rho_o*1000, 550.))
            sigma2_deuterium_cr = sp.integrate.odeint(dsigma2_dtD_upper, 0., rhos)[1]
            sigma2_18_cr = sp.integrate.odeint(dsigma2_dt18_upper, 0., rhos)[1]
            sigma2_17_cr = sp.integrate.odeint(dsigma2_dt17_upper, 0., rhos)[1]

            rhos = np.arange(550., rho, drho)
            sigma2_deuterium = sp.integrate.odeint(dsigma2_dtD_lower, sigma2_deuterium_cr, rhos, args = (accum,))
            sigma2_18 = sp.integrate.odeint(dsigma2_dt18_lower, sigma2_18_cr, rhos, args = (accum,))
            sigma2_17 = sp.integrate.odeint(dsigma2_dt17_lower, sigma2_17_cr, rhos, args = (accum,))





        return rhos, np.sqrt(sigma2_deuterium), np.sqrt(sigma2_18), np.sqrt(sigma2_17)


def sigma_ice(depth, age, temp, thinning, dt, ice_diffusivity = "ramseier"):
        """
        depth [m]
        age [yr]
        temp [K]
        dt [years]

        returns sigma_ice in meters
        """
        age_dt = np.arange(age[0], age[-1], dt)
        depth_dt = np.interp(age_dt, age, depth)
        temp_dt = np.interp(age_dt, age, temp)
        thinning_dt = np.interp(age_dt, age, thinning)
        rho_ice_dt = 0.9165*(1-1.53e-4*(temp_dt - 273.15)) ##Ice density Bader 1964

        if ice_diffusivity != "ramseier":
            try:
                ice_diffusivity_dt = eval("diffusivity.IceDiffusivity(temp_dt)."\
                    + ice_diffusivity + "()*3600*24*365.25 #m2yr-1")
                print(ice_diffusivity)
            except:
                print("wrong diffusivity code")

        else:
            ice_diffusivity_dt = diffusivity.IceDiffusivity(temp_dt).ramseier()*3600*24*365.25 #m2yr-1


        sigma_ice_dt = \
            np.sqrt(thinning_dt[0:-1]**2*sp.integrate.cumtrapz(2*ice_diffusivity_dt*thinning_dt**(-2), age_dt))

        sigma_ice = np.interp(age, age_dt[0:-1], sigma_ice_dt)

        return sigma_ice



def sampling_sigma(dx = 0.05):
    """
    def sampling_sigma(dx = 0.05):
    Returns the diffusion length due to discrete sampling at a resolution dx in cm"
    """

    samp_sigma = np.sqrt(2*np.log(np.pi/(2))/(np.pi**2)*dx**2)

    return samp_sigma








class SigmaToolbox():

    def __init__(self):
        return

    def experiment1(self, fileout):
        """
        Experiment 2\n
        Generates contour plots accum-temp-delsigmasq at close off in iceeq.\n
        """
        plt.close("all")
        P = 0.7
        rho = 804.3
        fo = 1.
        f1 = 1.
        temps = np.linspace(-60, -20, 800) + 273.15
        accums_ice = np.logspace(-2, -0.3, 800) #Accumulation in ice equivalent!!

        plt.close("all")
        sigma_o18_analyt = np.array(())
        sigma_D_analyt = np.array(())
        delsigma_sq_analyt = np.array(())


        fileout_nlt = fileout + ".nlt"
        fileout_inp = fileout + ".inp"
        f_nlt = open(fileout_nlt, "w")
        f_inp = open(fileout_inp, "w")

        inputs = np.transpose(np.vstack((temps, accums_ice)))
        np.savetxt(f_inp, inputs, fmt = ("%0.3f\t%0.5f"))
        f_inp.close()


        accums = accums_ice*917./1000. #accumulation in water eq. !!!

        for accum in accums:
            print(("Accum: %0.2e" %(accum)))
            delsigma_sq_analyt = np.array(())
            sigma_o18_analyt = np.array(())
            for temp in temps:
                sigma_instance = Sigma(P, fo = fo, f1 = f1, rho_o=350.1)

                sigma_o18_analyt_i = sigma_instance.analytical_HL(rho, T = temp, accum = accum)[1]
                # sigma_o18_analyt_i = sigma_o18_analyt_i*rho/917.
                sigma_o18_analyt = np.append(sigma_o18_analyt, sigma_o18_analyt_i)

                sigma_D_analyt_i = sigma_instance.analytical_HL(rho, T = temp, accum = accum)[0]
                sigma_D_analyt_i = sigma_D_analyt_i*rho/917.
                sigma_D_analyt = np.append(sigma_D_analyt, sigma_D_analyt_i)

                delsigma_sq_analyt_i = sigma_o18_analyt_i**2 - sigma_D_analyt_i**2

                delsigma_sq_analyt = np.append(delsigma_sq_analyt, delsigma_sq_analyt_i)
                print("\t\tTemp: %0.2f" %(temp))

            np.savetxt(f_nlt, np.reshape(sigma_o18_analyt, (1,-1)), delimiter = "\t", fmt = ("%0.3f"))

        f_nlt.close()

        sigma_o18_nlt = np.loadtxt(fileout_nlt)*100


        plt.figure(1, tight_layout = True, figsize = (8,6))
        axis = plt.gca()
        axis.set_yscale('log')
        axis.set_xscale('linear')
        evalu = np.array((2,4,5,6,7,8,9,10,11,12,13, 14,16,18,22,26,30,40))
        CPf = plt.contourf(temps, accums_ice, sigma_o18_nlt, evalu, cmap = cm.GnBu)
        plt.colorbar()
        CP = plt.contour(temps, accums_ice, sigma_o18_nlt, evalu, colors = "k", linewidths = 0.5, linestyles = 'dashed')
        plt.clabel(CP, v = evalu[::], colors = "k", inline = True, fontsize = 12, fmt = "%0.0f")

        plt.plot(-55.5 + 273.15, 0.024, "kx")
        plt.text(-55.5+ 273.15 + 0.45, 0.024 - 0.003, "Vostok")

        plt.plot(-54.5+ 273.15, 0.027, "kx")
        plt.text(-54.5+ 273.15 + 0.4, 0.027 + 0.001, "DomeC")

        plt.plot(-25+ 273.15, 0.087, "kx")
        plt.text(-25+ 273.15-2, 0.087 - 0.013, "SipleDome")

        plt.plot(-51+ 273.15, 0.076, "kx")
        plt.text(-51+ 273.15 - 3.2, 0.076 - 0.01, "South Pole")

        plt.plot(-32+ 273.15, 0.207, "kx")
        plt.text(-32+ 273.15 - 4, 0.207 - 0.03, "NGRIP")

        plt.plot(-30+ 273.15, 0.2, "kx")
        plt.text(-30+ 273.15 + 0.5, 0.2 - 0.025, "NEEM")

        plt.plot(-31.4+ 273.15, 0.248, "kx")
        plt.text(-31.4+ 273.15 + 0.2, 0.248 + 0.02, "GISP2")

        plt.plot(-31.7+ 273.15, 0.23, "kx")
        plt.text(-31.7+ 273.15 - 3.5, 0.23 - 0.005, "GRIP")

        plt.plot(-31.+ 273.15, 0.23, "kx")
        plt.text(-31.+ 273.15 + .5, 0.22 - 0.005, "WAIS-D")



        majorxLocator   = MultipleLocator(5)
        majorxFormatter = FormatStrFormatter('%d')
        minorxLocator   = MultipleLocator(1)

        majoryLocator   = FixedLocator((0.01, 0.1, 1))
        majoryFormatter = FixedFormatter(['0.01', "0.1", "1"])
        # minoryLocator   = FixedLocator(np.arange(0.01, 0.09, 0.01))

        axis.xaxis.set_major_locator(majorxLocator)
        axis.xaxis.set_major_formatter(majorxFormatter)
        axis.xaxis.set_minor_locator(minorxLocator)
        axis.tick_params(axis = "both", which = 'major', labelsize = '14')

        axis.yaxis.set_major_locator(majoryLocator)
        axis.yaxis.set_major_formatter(majoryFormatter)
        # axis.yaxis.set_minor_locator(minoryLocator)
        # axis.grid(True)
        plt.xlabel('Temperature [K]', fontsize = 14)
        plt.ylabel('Accumulation ice $[\mathrm{my}^{-1}]$', fontsize = 14)
        #plt.title('P = %0.2f, fo = %0.2f, f1 = %0.2f, rho: %0.1f' %(P, fo, f1, rho))
        plt.savefig("/Users/vasilis/Desktop/diffusion_contour_1.eps")

        plt.figure(2)
        axis = plt.gca()
        axis.set_yscale('log')
        axis.set_xscale('linear')
        CP = plt.contour(temps-273.15, accums_ice, sigma_o18_nlt, 20)
        plt.clabel(CP, color = "k", inline = True, fontsize = 12, fmt = "%0.0f")
        majorxLocator   = MultipleLocator(5)
        majorxFormatter = FormatStrFormatter('%d')
        minorxLocator   = MultipleLocator(1)

        majoryLocator   = MultipleLocator(0.1)
        majoryFormatter = FormatStrFormatter('%1.2f')
        minoryLocator   = MultipleLocator(0.01)

        axis.xaxis.set_major_locator(majorxLocator)
        axis.xaxis.set_major_formatter(majorxFormatter)
        axis.xaxis.set_minor_locator(minorxLocator)

        axis.yaxis.set_major_locator(majoryLocator)
        axis.yaxis.set_major_formatter(majoryFormatter)
        # axis.yaxis.set_minor_locator(minoryLocator)
        axis.grid(True)
        plt.xlabel('Temperature [C]', fontsize = 16)
        plt.ylabel('Accumulation ice [m/yr]', fontsize = 16)
        plt.title('P = %0.2f, fo = %0.2f, f1 = %0.2f, rho: %0.1f' %(P, fo, f1, rho))

        plt.show()
        plt.savefig("/Users/vasilis/Desktop/diffusion_contour.eps")

        return


    def contour_plot_cfm_paper(self, fileout):
        """
        Experiment 2\n
        Generates contour plots accum-temp-delsigmasq at close off in iceeq.\n
        """
        plt.close("all")
        P = 0.7
        rho = 804.3
        fo = 1.
        f1 = 1.
        temps = np.linspace(-60, -20, 400) + 273.15
        accums_ice = np.logspace(-2, -0.3, 400) #Accumulation in ice equivalent!!

        plt.close("all")
        sigma_o18_analyt = np.array(())
        sigma_D_analyt = np.array(())
        delsigma_sq_analyt = np.array(())


        fileout_nlt = fileout + ".nlt"
        fileout_inp = fileout + ".inp"
        f_nlt = open(fileout_nlt, "w")
        f_inp = open(fileout_inp, "w")

        inputs = np.transpose(np.vstack((temps, accums_ice)))
        np.savetxt(f_inp, inputs, fmt = ("%0.3f\t%0.5f"))
        f_inp.close()


        accums = accums_ice*917./1000. #accumulation in water eq. !!!

        for accum in accums:
            print(("Accum: %0.2e" %(accum)))
            delsigma_sq_analyt = np.array(())
            sigma_o18_analyt = np.array(())
            for temp in temps:
                sigma_instance = Sigma(P, fo = fo, f1 = f1, rho_o=350.1)

                sigma_o18_analyt_i = sigma_instance.analytical_HL(rho, T = temp, accum = accum)[1]
                # sigma_o18_analyt_i = sigma_o18_analyt_i*rho/917.
                sigma_o18_analyt = np.append(sigma_o18_analyt, sigma_o18_analyt_i)

                sigma_D_analyt_i = sigma_instance.analytical_HL(rho, T = temp, accum = accum)[0]
                sigma_D_analyt_i = sigma_D_analyt_i*rho/917.
                sigma_D_analyt = np.append(sigma_D_analyt, sigma_D_analyt_i)

                delsigma_sq_analyt_i = sigma_o18_analyt_i**2 - sigma_D_analyt_i**2

                delsigma_sq_analyt = np.append(delsigma_sq_analyt, delsigma_sq_analyt_i)
                print("\t\tTemp: %0.2f" %(temp))

            np.savetxt(f_nlt, np.reshape(sigma_o18_analyt, (1,-1)), delimiter = "\t", fmt = ("%0.3f"))

        f_nlt.close()

        sigma_o18_nlt = np.loadtxt(fileout_nlt)*100


        plt.figure(1, tight_layout = True, figsize = (8,6))
        axis = plt.gca()
        axis.set_yscale('log')
        axis.set_xscale('linear')
        evalu = np.array((1, 2,4,5,6,7,8,9,10,11,12,13, 14,16,18,22,26,30,40, 50))
        CPf = plt.contourf(temps, accums_ice, sigma_o18_nlt, evalu, cmap = cm.GnBu)
        plt.colorbar()
        CP = plt.contour(temps, accums_ice, sigma_o18_nlt, evalu, colors = "k", linewidths = 0.5, linestyles = 'dashed')
        plt.clabel(CP, evalu[::], colors = "k", inline = True, fontsize = 12, fmt = "%0.0f")

        plt.plot(-55.5 + 273.15, 0.024, "kx")
        plt.text(-55.5+ 273.15 + 0.45, 0.024 - 0.003, "Vostok")

        plt.plot(-54.5+ 273.15, 0.027, "kx")
        plt.text(-54.5+ 273.15 + 0.4, 0.027 + 0.001, "DomeC")

        plt.plot(-25+ 273.15, 0.087, "kx")
        plt.text(-25+ 273.15-1.5, 0.087 + 0.013, "SipleDome")

        plt.plot(-51+ 273.15, 0.076, "kx")
        plt.text(-51+ 273.15 - 3.2, 0.076 - 0.01, "South Pole")

        plt.plot(-32+ 273.15, 0.207, "kx")
        plt.text(-32+ 273.15 - 3.5, 0.207 - 0.02, "NGRIP")

        plt.plot(-30+ 273.15, 0.2, "kx")
        plt.text(-30+ 273.15 + 0.5, 0.2 - 0.025, "NEEM")

        plt.plot(-31.4+ 273.15, 0.248, "kx")
        plt.text(-31.4+ 273.15 -2.8, 0.248 + 0.02, "GISP2")

        plt.plot(-31.7+ 273.15, 0.23, "kx")
        plt.text(-31.7+ 273.15 - 3., 0.23 - 0.005, "GRIP")

        plt.plot(-31.+ 273.15, 0.23, "kx")
        plt.text(-31.+ 273.15 + .5, 0.22 - 0.005, "WAIS-D")



        plt.plot(temps, np.exp(-21.492+0.0811*temps), linewidth = 1.2, color = "k", linestyle = "dashed")



        majorxLocator   = MultipleLocator(5)
        majorxFormatter = FormatStrFormatter('%d')
        minorxLocator   = MultipleLocator(1)

        majoryLocator   = FixedLocator((0.01, 0.1, 1))
        majoryFormatter = FixedFormatter(['0.01', "0.1", "1"])
        # minoryLocator   = FixedLocator(np.arange(0.01, 0.09, 0.01))

        axis.xaxis.set_major_locator(majorxLocator)
        axis.xaxis.set_major_formatter(majorxFormatter)
        axis.xaxis.set_minor_locator(minorxLocator)
        axis.tick_params(axis = "both", which = 'major', labelsize = '14')

        axis.yaxis.set_major_locator(majoryLocator)
        axis.yaxis.set_major_formatter(majoryFormatter)
        # axis.yaxis.set_minor_locator(minoryLocator)
        # axis.grid(True)
        plt.xlabel('Temperature [K]', fontsize = 14)
        plt.ylabel('Accumulation ice equivalent $[\mathrm{my}^{-1}]$', fontsize = 14)
        #plt.title('P = %0.2f, fo = %0.2f, f1 = %0.2f, rho: %0.1f' %(P, fo, f1, rho))
        plt.savefig("/Users/vasilis/Desktop/diffusion_contour_1.eps")

        plt.figure(2)
        axis = plt.gca()
        axis.set_yscale('log')
        axis.set_xscale('linear')
        CP = plt.contour(temps-273.15, accums_ice, sigma_o18_nlt, 20)
        plt.clabel(CP, colors = "k", inline = True, fontsize = 12, fmt = "%0.0f")
        majorxLocator   = MultipleLocator(5)
        majorxFormatter = FormatStrFormatter('%d')
        minorxLocator   = MultipleLocator(1)

        majoryLocator   = MultipleLocator(0.1)
        majoryFormatter = FormatStrFormatter('%1.2f')
        minoryLocator   = MultipleLocator(0.01)

        axis.xaxis.set_major_locator(majorxLocator)
        axis.xaxis.set_major_formatter(majorxFormatter)
        axis.xaxis.set_minor_locator(minorxLocator)

        axis.yaxis.set_major_locator(majoryLocator)
        axis.yaxis.set_major_formatter(majoryFormatter)
        # axis.yaxis.set_minor_locator(minoryLocator)
        axis.grid(True)
        plt.xlabel('Temperature [C]', fontsize = 16)
        plt.ylabel('Accumulation ice [m/yr]', fontsize = 16)
        plt.title('P = %0.2f, fo = %0.2f, f1 = %0.2f, rho: %0.1f' %(P, fo, f1, rho))

        plt.show()
        plt.savefig("/Users/vasilis/Desktop/diffusion_contour.eps")

        return


    def experiment2(self, P = 0.75, temp = 244.15, accum = 0.22, rho_o = 350.0, \
                    fo = 1, f1 = 1, dz = 1., z_final = 100, fileout = False, plotFigs=False):
        """
        This experiment calculates diffusion length profiles as a function of
        depth and density. A text file is saved for D and O18 and O17 using the analytical HL
        approach. Diffusion lengths in m of firn at the current density.
        """

#        accum = accum*917./1000  #water equiv.

#        plt.close("all")
        sigma_o17_num = np.array(())
        sigma_o17_analyt = np.array(())
        sigma_o18_num = np.array(())
        sigma_o18_analyt = np.array(())
        sigma_D_num = np.array(())
        sigma_D_analyt = np.array(())
        delsigma_sq_num = np.array(())
        delsigma_sq_analyt = np.array(())

        herron_model = herron_lang.HL(temp = temp, accu = accum, rho_o = rho_o, fo = fo, f1 = f1).model(np.arange(0, z_final, dz))

        rhos = 1000*herron_model["rho_herron"]
        depths = herron_model["z"]

        sigma_instance = Sigma(P = P, rho_o = rho_o, fo = fo, f1 = f1)
        #for rho in rhos:
        #    print(depths[rhos==rho][0], rho, end=' ')
        #    sigma_o17_num = np.append(sigma_o18_num, sigma_instance.semi_analytical_HL(rho = rho, T = temp, accum = accum)[2])
        #    sigma_o18_num = np.append(sigma_o18_num, sigma_instance.semi_analytical_HL(rho = rho, T = temp, accum = accum)[1])
        #    sigma_D_num = np.append(sigma_D_num, sigma_instance.semi_analytical_HL(rho = rho, T = temp, accum = accum)[0])
        #    print(sigma_o18_num[-1])
        sigma_o17_analyt = sigma_instance.analytical_HL(rhos, T = temp, accum = accum)[2]
        sigma_o18_analyt = sigma_instance.analytical_HL(rhos, T = temp, accum = accum)[1]
        sigma_D_analyt = sigma_instance.analytical_HL(rhos, T = temp, accum = accum)[0]

#        print(sigma_o18_analyt)

        if plotFigs:
            plt.figure(10)
            plt.plot(sigma_o18_num, depths, "k", linewidth = 1)
            plt.plot(sigma_o18_analyt, depths, "b", linewidth = 1)
            plt.ylim(plt.ylim()[::-1])
            plt.figure(11)
            plt.grid(True, color = "b", linestyle = "-.")
            plt.plot(sigma_D_analyt, depths, "r", linewidth = 1)
            plt.plot(sigma_o18_analyt, depths, "b", linewidth = 1)
            plt.plot(sigma_o17_analyt, depths, "m", linewidth = 1)
            plt.ylim(plt.ylim()[::-1])
            plt.figure(12)
            plt.plot(rhos, depths, "r", linewidth = 1)
            plt.ylim(plt.ylim()[::-1])

        if not fileout == False:
            f = open(fileout, "w")
            f.write("Depth\tDensity\tsigma_D\tsigma_o18\tsigma_o17\n")
            data_out = np.transpose(np.vstack((depths, rhos, sigma_D_analyt, sigma_o18_analyt, sigma_o17_analyt)))
            np.savetxt(f, data_out, delimiter = "\t", fmt = ("%0.2f", "%0.3f", "%0.6f", "%0.6f", "%0.6f"))
            f.close()
#        plt.show()
        return sigma_o18_analyt


    def sigma_herron_profile(self, depth, age, temp, accu, thinning, dt = 1, p_atm = 0.75, fo = 1, f1 = 1, view = False):
        """
        Calculates diffusion lengths for Deuterium and O18 based on a temperature
        and an accumulation history on a dt time step with and without strain
        It saves a text file sigma_herron_profile.out
        def sigma_herron_profile(depth, age, temp, accu, thinning, dt, p_atm, fo, f1):

            depth: [m]
            age: [yr]
            temp: [K]
            accu: [m water eq.]
            thinning:
            dt: [yr]
            p_atm: [Atm]
            fo:
            f1:
        """

        age1 = np.arange(age[0], age[-1], dt)
        temp1 = sp.interpolate.interp1d(age, temp)(age1)
        accu1 = sp.interpolate.interp1d(age, accu)(age1)
        depth1 = sp.interpolate.interp1d(age, depth)(age1)
        thinning1 = sp.interpolate.interp1d(age, thinning)(age1)

        age = copy.deepcopy(age1)
        depth = copy.deepcopy(depth1)
        temp = copy.deepcopy(temp1)
        accu = copy.deepcopy(accu1)
        thinning = copy.deepcopy(thinning1)

        N = np.size(age)
        results = np.zeros((N,11))
        for i in np.arange(N):
            sigma_instance = Sigma(P = p_atm)
            results_HL = sigma_instance.analytical_HL(T = temp[i], accum = accu[i])
            sigma_D, sigma_o18, sigma_o18 = results_HL
            print(age[i], temp[i], accu[i], sigma_D)
            del_sigma_sq = sigma_o18**2 - sigma_D**2
            hl_instance = herron_lang.HL(temp = temp[i], accu = accu[i], rho_co = 804.3, fo = fo, f1 = f1)
            z_co = hl_instance.get_z_close_off()
            print(("Close off Depth at: %0.2e" %(z_co)))
            thinning_at_co = thinning[depth>=z_co][0]
            thinning_at_depth_i = thinning[i]
            sigma_D_strain_on = sigma_D*804.3/917.*thinning_at_depth_i
            sigma_o18_strain_on = sigma_o18*804.3/917.*thinning_at_depth_i
            del_sigma_sq_strain_on = sigma_o18_strain_on**2 - sigma_D_strain_on**2

            results_row = np.array((age[i], depth[i], temp[i], accu[i]*1000./917., \
                thinning[i], sigma_D, sigma_o18, del_sigma_sq, sigma_D_strain_on, \
                    sigma_o18_strain_on, del_sigma_sq_strain_on))
            results[i] = results_row

        print(results)
        f = open("./sigma_herron_profile.out", "w")
        f.write("Age [yr]\tDepth [m]\tT_site [K]\taccu [ma-1 ice]\tthinning\tsigmaD_raw[m]\tsigmao18_raw [m]\tdel_sigma_sq_raw [m2]\t")
        f.write("sigma_D_strain_on [m]\tsigma_o18_strain_on [m]\tdel_sigma_sq_strain_on [m]\n")

        np.savetxt(f, results, delimiter = "\t", fmt = ("%0.8f"))
        f.close()

        if view == True:
            plt.close("all")
            plt.figure(1)
            plt.subplot(311)
            plt.plot(results[:,1], results[:,5], "r")
            plt.subplot(312)
            plt.plot(results[:,1], results[:,6], "b")
            plt.subplot(313)
            plt.plot(results[:,1], results[:,7], "k")

            plt.figure(2)
            plt.subplot(311)
            plt.plot(results[:,1], results[:,8], "r")
            plt.subplot(312)
            plt.plot(results[:,1], results[:,9], "b")
            plt.subplot(313)
            plt.plot(results[:,1], results[:,10], "k")

        return results




    def transition(self, age_i, t_st, t_tr, dt, temp_i, temp_f,
                   accum_i, accum_f):
        """
        Ceates a diffusion length history based on a transition
        in temperature and accumulation modelled by a cdf.
        t_st refers to the time before and after the transition
        t_tr is the time the transition lasts (approx +- 2sigma)
        """


        t_final = 2*t_st + t_tr
        time = np.arange(0, t_final+dt, dt)
        t_loc = t_final/2

        rv = stats.norm(loc = t_loc, scale = float(t_tr/4))

        temp = temp_i + (temp_f - temp_i)*rv.cdf(time)
        accum = accum_i + (accum_f - accum_i)*rv.cdf(time)

        sigma_instance = Sigma(P = 0.75, rho_o = 330, rho_i = 917, rho_c = 550,\
                               rho_co = 804.3, fo = 1, f1 = 1)


        sigma_18 = np.zeros(np.size(temp))
        sigma_D = np.zeros(np.size(temp))
        sigma_sq = np.zeros(np.size(temp))

        for i in np.arange(np.size(temp)):
            sigma_18[i], sigma_D[i] = sigma_instance.semi_analytical_HL(T = temp[i], accum = accum[i])


        plt.clf()
        time = age_i - time
        sigma_18 = sigma_18*917./804.3
        sigma_D = sigma_D*917./804.3
        sigma_sq = sigma_D*2 - sigma_18**2
        plt.figure(845)
        plt.subplot(311)
        plt.plot(time, sigma_18)
        plt.subplot(312)
        plt.plot(time, sigma_D)
        plt.subplot(313)
        plt.plot(time, sigma_sq)

        return time, temp, accum, sigma_18, sigma_D, sigma_sq




    def power_spectrum(self, dt = 0.025, variance = 0.04, a1 = 0.1, P0 = 0.5, sigma = 0.05, view = False):

        """
        def power_spectrum(f, dt, variance, a1, P0, sigma, view = False)
        returns results = {"f": f, "k2": k2, "P": spectrum, "Ps": ps, "Nf": noise}

        """

        f_nyq = 1/(2*dt)
        f = np.linspace(0, f_nyq, 1000)

        def model_noise_markow(f, dt, variance, a1):
            Pn = variance*dt/(np.abs(1 - a1*np.exp(-1j*2*np.pi*f*dt))**2)

            return Pn

        def model_Ps(f, P0, diff_length):
            Ps = P0*np.exp(-(2*np.pi*f*diff_length)**2)

            return Ps


        ps = model_Ps(f, P0, sigma)
        noise = model_noise_markow(f, dt, variance, a1)
        spectrum = ps + noise

        k2 = (np.pi*2*f)**2


        if view == True:
            plt.clf()

            plt.figure(432)
            plt.semilogy(f, spectrum, "b", linewidth = 4)
            plt.semilogy(f, noise, "r", linewidth = 2)
            plt.semilogy(f, ps, "k", linewidth = 2)
            plt.ylim(ymin = 1e-7*P0)



        results = {"f": f, "k2": k2, "P": spectrum, "Ps": ps, "Nf": noise}
        return results






if __name__ == '__main__':
    run = SigmaToolbox()
    run.experiment2()

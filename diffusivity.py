"""
#######################################################################################
Changes history for diffusivity.py

08092014: Implemented O17fractionation factor in method FractionationFactor.o17()
08092014: Added diffusivity of O17 in air. Function AirDiffusivity.o17()
08092014: Modified FirnDiffusivity().o18() and FirnDiffusivity().deuterium()
          so they call all the air diffusivity, fractionation factors and
          sat. vap. press. from the relevant classes of the diffusivity module
08092014: Added FirnDiffusivity().o17() method.
09092014: Added small correction in FirnDiffusivity. Diffusivity returned is zero if
          where the density array is higher than rho_co
16102014: Added parametrizations for the Ellehoj fractionation factors.
          The fractionationFactor class methods now return a dictionary as
          for example alpha_D = {"Merlivat": alpha_D_Merlivat, "Ellehoj": alpha_D_Ellehoj}
          FirnDiffusivity methods get a parameter called f_factor_version that is
          "Merlivat"/"Ellehoj" for deuterium and "Majoube"/"Ellehoj" for o18 and o17.
29042015: Added various parametrizations for ice diffusivity based on 
          Ramseier, Delibaltas, Blicks, Sigfus. Everyting is in m2sec-1..!
20082016: Added ice diffusivity parametrization for Itagaki100 m2sec-1
#######################################################################################
"""
from __future__ import division
import numpy as np
import matplotlib.pylab as pl
import os.path
import sys



class FractionationFactor():
    """
    Fractionation Factors
    def __init__(self, T = 273.15)

    self.T: Temperature in K
    """
    def __init__(self, T = 273.15):
        self.T = T
        return

    def deuterium(self):
        alpha_D_Merlivat = 0.9098*np.exp(16288./(self.T**2)) #Merlivat Nief 1967
        alpha_D_Ellehoj = np.exp(0.2133 - 203.10/self.T + 48888./self.T**2)

        alpha_D = {"Merlivat": alpha_D_Merlivat, "Ellehoj": alpha_D_Ellehoj}
        return alpha_D

    def o18(self):
        alpha_o18_Majoube = 0.9722*np.exp(11.839/self.T) #Majoube 1970
        alpha_o18_Ellehoj = np.exp(0.0831 - 49.192/self.T + 8312.5/self.T**2)

        alpha_o18 = {"Majoube": alpha_o18_Majoube, "Ellehoj": alpha_o18_Ellehoj}
        return alpha_o18

    def o17(self):
        alpha_o17_Majoube = (self.o18()["Majoube"])**0.529  #Barkan2005, Majoube for o18
        alpha_o17_Ellehoj = (self.o18()["Ellehoj"])**0.529  #Barkan2005, Ellehoj 2012 for o18

        alpha_o17 = {"Majoube": alpha_o17_Majoube, "Ellehoj": alpha_o17_Ellehoj}
        return alpha_o17


    def deuterium_liquid(self):
        alpha_D_liquid = np.exp((52.612 - 76.248*(1000./self.T) + 24.844*(1e6/self.T**2))/1000)  #Majoube1971

        return alpha_D_liquid

    def o18_liquid(self):
        alpha_o18_liquid = np.exp((-2.0667 - 0.4156*(1e3/self.T) + 1.137*(1e6/self.T**2))/1000)  #Majoube 1971

        return alpha_o18_liquid


class P_Ice():
    """
    Various evaluations for saturation vapour pressure over ice
    Look in Murphy and Koop 2005
    """
    def __init__(self):
        return

    def clausius_clapeyron_simple(self, T):
        """
        def clausius_clapeyron(self, T):
        Simple evaluation using Clausius Clapeyron equation
        with constant latent heat of sublimation
        T: temp in K
        returns pressure in Pa
        """
        p_ice = np.exp(28.9074 - 6143.7/T)

        return p_ice

    def clausius_clapeyron_Lt(self, T):
        """
        def calusius_clapeyron_Lt(self, T):
        Evaluation using temp dependance of latent heat plus a numerical
        fit to experimental data
        T: temp in K
        returns pressure in Pa
        """
        p_ice = np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)
        return p_ice

    def sigfus_2000(self, T):
        """
        def sigfus_2000(self, T):
        expression used in Johnsen2000 paper on diffusion
        T: temp in K
        returns pressure in Pa
        """
        p_ice = 3.454e12*np.exp(-6133./T)
        return p_ice


class P_Water():
    """
    Various evaluations for saturation vapour pressure over water
    """
    def __init__(self):
        return

    def goff(self, T):
        """
        Goff Gratch equation 1946 pressure in Pa T in Kelvin
        """
        T_st = 373.15 #steaming temp
        p_steam = 101325
        log10_p_water = -7.90298*(T_st/T - 1) + 5.02808*np.log10(T_st/T)\
             - 1.3816e-7*(10**(11.344*(1 - T/T_st)) - 1) + \
                 8.1328e-3*(10**(-3.49149*(T_st/T - 1)) - 1) + np.log10(p_steam)

        p_water = 10**(log10_p_water)

        return p_water


class FirnDiffusivity():

    def __init__(self, rho, rho_co = 804.3, T = 218.5, P = 1):
        self.rho = rho
        self.rho_co = rho_co
        self.T = T
        self.P = P


        return


    def deuterium(self, f_factor_version = "Merlivat"):
        """
        Return Diffusivity in firn for HDO [m2sec-1]
        """
        if type(self.rho) in [float, int]:
            self.rho = np.array((self.rho, ))
        set_zero_at = np.where(self.rho>self.rho_co)[0]
        m = 18e-3 # kgr
        R = 8.314
        rho_ice = 917. #kgrm-3
        sat_vap_pres = P_Ice().sigfus_2000(self.T)  #sat. vapor pressure Pa
        Da = AirDiffusivity(T = self.T, P = self.P).h2o() # Air diffusivity m2sec-1
        Da_D = AirDiffusivity(T = self.T, P = self.P).deuterium() #Air diffusivity for D m2sec-1
        if f_factor_version == "Merlivat":
            alpha_D = FractionationFactor(T = self.T).deuterium()["Merlivat"] #fractionation factor D Merlivat
        elif f_factor_version == "Ellehoj":
            alpha_D = FractionationFactor(T = self.T).deuterium()["Ellehoj"] #fractionation factor D Ellehoj
        if f_factor_version not in ["Merlivat", "Ellehoj"]:
            print("f_factor_version not in [Merlivat, Ellehoj]")
            sys.exit()
        tau = 1./(1 - 1.3*(self.rho/rho_ice)**2)

        Df_D = (m*sat_vap_pres*Da_D)/(R*self.T*alpha_D*tau)*\
            (1./self.rho - 1./rho_ice)


        Df_D[set_zero_at] = 0




        return Df_D


    def o18(self, f_factor_version = "Majoube"):
        """
        Return Diffusivity in firn for H2O18 [m2sec-1]
        """
        if type(self.rho) in [float, int]:
            self.rho = np.array((self.rho, ))
        set_zero_at = np.where(self.rho>self.rho_co)[0]
        m = 18e-3 #kgr
        R = 8.314
        rho_ice = 917.  #kgrm-3
        sat_vap_pres = P_Ice().sigfus_2000(self.T)  #sat. vapor pressure Pa
        Da = AirDiffusivity(T = self.T, P = self.P).h2o() # Air diffusivity m2sec-1
        Da_o18 = AirDiffusivity(T = self.T, P = self.P).o18() #Air diffusivity for o18 m2sec-1

        if f_factor_version == "Majoube":
            alpha_o18 = FractionationFactor(T = self.T).o18()["Majoube"] #fractionation factor o18 Majoube
        elif f_factor_version == "Ellehoj":
            alpha_o18 = FractionationFactor(T = self.T).o18()["Ellehoj"] #fractionation factor o18 Ellehoj
        if f_factor_version not in ["Majoube", "Ellehoj"]:
            print("f_factor_version not in [Majoube, Ellehoj]")
            sys.exit()

        tau = 1./(1 - 1.3*(self.rho/rho_ice)**2)

        Df_o18 = (m*sat_vap_pres*Da_o18)/(R*self.T*alpha_o18*tau)*\
            (1./self.rho - 1./rho_ice)

        Df_o18[set_zero_at] = 0

        return Df_o18  #m2sec-1


    def o17(self, f_factor_version = "Majoube"):
        """
        Return Diffusivity in firn for H2O17 [m2sec-1]
        """
        if type(self.rho) in [float, int]:
            self.rho = np.array((self.rho, ))
        set_zero_at = np.where(self.rho>self.rho_co)[0]
        m = 18e-3 #kgr
        R = 8.314
        rho_ice = 917.  #kgrm-3
        sat_vap_pres = P_Ice().sigfus_2000(self.T)  #sat. vapor pressure Pa
        Da = AirDiffusivity(T = self.T, P = self.P).h2o() # Air diffusivity m2sec-1
        Da_o17 = AirDiffusivity(T = self.T, P = self.P).o17() #Air diffusivity for o18 m2sec-1

        if f_factor_version == "Majoube":
            alpha_o17 = FractionationFactor(T = self.T).o17()["Majoube"] #fractionation factor o17 Majoube - Barkan
        elif f_factor_version == "Ellehoj":
            alpha_o17 = FractionationFactor(T = self.T).o17()["Ellehoj"] #fractionation factor o17 Ellehoj - Barkan
        if f_factor_version not in ["Majoube", "Ellehoj"]:
            print("f_factor_version not in [Majoube, Ellehoj]")
            sys.exit()

        tau = 1./(1 - 1.3*(self.rho/rho_ice)**2)

        Df_o17 = (m*sat_vap_pres*Da_o17)/(R*self.T*alpha_o17*tau)*\
            (1./self.rho - 1./rho_ice)
        Df_o17[set_zero_at] = 0

        return Df_o17  #m2sec-1




class IceDiffusivity():
    """

    """

    def __init__(self, T = np.array((218.15,))):
        self.T = T  ##in Kelvin
        return


    def sigfus(self):
        """
        Uses parametrization from Sigfus 2000
        """
        D_ice = 1.255e-3*np.exp(-7273./self.T)  # m2sec-1
        return D_ice

    def ramseier(self):
        """
        Uses parametrization from Ramseier 1967
        """
        D_ice = 9.2e-4*np.exp(-7186/self.T) #m2sec-1
        return D_ice

    def blicks(self):
        """
        Uses parametrization from Blicks 1966
        """
        D_ice = 2.5e-3*np.exp(-7302/self.T)
        return D_ice

    def delibaltas(self):
        """
        Uses parametrization from Delibaltas 1966
        """
        D_ice = 0.0264*np.exp(-7881/self.T)
        return D_ice

    def itagaki100(self):
        """
        Uses parametrization from Itagaki 1964
        """
        D_ice = 0.014*np.exp(-7650/self.T)
        return D_ice



class AirDiffusivity():
    """
    Air Diffusivity
    def __init__(self, T = 218.5, P = 1):
    T: Temp in K
    P: Pressure in atm
    """
    def __init__(self, T = 218.5, P = 1):
        self.T = T
        self.P = P
        self.Da = 1e-4*0.211*(self.T/273.15)**1.94*(1./self.P) #Hall and Pruppacher1976
        return

    def deuterium(self):
        """
        Return Diffusivity in air for HDO [m2sec-1]
        """
        return self.Da*0.9755 #Merlivat1978 - Watch out Johnsen2001 has an error..!



    def o18(self):
        """
        Return Diffusivity in air for H2O18 [m2sec-1]
        """
        return self.Da*0.9723 #Merlivat1978 - Watch out Johnsen2001 has an error..!


    def o17(self):
        """
        Return Diffusivity in air for H2O17 [m2sec-1]
        """

        #return (self.o18())**0.518
        return self.Da*0.98555 #using 0.518 exponent and 0.9723 from Merlivat


    def h2o(self):
        return self.Da



class RayleighExperiment():
    """
    Rayleigh Fractionation Implementations
    def __init__(self):
    """

    def __init__(self):
        return



    def simple(self, To = np.array((283.9,)), T = np.array((263.15,)), out_full = False):
        """
        This is a simple Rayleigh experiment that calculates the
        isotopic composition of precipitation rp compared to the isotopic
        composition of an arbitrary initial vapor mass rvo.
        It uses the assumption that the relative change in the amount of vapor
        is equal to the relative change of the vapor pressure as described in
        IAEA isotopes vol. 2. For the vapor pressure we use Murphy2005 eq.2
        thus assuming lattent heat of sublimation is temperature independent
        and equal to 51059 Jmole-1

        def simple(To = np.array((283.25)), T = np.array((263.15))):
        To: temperature of the initial vapor mass [Kelvin]
        T: temperature of the cloud base at the sampling station [Kelvin]
        out_full: True or False for full output

        if out_full == True it returns the dictionary:
        {"fraction": fraction, "a_18": a_18, "a_D": a_D, "d18_v": d18_v,\
                       "dD_v": dD_v, "d18_pr": d18_pr, "dD_pr": dD_pr}

        else it returns a list:
        [d18_pr, dD_pr]

        """

        p_vo = P_Water().goff(To)
        p_v = P_Ice().clausius_clapeyron_Lt(T)

        fraction = p_v/p_vo

        a_18_site = FractionationFactor(T).o18()["Majoube"]
        a_D_site = FractionationFactor(T).deuterium()["Merlivat"]

        a_18_source = FractionationFactor(To).o18_liquid()
        a_D_source = FractionationFactor(To).deuterium_liquid()

        a_18_mean = FractionationFactor((T + To)/2).o18()["Majoube"]
        a_D_mean = FractionationFactor((T + To)/2).deuterium()["Merlivat"]



        d18_v = (1./a_18_source*fraction**(a_18_mean - 1) - 1)*1000
        dD_v = (1./a_D_source*fraction**(a_D_mean - 1) - 1)*1000

        d18_pr = a_18_site/a_18_source*fraction**(a_18_mean - 1) - 1
        dD_pr = a_D_site/a_D_source*fraction**(a_D_mean - 1) - 1

        d18_pr*=1000
        dD_pr*=1000

        full_output = {"fraction": fraction, "d18_v": d18_v,\
                       "dD_v": dD_v, "d18_pr": d18_pr, "dD_pr": dD_pr}

        simple_output = [d18_pr, dD_pr]

        if out_full == True:
            return full_output
        elif out_full == False:
            return simple_output




def interface():

    pl.ion()
    pl.clf()
    path = raw_input('Give full path of density data file (col1: depth, col2:dens \'tab\' delimited): ')
    T = np.float(raw_input('Give Temperature [K]: '))
    P = np.float(raw_input('Give ambient pressure [Atm]: '))

    depth, rho = np.loadtxt(path, delimiter = '\t', unpack = True)
    rho = rho*1000
    D = Diffusivity(rho, T, P)
    D_D = D.deuterium()
    D_o18 = D.o18()

    pl.subplot(121)
    pl.plot(D_D, depth, 'r-', D_o18, depth, 'b-')
    pl.grid(True)
    ymin, ymax = pl.ylim()
    pl.ylim(ymax, ymin)
    pl.ylabel(r'Depth [$m$]')
    pl.xlabel(r'Diffusivity [$m^{2} sec^{-1}$]')



    pl.subplot(122)
    pl.plot(D_D, rho, 'r-', D_o18, rho, 'b-')
    pl.grid(True)
    ymin, ymax = pl.ylim()
    pl.ylim(ymax, ymin)
    pl.ylabel(r'Density [$kgr m^{-3}$]')
    pl.xlabel(r'Diffusivity [$m^{2} sec^{-1}$]')



    data = np.transpose(np.vstack((depth, rho, D_o18, D_D)))
    filename = os.path.splitext(path)[0] + ".out"
    f = open(filename, "w")
    f.write("Depth\tDensity\tDiffusivity_Ox\tDiffusivity_D\n")
    np.savetxt(f, data, delimiter = "\t", fmt = ("%0.2f", "%0.3f", "%0.4e", "%0.4e"))
    f.close()





if __name__ == '__main__':
    interface()
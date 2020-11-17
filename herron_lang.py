
import numpy as np
import scipy as sp
from scipy.optimize import leastsq
import matplotlib.pyplot as pl
import os.path


class HL():
    """
    methods to fit firn density measurements with a Herron and Langway
        densification model [Herron and Langway 1980]


    Class HL version 1.2
        Added attributes self.rho_i, self.rho_c, self.rho_co
    Class HL version 1.1
        Added rho_to_depth method

    Class HL version 1.0
        methods to fit firn density measurements with a Herron and Langway
        densification model [Herron and Langway 1980]
    """

    def __init__(self, z_meas = np.zeros(10), rho_meas = np.zeros(10), \
        temp = 218.15, accu = 0.027, rho_o = 330., rho_i = 917.,\
        rho_c = 550., rho_co = 804., fudge = False, fo = 1, f1 = 1):

        """
        Arguments:
        -------------
            z_meas: 1d array with depth measurements

            rho_meas: 1d array with density measurements

            temp: Temperature in Kelvin

            accu: Accumulation in meter of !**water**! equivalent per year

            rho_o = 330: Initial surface snow density initiated at 0.33

            rho_i = 917.: Ice density

            rho_c = 550.: Density of densification transition

            rho_co = 804.: Close off density

            fudge = True If True the fo, f1 parameters are optimized for the
                    best fit to the data. If False fo, f1 are fixed and given by
                    the user. Default is 1 for both

            fo = 1.: Fudge factor

            f1 = 1.: Fudge factor


        The __call__ method:
        def __call__(self, z)
        HL densification model evaluated on z

        Arguments:
        ----------
        z [m]: depth profile to evaluate HL model


        Returns
        -------
        A tuple (z, rho)

        Other methods:

        Examples
        --------

        """
        self.z_meas = z_meas
        self.rho_meas = rho_meas/1000
        self.temp = temp
        self.accu = accu
        self.rho_i = rho_i/1000 #in Mgrm-3
        self.rho_c = rho_c/1000 #in Mgrm-3
        self.rho_co = rho_co/1000 #in Mgrm-3
        self.rho_o = rho_o/1000 #in Mgrm-3
        self.fudge = fudge
        self.R = 8.314
        self.fo_init = fo
        self.f1_init = f1

        return

    def model_top(self, z, p):
        fo = p[0]
        ko = fo*11*np.exp(-10160/(self.R*self.temp))
        Zo = np.exp(self.rho_i*ko*z + np.log(self.rho_o/\
            (self.rho_i - self.rho_o)))

        return self.rho_i*Zo/(1+Zo)

    def drho_dz_top(self, z, p):
        fo = p[0]
        ko = fo*11*np.exp(-10160/(self.R*self.temp))
        Zo = np.exp(self.rho_i*ko*z + np.log(self.rho_o/\
            (self.rho_i - self.rho_o)))

        drho_dz = (self.rho_i*(1+Zo) + self.rho_i*Zo)/((1+Zo)**2)*self.rho_i*ko*Zo

        return drho_dz

    def drho_dz_bottom(self, z, p):
        f1 = p[1]
        k1 = f1*575*np.exp(-21400/(self.R*self.temp))
        Z1 = np.exp(self.rho_i*k1*(z - self.get_hc(p))/np.sqrt(self.accu) + \
            np.log(self.rho_c/(self.rho_i - self.rho_c)))

        drho_dz = (self.rho_i*k1*Z1)/(np.sqrt(self.accu))*(self.rho_i*(1+Z1) + self.rho_i*Z1)/((1+Z1)**2)

        return drho_dz


    def drho_dt(self, z = np.arange(0, 200, 1), rho = None, p = [1., 1.]):
        fo = p[0]
        f1 = p[1]
        z_hc = self.get_hc(p)

        ko = fo*11*np.exp(-10160/(self.R*self.temp))
        k1 = f1*575*np.exp(-21400/(self.R*self.temp))

        drho_dt_array = np.zeros(np.size(z))
        drho_dt_array[z<=z_hc] = fo*11*np.exp(-10160/(self.R*self.temp))*self.accu*(self.rho_i - rho[z<=z_hc])
        drho_dt_array[z>z_hc] = f1*575*np.exp(-21400/(self.R*self.temp))*np.sqrt(self.accu)*(self.rho_i - rho[z>z_hc])

        return drho_dt_array




    def residual_top(self, p, z, d):
        res = d - self.model_top(z, p)

        return res

    def fit_fo(self):
        p = [self.fo_init, self.f1_init]
        p_f = leastsq(self.residual_top, p, args = \
            (self.z_meas[self.rho_meas < self.rho_c], self.rho_meas[self.rho_meas < self.rho_c]))[0]
        fo = p_f[0]
        return fo

    def get_hc(self, p):

        if self.fudge == True:
            fo = self.fit_fo()
        else:
            fo = p[0]
        ko = fo*11*np.exp(-10160/(self.R*self.temp))

        z_hc = 1./(self.rho_i*ko)*(np.log(self.rho_c/(self.rho_i - self.rho_c))\
             - np.log(self.rho_o/(self.rho_i - self.rho_o)))
        self.z_hc = z_hc

        return z_hc


    def model_bottom(self, z, p):
        f1 = p[1]
        k1 = f1*575*np.exp(-21400/(self.R*self.temp))
        Z1 = np.exp(self.rho_i*k1*(z - self.get_hc(p))/np.sqrt(self.accu) + \
            np.log(self.rho_c/(self.rho_i - self.rho_c)))

        return self.rho_i*Z1/(1+Z1)

    def residual_bottom(self, p, z, d):
        res = d - self.model_bottom(z, p)

        return res

    def fit_f1(self):
        p = [self.fo_init, self.f1_init]
        p_f = leastsq(self.residual_bottom, p, \
            args = (self.z_meas[self.rho_meas > self.rho_c], self.rho_meas[self.rho_meas > self.rho_c]))[0]

        return p_f[1]

    def time_scale_hl(self, z, rho_herron, p = [1.,1.]):
        """
        Calculates a time scale based on an HL model
        given depth, density and fudge coefs arrays
        returns z, age_herron
        """
        age_herron = np.array((0))
        fo = p[0]
        ko = fo*11*np.exp(-10160/(self.R*self.temp))
        f1 = p[1]
        k1 = f1*575*np.exp(-21400/(self.R*self.temp))
        t_55 = 1./(ko*self.accu)*np.log((self.rho_i - self.rho_o)/(self.rho_i - self.rho_c))

        for i in np.arange(np.size(z))[1:]:
            if rho_herron[i] <= self.rho_c:
                age_herron_i = 1./(self.accu*ko)*np.log((self.rho_i - self.rho_o)/(self.rho_i - rho_herron[i]))
                age_herron = np.append(age_herron, age_herron_i)
            elif rho_herron[i]> self.rho_c:
                if rho_herron[i] >= self.rho_i:
                    age_herron = np.append(age_herron, np.max(age_herron))
                else:
                    age_herron_i = 1./(np.sqrt(self.accu)*k1)*np.log((self.rho_i - self.rho_c)/(self.rho_i - rho_herron[i])) + t_55
                    age_herron = np.append(age_herron, age_herron_i)

        return z, age_herron


    def rho_to_depth(self, rho_array, p = [1., 1.]):
        """
        Calculates densities based on an HL model
        given a density array and a fudge coefs list
        returns depth, densities
        """
        fo = p[0]
        f1 = p[1]
        ko = fo*11*np.exp(-10160/(self.R*self.temp))
        k1 = f1*575*np.exp(-21400/(self.R*self.temp))

        depth_array = np.zeros(np.size(rho_array))
        i_up_part = np.where(rho_array <= self.rho_c)[0]
        i_bot_part = np.where(rho_array > self.rho_c)[0]

        if np.size(i_up_part) != 0:
            depth_array[i_up_part] = 1./(ko*self.rho_i)*\
                (np.log(rho_array[i_up_part]/(self.rho_i - rho_array[i_up_part]) - \
                    np.log(self.rho_o/(self.rho_i - self.rho_o))))
        else:
            pass

        if np.size(i_bot_part) != 0:
            depth_array[i_bot_part] = np.sqrt(self.accu)/(self.rho_i*k1)*\
                (np.log(rho_array[i_bot_part]/(self.rho_i - rho_array[i_bot_part])) - \
                    np.log(self.rho_c/(self.rho_i - self.rho_c))) + self.get_hc(p)
        else:
            pass


        return depth_array, rho_array

    def model(self, z = np.arange(0,200,1)):
        """
        Returns:\n
        model_results = {'z': z, 'rho_herron' : rho_herron, 'fo_final' : fo_final,\n
            'f1_final' : f1_final, 'z_hc' : z_hc, 'z_close_off' : z_close_off,\n
                'drho_dz' : drho_dz}
        """

        if self.fudge == True:
            fo_final = self.fit_fo()
            f1_final = self.fit_f1()

        else:
            fo_final = self.fo_init
            f1_final = self.f1_init

        p = [fo_final, f1_final]
        z_hc = self.get_hc(p)

        rho_herron_top = self.model_top(z[z <= z_hc], p)
        rho_herron_bottom = self.model_bottom(z[z >= z_hc], p)
        rho_herron = np.hstack((rho_herron_top, rho_herron_bottom))
        z_close_off = self.rho_to_depth(np.array((self.rho_co,)), p)[0][0]
        drho_dz = np.hstack((self.drho_dz_top(z[z<=z_hc], p), \
            self.drho_dz_bottom(z[z>z_hc], p)))
        drho_dt = self.drho_dt(z, rho_herron, p)
        timescale = self.time_scale_hl(z, rho_herron, p)[1]

        model_results = {'z': z, 'rho_herron' : rho_herron, 'fo_final' : fo_final,\
            'f1_final' : f1_final, 'z_hc' : z_hc, 'z_close_off' : z_close_off,\
                'drho_dz' : drho_dz, 'drho_dt' : drho_dt, 'time_scale_hl' : timescale}

        return model_results


    def get_fo(self):
        fo_final = self.fit_fo()

        return fo_final


    def get_f1(self):
        f1_final = self.fit_f1()

        return f1_final

    def get_z_close_off(self):
        return self.model(np.arange(0, 1500, 0.1))['z_close_off']

    def __call__(self, z = np.arange(0, 200, 1)):
        """
        The __call__ method:
        def __call__(self, z = np.arange(0, 200, 1))
        HL densification model evaluated from surface to z_final with
        dz resolution

        Arguments:
        ----------
        z: depth [m] to evaluate the HL model

        Returns
        -------
        A tuple (z, rho)


        Examples
        --------

        """
        return self.model(z)['z'], self.model(z)['rho_herron']


#hl_instance = HL()
#hl_model = hl_instance.model(np.arange(0,200,0.1))
#rhoH = hl_model['rho_herron']
#print(rhoH)

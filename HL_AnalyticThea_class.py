import numpy as np
import scipy as sp
from scipy.optimize import leastsq

class HL_Thea():
    """
        Class with methods to compute a densification model through analysis
        developed by [Herron and Langway 1980].

        Also available is a fit of real density measurements to the HL model
        through the parameters f1 and f0.
    """

    """
        Methods available:

        __call__ method:

    """
    def __init__(self, z_meas = np.zeros(10), rho_meas = np.zeros(10), \
                Temp_0 = 218.15, Acc_0 = 0.027, rho_0 = 330., rho_I = 917., \
                rho_Cr = 550., rho_CO = 804., opti = False, f0_init = 1, f1_init = 1):

        """
            Initializes the class with the following arguments needed:
                z_meas:         [1d array of floats] Depth measurements in [m]
                rho_meas:       [1d array of floats] density measurements in [kg m^-3]
                Temp_0:         [float] Temperature at surface in [K]
                Acc_0:          [float] Accumulation in meter of water equiv. per year
                rho_0 = 330.:   [float] Initial surface snow density. Initiated at 0.33
                rho_I = 917.:   [float] Ice density
                rho_Cr = 550.:  [float] Density at densification transition btw. stage 1 and 2
                rho_CO = 804.:  [float] Bubble close-off density
                opti = True:    [bool] If True, f0 and f1 are optimized to give the
                                       best fit to data. If False, f0 and f1 are fixed
                                       and given by user, default is 1.
                f0_init = 1:    [float] Optimization factor
                f1_init = 1:    [float] Optimization factor
        """



        self.z_meas = z_meas
        self.rho_meas = rho_meas/1000
        self.Temp_0 = Temp_0
        self.Acc_0 = Acc_0
        self.rho_0 = rho_0/1000
        self.rho_I = rho_I/1000
        self.rho_Cr = rho_Cr/1000
        self.rho_CO = rho_CO/1000
        self.opti = opti
        self.f0_init = f0_init
        self.f1_init = f1_init
        self.R = 8.314

        return

    def model_0(self, z, p):
        f0 = p[0]
        k0 = f0 * 11 * np.exp(-10160 / (self.Temp_0 * self.R))
        Z0 = np.exp(self.rho_I * k0 * z + np.log(self.rho_0 / (self.rho_I - self.rho_0)))
        rho0 = self.rho_I * Z0 / (1 + Z0)

        return rho0, k0, Z0

    def model_1(self, z, p):
        f1 = p[1]
        k1 = f1 * 575 * np.exp(-21400 / (self.Temp_0 * self.R))
        Z1 = np.exp(self.rho_I * k1 * (z - self.get_zCr(p)) / np.sqrt(self.Acc_0) \
            + np.log(self.rho_Cr / (self.rho_I - self.rho_Cr)))
        rho1 = self.rho_I * Z1 / (1 + Z1)

        return rho1, k1, Z1

    def drho_dz_0(self, z, p):
        rho0, k0, Z0 = self.model_0(z, p)

        return self.rho_I**2 * k0 * Z0 * ((1 + 2*Z0) / (1 + Z0)**2)

    def drho_dz_1(self, z, p):
        rho1, k1, Z1 = self.model_1(z, p)

        return self.rho_I**2 * k1 * Z1 * (1 / np.sqrt(self.Acc_0)) * ((1 + 2*Z1) / (1 + Z1)**2)

    def drho_dt(self, z, rho, p):
        rho0, k0, Z0 = self.model_0(z,p)
        rho1, k1, Z1 = self.model_1(z,p)

        drho_dt_arr = np.zeros_like(z)
        drho_dt_arr[z <= self.get_zCr(p)] = k0 * self.Acc_0 * (self.rho_I - self.rho_0)
        drho_dt_arr[z > self.get_zCr(p)] = k1 * np.sqrt(self.Acc_0)

        return drho_dt_arr

    def res_0(self, p, z, d):

        res = d - self.model_0(z,p)[0]
        return res

    def res_1(self, p, z, d):
        res = d - self.model_1(z,p)[0]
        return res

    def fit_f0(self):
        p = [self.f0_init, self.f1_init]
        p_fit = leastsq(self.res_0, p, \
                args = (self.z_meas[self.rho_meas < self.rho_Cr], self.rho_meas[self.rho_meas < self.rho_Cr]))[0]
        return p_fit[0]

    def fit_f1(self):
        p = [self.f0_init, self.f1_init]

        p_fit = leastsq(self.res_1, p, \
                args = (self.z_meas[self.rho_meas > self.rho_Cr], self.rho_meas[self.rho_meas > self.rho_Cr]))[0]
        return p_fit[1]

    def get_f0(self):
        return fit_f0()

    def get_f1(self):
        return fit_f1()

    def get_zCr(self, p):
        if self.opti == True:
            f0 = self.fit_f0()
        else:
            f0 = p[0]

        k0 = f0*11*np.exp(-10160/(self.R*self.Temp_0))

        return (1 / (self.rho_I * k0)) * \
                (np.log((self.rho_Cr) / (self.rho_I - self.rho_Cr)) \
                - np.log((self.rho_0) / (self.rho_I - self.rho_0)))

    def time_scale_HL(self, z, rho_H, p):
        age_H = []
        f0 = p[0]
        f1 = p[1]
        k0 = f0 * 11 * np.exp(-10160 / (self.Temp_0 * self.R))
        k1 = f1 * 575 * np.exp(-21400 / (self.Temp_0 * self.R))

        t055 = (1 / (k0 * self.Acc_0)) * np.log((self.rho_I - self.rho_0) \
                / (self.rho_I - self.rho_Cr))

        for i in range(len(z)):
            if rho_H[i] <= self.rho_Cr:
                age_H.append((1 / (k0 * self.Acc_0)) * np.log((self.rho_I - self.rho_0) \
                        / (self.rho_I - rho_H[i])))

            elif rho_H[i] > self.rho_Cr:
                if rho_H[i] >= self.rho_I:
                    age_H.append(np.max(age_H))

                else:
                     age_H.append((1 / (k1 * np.sqrt(self.Acc_0))) * np.log((self.rho_I - self.rho_Cr) \
                     / (self.rho_I - rho_H[i])) + t055)

        return z, age_H

    def dens_to_depth(self, rho_arr, p):
        f0 = p[0]
        f1 = p[1]
        k0 = f0 * 11 * np.exp(-10160 / (self.Temp_0 * self.R))
        k1 = f1 * 575 * np.exp(-21400 / (self.Temp_0 * self.R))

        depth = np.zeros_like(rho_arr)
        idx_0 = np.where(rho_arr <= self.rho_Cr)[0]
        idx_1 = np.where(rho_arr > self.rho_Cr)[0]

        if len(idx_0) != 0:
            depth[idx_0] = (1 / (self.rho_I * k0)) * (np.log(rho_arr[idx_0] / (self.rho_I - rho_arr[idx_0]))\
            - np.log(self.rho_0 / (self.rho_I - self.rho_0)))
        else:
            pass

        if len(idx_1) != 0:
            depth[idx_1] = (np.sqrt(self.Acc_0) / (self.rho_I * k1))\
                            * (np.log(rho_arr[idx_1] / (self.rho_I - rho_arr[idx_1])) \
                            - np.log(self.rho_Cr / (self.rho_I - self.rho_Cr))) + self.get_zCr(p)
        else:
            pass

        return rho_arr, depth

    def model(self, z):
        if self.opti == True:
            f0_fin = self.fit_f0()
            f1_fin = self.fit_f1()
        else:
            f0_fin = self.f0_init
            f1_fin = self.f1_init
        p = [f0_fin, f1_fin]
        zCr = self.get_zCr(p)

        rhoH = np.hstack((self.model_0(z[z <= zCr], p)[0], self.model_1(z[z > zCr], p)[0]))
        zCO = self.dens_to_depth(np.array((self.rho_CO,)), p)[0][0]
        drho_dz = np.hstack((self.drho_dz_0(z[z <= zCr], p), self.drho_dz_1(z[z > zCr], p)))
        drho_dt = self.drho_dt(z, rhoH, p)
        age = self.time_scale_HL(z, rhoH, p)[1]

        model_results = {'z': z, 'rhoHL': rhoH, 'f0_fin': f0_fin, 'f1_fin': f1_fin, 'zCr': zCr,\
                        'zCO': zCO, 'drho_dz': drho_dz, 'drho_dt': drho_dt, 'ageHL': age}

        return model_results

    def get_zCO(self):
        return self.model(np.arange(0,1500,0.1))['zCO']

    def __call__(self, z):
        return self.model(z)['z'], self.model(z)['rhoHL']
#
# hl_instance = HL_Thea(Acc_0=0.143, Temp_0=240.65)
# hl_model = hl_instance.model(np.arange(0,200,0.1))
# #ageH = hl_model['ageHL']
# rhoH = hl_model['rhoHL']
# print(rhoH)

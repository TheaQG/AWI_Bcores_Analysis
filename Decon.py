import numpy as np
import os
import pandas as pd
from pandas import ExcelWriter
import matplotlib.pyplot as plt
import openpyxl
import matplotlib as mpl
import scipy as sp
from scipy import stats
from scipy import signal
from scipy import fft
from scipy import io
from scipy import interpolate
from scipy import optimize
from scipy import linalg
from scipy import integrate
from scipy.fft import dct
import copy
import math

class SpectralDecon():

    def __init__(self, t, y, N_min):
        self.t = t
        self.y = y
        self.y = self.y - np.mean(self.y)
        self.N_min = N_min
        self.dt = self.t[1] - self.t[0]

        return



    def __call__(self):
        return



    def dct(self):
        data = copy.deepcopy(self.y)
        depth = copy.deepcopy(self.t)

        if data.size < self.N_min:
            N = math.ceil(self.N_min/data.size) * data.size
        else:
            N = data.size

        DCT = sp.fft.dct(data, 2, n = N, norm='ortho')
        freq = np.fft.fftfreq(2*N, self.dt)[:(2*N)//2]

        return freq, DCT



    def dct_psd(self):
        f, S = self.dct()
        P = abs(S)**2

        return f, P



    def idct(self):
        return



    def SpectralFit(self, printFitParams=True, **kwargs):
        f, P = self.dct_psd()

        def calc_res(params, x, y, dt, weights):
            P0, s_eta2, s_tot2, a1 = params

            Noise = self.func_Noise(x, s_eta2, a1, dt)
            Signal = self.func_Signal(x, P0, s_tot2)

            Pmod = Noise + Signal
            res = weights*(np.log10(y) - np.log10(np.copy(Pmod)))

            return res

        def sum2_res(params, x, y, dt, weights):

            return np.sum(calc_res(params, x, y, dt, weights)**2)

        #Define the default boundaries for the different parameters.
        boundas = {}
        boundas['P0_Min'] = 1e-5
        boundas['P0_Max'] = 10000
        boundas['s_eta2_Min'] = 1e-10
        boundas['s_eta2_Max'] = 10
        boundas['a1_Min'] = 1e-7
        boundas['a1_Max'] = 0.4
        boundas['s_tot2_Min'] = 1e-7
        boundas['s_tot2_Max'] = 1

        #If user has specified bounds for params, it is passed through here.
        if list(kwargs.keys()):
            print('Setting fit param boundaries to user specifics')
            for j in list(kwargs.keys()):
                if j in list(bounds.keys()):
                    bounds[j] = kwargs[j]
                    print(f'setting {j} as {kwargs[j]}')
        elif not list(kwargs.keys()):
            print('Using default boundaries for variance and a1')

        #Initial parameter guess.
        p0 = [0.005, 0.005, 0.01, 0.1]
        #Weights
        weights = np.ones_like(f)*1.

        #Optimization routine - minimizes residuals btw. data and model.
        params_fit, fit_func_val, fit_dict = sp.optimize.fmin_l_bfgs_b(sum2_res, p0, fprime=None, args = (f, P, self.dt, weights),\
                                        approx_grad=True, bounds = [(boundas['P0_Min'], boundas['P0_Max']), (boundas['s_eta2_Min'], boundas['s_eta2_Max']), \
                                        (boundas['s_tot2_Min'], boundas['s_tot2_Max']), (boundas['a1_Min'], boundas['a1_Max'])])

        P_fit = self.func_Noise(f, params_fit[1], params_fit[3], self.dt) + self.func_Signal(f, params_fit[0], params_fit[2])

        opt_fit_dict = {"P0_fit": params_fit[0], "s_eta2_fit": params_fit[1], "s_tot2_fit": params_fit[2], "a1_fit": params_fit[3]}

        P0_fit = opt_fit_dict['P0_fit']
        s_eta2_fit = opt_fit_dict['s_eta2_fit']
        s_tot2_fit = opt_fit_dict['s_tot2_fit']
        a1_fit = opt_fit_dict['a1_fit']

        if printFitParams:
            print('Fit Parameters:\n')
            print(f'P0 = {P0_fit}')
            print(f'Var = {s_eta2_fit}')
            print(f's_eta2 = {s_eta2_fit} m')
            print(f'Diff len = {s_tot2_fit*100} cm')
            print(f'a1 = {a1_fit}')


        w_PSD, P_PSD = self.dct_psd()


        Pnoise = self.func_Noise(w_PSD, opt_fit_dict['s_eta2_fit'],opt_fit_dict['a1_fit'], self.dt)
        Psignal = self.func_Signal(w_PSD, opt_fit_dict['P0_fit'], opt_fit_dict['s_tot2_fit'])

        s_eta2_fit = opt_fit_dict['s_eta2_fit']
        s_tot2_fit = opt_fit_dict['s_tot2_fit']

        print(f'Diff. len., fit [cm]: {s_tot2_fit*100:.3f}')

        return w_PSD, P_PSD, Pnoise, Psignal, P_fit



    def Filters(self, sigma):
        w_PSD, P_PSD, Pnoise, Psignal, P_fit = self.SpectralFit()

        OptFilter = Psignal / (Pnoise + Psignal)
#        sigma = 0.05#s_eta2_fit
        M = np.exp(-(2 * np.pi * w_PSD)**2 * sigma**2 / 2)

        R = OptFilter * M**(-1)

        return w_PSD, OptFilter, M, R



    def deconvolve(self, sigma):
        data = copy.deepcopy(self.y)
        depth = copy.deepcopy(self.t)

        w_PSD, OptF, M, R = self.Filters(sigma)

        if data.size < self.N_min:
            idx = math.ceil(self.N_min/data.size)
            R_short = R[0::idx]
            w_PSD_short = w_PSD[0::idx]
        else:
            R_short = R
            w_PSD_short = w_PSD

        data_f = sp.fft.dct(data, 2, norm='ortho')
        decon_f = data_f * R_short

        data_decon = sp.fft.dct(decon_f, 3, norm='ortho')

        return depth, data_decon



    def func_Noise(self, w, s_eta2, a1, dz):
        '''
            Arguments:
            ----------
                w:          [array of floats] Frequency array.
                s_eta2:     [float] Variance of noise
                a1:         [float] AR-1 process coefficient (red noise)
                dz:         [float] Frequency spacing.

            returns:
            --------
                Noise function corresponding to the given params.
        '''
        return (s_eta2**2 * dz) / (np.abs(1 - a1 * np.exp(- 2 * np.pi * 1j * w * dz))**2)



    def func_Signal(self, w, p0, s_tot2):
        '''
            Arguments:
            ----------
                w:          [array of floats] Frequency array.
                p0:         [float] Signal amplification.
                s_tot2:     [float] Diffusion length.

            returns:
            --------
                Signal function corresponding to passed params.
        '''
        return p0 * np.exp(- (2 * np.pi * w * s_tot2)**2)

    def plotSpectrum(self):
        w_PSD, P_PSD, Pnoise, Psignal, P_fit = self.SpectralFit()

        figPSDfit, axPSDfit = plt.subplots(figsize=(10,8))
        axPSDfit.set(ylim=(min(P_PSD)-min(P_PSD)*0.8, max(P_PSD)+max(P_PSD)*0.8), xlim=(min(w_PSD), max(w_PSD)),\
                    ylabel='$\delta^{18}$O [\permil]', xlabel='Frequency')
        axPSDfit.semilogy(w_PSD, P_PSD, label='$P_{data}$')
        axPSDfit.semilogy(w_PSD, Pnoise, label='$|\hat{\eta}|^2$')
        axPSDfit.semilogy(w_PSD, Psignal, label='$P_{signal}$')
        axPSDfit.semilogy(w_PSD, P_fit, color='k', label='$P_{fit}$')
        axPSDfit.legend()
        axPSDfit.set_ylim((1e-12,10))

        return

    def plotFilters(self, sigma):
        w_PSD, OptFilter, M, R = self.Filters(sigma)

        figFilters, axFilters = plt.subplots(figsize=(8,8))
        axFilters.loglog(w_PSD, OptFilter, label='$\phi$')
        axFilters.loglog(w_PSD, M, label='$M$')
        axFilters.loglog(w_PSD, M**(-1), label='$M^{-1}$')
        axFilters.loglog(w_PSD, OptFilter * M**(-1), label='$R = \phi\cdot M^{-1}$')
        axFilters.set_ylim((10**(-2),100));
        axFilters.legend()

        return figFilters, axFilters



    def plotDecon(self, sigma):
        depth, data_decon = self.deconvolve(sigma)

        figDecon, axDecon = plt.subplots(figsize=(14,7))

        axDecon.plot(depth,data_decon, label='decon')
        axDecon.plot(self.t,self.y, label='data')
        axDecon.legend()
        return figDecon, axDecon

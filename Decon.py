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
from mpmath import *
import multiprocessing

class SpectralDecon():
    '''
        Methods available:
            __init__(self, t, y, N_min):
                    Initializes class, given three arguments, t, y, N_min.
                    Detrends y by subtracting the mean, to center around zero.
                    Defines the sampling, Delta = dt, as the difference between t1 and t0.

            __call__():
                    Call function,

            dct(self):
                    Method to compute the DCT, Discrete Cosine Transform, of the data.
                    Computes for N (minimum N_min) spectral frequencies - if y.size is
                    smaller than N_min, then N = ceil(N_min/y.size) * y.size.
                    Based on

            dct_psd(self, N):
                    From the computed DCT signal, S, this method computes the PSD
                    of the DCT by PSD = |S|^2.

            SpectralFit(self, printFitParams=True, **kwargs):
                    Given the computed PSD signal, (f,P), estimates parameters of noise and
                    signal by minimizing the residuals between the PSD and a fit made to it
                    (fit based on empirical noise+signal models).

            Filters(self, sigma):
                    Computes the spectral filters to amplify signal and minimize noise.

            deconvolve(self, sigma):
                    Computes the deconvolution of data multiplied with restoration filter.


            func_Noise(self, w, s_eta2, a1, dz):
                    Emipirical (red) noise function to minimize according to data.

            func_Signal(self, w, p0, s_tot2):
                    Empirical signal function to minimize according to data.


            plotSpectrum(self):
                    Plots the spectrum with fitted function, and visual of signal and noise.

            plotFilters(self, sigma):
                    Method to plot transfer, optimal and restoration filter.

            plotDecon(self, sigma):
                    Plot the raw and the deconvoluted data.
    '''

    def __init__(self, t, y, N_min, transType = 'DCT', resampling = 'uniform'):
        '''
            Initialize the class instance.


            Arguments:
            ----------
                t:              [array of floats] Input time series x data.
                y:              [array of floats] Input time series y data.
                N_min:          [int] Minimum number of points to use for spectral transform.

            returns:
            --------
                None

        '''
        self.t = t
        self.y = y
        self.y = self.y - np.mean(self.y) # Detrending
        self.N_min = N_min
        self.transType = transType
        self.resampling = resampling
        self.dt = self.t[1] - self.t[0] # Setting sampling size

        return



    def __call__(self):
        '''
            Call function.


            Arguments:
            ----------

            returns:
            --------
        '''
        return

    def interpData(self, DeltaInput=False, DeltaIn=0.):

        d = self.t
        x = self.y

        if DeltaInput:
            Delta = DeltaIn
        else:
            diff = np.diff(d)
            Delta = round(min(diff), 3)

        d_min = Delta * np.ceil(d[0]/Delta)
        d_max = Delta * np.floor(d[-1]/Delta)

        n = int(1 + (d_max - d_min)/Delta)

        j_arr = np.linspace(0,n,n)
        dhat0 = d_min + (j_arr - 1)*Delta

        f = interpolate.CubicSpline(d,x)

        xhat0 = f(dhat0)

        dhat = dhat0[(dhat0 >= min(d)) & (dhat0 <= max(d))]
        xhat = xhat0[(dhat0 >= min(d)) & (dhat0 <= max(d))]

        return dhat, xhat, Delta

    def Ndct(self):
        '''
            Directly computes the cosine transform for non-uniformly distributed input data.
            Uses DCT-II with orthonormal normalization.
            Frequencies are uniformly distributed and determined through FFT.

            Arguments:
            ----------
                None

            Returns:
            --------
                freq:           [array of floats] Positive spectral frequencies (transform of t data).
                NDCT:           [array of floats] Amplitude array of transformed y data.
        '''

        data = copy.deepcopy(self.y)
        depth = copy.deepcopy(self.t)

        if data.size < self.N_min:
            N = math.ceil(self.N_min/data.size) * data.size
        else:
            N = data.size

        dt = np.mean(np.diff(depth))
        freq = np.fft.fftfreq(2*N, dt)[:(2*N)//2]

        D = np.cos(np.outer(2*np.pi*freq,(depth + 1/(2*N))))

        NDCT = 2*D.dot(data)
        NDCT[0] = NDCT[0]*np.sqrt(1/(4*N))
        NDCT[1:] = NDCT[1:]*np.sqrt(1/(2*N))

        return freq, NDCT

    def INdct(self, freq, Amp):
        '''
            Directly computes the inverse cosine transform for non-uniformly distributed input data.
            Uses DCT-III with orthonormal normalization.

            Arguments:
            ----------
                None

            Returns:
            --------
                s:             [array of floats] Signal (transform of w data).
        '''
        depth = copy.deepcopy(self.t)

        N = freq.size
        Di = np.cos(np.outer((depth + 1/(2*N)), 2*np.pi*freq[1:]))

        s = Amp[0]/np.sqrt(N) + np.sqrt(2/N) * Di.dot(Amp[1:])
        return s

    def dct(self, N_in = 0):
        '''
            Computes the DCT through SciPy FFT package. Uses 2nd DCT with orthonormal normalization.
            Frequencies are computed through a regular FFT, taking only the positive frequencies, as
            this is in the definition of the DCT.
            Computes for N (minimum N_min) spectral frequencies - if y.size is
            smaller than N_min, then N = ceil(N_min/y.size) * y.size.


            Arguments:
            ----------
                None

            returns:
            --------
                freq:           [array of floats] Positive spectral frequencies, the transformed t data.
                DCT:            [array of floats] Amplitude array of the transformed y data.
        '''
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
        '''
            From the DCT/NDCT/FFT method, computes the Power Spectral Density (PSD).


            Arguments:
            ----------
                None
            returns:
            --------
                f:              [array of floats] Positive spectral frequencies, same as in self.dct.
                P:              [array of floats] Computed Power Spectral Density, only positive values.
        '''

        if self.transType == 'DCT':
            f, S = self.dct()
            P = abs(S)**2
        elif self.transType == 'NDCT':
            f, S = self.Ndct()
            P = abs(S)**2
        elif self.transType == 'FFT':
            f, P = self.fft_psd()

        return f, P



    def fft(self, N_in = 2000):
        '''
            Arguments:
            ----------
                N:              [int] Number of points(must be as power of twos)

            returns:
            --------
                w:              [array of floats] Frequencies, both negative and positive.
                A:              [array of floats] Amplitude array containing real+complex values.
        '''
        data = copy.deepcopy(self.y)
        depth = copy.deepcopy(self.t)

        if data.size < self.N_min:
            N = math.ceil(self.N_min/data.size) * data.size
        else:
            N = data.size

        A = sp.fft.fft(data, n = 2*N, norm='ortho')
        w = np.fft.fftfreq(2*N, self.dt)[:(2*N)//2]

        return w, A, N

    def fft_psd(self, N_in = 2000):
        '''
            Arguments:
            ----------
                N:              [int] Number of points(must be as power of twos).

            returns:
            --------
                w_pos:          [array of floats] Positive frequencies of the spectrum.
                P:              [array of floats] The FFT-generated PSD of the time series (real and positive)
        '''
        w, s, N = self.fft()

        w_pos = w[np.where(w>=0)]
        s_pos = s[np.where(w>=0)]

        P_pos = np.abs(s_pos)**2

        return w_pos, P_pos


    def SpectralFit(self, printFitParams=True, printDiffLen=True, printParamBounds=False,**kwargs):
        '''

            Arguments:
            ----------
                printFitParams: [bool] To print fitted parameters or not.
                **kwargs:       [tuple] Contains user specified boundaries for fit parameters

            returns:
            --------
                w_PSD:          [array of floats] Spectral frequencies.
                P_PSD:          [array of floats] Power Spectral Density.
                Pnoise:         [array of floats] Estimated noise function PSD.
                Psignal:        [array of floats] Estimated signal function, PSD.
                P_fit:          [array of floats] Estimated fit to PSD data.
                opt_fit_dict:   [tuple] Dictionary containing estimated fit parameters
                params_fit:     [array of floats] Estimated fit parameters, in array, from scipy.optimize.
                fit_func_val:   [array of floats] Estimated fit, from scipy.optimize.
                fit_dict:       [tuple] Dictionary from scipy.optimize.
        '''
        f, P = self.dct_psd()


        def calc_res(params, x, y, dt, weights):
            '''
                Calculates the log residuals between data, y, and model estimated from
                x and a given set of fit parameters.


                Arguments:
                ----------
                    params:     [array of floats] Parameters to compute the model estimate from.
                    x:          [array of floats] Data, x values.
                    y:          [array of floats] Data, y values.
                    dt:         [float] Spacing of x data.
                    weights:    [array of floats] Weights to fudge residuals.

                returns:
                --------
                    res:        [array of floats] Residuals of data vs model.
            '''
            P0, s_eta2, s_tot2, a1 = params

            Noise = self.func_Noise(x, s_eta2, a1, dt)
            Signal = self.func_Signal(x, P0, s_tot2)

            Pmod = Noise + Signal
            res = weights*(np.log10(y) - np.log10(np.copy(Pmod)))

            return res

        def sum2_res(params, x, y, dt, weights):
            '''
                Calculates the squared sum of residuals.

                Arguments:
                ----------
                    params:     [array of floats] Parameters to compute the model estimate from.
                    x:          [array of floats] Data, x values.
                    y:          [array of floats] Data, y values.
                    dt:         [float] Spacing of x data.
                    weights:    [array of floats] Weights to fudge residuals.

                returns:
                --------
                    sum2_res:   [float] Value of the sum of the squarred residuals.
                                        (We seek to minimize this).
            '''
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
                    if printParamBounds:
                        print(f'setting {j} as {kwargs[j]}')
        elif not list(kwargs.keys()):
            if printParamBounds:
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

        if printDiffLen:
            print(f'Diff. len., fit [cm]: {s_tot2_fit*100:.3f}')

        return w_PSD, P_PSD, Pnoise, Psignal, P_fit, opt_fit_dict, params_fit, fit_func_val, fit_dict



    def Filters(self, sigma):
        '''
            Computes spectral filters. Computes optimalfilter, OptFilter = Psignal / (Psignal + Pnoise),
            exponential transfer function, M, and restoration filter, R = OptFilter * M^(-1).


            Arguments:
            ----------
                sigma:              [float] Theoretical diffusion length to be used in transfer function.

            returns:
            --------
                w_PSD:              [array of floats] Spectral frequencies
                OptFilter:          [array of floats] Optimal filter, as a function of frequencies.
                M:                  [array of floats] Transfer function, filtering due to diffusion.
                R:                  [array of floats] Total restoration filter.

        '''
        w_PSD, P_PSD, Pnoise, Psignal, P_fit, _, _ , _, _ = self.SpectralFit(printFitParams=False, printDiffLen=False, printParamBounds=False)

        OptFilter = Psignal / (Pnoise + Psignal)

        M = np.exp(-(2 * np.pi * w_PSD)**2 * sigma**2 / 2)

        R = OptFilter * M**(-1)

        return w_PSD, OptFilter, M, R



    def Filters2(self, sigma_arr, z_arr):
        '''
            !!!!!!!!!!!! NOT WORKING !!!!!!!!!!!!

            Computes spectral filters. Computes optimalfilter, OptFilter = Psignal / (Psignal + Pnoise),
            exponential transfer function, M, and restoration filter, R = OptFilter * M^(-1).


            Arguments:
            ----------
                sigma:              [float] Theoretical diffusion length to be used in transfer function.

            returns:
            --------
                w_PSD:              [array of floats] Spectral frequencies
                OptFilter:          [array of floats] Optimal filter, as a function of frequencies.
                M:                  [array of floats] Transfer function, filtering due to diffusion.
                R:                  [array of floats] Total restoration filter.

        '''
        w_PSD, P_PSD, Pnoise, Psignal, P_fit, _, _ , _, _ = self.SpectralFit(printFitParams=False, printDiffLen=False, printParamBounds=False)

        OptFilter = Psignal / (Pnoise + Psignal)


        x_arr = z_arr
        y_arr = sigma_arr

        # Define objective function
        def fct(x,a,b,c):
            return a*(x**2) + b*x + c

        # Fit curve to data
        from scipy.optimize import curve_fit
        popt, test = curve_fit(fct, x_arr, y_arr)
        a, b, c = popt

        Npoints = len(OptFilter)
        delta = (max(x_arr) - min(x_arr))/Npoints

        xInt = np.arange(min(x_arr), max(x_arr), delta)

        yInt = fct(xInt, a, b, c)

        sigma = yInt
        zInt = xInt
        #m = np.exp(-(depthDiff**2/(2*sig_arr**2)))
        #mNorm = 1/sum(m) * m
        #M = sp.fft.dct(mNorm, 2, norm='ortho')

        expo = -zInt**2/(2*sigma**2)


        mp.dps = 10
        test_mpf = [exp(x) for x in expo]

        #test_mpfN = [(x - minV)/(maxV - minV) for x in test_mpf]


        #sumtest = 1/(10*np.pi)*sum(test_mpf)#/(2*np.pi)
        minT = min(test_mpf)
        maxT = max(test_mpf)
        test_mpfN = [(x - minT)/(maxT - minT) for x in test_mpf]

        m = test_mpf
        M = sp.fft.dct(test_mpfN, 2, norm='ortho')#[1:]
#        sumM = sum(M)
        Mfin = (M - min(M))/(max(M) - min(M))
#        Mfin = np.exp(-(2*np.pi*w_PSD)**2 * sigma**2 / 2)
        #M = np.exp(-(2 * np.pi * w_PSD)**2 * sigma**2 / 2)

        R = OptFilter * (Mfin)**(-1)

        return w_PSD, OptFilter, Mfin, R, m

    def deconvolve(self, sigma, invTransType='DCT'):#, depthDiff = np.array([])):
        '''
            Deconvolution of the restored spectral data, DCT(data) * R.
            Takes in to account that data and R are of different lengths, so R is discretized
            to a lower resolution to be able to be multiplied with the data.

            Arguments:
            ----------
                sigma:              [float] Theoretical diffusion length to be used in transfer function.

            returns:
            --------
                depth:              [array of floats] Original x data of time series.
                data_decon:         [array of floats] Deconvolution of data multiplied with restoration filter.
        '''
        data = copy.deepcopy(self.y)
        depth = copy.deepcopy(self.t)

        #if hasattr(sigma, "__len__"):
        #    w_PSD, OptF, M, R, m = self.Filters2(sigma, depthDiff)
        #else:
        sigma_use=sigma
        w_PSD, OptF, M, R = self.Filters(sigma_use)



        if data.size < self.N_min:
            idx = math.ceil(self.N_min/data.size)
            R_short = R[0::idx]
            w_PSD_short = w_PSD[0::idx]

        else:
            R_short = R
            w_PSD_short = w_PSD



        if self.transType == 'DCT':
            data_f = sp.fft.dct(data, 2, norm = 'ortho')
            decon_f = data_f * R_short
            data_decon = sp.fft.dct(decon_f, 3, norm='ortho')
            depth_decon = depth

        elif self.transType == 'NDCT':
            if self.resampling == 'uniform':
                data_Nuni = data
                depth_Nuni = depth
                N = len(R_short)+1
                delta = (max(depth_Nuni) - min(depth_Nuni))/N#np.mean(np.diff(depth_Nuni))

                depth_uni, data_uni, Delta = self.interpData(DeltaInput=True, DeltaIn = delta)
                data_uni_f = sp.fft.dct(data_uni, 2, norm='ortho')
                decon_f = data_uni_f  * R_short
                data_decon = sp.fft.dct(decon_f, 3, norm='ortho')
                depth_decon = depth_uni
            else:
                w_f, data_f = self.Ndct()
                decon_f = data_f * R
                data_decon = self.INdct(w_PSD, decon_f)
                depth_decon = depth

        elif self.transType == 'FFT':
            data_f = np.real(sp.fft.fft(data, n = 2*len(data), norm='ortho')[:len(data)]*2)
            decon_f = data_f * R_short
            data_decon = np.real(sp.fft.ifft(decon_f, n=2*len(data), norm='ortho')[:len(data)]*2)#[:Ns]
            depth_decon = depth

        return depth_decon, data_decon

    def deconvolve2(self, z_arr, sigma_arr):
        '''
            Deconvolution when taking an entire array of diffusion lengths (diff. lens. at the given depth).
            Deconvolution of the restored spectral data, DCT(data) * R.
            Takes in to account that data and R are of different lengths, so R is discretized
            to a lower resolution to be able to be multiplied with the data.

            Arguments:
            ----------
                z_arr:                  [array of floats] Depth array from section defining what diffusion lengths to use.
                sigma_arr:              [array of floats] Theoretical diffusion length array to be used in transfer function.

            returns:
            --------
                depth:                  [array of floats] Original x data of time series.
                data_decon:             [array of floats] Deconvolution of data multiplied with restoration filter.
        '''
        data = copy.deepcopy(self.y)
        depth = copy.deepcopy(self.t)
        N = self.N_min

        if len(sigma_arr) != len(self.y):
            x_arr = z_arr
            y_arr = sigma_arr

            # Define objective function
            def fct(x,a,b,c):
                return a*(x**2) + b*x + c

            # Fit curve to data
            from scipy.optimize import curve_fit
            popt, test = curve_fit(fct, x_arr, y_arr)
            # Define fit parameters
            a, b, c = popt

            # Define new x values
            Npoints = len(data)
            delta = (max(x_arr) - min(x_arr))/Npoints
            xInt = np.arange(min(x_arr), max(x_arr), delta)

            # Calculate new y values
            yInt = fct(xInt, a, b, c)

            sigmas = yInt
            zInt = xInt
        elif len(sigma_arr) == len(self.y):
            sigmas = sigma_arr
            zInt = z_arr

        # Define a multiprocessing pool-instance with 4 threads (N_CPU = 4)
        a_pool = multiprocessing.Pool(4)

        # Calculate the new deconvoluted individual points of the data by folding N = len(data) times with
        # N = len(data) transfer functions and pull out one individual data point each time.
        # Run it on multiple processors to speed it up.

        result = a_pool.starmap(self.funcDecon, zip(range(len(data)), sigmas))

        # Save the results to numpy arrays
        r = np.array(result,dtype=np.dtype('float,float,float,float'))
        depthBD = r['f0']
        dataBD = r['f1']
        sigs = r['f3']


        return depthBD, dataBD

    def funcDecon(self, idx, sigma):
        sigma_in = sigma
        depth, data = self.deconvolve(sigma=sigma_in)

        return depth[idx], data[idx], idx, sigma_in

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

    def plotSpectrum(self, showFig=True):
        '''

            Arguments:
            ----------
                showFig:            [bool] Show figure or not. Default: True.

            returns:
            --------
                figPSDfit:          [mpl.figure.Figure] Matplotlib figure of PSD with signal, noise and signal+noise fit.
                axPSDfit:           [mpl.axes._subplots.AxesSubplot] MPL axes handle.
        '''
        w_PSD, P_PSD, Pnoise, Psignal, P_fit, _, _, _, _ = self.SpectralFit(printFitParams=False)

        figPSDfit, axPSDfit = plt.subplots(figsize=(10,8))
        axPSDfit.grid(linestyle='--',lw=1.3, which='both')
        axPSDfit.set(ylim=(min(P_PSD)-min(P_PSD)*0.8, max(P_PSD)+max(P_PSD)*0.8), xlim=(min(w_PSD), max(w_PSD)),\
                    ylabel='$\delta^{18}$O [\permil]', xlabel='Frequency')
        axPSDfit.semilogy(w_PSD, P_PSD, label='$P_{data}$')
        axPSDfit.semilogy(w_PSD, Pnoise, label='$|\hat{\eta}|^2$')
        axPSDfit.semilogy(w_PSD, Psignal, label='$P_{signal}$')
        axPSDfit.semilogy(w_PSD, P_fit, color='k', label='$P_{fit}$')
        axPSDfit.legend()
        axPSDfit.set_ylim((1e-12,10))
        if not showFig:
            plt.close()

        return figPSDfit, axPSDfit

    def plotFilters(self, sigma, showFig=True):
        '''

            Arguments:
            ----------
                sigma:              [float] Theoretical diffusion length to be used in transfer function.
                showFig:            [bool] Show figure or not. Default: True.

            returns:
            --------
                figFilters:          [mpl.figure.Figure] Matplotlib figure of OptFilter, M, M^(-1) and R filters.
                axFilters:           [mpl.axes._subplots.AxesSubplot] MPL axes handle.
        '''
        w_PSD, OptFilter, M, R = self.Filters(sigma)

        figFilters, axFilters = plt.subplots(figsize=(8,8))
        axFilters.loglog(w_PSD, OptFilter, label='$\phi$')
        axFilters.loglog(w_PSD, M, label='$M$')
        axFilters.loglog(w_PSD, M**(-1), label='$M^{-1}$')
        axFilters.loglog(w_PSD, OptFilter * M**(-1), label='$R = \phi\cdot M^{-1}$')
        axFilters.set_ylim((10**(-2),100));
        axFilters.legend()

        if not showFig:
            plt.close()

        return figFilters, axFilters



    def plotDecon(self, sigma, showFig=True):
        '''

            Arguments:
            ----------
                sigma:              [float] Theoretical diffusion length to be used in transfer function.
                showFig:            [bool] Show figure or not. Default: True.

            returns:
            --------
                figDecon:          [mpl.figure.Figure] Matplotlib figure of deconvoluted and raw data.
                axDecon:           [mpl.axes._subplots.AxesSubplot] MPL axes handle.
        '''
        depth, data_decon = self.deconvolve(sigma)

        figDecon, axDecon = plt.subplots(figsize=(14,7))

        axDecon.plot(depth,data_decon, label='decon')
        axDecon.plot(self.t,self.y, label='data')
        axDecon.legend()
        if not showFig:
            plt.close()

        return figDecon, axDecon

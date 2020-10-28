import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from MEM_class import MEM
import copy

class SignalToF():
    '''

    '''
    '''
        Methods available:
            __init__():
                    Initializes class, given three arguments, t, y, psdType.
                    Detrends y by subtracting the mean, to center around zero.
                    Gives an error if len(y) is not an even number.

            __call__():
                    Call function,

            fft(self, N):
                    Method to compute the fast fourier transform of the time
                    series passed. Computes it for N points, where N must be
                    a power of two (FFT constraint).

            fft_psd(self, N):
                    From the computed FFT signal, this method computes the PSD
                    of the FFT by conversion to real and positive values.

            mem(self, M, N):
                    On the basis of the MEM_class.py class, method computes the
                    PSD through Burg's MEM. Based on M nodes and N points.

            OptFilterFit(self, **kwargs):
                    Given the psdType and from here the computed PSD signal, (f,P),
                    estimates the optimal filter(to amplify signal and kill noise)
                    by minimizing the residuals between the PSD and a fit made to it
                    (fit based on empirical noise+signal models).

            func_Noise(self, w, s_eta2, a1, dz):
                    Emipirical (red) noise function to minimize according to data.

            func_Signal(self, w, p0, s_tot2):
                    Empirical signal function to minimize according to data.

            convolve(self):
                    Deconvolution of data given the found optimal filter, transfer
                    function and data.

            plotFilters(self):
                    Method to plot transfer, optimal and restoration filter.

    '''
    def __init__(self, t, y, psdType):
        '''
            Arguments:
            ----------
                t:              [array of floats] Input time series x data.
                y:              [array of floats] Imput time series y data.
                psdType:        [str] {'FFT', 'MEM'}. What spectral transfor-
                                mation method to use.

            returns:
            --------
                None
        '''

        self.t = t
        self.y = y
        self.y = self.y - np.mean(self.y)
        self.psdType = psdType

        if np.size(self.y)%2 != 0:
            print('ERROR: Give signal with even number of points.')

        return

    def __call__(self):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        return

    def fft(self, N = 8192):
        '''
            Arguments:
            ----------
                N:              [int] Number of points(must be as power of twos)

            returns:
            --------
                w:              [array of floats] Frequencies, both negative and positive.
                A:              [array of floats] Amplitude array containing real+complex values.
        '''
        dt = self.t[1] - self.t[0]

        w = np.fft.fftfreq(N, dt)
        A = np.fft.fft(self.y, n = N) * dt

        return w, A

    def fft_psd(self, N = 8192):
        '''
            Arguments:
            ----------
                N:              [int] Number of points(must be as power of twos).

            returns:
            --------
                w_pos:          [array of floats] Positive frequencies of the spectrum.
                P:              [array of floats] The FFT-generated PSD of the time series (real and positive)
        '''
        w, s = self.fft(N)
        aS2 = np.abs(s)**2
        P = (aS2[np.where(w >= 0)] + aS2[np.where(w < 0)][::-1]) / N**2
        w_pos = w[np.where(w >= 0)]

        return w_pos, P

    def mem(self, M, N):
        '''
            Arguments:
            ----------
                M:              [int] Number of nodes to compute PSD via Burg's MEM.
                N:              [int] Number of points to generate PSD from.

            returns:
            --------
                w:              [array of floats] (Positive) frequencies.
                P:              [array of floats] (Real and positive) values of the PSD.
        '''

        MEM_inst = MEM(t_data = self.t, y_data = self.y, M = M)
        w, P = MEM_inst(t_data = self.t, y_data = self.y, M = M, N = N, view = False)

        return w, P

    def OptFilterFit(self, **kwargs):
        '''
            Arguments:
            ----------
                **kwargs:       [tuple] Contains user specified boundaries for fit parameters

            returns:
            --------
                opt_fit_dict:   [tuple] Dictionary containing estimated fit parameters
                P_fit:          [array of floats] Estimated fit to PSD data.
                params_fit:     [array of floats] Estimated fit parameters, in array, from scipy.optimize.
                fit_func_val:   [array of floats] Estimated fit, from scipy.optimize.
                fit_dict:       [tuple] Dictionary from scipy.optimize.
        '''

        # Calculate the PSD throuhg either FFT or MEM
        if self.psdType == 'FFT':
            f, P = self.fft_psd(N=8192)#fft_psd()#
        elif self.psdType == 'MEM':
            f, P = self.mem(M=100, N=8192)

        dt = self.t[1] - self.t[0]


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

            if self.psdType == 'FFT':
                res = weights*(np.log10(y) - np.log10(np.copy(Pmod)))
            elif self.psdType == 'MEM':
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
        boundas['P0_Min'] = 1e-15
        boundas['P0_Max'] = 10
        boundas['s_eta2_Min'] = 1e-10
        boundas['s_eta2_Max'] = 0.1
        boundas['a1_Min'] = 1e-7
        boundas['a1_Max'] = 0.3
        boundas['s_tot2_Min'] = 1e-7
        boundas['s_tot2_Max'] = 0.5

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
        params_fit, fit_func_val, fit_dict = sp.optimize.fmin_l_bfgs_b(sum2_res, p0, fprime=None, args = (f, P, dt, weights),\
                                        approx_grad=True, bounds = [(boundas['P0_Min'], boundas['P0_Max']), (boundas['s_eta2_Min'], boundas['s_eta2_Max']), \
                                        (boundas['s_tot2_Min'], boundas['s_tot2_Max']), (boundas['a1_Min'], boundas['a1_Max'])])

        P_fit = self.func_Noise(f, params_fit[1], params_fit[3], dt) + self.func_Signal(f, params_fit[0], params_fit[2])

        self.P_fit = P_fit
        self.params_fit = params_fit
        self.fit_func_val = fit_func_val
        self.fit_dict = fit_dict

        opt_fit_dict = {"P0_fit": params_fit[0], "s_eta2_fit": params_fit[1], "s_tot2_fit": params_fit[2], "a1_fit": params_fit[3]}

        return opt_fit_dict, P_fit, params_fit, fit_func_val, fit_dict



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



    def convolve(self):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''

        if self.psdType == 'FFT':
            w, P = self.fft_psd(N=8192)
        elif self.psdType == 'MEM':
            w, P = self.mem(M=100, N=8192)

        dt = self.t[1] - self.t[0]


        params_fit = self.OptFilterFit()[2]

        p0 = params_fit[0]
        s_eta2 = params_fit[1]
        s_tot2 = params_fit[2]
        a1 = params_fit[3]


        wConv = np.linspace(min(w), max(w), np.size(self.y))

        F_Signal = self.func_Signal(wConv, p0, s_tot2)
        F_Noise = self.func_Noise(wConv, s_eta2, a1, dt)

        F = F_Signal / (F_Signal + F_Noise)
        T = np.exp(- ((2 * np.pi * wConv * s_eta2)**2)/2)
        T = ((np.sin(np.pi * wConv * dt)) / (np.pi * wConv * dt)) * T
        R = F * T**(-1)

        s = copy.deepcopy(self.y)


        Rf = np.fft.fft(R) * dt
        Sf = np.fft.fft(s) * dt

        conv = np.fft.ifft(Sf * Rf) / dt


        return wConv, F, T, R, F_Signal, F_Noise, conv

    def plotFilters(self):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        w, F, T, Fs, Fn, Dest = self.convolve()

        figOptFilt, axOptFilt = plt.subplots(figsize = (10,8))

        axOptFilt.loglog(self.t, F)

        figTrans, axTrans = plt.subplots(figsize = (10,8))

        axTrans.loglog(self.t, T)

        return

    def invFft(self):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        return

    def invMem(self):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        return

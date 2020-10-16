import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from MEM_class import MEM
import copy

class SignalToF():
    '''

    '''
    def __init__(self, t, y, psdType):
        '''
            Arguments:
            ----------

            returns:
            --------
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

            returns:
            --------
        '''
        dt = self.t[1] - self.t[0]

        w = np.fft.fftfreq(N, dt)
        A = np.fft.fft(self.y, n = N) #* dt

        return w, A

    def fft_psd(self, N = 8192):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        w, s = self.fft(N)
        aS2 = np.abs(s)**2
        P = (aS2[np.where(w >= 0)] + aS2[np.where(w < 0)][::-1]) / N**2
        w_pos = w[np.where(w > 0)]

        return w_pos, P

    def mem(self, M, N):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        MEM_inst = MEM(t_data = self.t, y_data = self.y, M = M)
        w, P = MEM_inst(t_data = self.t, y_data = self.y, M = M, N = N, view = False)

        return w, P

    def OptFilterFit(self, **kwargs):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        if self.psdType == 'FFT':
            f, P = self.fft_psd(N=8192)#fft_psd()#
        elif self.psdType == 'MEM':
            f, P = self.mem(M=100, N=8192)

        dt = self.t[1] - self.t[0]


        def calc_res(params, x, y, dt, weights):

            P0, s_eta2, s_tot2, a1 = params

            Noise = self.func_Noise(x, s_eta2, a1, dt)
            Signal = self.func_Signal(x, P0, s_tot2)

            Pmod = Noise + Signal

            if self.psdType == 'FFT':
                res = weights*(np.log10(y[:-1]) - np.log10(np.copy(Pmod)))
            elif self.psdType == 'MEM':
                res = weights*(np.log10(y) - np.log10(np.copy(Pmod)))

            return res

        def sum2_res(params, x, y, dt, weights):
            return np.sum(calc_res(params, x, y, dt, weights)**2)

        boundas = {}
        boundas['P0_Min'] = 1e-15
        boundas['P0_Max'] = 10
        boundas['s_eta2_Min'] = 1e-10
        boundas['s_eta2_Max'] = 0.1
        boundas['a1_Min'] = 1e-7
        boundas['a1_Max'] = 0.3
        boundas['s_tot2_Min'] = 1e-7
        boundas['s_tot2_Max'] = 0.5

        if list(kwargs.keys()):
            print('Setting fit param boundaries to user specifics')
            for j in list(kwargs.keys()):
                if j in list(bounds.keys()):
                    bounds[j] = kwargs[j]
                    print(f'setting {j} as {kwargs[j]}')
        elif not list(kwargs.keys()):
            print('Using default boundaries for variance and a1')

        p0 = [0.005, 0.005, 0.01, 0.1]

        weights = np.ones_like(f)*1.

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

            returns:
            --------
        '''
        return (s_eta2**2 * dz) / (np.abs(1 - a1 * np.exp(- 2 * np.pi * 1j * w * dz))**2)

    def func_Signal(self, w, p0, s_tot2):
        '''
            Arguments:
            ----------

            returns:
            --------
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
        T = np.exp(- ((2 * np.pi * wConv )**2 * (s_eta2)**2)/2)
        R = F * T**(-1)

        s = copy.deepcopy(self.y)


        Rf = R
        Sf = np.fft.fft(s) * dt

        conv = np.fft.ifft(Sf * Rf) / dt


        return wConv, F, T, R, F_Signal, F_Noise, conv

    def plotFilters(self):
        w, F, T, Fs, Fn, Dest = self.convolve()

        figOptFilt, axOptFilt = plt.subplots(figsize = (10,8))

        axOptFilt.loglog(self.t, F)

        figTrans, axTrans = plt.subplots(figsize = (10,8))

        axTrans.loglog(self.t, T)


        '''
            Arguments:
            ----------

            returns:
            --------
        '''
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

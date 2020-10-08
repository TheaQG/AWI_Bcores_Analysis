import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from MEM_class import MEM

class SignalToF():
    '''

    '''
    def __init__(self, t, y):
        '''
            Arguments:
            ----------

            returns:
            --------
        '''
        self.t = t
        self.y = y

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
        A = np.fft.fft(self.y, n = N) * dt

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
        def func_Noise(w, s_eta2, a1, dz):
            return (s_eta2 * dz) / (abs(1 + a1 * np.exp(- 2 * np.pi * 1j * w * dz))**2)



        def func_Signal(w, p0, s_tot2):
            return p0 * np.exp(- (2 * np.pi * w)**2 * s_tot2)

        def calc_res(x, y, params, dt, weights):

            P0, s_tot2, s_eta2, a1 = params

            Noise = func_Noise(x, s_eta2, a1, dt)
            Signal = func_Signal(x, P0, s_tot2)

            P = Noise + Signal
            res = weights*(np.log10(y) - np.log10(model))

            return res

        def sum2_res(x, y, params, dt, weights):
            return np.sum(calc_res(x, y, params, dt, weights)**2)

        bounds = {}
        bounds['s_eta2_Min'] = 0.00001
        bounds['s_eta2_Max'] = 2.
        bounds['a1_Min'] = 0.00001
        bounds['a1_Max'] = 0.2
        bounds['s_tot2_Min'] = 0.00001
        bounds['s_tot2_Max'] = 0.5

        if list(kwargs.keys()):
            print('Setting fit param boundaries to user specifics')
            for j in list(kwargs.keys()):
                if j in list(bounds.keys()):
                    bounds[j] = kwargs[j]
                    print(f'setting {j} as {kwargs[j]}')
        elif not list(kwargs.keys()):
            print('Using default boundaries for variance and a1')

        p0 = [0.5, 0.005, 0.01, 0.1]
        


        return

    def convolve(self):
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

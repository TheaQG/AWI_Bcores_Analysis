import numpy as np
import matplotlib.pyplot as plt
from MEM_class import MEM
from scipy.optimize import curve_fit

class WienerFilter():
    def __init__(self, name, depth_data, d18O_data, view_PSD = True, view_fit = True, print_DiffLen = True):
        self.depth_data = depth_data
        self.d18O_data = d18O_data
        self.view_PSD = view_PSD
        self.view_fit = view_fit
        self.print_DiffLen = print_DiffLen
        self.name = name
        return

    def __call__(self):
        return

    def calc_PSD(self):
        MEM_instance = MEM(t_data = self.depth_data, y_data = self.d18O_data, M = 70)
        MEM_power = MEM_instance(t_data = self.depth_data, y_data = self.d18O_data, M = 100, N = np.size(d18O_data), view = False)

        if self.view_PSD:
            figPSD, axPSD = plt.subplots(figsize = (10,6))
            axPSD.set(xlabel = 'Frequency [cycles/m]', ylabel = r'Power [$\permil ^2$ m]')
            axPSD.semilogy(MEM_power[0], MEM_power[1], label = self.name)
            axPSD.legend()

        return MEM_power

    def fit_PSD(self):
        def func_Noise(w, s_eta2, a1):
            dz = 0.02
            return (s_eta2 * dz) / (abs(1 + a1 * np.exp(- 2 * np.pi * 1j * w * dz))**2)



        def func_Signal(w, p0, s_tot2):
            return p0 * np.exp(- (2 * np.pi * w)**2 * s_tot2)


        x_All = self.calc_PSD()[0]
        y_All = self.calc_PSD()[1]

        x_Signal = x_all[x_all < 5]
        y_Signal = y_all[x_all < 5]

        x_Noise = x_all[x_all >= 5]
        y_Noise = y_all[x_all >= 5]

        popt_Signal, pvoc_Signal = curve_fit(func_Signal, x, y)
        popt_noise, pvoc_noise = curve_fit(func_Noise, x, y)

        s_eta2_est = popt_Noise[0]
        a1_est = popt_Noise[1]
        p0_est = popt_Signal[0]
        s_tot2_est = popt_Signal[1]

        if self.print_DiffLen:
            print(f'Diffusion length estimate: {s_tot2_est}')

        yNew_Signal = func_Signal(x_Signal, p0_est, s_tot2_est)
        yNew_Signal_allF = func_Signal(x_all, p0_est, s_tot2_est)

        yNew_Noise = func_Noise(x_Noise, s_eta2_est, a1_est)
        yNew_Noise_allF = func_Noise(x_all, s_eta2_est, a1_est)

        y_Combo = yNew_Signal_allF + yNew_Noise_allF

        if view_fit:
            figFit, axFit = plt.subplots(figsize=(7,4))
            axFit.semilogy(x_All, y_All, color='k', label='Data')
            axFit.semilogy(x_Signal, yNew_Signal, label='Fit to signal')
            axFit.semilogy(x_all, yNew_Noise_allF, label='Fit to noise')

            axFit.semilogy(x_all, y_Combo, label='Fit to all')
            axFit.legend()
#            axFit.set(ylim = (10**-4.5,10**0.5) )

        return [x_Signal, y_Signal, yNew_Signal], [x_Noise, y_Noise]

    def deconvolce(self):
        filt = y_testAll/y_ALL
        trans = np.exp(- ((2 * np.pi * x_all)**2 * s_tot2_est) / 2)
        transInv = 1/trans

        d18O_mean = np.mean(self.d18O_data)
        F = np.fft.fft(d18O_data_2[1:-1:2]) * filt[:-1] * transInv[:-1]
        f = np.fft.ifft(F)


        return

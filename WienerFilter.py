import numpy as np
import matplotlib.pyplot as plt
from MEM_class import MEM
from scipy.optimize import curve_fit

class WienerFilter():
    def __init__(self, name, depth_data, d18O_data, N, M, bar, view_PSD = True, view_fit = True, view_decon = True, print_DiffLen = True):
        self.depth_data = depth_data
        self.d18O_data = d18O_data
        self.N = N
        self.M = np.copy(M)
        self.bar = bar
        self.view_PSD = view_PSD
        self.view_fit = view_fit
        self.view_decon = view_decon
        self.print_DiffLen = print_DiffLen
        self.name = name
        return

    def __call__(self):
        f_data = self.deconvolve()
        f_depth = self.depth_data[1:-1:2]

        _,_,_, fit_params = self.fit_PSD(view=False, print_DiffLen=False)
        diff_len = np.sqrt(fit_params[3])*1000

        if self.view_decon:
            fig_decon, ax_decon = plt.subplots(figsize=(7,4))
            ax_decon.plot(self.depth_data, self.d18O_data+np.mean(self.d18O_data))
            ax_decon.plot(f_depth, f_data+np.mean(self.d18O_data))

        return f_depth, f_data+np.mean(self.d18O_data), fit_params

    def calc_PSD(self):
        MEM_instance = MEM(t_data = self.depth_data, y_data = self.d18O_data, M = self.M)
        MEM_power = MEM_instance(t_data = self.depth_data, y_data = self.d18O_data, M = self.M, N = self.N, view = False)

        if self.view_PSD:
            figPSD, axPSD = plt.subplots(figsize = (7,4))
            axPSD.set(xlabel = 'Frequency [cycles/m]', ylabel = r'Power [$\permil ^2$ m]')
            axPSD.semilogy(MEM_power[0], MEM_power[1], label = self.name)
            axPSD.legend()

        return MEM_power

    def fit_PSD(self, view, print_DiffLen):
        def func_Noise(w, s_eta2, a1):
            dz = 0.02
            return (s_eta2 * dz) / (abs(1 + a1 * np.exp(- 2 * np.pi * 1j * w * dz))**2)



        def func_Signal(w, p0, s_tot2):
            return p0 * np.exp(- (2 * np.pi * w)**2 * s_tot2)


        x_All = self.calc_PSD()[0]
        y_All = self.calc_PSD()[1]

        x_Signal = x_All[x_All < self.bar]
        y_Signal = y_All[x_All < self.bar]

        x_Noise = x_All[x_All >= self.bar]
        y_Noise = y_All[x_All >= self.bar]

        popt_Signal, pvoc_Signal = curve_fit(func_Signal, x_Signal, y_Signal)
        popt_Noise, pvoc_Noise = curve_fit(func_Noise, x_Noise, y_Noise, bounds=([-100,-100],[100,100]))

        s_eta2_est = popt_Noise[0]
        a1_est = popt_Noise[1]
        p0_est = popt_Signal[0]
        s_tot2_est = popt_Signal[1]

        fit_params = [s_eta2_est, a1_est, p0_est, s_tot2_est]

        if print_DiffLen:
            print(f'Diffusion length estimate (mm): {s_tot2_est*1000:.4f}\n')


        yNew_Signal = func_Signal(x_Signal, p0_est, s_tot2_est)
        yNew_Signal_allF = func_Signal(x_All, p0_est, s_tot2_est)

        yNew_Noise = func_Noise(x_Noise, s_eta2_est, a1_est)
        yNew_Noise_allF = func_Noise(x_All, s_eta2_est, a1_est)

        y_Combo = yNew_Signal_allF + yNew_Noise_allF

        if view:
            figFit, axFit = plt.subplots(figsize=(7,4))
            axFit.semilogy(x_All, y_All, color='k', label='Data')
            axFit.semilogy(x_Signal, yNew_Signal, label='Fit to signal')
            axFit.semilogy(x_All, yNew_Noise_allF, label='Fit to noise')

            axFit.semilogy(x_All, y_Combo, label='Fit to all')
            axFit.legend()
            axFit.set_ylim((np.min(y_All)-np.min(y_All)/20,np.max(y_All)+np.max(y_All)/20))

        Signal_data = [x_Signal, y_Signal, yNew_Signal_allF]
        Noise_data = [x_Noise, y_Noise, yNew_Noise_allF]
        Combo_data = [x_All, y_All, y_Combo]

        return Signal_data, Noise_data, Combo_data, fit_params

    def deconvolve(self):
        Signal_data, Noise_data, All_data, fit_Params = self.fit_PSD(view=self.view_fit, print_DiffLen=True)
        filt = Signal_data[2]/All_data[2]

        trans = np.exp(- ((2 * np.pi * All_data[0])**2 * 0.001) / 2)
        transInv = 1/trans

#        d18O_mean = np.mean(self.d18O_data)
        d18O_data_Filt = self.d18O_data[1:-1:2]

        if len(d18O_data_Filt) == len(filt):
            F = np.fft.fft(d18O_data_Filt) * filt * transInv

        elif (len(d18O_data_Filt) + 1) == len(filt):
            F = np.fft.fft(d18O_data_Filt) * filt[:-1] * transInv[:-1]

        else:
            F = np.zeros_like(d18O_data_Filt)
            print('Something is wrong')

        f = np.fft.ifft(F)


        return f
    # def plotFilters(self):
    #     Signal_data, Noise_data, Combo_data, fit_params = self.fit_PSD()
    #
    #     filt = Signal_data[2]/All_data[2]
    #
    #     trans = np.exp(- ((2 * np.pi * All_data[0])**2 * fit_Params[3]) / 2)
    #     transInv = 1/trans
    #
    #     R = filt * transInv
    #
    #     figFilt, axFilt = plt.subplots(figsize=(7,4))
    #     axFilt.loglog()
    #
    #     return

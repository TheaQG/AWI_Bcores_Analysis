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
from scipy import signal

from Decon import SpectralDecon


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{textcomp}',
    r'\usepackage{wasysym}']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'STIXGeneral'

class BackDiffuse():

    def __init__(self, coreName, d18OData, coreSpecs, depthMin, depthMax, ysInSec, interpAll = False, diffDensData_in = True, diffLenData = None, densData = None):
        self.coreName = coreName
        self.d18OData = d18OData
        self.coreSpecs = coreSpecs
        self.depthMin = depthMin
        self.depthMax = depthMax
        self.ysInSec = ysInSec
        self.interpAll = interpAll
        self.diffDensData_in = diffDensData_in

        self.densData = densData
        self.diffLenData = diffLenData

        return




    def __call__():
        return


    def interpCores(self, pad = 1):
        isoData = self.d18OData
        d_in = isoData['depth']
        x_in = isoData['d18O']


        if self.interpAll:
            valMin = d_in.min()
            valmax = d_in.max()
        else:
            valMin = self.depthMin - pad
            valMax = self.depthMax + pad


        d = d_in[(d_in >= valMin) & (d_in <= valMax)]
        x = x_in[(d_in >= valMin) & (d_in <= valMax)]

        diff = np.diff(d)
        Delta = round(min(diff), 3)

        d_min = Delta * np.ceil(d.values[0]/Delta)
        d_max = Delta * np.floor(d.values[-1]/Delta)

        n = int(1 + (d_max - d_min)/Delta)

        j_arr = np.linspace(0,n,n)
        dhat0 = d_min + (j_arr - 1)*Delta

        f = interpolate.CubicSpline(d,x)

        xhat0 = f(dhat0)

        dhat = dhat0[(dhat0 >= self.depthMin) & (dhat0 <= self.depthMax)]
        xhat = xhat0[(dhat0 >= self.depthMin) & (dhat0 <= self.depthMax)]

        return dhat, xhat, Delta




    def diffProfile(self):
        if self.diffDensData_in:
            diffDepth = self.diffLenData['Depth']
            diffData = self.diffLenData['sigma_o18']
            return diffDepth, diffData
        else:
            print('Compute diff len profile first!')

            return [], []

    def densProfile(self):
        if self.diffDensData_in:
            densDepth = self.densData['depth']
            densHL = self.densData['HLmodel']

            return densDepth, densHL

        else:
            print('Compute density profile first"')

            return [], []

    def spectralEstimate(self,N=2000):
        dInt, d18OInt, Delta = self.interpCores()

        decon_inst = SpectralDecon(dInt, d18OInt, N)
        w_PSD, P_PSD, Pnoise, Psignal, P_fit, opt_fit_dict, params_fit, fit_func_val, fit_dict = decon_inst.SpectralFit(printFitParams=False, printDiffLen=False)
        diffLen_FitEst = opt_fit_dict['s_tot2_fit']


        return diffLen_FitEst

    def diffLenEstimateHL(self):
        diffDepth, diffData = self.diffProfile()

        zSec = diffDepth[(diffDepth >= self.depthMin) & (diffDepth <= self.depthMax)]
        diffSec = diffData[(diffDepth >= self.depthMin) & (diffDepth <= self.depthMax)]

        sigma_range = [diffSec.min(), diffSec.max()]

        return sigma_range

    def backDiffused(self, N=2000):
        sigma_rangeHL = self.diffLenEstimateHL()
        sigma_FitEst = self.spectralEstimate()

        dInt, d18OInt, Delta = self.interpCores()

        diffLen0 = min(min(sigma_rangeHL), sigma_FitEst) - 0.01
        print(f'Starting sigma: {diffLen0*100:.2f} [cm]')

        decon_inst = SpectralDecon(dInt, d18OInt, N)

        depth0, dataD0 = decon_inst.deconvolve(diffLen0)

        from scipy import signal
        N_peaks = 0

        depth = depth0
        data = dataD0
        diffLen = diffLen0
        while N_peaks != self.ysInSec:
            depth, data = decon_inst.deconvolve(diffLen)
            idxPeak = signal.find_peaks(data, distance=3)[0]
            N_peaks = len(idxPeak)
            print(len(idxPeak))

            if N_peaks > self.ysInSec:
                diffLen -= 0.0005
            if N_peaks < self.ysInSec:

                diffLen += 0.0005
        print(f'Final sigma: {diffLen*100:.2f} [cm]')
        depthEst = depth
        dataEst = data
        diffLenFin = diffLen

        return depthEst, dataEst, diffLenFin, idxPeak

#temp_holo = 213.15, delta_holo = -51, slope = 0.69

    def DeltaToTemp(self, temp_holo = 213.15, delta_holo = -51, slope = 0.69):
        depth, data,_,_ = self.backDiffused()
        temp = np.zeros(np.size(data)) + temp_holo + slope*(data - delta_holo)

        return depth, temp


    def interpDiffData():
        return

    def countPeaks():
        return
#['Crete', 'SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
coreName = 'SiteG'

d18OData = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+coreName+'_det.txt',',')
densities = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/'+coreName+'DepthDens_w_Models.txt','\t')
diffLens = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/'+coreName+'_DepthDiff.txt','\t')
specsCores = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt',',')

dTamb = 60.50
dLaki = 69.40
depth_LT = d18OData['depth'][(d18OData['depth'] >= dTamb) & (d18OData['depth'] <= dLaki)]
d18O_LT = d18OData['d18O'][(d18OData['depth'] >= dTamb) & (d18OData['depth'] <= dLaki)]

inst = BackDiffuse(coreName, d18OData, specsCores, dTamb, dLaki, 32, diffLenData=diffLens[['Depth','sigma_o18']], densData=densities)
test, test2, test3 = inst.interpCores()

diffLen = inst.spectralEstimate()

difflenEstHL = inst.diffLenEstimateHL()

depth, data, diffLen, peaks = inst.backDiffused()

fig, ax = plt.subplots(figsize=(10,7))
ax.plot(depth_LT, d18O_LT-np.mean(d18O_LT),color='k', lw=1, label = 'Data')
ax.plot(depth, data, lw=1, label='Back diffused')
ax.plot(depth[peaks], data[peaks],'.',lw=1, label='Estimated peaks')
ax.set(xlabel = 'Depth [m]', ylabel = '$\delta^{18}$O [\permil]', title=coreName)
ax.legend(fontsize=16)
fig.tight_layout()
fig.savefig(coreName + '_peaks.jpg')



depthT, dataT = inst.DeltaToTemp()

fig2, ax2 = plt.subplots(figsize=(10,7))
ax2.plot(depthT, dataT-273.15,color='k', lw=1)
ax2.set(xlabel = 'Depth [m]', ylabel = 'Temperature [C]', title=coreName)
fig2.savefig(coreName + '_Temp.jpg')

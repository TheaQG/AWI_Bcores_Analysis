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


class BackDiffuse():
    '''
        Methods available:
            __init__(self, coreName, d18OData, coreSpecs, depthMin, depthMax, ysInSec, interpAll = False, diffDensData_in = True, diffLenData = None, densData = None):
                    Initializes class given six required arguments and four optional.

            __call__():

            interpCores(self, pad = 1):
                    Interpolates core d18O data through a cubic spline with new
                    sampling size of the size of the minimum sampling size available.
                    Interpolates between depthMin and depthMax with a padding of
                    size pad, default pad = 1 [m].

            diffProfile(self):
                    Sets the diffusion profile (depth, diffusion length data) as
                    the given input, diffLenData. If None passed, exits class.
                    Need to build extra code to compute diff len if None is passed...

            densProfile(self):
                    Sets the density profile (depth, density data) as
                    the given input, densData. If None passed, exits class.
                    Need to build extra code to compute dens if None is passed...

            spectralEstimate(self,N=2000):
                    Based on SpectralDecon. Computes a diff len estimate from
                    the PSD. This diff len is then used together with theoretical
                    est. diff len to choose minimal starting diff len.

            diffLenEstimateHL(self):
                    Gives the minimum and maximum diffusion length based on HL
                    diffusion profile between depthMin and depthMax. This range
                    is used together with the spectral estimate to choose min
                    diff len 0.

            backDiffused(self, N=2000):
                    Estimates the highest diffusion length to deconvolve with that
                    will still give ysInSec peaks in depth section [depthMin, depthMax].
                    Need to estimate peaks in a better eay than SciPy...

            DeltaToTemp(self, temp_holo = 213.15, delta_holo = -51, slope = 0.69):
                    Calculates the temperature estimate in [K] of a given depth, based
                    on the measured delta values and a slope, a holocene temperature
                    estimate and a holocene delta estimate.
    '''
    def __init__(self, coreName, d18OData, coreSpecs, depthMin, depthMax, ysInSec, interpAll = False, diffDensData_in = True, diffLenData = None, densData = None):
        '''
            Initialize the class instance.


            Arguments:
            ----------
                coreName:       [str] Name of examined core. Needs to match name in corespecs and data files.
                d18OData:       [pd.DataFrame] d18O and depth data of core. May need to change to np array...
                coreSpecs:      [pd.DataFrame] Dataframe containing all available cores specs.
                depthMin:       [float] Minimum depth value to examine.
                depthMax:       [float] Maximum depth value to examine
                ysInSec:        [int] How many years in section [depthMin, depthMax].
                interpAll:      [bool] Default = False. If wanting to interpolate whole core. Not advised.
                diffDensData_in:    [bool] Default = True
                diffLenData:    [pd.DataFrame] Default None. Diff len and depth data of core. May need to change to np array...
                densData        [pd.DataFrame] Default None. Dens and depth data of core. May need to change to np array...

            returns:
            --------
                None

        '''

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
        '''


            Arguments:
            ----------

            returns:
            --------

        '''
        return


    def interpCores(self, pad = 1):
        '''
            Computes interpolation btw. depthMin-pad and depthMax+pad, or for whole core.
            Padding is added to avoid undesired border effects.
            Uses cubic spline interpolation with a new sample size equal to the
            minimal sample size in the original data set.


            Arguments:
            ----------
                pad:            [float] Padding to avoid werid borders. Default = 1 [m].

            returns:
            --------
                dhat:           [array of floats] Depth data corresponding to interpolated data.
                xhat:           [array of floats] d18O data corresponding to interpolated data.
                Delta:          [float] New sample size.

        '''
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
        '''
            Gets the diffusion length profile data (diff len v. depth) from passed files.
            If no passed file tells you to create such file.
            (In future: make method compute theoretical/empirical diffLenData file if possible)

            Arguments:
            ----------
                None

            returns:
            --------
                diffDepth:      [arr of floats] Depth data.
                diffData:       [arr of floats] Diffusion length data

        '''
        if self.diffDensData_in:
            diffDepth = self.diffLenData['Depth']
            diffData = self.diffLenData['sigma_o18']
            return diffDepth, diffData
        else:
            print('Compute diff len profile first!')

            return [], []

    def densProfile(self):
        '''
            Gets the density profile data (dens v. depth) from passed files.
            If no passed file tells you to create such file.
            (In future: make method compute theoretical/empirical densData file if possible)
            Needs also to incorporate HLmodelOpti and densMeas if existing.

            Arguments:
            ----------
                None

            returns:
            --------
                densDepth:      [arr of floats] Depth data.
                densHL:         [arr of floats] Density data, HL model.

        '''
        if self.diffDensData_in:
            densDepth = self.densData['depth']
            densHL = self.densData['HLmodel']

            return densDepth, densHL

        else:
            print('Compute density profile first"')

            return [], []

    def spectralEstimate(self,N=2000):
        '''
            Based on SpectralDecon instance, computes all parts for a spectral transform
            and fit given a depth series. Returns only estimated diff len.


            Arguments:
            ----------
                N:              [int] Number of points to generate spectral data with.

            returns:
            --------
                diffLen_FitEst: [float] From fit estimated diffusion length.

        '''
        dInt, d18OInt, Delta = self.interpCores()

        decon_inst = SpectralDecon(dInt, d18OInt, N)
        w_PSD, P_PSD, Pnoise, Psignal, P_fit, opt_fit_dict, params_fit, fit_func_val, fit_dict = decon_inst.SpectralFit(printFitParams=False, printDiffLen=False)
        diffLen_FitEst = opt_fit_dict['s_tot2_fit']


        return diffLen_FitEst

    def diffLenEstimateHL(self):
        '''
            Based on HL model, from self.diffProfile determines min and max values of
            diff len in interval [depthMin, depthMax]. Returns diff len range.


            Arguments:
            ----------
                None

            returns:
            --------
                sigma_rangeHL:  [arr of 2 floats] Max and min diff len estimated from theoretical/empirical data.

        '''
        diffDepth, diffData = self.diffProfile()

        zSec = diffDepth[(diffDepth >= self.depthMin) & (diffDepth <= self.depthMax)]
        diffSec = diffData[(diffDepth >= self.depthMin) & (diffDepth <= self.depthMax)]

        sigma_range = [diffSec.min(), diffSec.max()]

        return sigma_range

    def backDiffused(self, N=2000):
        '''
            Method to compute the maximal diffusion length that still give ysInSec
            peaks. Computes first any value that returns ysInSec peaks, and computes
            then the maximum value that still returns that exact number of peaks.

            Arguments:
            ----------
                N:              [int] Number of points to generate spectral data with.

            returns:
            --------
                depthEst:       [arr of floats] Depth to estimated data.
                dataEst:        [arr of floats] Backdiffused d18O data.
                diffLenFin:     [float] Final diffusion length estimate.
                idxPeak:        [arr of idxs] Number of peaks in the final data set.

        '''
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

        while N_peaks == self.ysInSec:
            depth, data = decon_inst.deconvolve(diffLen)
            idxPeak = signal.find_peaks(data, distance=3)[0]
            N_peaks = len(idxPeak)
            print(len(idxPeak))
            diffLen += 0.0001

        diffLen -= 0.0002
        depth, data = decon_inst.deconvolve(diffLen)
        idxPeak = signal.find_peaks(data, distance=3)[0]
        N_peaks = len(idxPeak)

        print(f'Final sigma: {diffLen*100:.2f} [cm]')
        print(f'Final # of peaks: {N_peaks}')
        depthEst = depth
        dataEst = data
        diffLenFin = diffLen

        return depthEst, dataEst, diffLenFin, idxPeak

#temp_holo = 213.15, delta_holo = -51, slope = 0.69

    def DeltaToTemp(self, temp_holo = 213.15, delta_holo = -51, slope = 0.69):
        '''
            Method to estimate temperature in [k] given a d18O depth series, along with specific
            temp and d18O estimates for holocene.


            Arguments:
            ----------
                temp_holo:      [float] Estimated holocene mean temperature.
                delta_holo:     [float] Estimated holocene mean d18O .
                slope:          [float] Estimated slope.

            returns:
            --------
                depth:          [arr of floats] Depth data.
                temp:           [arr of floats] Estimated temperature data.
        '''
        depth, data,_,_ = self.backDiffused()
        temp = np.zeros(np.size(data)) + temp_holo + slope*(data - delta_holo)

        return depth, temp


    def interpDiffData():
        '''


            Arguments:
            ----------

            returns:
            --------

        '''
        return

    def countPeaks():
        '''


            Arguments:
            ----------

            returns:
            --------

        '''
        return









mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{textcomp}',
    r'\usepackage{wasysym}']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'STIXGeneral'


    # Core names of cores available
coreNames = ['Crete', 'SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']

    # Selecting core name
coreName = 'SiteB'

    # Reading datafiles for specific core
d18OData = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+coreName+'_det.txt',',')
densities = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/'+coreName+'DepthDens_w_Models.txt','\t')
diffLens = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/'+coreName+'_DepthDiff.txt','\t')
specsCores = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt',',')
specIdx = specsCores['CoreName'][specsCores['CoreName'] == coreName].index[0]

    # Set the specs for depth of Laki and Tambora eruptions for core
specsCore = specsCores.iloc[specIdx]
dTamb = np.float64(specsCore['dTamb'])
dLaki = np.float64(specsCore['dLaki'])

    # (FOR PLOTTING) Make array of only d18O data between Laki and Tamb
depth_LT = d18OData['depth'][(d18OData['depth'] >= dTamb) & (d18OData['depth'] <= dLaki)]
d18O_LT = d18OData['d18O'][(d18OData['depth'] >= dTamb) & (d18OData['depth'] <= dLaki)]

    # Create instance of back diffusion
inst = BackDiffuse(coreName, d18OData, specsCores, dTamb, dLaki, 34, diffLenData=diffLens[['Depth','sigma_o18']], densData=densities)

    # Make spectral estimate of diff len
diffLen = inst.spectralEstimate()

    # Make model/empiric estimate of diff len
difflenEstHL = inst.diffLenEstimateHL()

    # Compute final depth/d18O back diffused data w. final diff len and No. peaks
depth, data, diffLen, peaks = inst.backDiffused()

    # Plot original data, back diffused data and peak estimations
fig, ax = plt.subplots(figsize=(10,7))
ax.plot(depth, data, lw=1, label='Back diffused')
ax.plot(depth_LT, d18O_LT-np.mean(d18O_LT),color='k', lw=1, label = 'Data')
ax.plot(depth[peaks], data[peaks],'.',lw=1, label='Estimated peaks')
ax.set(xlabel = 'Depth [m]', ylabel = '$\delta^{18}$O [\permil]', title=coreName+'$, \sigma_{fin} =$ ' + f'{diffLen*100:.2f} [cm]')
ax.legend(fontsize=16)
fig.tight_layout()
fig.savefig(coreName + '_peaks.jpg')


    # Compute temperature estimate in [K]
depthT, dataT = inst.DeltaToTemp()

    # Plot Temperature data in [C]
fig2, ax2 = plt.subplots(figsize=(10,7))
ax2.plot(depthT, dataT-273.15,color='k', lw=1)
ax2.set(xlabel = 'Depth [m]', ylabel = 'Temperature [C]', title=coreName)
fig2.savefig(coreName + '_Temp.jpg')

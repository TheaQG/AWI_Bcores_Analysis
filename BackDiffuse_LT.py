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
from SignalAttenuation import Attenuation, AnnualLayerThick
from Decon import SpectralDecon

'''
        *********** TODO************
        - (12/02/21) Create method w. diffusion length as fct., sigma(z), instead of fixed.
'''

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
    def __init__(self, coreName, d18OData, coreSpecs, depthMin, depthMax, ysInSec, interpAll = False, diffDensData_in = True, diffLenData = None, densData = None, Dist=3, transType='DCT'):
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
        self.Dist = Dist
        self.densData = densData
        self.diffLenData = diffLenData
        self.transType = transType
        return




    def __call__():
        '''


            Arguments:
            ----------

            returns:
            --------

        '''
        return


    def interpCores(self, pad = 1, DeltaInput = False, DeltaIn = 0):
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

        if DeltaInput:
            Delta = DeltaIn
        else:
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

    def spectralEstimate(self,N=4000):
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

        decon_inst = SpectralDecon(dInt, d18OInt, N, self.transType)
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

    def backDiffused(self, N=4000, print_Npeaks=True, theoDiffLen=True, diffLenStart_In=0., diffLenEnd_In=0.1, interpAfterDecon=True, newDelta=0, interpBFDecon=True, diffLenMin = 0.):
        '''
            Method to compute the maximal diffusion length that still give ysInSec
            peaks. Computes first any value that returns ysInSec peaks, and computes
            then the maximum value that still returns that exact number of peaks.

            Arguments:
            ----------
                N:                  [int] Default = 4000. Number of points to generate spectral data with.
                print_Npeaks:       [bool] Default = True.
                theoDiffLen:        [bool] Default = True.
                diffLenStart_In:    [float] Default = 0.
                diffLenEnd_In:      [float] Default = 0.1.
                interpAfterDecon:   [bool] Default = True.
                newDelta:           [float] Default = 0.
                interpBFDecon:      [bool] Default = True.
                diffLenMin:         [float] Default = 0.

            returns:
            --------
                depthEst:           [arr of floats] Depth to estimated data.
                dataEst:            [arr of floats] Backdiffused d18O data.
                diffLenFin:         [float] Final diffusion length estimate.
                idxPeak:            [arr of idxs] Number of peaks in the final data set.
                arr_diffLens:       [list of floats]
                arr_Npeaks:         [list of ints]
                arr_depth:          [list of arrs]
                arr_data:           [list of arrs]
        '''


            ###############
            ## INITIALIZATION AND COMPUTATION OF FIRST ESTIMATES
            ###############


            # Compute theoretical diffusion length estimates
        sigma_rangeHL = self.diffLenEstimateHL()

            # Compute diffusion length estimate of spectral fit
        sigma_FitEst = self.spectralEstimate()

            # If interpolation before deconvolution is wanted, then interpolate the isotope data
        if interpBFDecon:
            dInt, d18OInt, Delta = self.interpCores()
        else:
            isoData = self.d18OData
            dInt = np.asarray(isoData['depth'])
            d18OInt = np.asarray(isoData['d18O'])
            Delta = dInt[1] - dInt[0]

            # Decide if the first diffusion length estimate should be given from the theoretical estimate.
        if theoDiffLen:
            diffLen0 = min(min(sigma_rangeHL), sigma_FitEst) - 0.02
            print(f'Sigma fit: {sigma_FitEst*100:.2f}')
            print(f'Min sigma analyt: {min(sigma_rangeHL)*100:.2f}')

            # If not, then use the input diffusion length
        else:
            diffLen0 = diffLenStart_In
        print(f'Starting sigma: {diffLen0*100:.2f} [cm]')

            # Create an instance of the spectralDecon class from Decon.py, with the (non-) interpolated data.
        decon_inst = SpectralDecon(dInt, d18OInt, N, self.transType)

            # Use the wanted first sigma estimate to compute the first estimate of the back diffused data.
        depth0, dataD0 = decon_inst.deconvolve(diffLen0)



            ###############
            ## INITIATE THE PEAK COUNTING/CONSTRAINT MODULE OF ALGORITHM
            ###############

        from scipy import signal
            # Set the initial number of counted peaks to 0
        N_peaks = 0
            # Set the data to the initial values
        depth = depth0
        data = dataD0
        diffLen = diffLen0
            # Create empty lists for storing data for each run
        arr_diffLens = []
        arr_Npeaks = []
        arr_depth = []
        arr_data = []
            # Start the count at 0
        i = 0

            ###############
            ## SIGMA OPTIMIZATION ROUTINE UNDER CONSTRAINTS
            ## MODULE 1: N != ysInSec
            ###############

            # Run this routine while the number of counted peaks is different from the number of wanted peaks (ysInSec: years in section)
        while N_peaks != self.ysInSec:
                # Compute back diffused signal data from current sigma (diffLen)
            depth, data = decon_inst.deconvolve(diffLen)

                # Interpolate the BD signal, if wanted
            if interpAfterDecon:
                    # If no interpolation size (newDelta) is inputted, choose half the first sample size
                if newDelta == 0:
                    newDeltaUse = (depth[1] - depth[0])/2
                    # Otherwise, use the inputted new sample size (newDelta)
                else:
                    newDeltaUse = newDelta
                    # Interpolate and define new data to use
                newDepth, newData, _ = interpCores2(depth[0], depth[-1], pd.Series(depth), pd.Series(data), DeltaInput=True, DeltaIn=newDeltaUse)

                    # Calculate the expected average distance between peaks (i.e. length of section divided by years in section)
                ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
                    # Calculate the expected average number of points between peaks
                ave_Npoints = ave_dist/newDeltaUse
                    # Determine a minimum distance between peaks. Based on the self.Dist variable.
                min_peakDist = int(ave_Npoints/self.Dist)

                # If no interpolation wanted, use signal data as new data to use
            else:
                newDepth = depth
                newData = data
                    # Calculate the expected average distance between peaks (i.e. length of section divided by years in section)
                ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
                    # Calculate the expected average number of points between peaks
                ave_Npoints = ave_dist/(newDepth[-1] - newDepth[0])
                    # Determine a minimum distance between peaks. Based on the self.Dist variable.
                min_peakDist = int(ave_Npoints/self.Dist)

                # If the computed minimum peak distance is zero, then add no constraint when using find_peaks.
            if min_peakDist==0:
                idxPeak = signal.find_peaks(newData)[0]
                # Otherwise, use the computed minimum peak distance as constraint.
            else:
                idxPeak = signal.find_peaks(newData, distance=min_peakDist)[0]

                # Calculate the number of counted peaks as the length of the idxPeak array.
            N_peaks = len(idxPeak)

                # Append the storing lists with their respective data
            arr_diffLens.append(diffLen)
            arr_Npeaks.append(N_peaks)
            arr_depth.append(newDepth)
            arr_data.append(newData)

                # If wanted, print the current diffusion length and the corresponding peaks counted
            if print_Npeaks:
                print(f'N peaks: {len(idxPeak)}')
                print(f'sigma: {diffLen*100}:.2f')

                # If the number of counted peaks is larger than the number of years in the section,
                # then subtract a small amount from the diffusion length estimate.
            if N_peaks > self.ysInSec:
                diffLen -= 0.0001
                i += 1
                if i%100 == 0:
                    print(f'{i}. Npeaks: {N_peaks}, diffLen: {diffLen*100:.3f} cm')

                # If the number of counted peaks is smaller than the number of years in the section,
                # then add a small amount from the diffusion length estimate.
            if N_peaks < self.ysInSec:
                diffLen += 0.0001005
                i += 1
                if i%100 == 0:
                    print(f'{i}. Npeaks: {N_peaks}, diffLen: {diffLen*100:.3f} cm')

                # If the diffusion length is larger than or equal to the maximal wanted diffusion length,
                # then stop the algorithm.
            if diffLen >= diffLenEnd_In:
                break
            if diffLen <= diffLenMin:
                break

            ###############
            ## SIGMA OPTIMIZATION ROUTINE UNDER CONSTRAINTS
            ## MODULE 2: N == ysInSec
            ###############

            # Run this routine while the number of counted peaks is equal to the number of wanted peaks (ysInSec: years in section)
        while N_peaks == self.ysInSec:
                # Compute back diffused signal data from current sigma (diffLen)
            depth, data = decon_inst.deconvolve(diffLen)
                # Interpolate the BD signal, if wanted
            if interpAfterDecon:
                    # If no interpolation size (newDelta) is inputted, choose half the first sample size
                if newDelta == 0:
                    newDeltaUse = (depth[1] - depth[0])/2
                    # Otherwise, use the inputted new sample size (newDelta)
                else:
                    newDeltaUse = newDelta
                    # Interpolate and define new data to use
                newDepth, newData, _ = interpCores2(depth[0], depth[-1], pd.Series(depth), pd.Series(data), DeltaInput=True, DeltaIn=newDeltaUse)
                    # Calculate the expected average distance between peaks (i.e. length of section divided by years in section)
                ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
                    # Calculate the expected average number of points between peaks
                ave_Npoints = ave_dist/newDeltaUse
                    # Determine a minimum distance between peaks. Based on the self.Dist variable.
                min_peakDist = int(ave_Npoints/self.Dist)

                # If no interpolation wanted, use signal data as new data to use
            else:
                newDepth = depth
                newData = data
                Delta = depth[1]-depth[0]
                    # Calculate the expected average distance between peaks (i.e. length of section divided by years in section)
                ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
                    # Calculate the expected average number of points between peaks
                ave_Npoints = ave_dist/Delta
                    # Determine a minimum distance between peaks. Based on the self.Dist variable.
                min_peakDist = int(ave_Npoints/self.Dist)

                # If the computed minimum peak distance is zero, then add no constraint when using find_peaks.
            if min_peakDist==0:
                idxPeak = signal.find_peaks(newData)[0]
                # Otherwise, use the computed minimum peak distance as constraint.
            else:
                idxPeak = signal.find_peaks(newData, distance=min_peakDist)[0]

                # Calculate the number of counted peaks as the length of the idxPeak array.
            N_peaks = len(idxPeak)

                # Append the storing lists with their respective data
            arr_diffLens.append(diffLen)
            arr_Npeaks.append(N_peaks)
            arr_depth.append(newDepth)
            arr_data.append(newData)

                # If wanted, print the current diffusion length and the corresponding peaks counted
            if print_Npeaks:
                print(f'N peaks: {len(idxPeak)}')
                print(f'sigma: {diffLen*100}:.2f')

                # Add a small positive perturbation to the diffusion length estimate.
            diffLen += 0.0001
            i += 1
            if i%100 == 0:
                print(f'{i}. Npeaks: {N_peaks}, diffLen: {diffLen*100:.3f} cm')


            ###############
            ## FINAL PART OF ALGORITHM. DIFF LEN MAX HAS BEEN FOUND AS THE
            ## THE LARGEST TO STILL GENERATE N = YSINSEC PEAKS
            ###############


            # Define the final sigma as the last sigma found minus twice the perturbation size (to get below ysInSec + 1 peaks)
        diffLen -= 0.0002

            # Use this final diffusion length estimate to compute the back diffused data
        depth, data = decon_inst.deconvolve(diffLen)

            # If interpolation after deconvolution is wanted, then interpolate
        if interpAfterDecon:
                # If the input new sample size is 0, use half of the original sampling size
            if newDelta == 0:
                newDeltaUse = (depth[1] - depth[0])/2
                # Otherwise, use input new sample size
            else:
                newDeltaUse = newDelta

                # Interpolate data with the given new sample size and create new signal data.
            newDepth, newData, _ = interpCores2(depth[0], depth[-1], pd.Series(depth), pd.Series(data), DeltaInput=True, DeltaIn=newDeltaUse)
            ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
            ave_Npoints = ave_dist/newDeltaUse
            min_peakDist = int(ave_Npoints/self.Dist)
        else:
            newDepth = depth
            newData = data
            Delta = newDepth[1] - newDepth[0]
                # Calculate the expected average distance between peaks (i.e. length of section divided by years in section)
            ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
                # Calculate the expected average number of points between peaks
            ave_Npoints = ave_dist/Delta
                # Determine a minimum distance between peaks. Based on the self.Dist variable.
            min_peakDist = int(ave_Npoints/self.Dist)

            # If the computed minimum peak distance is zero, then add no constraint when using find_peaks.
        if min_peakDist==0:
            idxPeak = signal.find_peaks(newData)[0]
            # Otherwise, use the computed minimum peak distance as constraint.
        else:
            idxPeak = signal.find_peaks(newData, distance=min_peakDist)[0]

            # Calculate the number of counted peaks as the length of the idxPeak array.
        N_peaks = len(idxPeak)

            # Append the storing lists with their respective data
        arr_diffLens.append(diffLen)
        arr_Npeaks.append(N_peaks)
        arr_depth.append(newDepth)
        arr_data.append(newData)

            # Print the final optimized variables: diffusion length, number of peaks, new (interpolated) sample size
        print(f'Final sigma: {diffLen*100:.2f} [cm]')
        print(f'Final # of peaks: {N_peaks}')
        print(f'Delta: {depth[1]-depth[0]:.3f}')
        print(f'Delta new: {newDepth[1]-newDepth[0]:.3f}')

            # Set the final estimates as the values found in the last iteration.
        depthEst = newDepth
        dataEst = newData
        diffLenFin = diffLen


        return depthEst, dataEst, diffLenFin, idxPeak, arr_diffLens, arr_Npeaks, arr_depth, arr_data

    def BackDiffused_constraints(self, sigDel0_in = 0.03, lSecs_in = 5, acceptPct_dist_in = 2/4, acceptPct_prom_in = 2/4, epsilon_in = 1e-10, kmax_in = 50, N_sigs_in = 5, N=2000, print_Npeaks=True, theoDiffLen=True, diffLenStart_In=0, diffLenEnd_In=0.1, interpAfterDecon=True, newDelta=0, interpBFDecon=True):
        '''
            Method to compute the maximal diffusion length that still give ysInSec
            peaks.

            Arguments:
            ----------
                sigDel0_in:         [float] Default = 0.03.
                lSecs_in:           [int] Default = 5.
                acceptPct_dist_in:  [float] Default = 2/4.
                acceptPct_prom_in:  [float] Default = 2/4.
                epsilon_in:         [float] Default = 1e-10.
                kmax_in:            [int] Default = 50.
                N_sigs_in:          [int] Default = 5.
                N:                  [int] Default = 4000. Number of points to generate spectral data with.
                print_Npeaks:       [bool] Default = True.
                theoDiffLen:        [bool] Default = True.
                diffLenStart_In:    [float] Default = 0.
                diffLenEnd_In:      [float] Default = 0.1.
                interpAfterDecon:   [bool] Default = True.
                newDelta:           [float] Default = 0.
                interpBFDecon:      [bool] Default = True.
                diffLenMin:         [float] Default = 0.

            returns:
            --------
                newDepth:           [arr of floats] Depth to estimated data.
                newData:            [arr of floats] Backdiffused d18O data.
                diffLenFin:         [float] Final diffusion length estimate.
                newPs:              [arr of idxs] Peak positions in the final data set.
                newTs:              [arr of idxs] Trough positions in the final data set.
                newPattern:         [bool] Is pattern present?
        '''


            ###############
            ## INITIALIZATION AND ALT ESTIMATE COMPUTATION
            ###############

            # Define length [m] of sections to compute ALT from
        lSecs = lSecs_in

            # Set entire core data (depth and d18O) as separate np-arrays
        isoData = self.d18OData
        depth_ALT = np.asarray(isoData['depth'])
        d18O_ALT = np.asarray(isoData['d18O'])

            # Create annual layer thickness instance
        inst_ALT = AnnualLayerThick(depth_ALT, d18O_ALT, lSecs)
            # Compute ALT for entire core.
        fksMax, ls, lMean, lStd, vals = inst_ALT.ALT_fullCore()
            # Compute an estimate for ALT at LT depth
        vals_use = vals[:-1]
        l_LT = np.mean(lMean[(vals_use > self.depthMin) & (vals_use < self.depthMax)])
        ALT_LT = l_LT

            # If interpolation before deconvolution is wanted, then interpolate the isotope data
        if interpBFDecon:
            dInt, d18OInt, Delta = self.interpCores()
        else:
            isoData = self.d18OData
            dInt = np.asarray(isoData['depth'])
            d18OInt = np.asarray(isoData['d18O'])
            Delta = dInt[1] - dInt[0]


            # Estimate # of points in a layer cycle, based on newly interpolated data
        points_ALT = ALT_LT * len(dInt)/(max(dInt) - min(dInt))
            # Set the accepted percentage of this distance estimate
        acceptPct_dist = acceptPct_dist_in
        dist = np.floor(points_ALT * acceptPct_dist)

            # Set the accepted percentage of the prominence
        acceptPct_prom = acceptPct_prom_in

            # Compute theoretical diffusion length estimates
        sigma_rangeHL = self.diffLenEstimateHL()
            # Compute diffusion length estimate of spectral fit
        sigma_FitEst = self.spectralEstimate()


            # Decide if the first diffusion length estimate should be given from the theoretical estimate.
        if theoDiffLen:
            diffLen0 = min(min(sigma_rangeHL), sigma_FitEst) - 0.02
            print(f'Sigma fit: {sigma_FitEst*100:.2f}')
            print(f'Min sigma analyt: {min(sigma_rangeHL)*100:.2f}')

            # If not, then use the input diffusion length
        else:
            diffLen0 = diffLenStart_In
        print(f'Starting sigma: {diffLen0*100:.2f} [cm]')


            # Create an instance of the spectralDecon class from Decon.py, with the (non-) interpolated data.
        decon_inst = SpectralDecon(dInt, d18OInt, N, self.transType)

            # Set the minimum peak distance as the accepted percentage of the ALT estimate
        min_peakDist = dist

            # Define the initial diff. len. perturbation as the inputted
        sigDel0 = sigDel0_in

            # Define initial diff len boundaries as twice the perturbation to each site of the diff len init estimate
        sigMin0 = diffLen0 - 2*sigDel0
        sigMax0 = diffLen0 + 2*sigDel0

            # Set the minimal difference as the inputted
        epsilon = epsilon_in
            # Set the maximal and initial number of itarations
        kmax = kmax_in
        k = 0

            # Set the sigma boundaries as the initial
        sigMin = sigMin0
        sigMax = sigMax0
            # Set the number of section boundaries to search through to the inputted
        N_sigs = N_sigs_in
            # Define an array with the initial diff len bounds/sections
        diffLens = np.linspace(sigMin, sigMax, N_sigs)

            # Set the number of peaks and troughs as the inputted from initialization of the class
        Npeaks = self.ysInSec
        Ntroughs = self.ysInSec


            ###############
            ## INITIATE THE PEAK COUNTING/CONSTRAINT MODULE OF ALGORITHM
            ###############


            # Continue to run this routine, while the difference in diff lens is larger than the inputted min diff epsilon
        while (diffLens[1] - diffLens[0]) > epsilon:

                # Define an array with the diff len bounds/sections
            diffLens = np.linspace(sigMin, sigMax, N_sigs)

                # Initialize empty arrays and lists to store final results in
                    # Number of peaks
            lenPs = np.zeros(len(diffLens))
                    # Number of troughs
            lenTs = np.zeros(len(diffLens))
                    # Pattern? (True/False = 1.0/0.0)
            patterns = np.zeros(len(diffLens))
                    # Peak positions
            Pss = []
                    # Trough positions
            Tss = []

                # Now, run through the diff len bounds under consideration
            for i in range(len(diffLens)):
                    # Set diff len as the i'th element in the array
                diffLen0 = diffLens[i]

                    # Compute the deconvolved data with the current diff len
                depth0, dataD0 = decon_inst.deconvolve(diffLen0)

                    # Set the minimum prominence as the std of the data times the acceptance percentage defined above
                promT = np.std(dataD0) * acceptPct_prom
                promP = np.std(dataD0) * acceptPct_prom

                    # Interpolate the BD signal, if wanted
                if interpAfterDecon:
                        # If no interpolation size (newDelta) is inputted, choose half the first sample size
                    if newDelta == 0:
                        newDeltaUse = (depth0[1] - depth0[0])/2
                        # Otherwise, use the inputted new sample size (newDelta)
                    else:
                        newDeltaUse = newDelta
                        # Interpolate and define new data to use
                    newDepth, newData, _ = interpCores2(depth0[0], depth0[-1], pd.Series(depth0), pd.Series(dataD0), DeltaInput=True, DeltaIn=newDeltaUse)

                    # If no interpolation wanted, use signal data as new data to use
                else:
                    newDepth = depth0
                    newData = dataD0

                    # Use the method find_constrainedPeaks to estimate peaks, troughs, patterns
                pattern, start, end, Ps, Ts = self.find_constrainedPeaks(newData, min_peakDist, promP)

                # if k%10 == 0:
                #     print(f'Pattern?\t {pattern}')
                #     print(f'N peaks: {len(Ps)}')
                #     print(f'N troughs: {len(Ts)}\n')

                    # Save the results of interest
                Pss.append(Ps)
                Tss.append(Ts)
                lenPs[i] = len(Ps)
                lenTs[i] = len(Ts)
                patterns[i] = pattern

                # Find the first boundary (i.e. diff len) to exceed the expected number of peaks
            firstTrue = np.where(lenPs > Npeaks)[0][0]

                # Set the new diff len min/max boundaries to the [first - 1]/first to exceed Npeaks
            sigMin = diffLens[firstTrue - 1]
            sigMax = diffLens[firstTrue]

                # Add one to the iterator k
            k += 1




            # Save the positions of the final peaks and troughs and pattern
        newPs = Pss[firstTrue - 1]
        newTs = Tss[firstTrue - 1]
        newPattern = patterns[firstTrue - 1]

            # Define the final diff len estimate as the final minimum bound found
        diffLenFin = sigMin

            # Compute the depth and data corresponding to the final diff len estimate
        depth, data = decon_inst.deconvolve(diffLenFin)

            # Interpolate the BD signal, if wanted
        if interpAfterDecon:
                # If no interpolation size (newDelta) is inputted, choose half the first sample size
            if newDelta == 0:
                newDeltaUse = (depth[1] - depth[0])/2
                # Otherwise, use the inputted new sample size (newDelta)
            else:
                newDeltaUse = newDelta
                # Interpolate and define new data to use
            newDepth, newData, _ = interpCores2(depth[0], depth[-1], pd.Series(depth), pd.Series(data), DeltaInput=True, DeltaIn=newDeltaUse)

            # If no interpolation wanted, use signal data as new data to use
        else:
            newDepth = depth
            newData = data

        print(f'Final sigma: {diffLenFin*100:.2f} [cm]')
        print(f'Final # of peaks: {Npeaks}')
        print(f'Delta: {depth[1]-depth[0]:.3f}')
        print(f'Delta new: {newDepth[1]-newDepth[0]:.3f}')

        return newDepth, newData, diffLenFin, newPs, newTs, newPattern#, idxPeak, arr_diffLens, arr_Npeaks, arr_depth, arr_data

    def BackDiffuse_manuel(self, sigma, N = 2000, newDelta = 0.01, interpAfterDecon=True):
        diffLen = sigma
        dInt, d18OInt, Delta = self.interpCores()

        print(f'Sigma input: {diffLen*100:.2f} [cm]')

        decon_inst = SpectralDecon(dInt, d18OInt, N, self.transType)

        depthBD, dataBD = decon_inst.deconvolve(diffLen)

        if interpAfterDecon:
            if newDelta == 0:
                newDeltaUse = (depthBD[1] - depthBD[0])/2
            else:
                newDeltaUse = newDelta
            newDepth, newData, _ = interpCores2(depthBD[0], depthBD[-1], pd.Series(depthBD), pd.Series(dataBD), DeltaInput=True, DeltaIn=newDeltaUse)
            ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
            ave_Npoints = ave_dist/newDeltaUse
            min_peakDist = int(ave_Npoints/self.Dist)

        else:
            newDepth = depthBD; newData = dataBD
            ave_dist = (newDepth[-1] - newDepth[0])/self.ysInSec
            ave_Npoints = ave_dist/newDelta
            min_peakDist = int(ave_Npoints/self.Dist)

        from scipy import signal

        if min_peakDist==0:
            idxPeak = signal.find_peaks(newData)[0]
            idxTrough = signal.find_peaks(-newData)[0]
        else:
            idxPeak = signal.find_peaks(newData, distance=min_peakDist)[0]
            idxTrough = signal.find_peaks(-newData, distance=min_peakDist)[0]

        N_peaks = len(idxPeak)
        N_troughs = len(idxTrough)

        allData = pd.DataFrame({'depth': newDepth, 'd18O': newData})#, 'idxPeaks': idxPeak0, 'Npeaks': N_peaks0})

        return allData, idxPeak, N_peaks, idxTrough, N_troughs, diffLen




    def BackDiffuse_Theo(self, N=2000, newDelta=0.01, interpAfterDecon=True):
        '''


            Arguments:
            ----------

            returns:
            --------

        '''
        sigma_rangeHL = self.diffLenEstimateHL()
        sigma_FitEst = self.spectralEstimate()
        diffLenTheo0 = sigma_rangeHL[0]
        diffLenTheo1 = sigma_rangeHL[1]

        dInt, d18OInt, Delta = self.interpCores()

        print(f'Theo. sigma Min: {sigma_rangeHL[0]*100:.2f} [cm]')
        print(f'Theo. sigma Max: {sigma_rangeHL[1]*100:.2f} [cm]')

        decon_inst = SpectralDecon(dInt, d18OInt, N, self.transType)

        depth0, dataD0 = decon_inst.deconvolve(diffLenTheo0)
        depth1, dataD1 = decon_inst.deconvolve(diffLenTheo1)

        if interpAfterDecon:
            if newDelta == 0:
                newDeltaUse = (depth0[1] - depth0[0])/2
            else:
                newDeltaUse = newDelta
            newDepth0, newData0, _ = interpCores2(depth0[0], depth0[-1], pd.Series(depth0), pd.Series(dataD0), DeltaInput=True, DeltaIn=newDeltaUse)
            ave_dist0 = (newDepth0[-1] - newDepth0[0])/self.ysInSec
            ave_Npoints0 = ave_dist0/newDeltaUse
            min_peakDist0 = int(ave_Npoints0/self.Dist)

            newDepth1, newData1, _ = interpCores2(depth1[0], depth1[-1], pd.Series(depth1), pd.Series(dataD1), DeltaInput=True, DeltaIn=newDeltaUse)
            ave_dist1 = (newDepth1[-1] - newDepth1[0])/self.ysInSec
            ave_Npoints1 = ave_dist1/newDeltaUse
            min_peakDist1 = int(ave_Npoints1/self.Dist)

        else:
            newDepth0 = depth0; newData0 = dataD0
            ave_dist0 = (newDepth0[-1] - newDepth0[0])/self.ysInSec
            ave_Npoints0 = ave_dist0/newDelta0
            min_peakDist0 = int(ave_Npoints0/self.Dist)

            newDepth1 = depth1; newData1 = dataD1
            ave_dist1 = (newDepth1[-1] - newDepth1[0])/self.ysInSec
            ave_Npoints1 = ave_dist1/newDelta
            min_peakDist1 = int(ave_Npoints1/self.Dist)

        from scipy import signal

        if min_peakDist0==0:
            idxPeak0 = signal.find_peaks(newData0)[0]
            idxTrough0 = signal.find_peaks(-newData0)[0]
        else:
            idxPeak0 = signal.find_peaks(newData0, distance=min_peakDist0)[0]
            idxTrough0 = signal.find_peaks(-newData0, distance=min_peakDist0)[0]

        if min_peakDist1==0:
            idxPeak1 = signal.find_peaks(newData1)[0]
            idxTrough1 = signal.find_peaks(-newData1)[0]
        else:
            idxPeak1 = signal.find_peaks(newData1, distance=min_peakDist1)[0]
            idxTrough1 = signal.find_peaks(-newData1, distance=min_peakDist1)[0]

        N_peaks0 = len(idxPeak0)
        N_peaks1 = len(idxPeak1)

        N_troughs0 = len(idxTrough0)
        N_troughs1 = len(idxTrough1)

        allData_TheoMin = pd.DataFrame({'depth': newDepth0, 'd18O': newData0})#, 'idxPeaks': idxPeak0, 'Npeaks': N_peaks0})
        allData_TheoMax = pd.DataFrame({'depth': newDepth1, 'd18O': newData1})#, 'idxPeaks': idxPeak1, 'Npeaks': N_peaks1})

        return allData_TheoMin, idxPeak0, N_peaks0, idxTrough0, N_troughs0, diffLenTheo0, allData_TheoMax, idxPeak1, N_peaks1, idxTrough1, N_troughs1, diffLenTheo1
# #temp_holo = 213.15, delta_holo = -51, slope = 0.69
#
#     def DeltaToTemp(self, temp_holo = 213.15, delta_holo = -51, slope = 0.69):
#         '''
#             Method to estimate temperature in [k] given a d18O depth series, along with specific
#             temp and d18O estimates for holocene.
#
#
#             Arguments:
#             ----------
#                 temp_holo:      [float] Estimated holocene mean temperature.
#                 delta_holo:     [float] Estimated holocene mean d18O .
#                 slope:          [float] Estimated slope.
#
#             returns:
#             --------
#                 depth:          [arr of floats] Depth data.
#                 temp:           [arr of floats] Estimated temperature data.
#         '''
#         depth, data,_,_ = self.backDiffused()
#         temp = np.zeros(np.size(data)) + temp_holo + slope*(data - delta_holo)
#
#         return depth, temp
#

    def interpDiffData():
        '''


            Arguments:
            ----------

            returns:
            --------

        '''
        return

    def find_constrainedPeaks(self, data, dist, prom):

            # Find peaks and troughs in data, given distance and prominence provided
        peaksBD = sp.signal.find_peaks(data, distance = dist, prominence=prom)[0]
        troughsBD = sp.signal.find_peaks(-data, distance = dist, prominence=prom)[0]


            # Create lists of ones (representing peaks) and zeros (representing troughs)
        peaksBD_lst = np.ones(len(peaksBD))
        troughsBD_lst = np.zeros(len(troughsBD))


            # Create array containing (unsorted) peaks and troughs positions and corresponding P (1) or T (0).
        exts = np.concatenate((peaksBD,troughsBD))
        exts_lst = np.concatenate((peaksBD_lst,troughsBD_lst))
            # Sort both lists at the same time, thus matching up idxs and P/T (1/0) value
        list1, list2 = (np.array(t) for t in zip(*sorted(zip(exts, exts_lst))))


            # Check if list is divisible by two. If not, append P/T array with -1 and idx array with 0
            # at either first or last position depending on starting with P (1) or T (0)
        if len(list1)%2 != 0:
            if list2[0] == 1:
                listNew_lst = np.append(list2,-1)
                listNew = np.append(list1,0)
            elif list2[0] == 0:
                listNew_lst = np.insert(list2,0,-1)
                listNew = np.insert(list1,0,0)
        else:
            listNew_lst = list2
            listNew = list1


            # Reshape arrays into 2 x len(Ps/Ts)
        PTs = listNew_lst.reshape((int(len(listNew_lst)/2)),2)
        PTs_idx = listNew.reshape((int(len(listNew)/2)),2)


            # Get only positive values from list (remove the appended values from earlier)
        listNew_lstPos = listNew_lst[listNew_lst>=0]
        listNew_Pos = listNew[listNew_lst>=0]


            # If the first value in the array is a peak (1), make lists with all (estimated) peaks (every second
            # element) and all troughs (every second + 1 element)
        if listNew_lstPos[0] == 1:

            listNew_lstPs = listNew_lstPos[::2]
            listNew_Ps = listNew_Pos[::2]

            listNew_lstTs = listNew_lstPos[1::2]
            listNew_Ts = listNew_Pos[1::2]

            # If the first value in the array is a trough (0), make list with all (estimated) peaks (every second
            # + 1 element) and all troughs (every second element)
        elif listNew_lstPos[0] == 0:
            listNew_lstPs = listNew_lstPos[1::2]
            listNew_Ps = listNew_Pos[1::2]

            listNew_lstTs = listNew_lstPos[::2]
            listNew_Ts = listNew_Pos[::2]


            # Check if [..PTPTPTPTPT..] pattern exists by summing all estimated peaks and checking if sum(P) == len(Ps)
            # and summing all estimated troughs and checking if sum(T) == 0
            # If so, return pattern = True and what the pattern starts and ends with (P=1 and T=0)

        if (sum(listNew_lstPs) == len(peaksBD)) & (sum(listNew_lstTs) == 0):
            pattern = True
            patternStart = listNew_lstPos[0]
            patternEnd = listNew_lstPos[-1]
            # If not, return pattern = False and start and end as empty lists.
        else:
            pattern = False
            patternStart = []
            patternEnd = []

        return pattern, patternStart, patternEnd, listNew_Ps, listNew_Ts



def interpCores2(valMin, valMax, d_In, x_In, pad = 1, DeltaInput = False, DeltaIn = 0, interpAll=True):
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

    d_in = d_In
    x_in = x_In


    if interpAll:
        valMin = d_in.min()
        valmax = d_in.max()
    else:
        valMin = d_in.min - pad
        valMax = d_in.max + pad

    d = d_in[(d_in >= valMin) & (d_in <= valMax)]
    x = x_in[(d_in >= valMin) & (d_in <= valMax)]

    if DeltaInput:
        Delta = DeltaIn
    else:
        diff = np.diff(d)
        Delta = round(min(diff), 3)

    d_min = Delta * np.ceil(d.values[0]/Delta)
    d_max = Delta * np.floor(d.values[-1]/Delta)

    n = int(1 + (d_max - d_min)/Delta)

    j_arr = np.linspace(0,n,n)
    dhat0 = d_min + (j_arr - 1)*Delta

    f = interpolate.CubicSpline(d,x)

    xhat0 = f(dhat0)

    dhat = dhat0[(dhat0 >= d_min) & (dhat0 <= d_max)]
    xhat = xhat0[(dhat0 >= d_min) & (dhat0 <= d_max)]

    return dhat, xhat, Delta






mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{textcomp}',
    r'\usepackage{wasysym}']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'STIXGeneral'

#
#     # Core names of cores available
# coreNames = ['Crete', 'SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
#
#     # Selecting core name
# coreName = 'SiteB'
#
#     # Reading datafiles for specific core
# d18OData = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+coreName+'_det.txt',',')
# densities = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/'+coreName+'DepthDens_w_Models.txt','\t')
# diffLens = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/'+coreName+'_DepthDiff.txt','\t')
# specsCores = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt',',')
# specIdx = specsCores['CoreName'][specsCores['CoreName'] == coreName].index[0]
#
#     # Set the specs for depth of Laki and Tambora eruptions for core
# specsCore = specsCores.iloc[specIdx]
# dTamb = np.float64(specsCore['dTamb'])
# dLaki = np.float64(specsCore['dLaki'])
#
#     # (FOR PLOTTING) Make array of only d18O data between Laki and Tamb
# depth_LT = d18OData['depth'][(d18OData['depth'] >= dTamb) & (d18OData['depth'] <= dLaki)]
# d18O_LT = d18OData['d18O'][(d18OData['depth'] >= dTamb) & (d18OData['depth'] <= dLaki)]
#
#     # Create instance of back diffusion
# inst = BackDiffuse(coreName, d18OData, specsCores, dTamb, dLaki, 34, diffLenData=diffLens[['Depth','sigma_o18']], densData=densities)
#
#     # Make spectral estimate of diff len
# diffLen = inst.spectralEstimate()
#
#     # Make model/empiric estimate of diff len
# difflenEstHL = inst.diffLenEstimateHL()
#
#     # Compute final depth/d18O back diffused data w. final diff len and No. peaks
# depth, data, diffLen, peaks = inst.backDiffused()
#
#     # Plot original data, back diffused data and peak estimations
# fig, ax = plt.subplots(figsize=(10,7))
# ax.plot(depth, data, lw=1, label='Back diffused')
# ax.plot(depth_LT, d18O_LT-np.mean(d18O_LT),color='k', lw=1, label = 'Data')
# ax.plot(depth[peaks], data[peaks],'.',lw=1, label='Estimated peaks')
# ax.set(xlabel = 'Depth [m]', ylabel = '$\delta^{18}$O [\permil]', title=coreName+'$, \sigma_{fin} =$ ' + f'{diffLen*100:.2f} [cm]')
# ax.legend(fontsize=16)
# fig.tight_layout()
# fig.savefig(coreName + '_peaks.jpg')
#
#
#     # Compute temperature estimate in [K]
# depthT, dataT = inst.DeltaToTemp()
#
#     # Plot Temperature data in [C]
# fig2, ax2 = plt.subplots(figsize=(10,7))
# ax2.plot(depthT, dataT-273.15,color='k', lw=1)
# ax2.set(xlabel = 'Depth [m]', ylabel = 'Temperature [C]', title=coreName)
# fig2.savefig(coreName + '_Temp.jpg')

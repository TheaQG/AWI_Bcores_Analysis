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
from scipy.interpolate import interp1d, CubicSpline
from Decon import SpectralDecon

'''
        *********** TODO************

        -
'''

class Interpolation():
    '''
        Methods available:
            __init__(self, x_data, y_data, interval, interpType, DeltaInput=False, samplingSize=1, GapInput=False, interpTypeGap='None',  gapInterval=np.array([1,2])):
                    Initializes class with 4 required arguments and 5 optional arguments.

            __call__(self):
                    Call function interpolates data as-is with the given interpolation type passed.

            interpolateData(self, xMin_in, xMax_in, x_data, y_data, DeltaIn, pad=False, valPad=0):
                    Interpolate all data (datax, datay) passed to method in interval [xMin:xMax].
                    Uses passed samplingsize for interpolated resampling.

            interpolateGap(self, getData_BelAbo=False):
                    Uses interpolateData method to interpolate passed gap in data.

            interpolateDataAndGap(self, xlabel='depth', ylabel='d18O'):
                    Uses above methods to give a resampled signal of both data and gap.
                    Uses the passed resampling interval.
    '''

    def __init__(self, x_data, y_data, interval, interpType, DeltaInput=False, samplingSize=1, GapInput=False, interpTypeGap='None',  gapInterval=np.array([1,2])):
        '''
            Initialize the class instance.


            Arguments:
            ----------
                x_data:         [pd.Series] Measured (uneven) x data. Must match y_data.
                y_data:         [pd.Series] Measured (uneven) x data. Must match x_data.
                interval:       [np.array] Interval to interpolate/resample between.
                interpType:     [str] 'Linear' or 'CubicSpline'. Defines interpolation type.
                DeltaInput:     [bool] Default: False. Per default uses resampling size of
                                    min. sampling size from data. Else uses passed samplingSize.
                samplingSize    [float] Default: 1. If DeltaInput = True, use this resampling size.
                GapInput:       [bool] Default: False. If True, then data has a gap.
                interpTypeGap   [str] Default: 'None'. 'Linear' or 'CubicSpline'. If GapInput=True
                                    what interpolation type to use.
                gapInterval     [np.array] Default: np.array([1,2]). Gap placement in x data.

            returns:
            --------
                None

        '''
        self.x_data = x_data
        self.y_data = y_data
        self.xMin = interval[0]
        self.xMax = interval[1]
        self.interpType = interpType
        self.DeltaInput = DeltaInput
#        self.samplingSize = samplingSize
        self.GapInput = GapInput
        self.interpTypeGap = interpTypeGap
        self.gapMin = gapInterval[0]
        self.gapMax = gapInterval[1]


        if self.DeltaInput:
            self.samplingSize = samplingSize

        else:
            diff = np.diff(self.x_data)
            Delta = round(min(diff), 3)
            self.samplingSize = Delta


        return

    def __call__(self):
        '''


            Arguments:
            ----------
                None

            returns:
            --------
                xHat:       [np.array] Interpolated/resampled x data.
                yHat:       [np.array] Interpolated/resampled y data.
                Delta:      [float] Resampling interval used for interoaltion.
        '''

        xHat, yHat, Delta = self.interpolateData(xMin_in=self.xMin, xMax_in=self.xMax, x_data=self.x_data, y_data=self.y_data, DeltaIn = self.samplingSize, pad=False, valPad=0)

        return xHat, yHat, Delta


    def interpolateData(self, xMin_in, xMax_in, x_data, y_data, DeltaIn, pad=False, valPad=0):
        '''
            Computes interpolation btw. xMin and xMax, with optional padding.
            Padding is added to avoid undesired border effects.
            Uses cubic spline or linear interpolation with a new sample size given
            as DeltaIn, passed when initializing class.


            Arguments:
            ----------
                pad:            [float] Padding to avoid werid borders. Default = 1 [m].

            returns:
            --------
                dhat:           [array of floats] Depth data corresponding to interpolated data.
                xhat:           x_dataBelGap_int, y_dataBelGap_int, DeltaBelGap = [array of floats] d18O data corresponding to interpolated data.
                Delta:          [float] New sample size.

        '''

            # Create np.arrays instead of pd.Series
        x_arr = np.asarray(x_data)
        y_arr = np.asarray(y_data)

            # Define x interval to interpolate over - with or without padding.
        if not pad:
            xMin = xMin_in
            xMax = xMax_in
        else:
            xMin = xMin_in - valPad
            xMax = xMax_in + valPad

            # Define actual data to use for interpolation as data in interval.
        x = x_arr[(x_arr >= xMin) & (x_arr <= xMax)]
        y = y_arr[(x_arr >= xMin) & (x_arr <= xMax)]

            # Define sampling size as Delta, given in initialization.
        Delta = self.samplingSize

            # Define new interval, given the resampling size.
        xHat_min = Delta * np.ceil(x[0]/Delta)
        xHat_max = Delta * np.floor(x[-1]/Delta)

            # Calculate the new number of points in resampled signal.
        n = int(1 + (xHat_max - xHat_min)/Delta)

            # Compute the placement of the evenly spaced new data points.
        j_arr = np.linspace(1,n+1,n+1)
        xHat0 = xHat_min + (j_arr - 1)*Delta

            # Define interpolation function
        if self.interpType == 'CubicSpline':
            f = interpolate.CubicSpline(x,y)
        elif self.interpType == 'Linear':
            f = interpolate.interp1d(x,y, kind='linear')

            # Compute new, resampled y-data from interpolation.
        yHat0 = f(xHat0)

            # Define final interpolated/resampled data to be inside new interval.
        xHat = xHat0[(xHat0 >= xHat_min) & (xHat0 <= xHat_max)]
        yHat = yHat0[(xHat0 >= xHat_min) & (xHat0 <= xHat_max)]


        return xHat, yHat, Delta

    def interpolateGap(self, getData_BelAbo=False):
        '''
            Uses method interpolateData to interpolate across the defined gap in data.

            (For use in interpolateDataAndGap: If getData_BelAbo = True, the method
            saves the resampled data from above and below gap in self.dataAboGap_int
            and self.dataBelGap_int.)


            Arguments:
            ----------
                getData_BelAbo:     [bool] Default: False. Save interpolated data above
                                        and below gap in class self.?

            returns:
            --------
                xGap:               [np.array] Interpolated/resampled x data in gap.
                yGap:               [np.array] Interpolated/resampled y data in gap.
                DeltaAll:           [float] Resampling interval.

        '''

            # Define length and placement of gap (and data)
        gapMin = self.gapMin
        gapMax = self.gapMax
        gapL = gapMax - gapMin
        dataL = self.xMax - self.xMin

            # Define data as data in interval passed.
        x_data = self.x_data[(self.x_data >= self.xMin) & (self.x_data <= self.xMax)]
        y_data = self.y_data[(self.x_data >= self.xMin) & (self.x_data <= self.xMax)]

            # Define what data is 'below'(left of) gap and what is 'above'(right of) gap
        x_dataBelGap = x_data[x_data <= gapMin]
        y_dataBelGap = y_data[x_data <= gapMin]

        x_dataAboGap = x_data[x_data >= gapMax]
        y_dataAboGap = y_data[x_data >= gapMax]

            # Find sampling intervals in data above and below gap.
        diffsBel = np.diff(x_dataBelGap)
        diffsAbo = np.diff(x_dataAboGap)

            # Define the sampling interval to use as the minimum of all samplings used.
        diffMin = np.min([np.min(diffsBel), np.min(diffsAbo)])

            # Compute how many points to resample in entire data set.
            # Calculate resampling interval based on this.
        N_pointsAll = np.floor(dataL/diffMin)
        DeltaAll = dataL/N_pointsAll

            # Use above method to interpolate data below and above gap.
        x_dataBelGap_int, y_dataBelGap_int, DeltaBelGap = self.interpolateData(xMin_in = min(x_dataBelGap), xMax_in = max(x_dataBelGap), x_data=x_dataBelGap, y_data=y_dataBelGap, DeltaIn = DeltaAll)
        x_dataAboGap_int, y_dataAboGap_int, DeltaAboGap = self.interpolateData(xMin_in = min(x_dataAboGap), xMax_in = max(x_dataAboGap), x_data=x_dataAboGap, y_data=y_dataAboGap, DeltaIn = DeltaAll)

            # If wanted, save interpolated data above and below gap.
        if getData_BelAbo:
            self.x_dataBelGap_int = x_dataBelGap_int
            self.y_dataBelGap_int = y_dataBelGap_int

            self.x_dataAboGap_int = x_dataAboGap_int
            self.y_dataAboGap_int = y_dataAboGap_int

            # Calculate placement of data points to be in gap.
        xGap = np.arange(x_dataBelGap_int[-1] + DeltaAll, x_dataAboGap_int[0], DeltaAll)

            # Define what interpolation type to use.
        interpTypeGap = self.interpTypeGap

            # Define interpolation function and compute y values for data points.
        if interpTypeGap == 'Linear':
                # If 'Linear', use only one data point on each side of the gap.
            f_interp = interp1d([x_dataBelGap_int[-1],x_dataAboGap_int[0]], [y_dataBelGap_int[-1], y_dataAboGap_int[0]])
            yGap = f_interp(xGap)
        elif interpTypeGap == 'CubicSpline':
                # If 'CubicSpline' use ?? data points on each side of the gap.
            xDataCS = np.concatenate((x_dataBelGap_int,x_dataAboGap_int), axis=0)
            yDataCS = np.concatenate((y_dataBelGap_int,y_dataAboGap_int), axis=0)

            f_interp = CubicSpline(xDataCS, yDataCS)
            yGap = f_interp(xGap)


        return xGap, yGap, DeltaAll


    def interpolateDataAndGap(self, xlabel='depth', ylabel='d18O'):
        '''


            Arguments:
            ----------
                xlabel:         [str] Default: 'depth'. Column name for final dataframe.
                ylabel:         [str] Default: 'd18O'. Column name for final dataframe.

            returns:
            --------
                dataComb:       [pd.DataFrame] All interpolated data, gap and data combined .
                dataGap:        [pd.DataFrame] Only data in gap.
                dataBelGap:     [pd.DataFrame] Only data below gap.
                dataAboGap:     [pd.DataFrame] Only data above gap.

        '''
            # Compute interolated/resampled data in gap. Save computed data
            # above and below gap.
        xGap, yGap, DeltaAll = self.interpolateGap(getData_BelAbo=True)

            # Define x and y data from below and above gap.
        xBelGap_int = self.x_dataBelGap_int
        yBelGap_int = self.y_dataBelGap_int

        xAboGap_int = self.x_dataAboGap_int
        yAboGap_int = self.y_dataAboGap_int

            # Combine below gap, gap, and above gap to one array.
        xAll = np.concatenate((xBelGap_int, xGap, xAboGap_int))
        yAll = np.concatenate((yBelGap_int, yGap, yAboGap_int))

            # Save as dataframes with passed x- and y-labels.
        dataBelGap = pd.DataFrame({xlabel:xBelGap_int, ylabel:yBelGap_int}, index=None)
        dataAboGap = pd.DataFrame({xlabel:xAboGap_int, ylabel:yAboGap_int}, index=None)
        dataGap = pd.DataFrame({xlabel:xGap, ylabel:yGap}, index=None)
        dataComb = pd.DataFrame({xlabel:xAll, ylabel:yAll}, index=None)

        return dataComb, dataGap, dataBelGap, dataAboGap

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


class Interpolation():
    '''
        Methods available:
            __init__(self, x_data, y_data, interval, interpType, DeltaInput=False, samplingSize=1, GapInput=False, interpTypeGap='None',  gapInterval=np.array([1,2])):

            __call__(self):

            interpolateData(self, xMin_in, xMax_in, x_data, y_data, DeltaIn, pad=False, valPad=0):

            interpolateGap(self, getData_BelAbo=False):

            interpolateDataAndGap(self, xlabel='depth', ylabel='d18O'):
    '''

    def __init__(self, x_data, y_data, interval, interpType, DeltaInput=False, samplingSize=1, GapInput=False, interpTypeGap='None',  gapInterval=np.array([1,2])):
        '''
            Initialize the class instance.


            Arguments:
            ----------
                x_data:         []
                y_data:         []
                interval:       []
                interpType:     []
                DeltaInput:     [] False
                samplingSize    [] =1
                GapInput:       []=False
                interpTypeGap   []='None'
                gapInterval     []=np.array([1,2])

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
        self.samplingSize = samplingSize
        self.GapInput = GapInput
        self.interpTypeGap = interpTypeGap
        self.gapMin = gapInterval[0]
        self.gapMax = gapInterval[1]


        if self.DeltaInput:
            Delta = self.samplingSize
            self.Delta = Delta
        else:
            diff = np.diff(self.x_data)
            Delta = round(min(diff), 3)
            self.Delta = Delta


        return

    def __call__(self):
        '''


            Arguments:
            ----------
                None

            returns:
            --------
                xHat:       []
                yHat:       []
                Delta:      []
        '''

        xHat, yHat, Delta = self.interpolateData(xMin_in=self.xMin, xMax_in=self.xMax, x_data=self.x_data, y_data=self.y_data, DeltaIn = self.samplingSize, pad=False, valPad=0)

        return xHat, yHat, Delta


    def interpolateData(self, xMin_in, xMax_in, x_data, y_data, DeltaIn, pad=False, valPad=0):
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
                xhat:           x_dataBelGap_int, y_dataBelGap_int, DeltaBelGap = [array of floats] d18O data corresponding to interpolated data.
                Delta:          [float] New sample size.

        '''

        x_arr = np.asarray(x_data)
        y_arr = np.asarray(y_data)

        if not pad:
            xMin = xMin_in
            xMax = xMax_in
        else:
            xMin = xMin_in - valPad
            xMax = xMax_in + valPad

        x = x_arr[(x_arr >= xMin) & (x_arr <= xMax)]
        y = y_arr[(x_arr >= xMin) & (x_arr <= xMax)]

        Delta = self.Delta

        xHat_min = Delta * np.ceil(x[0]/Delta)
        xHat_max = Delta * np.floor(x[-1]/Delta)

        n = int(1 + (xHat_max - xHat_min)/Delta)

        j_arr = np.linspace(1,n+1,n+1)
        xHat0 = xHat_min + (j_arr - 1)*Delta


        if self.interpType == 'CubicSpline':
            f = interpolate.CubicSpline(x,y)
        elif self.interpType == 'Linear':
            f = interpolate.interp1d(x,y, kind='linear')

        yHat0 = f(xHat0)

        xHat = xHat0[(xHat0 >= xHat_min) & (xHat0 <= xHat_max)]
        yHat = yHat0[(xHat0 >= xHat_min) & (xHat0 <= xHat_max)]


        return xHat, yHat, Delta

    def interpolateGap(self, getData_BelAbo=False):
        '''


            Arguments:
            ----------
                getData_BelAbo:     [] = False

            returns:
            --------
                xGap:               []
                yGap:               []
                DeltaAll:           []

        '''
        gapMin = self.gapMin
        gapMax = self.gapMax
        gapL = gapMax - gapMin
        dataL = self.xMax - self.xMin

        x_data = self.x_data[(self.x_data >= self.xMin) & (self.x_data <= self.xMax)]
        y_data = self.y_data[(self.x_data >= self.xMin) & (self.x_data <= self.xMax)]


        x_dataBelGap = x_data[x_data <= gapMin]
        y_dataBelGap = y_data[x_data <= gapMin]

        x_dataAboGap = x_data[x_data >= gapMax]
        y_dataAboGap = y_data[x_data >= gapMax]

        diffsBel = np.diff(x_dataBelGap)
        diffsAbo = np.diff(x_dataAboGap)

        diffMin = np.min([np.min(diffsBel), np.min(diffsAbo)])

        N_pointsAll = np.floor(dataL/diffMin)
        DeltaAll = dataL/N_pointsAll


        x_dataBelGap_int, y_dataBelGap_int, DeltaBelGap = self.interpolateData(xMin_in = min(x_dataBelGap), xMax_in = max(x_dataBelGap), x_data=x_dataBelGap, y_data=y_dataBelGap, DeltaIn = DeltaAll)
        x_dataAboGap_int, y_dataAboGap_int, DeltaAboGap = self.interpolateData(xMin_in = min(x_dataAboGap), xMax_in = max(x_dataAboGap), x_data=x_dataAboGap, y_data=y_dataAboGap, DeltaIn = DeltaAll)

        if getData_BelAbo:
            self.x_dataBelGap_int = x_dataBelGap_int
            self.y_dataBelGap_int = y_dataBelGap_int

            self.x_dataAboGap_int = x_dataAboGap_int
            self.y_dataAboGap_int = y_dataAboGap_int

        xGap = np.arange(x_dataBelGap_int[-1] + DeltaAll, x_dataAboGap_int[0], DeltaAll)

        interpTypeGap = self.interpTypeGap

        if interpTypeGap == 'Linear':
            f_interp = interp1d([x_dataBelGap_int[-1],x_dataAboGap_int[0]], [y_dataBelGap_int[-1], y_dataAboGap_int[0]])
            yGap = f_interp(xGap)
        elif interpTypeGap == 'CubicSpline':
            xDataCS = np.concatenate((x_dataBelGap_int[-11:],x_dataAboGap_int[:10]), axis=0)
            yDataCS = np.concatenate((y_dataBelGap_int[-11:],y_dataAboGap_int[:10]), axis=0)

            f_interp = CubicSpline(xDataCS, yDataCS)
            yGap = f_interp(xGap)


        return xGap, yGap, DeltaAll


    def interpolateDataAndGap(self, xlabel='depth', ylabel='d18O'):
        '''


            Arguments:
            ----------
                xlabel:         []='depth'
                ylabel:         []='d18O'

            returns:
            --------
                dataComb:       []
                dataGap:        []
                dataBelGap:     []
                dataAboGap:     []
                
        '''
        xGap, yGap, DeltaAll = self.interpolateGap(getData_BelAbo=True)

        xBelGap_int = self.x_dataBelGap_int
        yBelGap_int = self.y_dataBelGap_int

        xAboGap_int = self.x_dataAboGap_int
        yAboGap_int = self.y_dataAboGap_int

        xAll = np.concatenate((xBelGap_int, xGap, xAboGap_int))
        yAll = np.concatenate((yBelGap_int, yGap, yAboGap_int))

        dataBelGap = pd.DataFrame({xlabel:xBelGap_int, ylabel:yBelGap_int}, index=None)
        dataAboGap = pd.DataFrame({xlabel:xAboGap_int, ylabel:yAboGap_int}, index=None)
        dataGap = pd.DataFrame({xlabel:xGap, ylabel:yGap}, index=None)
        dataComb = pd.DataFrame({xlabel:xAll, ylabel:yAll}, index=None)

        return dataComb, dataGap, dataBelGap, dataAboGap

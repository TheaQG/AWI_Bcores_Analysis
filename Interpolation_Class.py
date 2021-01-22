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


class Interpolations():
    '''
        Methods available:
            __init__():

    '''

    def __init__(self, x_data, y_data, interval, interpType, DeltaInput=False, samplingSize=1, GapInput=False, interpTypeGap='None',  gapInterval=np.array([1,2]), interpWithGap=False):

        self.x_data = x_data
        self.y_data = y_data
        self.xMin = interval[0]
        self.xMax = interval[1]
        self.interpType = interpType
        self.DeltaInput = DeltaInput
        self.samplingSize = samplingSize
        self.GapInput = GapInput
        self.gapMin = gapInterval[0]
        self.gapMax = gapInterval[1]
        self.interpWithGap = interpWithGap

        return

    def interpolateData(self, pad=False, valPad=0):
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

        x_arr = self.x_data
        y_arr = self.y_data

        if not pad:
            xMin = self.xMin
            xMax = self.xMax
        else:
            xMin = self.xMin - valPad
            xMax = self.xMax + valPad

        x = x_arr[(x_arr >= xMin) & (x_arr <= xMax)]
        y = y_arr[(x_arr >= xMin) & (x_arr <= xMax)]

        if self.DeltaInput:
            Delta = self.samplingSize
        else:
            diff = np.diff(x)
            Delta = round(min(diff), 3)

        xHat_min = Delta * np.ceil(x.values[0]/Delta)
        xHat_max = Delta * np.floor(x.values[-1]/Delta)

        n = int(1 + (xHat_max - xHat_min)/Delta)

        j_arr = np.linspace(0,n,n)
        xhat0 = xHat_min + (j_arr - 1)*Delta

        if self.interpType == 'CubicSpline':
            f = interpolate.CubicSpline(x,y)
        elif self.interpType == 'Linear':
            f = interpolate.interp1d(x,y, kind='linear')

        yHat0 = f(xHat0)

        xHat = xHat0[(xHat0 >= xHat_min) & (xHat0 <= xHat_max)]
        yHat = yHat0[(xHat0 >= xHat_min) & (xHat0 <= xH_max)]


        return xHat, yHat, Delta

    def interpolateGap(self):
        gapMin = self.gapMin
        gapMax = self.gapMax
        gapLen = gapMax - gapMin
        return

    def interpolateDataAndGap(self):

        return

    def

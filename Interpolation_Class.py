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

    def interpolateData(self):

        return x_int, y_int

    def interpolateGap(self):

        return

    def interpolateDataAndGap(self):

        return

    def

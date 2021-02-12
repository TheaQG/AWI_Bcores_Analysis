import numpy as np
import pandas as pd
from pandas import ExcelWriter
import matplotlib.pyplot as plt
import openpyxl
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



import sys
import os
sys.path.append('../')

from BackDiffuse_LT import BackDiffuse
from GetCoreData_fct import GetCoreData
from Interpolation_Class import Interpolation



'''
        *********** TODO************
        - (12/02/21) Peak Detection Methods
        - (12/02/21) Linear Timescale
        - (12/02/21) Detrending/Standardization

'''

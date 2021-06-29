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
from BackDiffuse_LT import BackDiffuse


'''
    ESTIMATING ANNUAL LAYER THICKNESS
'''
"""
    - Amplitude attenuation
"""
"""
    - Annual Layer Thicknes
"""


################################################################################

'''
    FINAL SIGMA ESTIMATES
'''
"""
    - "First" guess ...
"""
sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']

for i in range(len(sites)):
    site = sites[i]
    print('\n##########'+site+'##########\n')
    N_InInt = 33

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    N = 500
    diffLens = np.zeros(N)
    dTambs = np.zeros(N)
    dLakis = np.zeros(N)

    for i in range(N):
        print(i)
        diffLens[i], dTambs[i], dLakis[i] = Calc_diffLen_Gauss(site, 33, CoresSpecs)

    np.savetxt(site+'diffLens_GaussDistwDepths.csv', np.array([diffLens,dTambs,dLakis]))



"""
    - Spectral transform effects (FFT/DCT/NDCT)
        - Sigma estimate
        - Speed
"""
"""
    - Sig_const v. sig(z)
        - Sigma estimate
        - Speed
"""
"""
    - No constraints/Constraints
        - Sigma estimate
        - Speed
"""
"""
    - Location distributions
"""

################################################################################

'''
    TEMPERATURE ESTIMATES
'''
"""
    - Location distributions
"""
"""
    - Steady state solutions
        - Accumulation distribution
"""
"""
    - Non steady-state??
"""

################################################################################

'''
    SENSITIVITY TESTING
'''
"""
    - Interpolation effect
        - Before deconvolution
        - After deconvolution
"""
"""
    - N peaks v. sigma
        - Pattern/No pattern
"""

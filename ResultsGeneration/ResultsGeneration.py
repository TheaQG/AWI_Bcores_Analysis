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

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{textcomp}',
    r'\usepackage{wasysym}']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'STIXGeneral'

saveFigs = False

import sys
import os

sys.path.append('../../')
sys.path.append('../')

from GetCoreData_fct import GetCoreData
from BackDiffuse_LT import BackDiffuse
from Interpolation_Class import Interpolation
from HL_AnalyticThea_class import HL_Thea
from DiffusionProfiles_calculations import DiffusionLength
from transforms import transforms
from Decon import SpectralDecon
from sigmaSolver import sigma_Solver
from SignalAttenuation import Attenuation, AnnualLayerThick


import Results_LTLocationDistribution
import Results_diffLenVPeaks
import Results_EffectsSigma_SpectralTransforms
import Results_EffectsSigma_Interpolation


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

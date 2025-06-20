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

from src.GetCoreData_fct import GetCoreData
from src.BackDiffuse_LT import BackDiffuse
from src.Interpolation_Class import Interpolation
from src.HL_AnalyticThea_class import HL_Thea
from src.DiffusionProfiles_calculations import DiffusionLength
from src.transforms import transforms
from src.Decon import SpectralDecon
from src.sigmaSolver import sigma_Solver
from src.SignalAttenuation import Attenuation, AnnualLayerThick


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

sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']

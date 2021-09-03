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
from DiffLen_UncertaintyEst import Calc_diffLen_Gauss, Calc_diffLen_Gauss_1const, Calc_diffLen_Gauss_MonthVar

'''
    SENSITIVITY TESTING
'''
"""
    - Location distributions
"""
sites = ['Crete']#['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'

print('\t\t####################################')
print('\t\t### ESTIMATING DIFF LEN (GAUSS) ####')
print('\t\t##### VARIATION IN L & T BOTH ######')
print('\t\t####################################')
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

    np.savetxt(pathResults + site+'_diffLens_GaussDistwDepths.csv', np.array([diffLens,dTambs,dLakis]))

print('\t\t####################################')
print('\t\t### ESTIMATING DIFF LEN (GAUSS) ####')
print('\t\t###### VARIATION IN LAKI ONLY ######')
print('\t\t####################################')
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
        diffLens[i], dTambs[i], dLakis[i] = Calc_diffLen_Gauss_1const(site, 33, CoresSpecs, eruption='Laki')

    np.savetxt(pathResults + site+'_diffLens_GaussDistwDepths_Laki.csv', np.array([diffLens,dTambs,dLakis]))






print('\t\t####################################')
print('\t\t### ESTIMATING DIFF LEN (GAUSS) ####')
print('\t\t###### VARIATION IN TAMB ONLY ######')
print('\t\t####################################')
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
        diffLens[i], dTambs[i], dLakis[i] = Calc_diffLen_Gauss_1const(site, N_InInt, CoresSpecs, eruption='Tambora')

    np.savetxt(pathResults + site+'_diffLens_GaussDistwDepths_Tamb.csv', np.array([diffLens,dTambs,dLakis]))









print('\t\t####################################')
print('\t\t### ESTIMATING DIFF LEN (GAUSS) ####')
print('\t\t###### 2 MONTH VARIATION ONLY ######')
print('\t\t####################################')
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
        diffLens[i], dTambs[i], dLakis[i] = Calc_diffLen_Gauss_MonthVar(site, N_InInt, CoresSpecs, Nmonths=2, transType_in='NDCT')

    np.savetxt(pathResults + site+'_diffLens_GaussDistwDepths_sigNmonths2_NDCT2.csv', np.array([diffLens,dTambs,dLakis]))

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
import csv
from itertools import zip_longest

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

'''
    EFFECTS ON SIGMA ESTIMATES
'''
"""
    - No constraints/Constraints
        - Sigma estimate
            - Qualitative analysis for visual inspection
"""
sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'

for i in range(len(sites)):

    site = sites[i]
    print('\n###### ' + site + '######\n')

    N_InInt = 33

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']


    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']
    accum0 = CoreSpecs['Accum0']
    accumIE = CoreSpecs['Accum1']
    Temp0 = CoreSpecs['T0']

    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
    data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
    data_dens = DataAll[4]; data_dens_LT = DataAll[5]
    data_diff = DataAll[6]; data_diff_LT = DataAll[7]


    depth = data_d18O['depth']
    d18O = data_d18O['d18O']

    depth_LT = data_d18O_LT['depth']
    d18O_LT = data_d18O_LT['d18O']

    dataAll = pd.DataFrame({'depth':depth,'d18O':d18O}, index=None)

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=5, transType = 'DCT')

    print('Generating results w/o constraints\n')
    depthEst, dataEst, diffLenFin, idxPeak, arr_diffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(print_Npeaks=False)


    header = ['depth','data','idxPeaks','sigma']
    rows = zip_longest(depthEst, dataEst, idxPeak, np.array([diffLenFin]), fillvalue='NaN')

    with open(pathResults+site+'_EffectsSigma_NoConstraints.csv', 'w', newline='') as f:
        csv.writer(f).writerow(header)
        csv.writer(f).writerows(rows)

    print('\nGenerating results w constraints:\n')
    dep, dat, diff, Ps, Ts, pats = inst.BackDiffused_constraints()

    headerC = ['depth','data','idxPeaks','sigma']
    rowsC = zip_longest(dep, dat, Ps, np.array([diff]), fillvalue='NaN')

    with open(pathResults+site+'_EffectsSigma_Constraints.csv', 'w', newline='') as f:
        csv.writer(f).writerow(headerC)
        csv.writer(f).writerows(rowsC)



"""
    - No constraints/Constraints
        - Sigma estimate
            - Quantitative analysis w. uncertainties based on ALT estimate
"""



"""
    - No constraints/Constraints
        - Speed
"""

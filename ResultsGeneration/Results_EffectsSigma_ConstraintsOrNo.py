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
import time

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
from src.DiffLen_UncertaintyEst import Calc_diffLen_Gauss_MonthVar

'''
    EFFECTS ON SIGMA ESTIMATES
'''
"""
    - No constraints/Constraints
        - Sigma estimate
            - Qualitative analysis for visual inspection
"""
#sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
sites = ['Crete']
pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'
#
# for i in range(len(sites)):
#
#     site = sites[i]
#     print('\n###### ' + site + '######\n')
#
#     N_InInt = 33
#
#     CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
#
#     coreNames = CoresSpecs['CoreName']
#
#
#     core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
#     CoreSpecs = CoresSpecs.iloc[core_idx]
#     dTamb = CoreSpecs['dTamb']
#     dLaki = CoreSpecs['dLaki']
#     accum0 = CoreSpecs['Accum0']
#     accumIE = CoreSpecs['Accum1']
#     Temp0 = CoreSpecs['T0']
#
#     DataAll = GetCoreData(site, 'Alphabet')
#
#     data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
#     data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
#     data_dens = DataAll[4]; data_dens_LT = DataAll[5]
#     data_diff = DataAll[6]; data_diff_LT = DataAll[7]
#
#
#     depth = data_d18O['depth']
#     d18O = data_d18O['d18O']
#
#     depth_LT = data_d18O_LT['depth']
#     d18O_LT = data_d18O_LT['d18O']
#
#     dataAll = pd.DataFrame({'depth':depth,'d18O':d18O}, index=None)
#
#     inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=5, transType = 'DCT')
#
#     print('Generating results w/o constraints\n')
#     depthEst, dataEst, diffLenFin, idxPeak, arr_diffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(print_Npeaks=False)
#
#
#     header = ['depth','data','idxPeaks','sigma']
#     rows = zip_longest(depthEst, dataEst, idxPeak, np.array([diffLenFin]), fillvalue='NaN')
#
#     with open(pathResults+site+'_EffectsSigma_NoConstraints.csv', 'w', newline='') as f:
#         csv.writer(f).writerow(header)
#         csv.writer(f).writerows(rows)
#
#     print('\nGenerating results w constraints:\n')
#     dep, dat, diff, Ps, Ts, pats = inst.BackDiffused_constraints()
#
#     headerC = ['depth','data','idxPeaks','sigma']
#     rowsC = zip_longest(dep, dat, Ps, np.array([diff]), fillvalue='NaN')
#
#     with open(pathResults+site+'_EffectsSigma_Constraints.csv', 'w', newline='') as f:
#         csv.writer(f).writerow(headerC)
#         csv.writer(f).writerows(rowsC)



"""
    - No constraints/Constraints
        - Sigma estimate
            - Quantitative analysis w. uncertainties based on ALT estimate
"""
#sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']




for i in range(len(sites)):
    site = sites[i]
    print('\n##########'+site+'##########\n')
    N_InInt = 33

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    N = 500
    diffLens_w = np.zeros(N)
    dTambs_w = np.zeros(N)
    dLakis_w = np.zeros(N)
    t_w = np.zeros(N)

    diffLens_wo = np.zeros(N)
    dTambs_wo = np.zeros(N)
    dLakis_wo = np.zeros(N)
    t_wo = np.zeros(N)

    for i in range(N):
        print('\n'+str(i))
        t_w[i], diffLens_w[i], dTambs_w[i], dLakis_w[i] = Calc_diffLen_Gauss_MonthVar(site, 33, CoresSpecs, constraints=True, timeIt = True)
        t_wo[i], diffLens_wo[i], dTambs_wo[i], dLakis_wo[i] = Calc_diffLen_Gauss_MonthVar(site, 33, CoresSpecs, constraints=False, timeIt=True)

        print(f't_w: {t_w[i]:.2f}')
        print(f't_wo: {t_wo[i]:.2f}\n')
    np.savetxt(pathResults + site+'_diffLens_ConstraintsGauss.csv', np.array([t_w, diffLens_w,dTambs_w,dLakis_w]))
    np.savetxt(pathResults + site+'_diffLens_NoConstraintsGauss.csv', np.array([t_wo, diffLens_wo,dTambs_wo,dLakis_wo]))



"""
    - No constraints/Constraints
        - Speed
"""

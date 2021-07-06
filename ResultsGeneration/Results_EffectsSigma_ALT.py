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


sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']



'''
    ALT ESTIMATION EFFECTS
'''
"""
    - Shift vs. ALT at LT
"""

pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'

#
# for i in range(len(sites)):
#     site = sites[i]
#
#     CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
#
#     coreNames = CoresSpecs['CoreName']
#
#     print('\n\n\t####'+site+'####\n')
#     core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
#     CoreSpecs = CoresSpecs.iloc[core_idx]
#     dTamb = CoreSpecs['dTamb']
#     dLaki = CoreSpecs['dLaki']
#
#     DataAll = GetCoreData(site, 'Alphabet')
#
#     data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
#
#     depth = np.asarray(data_d18O['depth'])
#     d18O = np.asarray(data_d18O['d18O'])
#
#     lsec = 5
#     shifts = np.linspace(0.1,lsec,70)
#
#     l_shift = np.zeros(len(shifts))
#     lStd_shift = np.zeros(len(shifts))
#     dataAll = []
#
#     ALT_inst = AnnualLayerThick(depth, d18O, lsec)
#
#     for j in range(len(shifts)):
#
#         fks, ls, lMean, lStd, secs = ALT_inst.ALT_fullCore_seq(shift=shifts[j], printItes=True)
#
#         dataAll.append(np.asarray(secs))
#         dataAll.append(np.asarray(lMean))
#         dataAll.append(np.asarray(lStd))
#
#         fks_LT = fks[(secs >= dTamb) & (secs <= dLaki)]
#         fks_LT_pos = fks_LT[fks_LT>0]
#
#         l_LT = np.mean(1/(fks_LT_pos))
#         lStd_LT = np.std(1/(fks_LT_pos))
#
#         l_shift[j] = l_LT
#         lStd_shift[j] = lStd_LT
#
#         print(f'ite {j}/{len(shifts)}. ALT at LT: {l_LT:.3f} +/- {lStd_LT:.3f}')
#
#     N = int(len(dataAll)/3)
#
#     elems = np.arange(0, N)
#     reps = np.repeat(elems, 3)
#     repsStr = np.asarray(list(map(str, reps)))
#
#     headerSolo = np.asarray(N*['secs','lMean','lStd'])
#     header = np.char.add(headerSolo,repsStr)
#
#     df_dataAll = pd.DataFrame(dataAll, header).T
#
#
#     df_dataAll.to_csv(pathResults+site+'_ALTvDepth_AllShifts.csv')
#
#
#     np.savetxt(pathResults+site+'_ALTvShifts.csv',np.array([shifts,l_shift, lStd_shift]).T, delimiter=",", header = 'shifts,lMeans,lStds',comments='')
#






"""
    - Len. of section vs. ALT at LT
"""
for i in range(len(sites)):
    site = sites[i]

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']

    print('\n\n\t####'+site+'####\n')
    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']

    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]

    depth = np.asarray(data_d18O['depth'])
    d18O = np.asarray(data_d18O['d18O'])

    lens = np.linspace(1,25, 100)
    l_lens = np.zeros(len(lens))
    lStd_lens = np.zeros(len(lens))
    dataAll_lens = []

    for j in range(len(lens)):
        ALT_inst = AnnualLayerThick(depth, d18O, lens[j])
        fks, ls, lMean, lStd, secs = ALT_inst.ALT_fullCore_seq(shift=1.0, printItes=False)

        dataAll_lens.append(np.asarray(secs))
        dataAll_lens.append(np.asarray(lMean))
        dataAll_lens.append(np.asarray(lStd))

        fks_LT = fks[(secs >= dTamb) & (secs <= dLaki)]
        fks_LT_pos = fks_LT[fks_LT>0]

        l_LT = np.mean(1/(fks_LT_pos))
        lStd_LT = np.std(1/(fks_LT_pos))
        print(f'Sec. len: {lens[j]:.2f}.\nALT at LT: {l_LT:.3f} +/- {lStd_LT:.3f}\n\n')
        l_lens[j] = l_LT
        lStd_lens[j] = lStd_LT

    N = int(len(dataAll_lens)/3)

    elems = np.arange(0, N)
    reps = np.repeat(elems, 3)
    repsStr = np.asarray(list(map(str, reps)))

    headerSolo = np.asarray(N*['secs','lMean','lStd'])
    header = np.char.add(headerSolo,repsStr)

    df_dataAll = pd.DataFrame(dataAll_lens, header).T


    df_dataAll.to_csv(pathResults+site+'_ALTvDepth_AllLens.csv')


    np.savetxt(pathResults+site+'_ALTvLens.csv',np.array([lens,l_lens, lStd_lens]).T, delimiter=",", header = 'lens,lMeans,lStds',comments='')

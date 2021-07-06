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

import sys
import os

sys.path.append('../../')
sys.path.append('../')

from Interpolation_Class import Interpolation


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{textcomp}',
    r'\usepackage{wasysym}']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'STIXGeneral'

saveFigs = False

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
    SENSITIVITY TESTING
'''
"""
    - Interpolation effect
        - Before deconvolution
        - After deconvolution
"""
sites = ['SiteA' , 'SiteB', 'SiteD', 'SiteE', 'SiteG']

shift_in = 1.5
lSecs_in = 7
pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'
#
# print('\t\t####################################')
# print('\t\t### INTERPOLATION AF/BF EFFECTS ####')
# print('\t\t###### INTERPOLATION BF DECON ######')
# print('\t\t####################################')
#
# deltaMins = [0.02,0.022,0.01,0.02,0.02]
# deltaMaxs = [0.13,0.12,0.14,0.12,0.11]
#
# interpType = 'CubicSpline'
#
#
# for j in range(len(sites)):
#     delta_arr = np.linspace(deltaMins[j],deltaMaxs[j],2)
#         # Load data
#     site = sites[j]
#     N_InInt = 33
#
#     print(f'\n {site}')
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
#     isoData = data_d18O
#     def avg(a):
#         return a[a > 0].mean()
#     def std(a):
#         return a[a>0].std()
#
#     try:
#         pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'
#         data = pd.read_csv(pathResults + site + '_ALT_FullCore_Pshift_'+str(int(shift_in))+'_lSecs_'+str(lSecs_in)+'.csv')
#         print('ALT file exists. Loading ALT data.')
#
#         lDCT = np.asarray(data['lDCT']);lNDCT = np.asarray(data['lNDCT']);lFFT = np.asarray(data['lFFT']);
#         vals_use = data['depth']
#
#         lks = np.c_[lDCT,lNDCT,lFFT]
#         lks_LT = lks[(vals_use>=dTamb)&(vals_use<=dLaki)]
#
#         l_LT = avg(lks_LT)
#         lStd_LT = std(lks_LT)
#
#         lMean = data['lMean']
#         lStd = data['lStd']
#         vals_use = data['depth']
#
#         # Otherwise compute ALTs
#     except:
#         print('ALT file does NOT exist. Computing ALT for core.')
#
#         depth_ALT = np.asarray(isoData['depth'])
#         d18O_ALT = np.asarray(isoData['d18O'])
#
#             # Create annual layer thickness instance
#         inst_ALT = AnnualLayerThick(depth_ALT, d18O_ALT, lSecs_in)
#             # Compute ALT for entire core.
#         fksMax, ls, lMean, lStd, vals_use = inst_ALT.ALT_fullCore_seq(shift=shift_in, printItes=False)
#         lks_LT = ls[(vals_use>=dTamb)&(vals_use<=dLaki)]
#
#         l_LT = avg(lks_LT)
#         lStd_LT = std(lks_LT)
#             # Compute an estimate for ALT at LT depth
#         #l_LT = np.mean(lMean[(vals_use > self.depthMin) & (vals_use < self.depthMax)])
#
#     ALT_LT = l_LT
#
#     interval = np.array([min(depth_LT), max(depth_LT)])
#     interpTypeAll = interpType
#
#     diffLens = np.zeros(len(delta_arr))
#     Npeaks = np.zeros(len(delta_arr))
#     #depths_BD = []
#     #datas_BD = []
# #    depth_ints = []
# #    data_ints = []
#
#     for i in range(len(delta_arr)):
#         print(f'\nRun {i}')
#         print(f'Delta: {delta_arr[i]:.3f}\n')
#         inst = Interpolation(depth_LT, pd.Series(d18O_LT), interval, interpTypeAll, DeltaInput=True, samplingSize=delta_arr[i])
#         depth_LT_int1, d18O_LT_int1, Delta = inst()
#
# #        depth_ints.append(depth_LT_int1)
# #        data_ints.append(d18O_LT_int1)
#
#         dataAll = pd.DataFrame({'depth':depth_LT_int1,'d18O':d18O_LT_int1}, index=None)
#
#         inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)
#         depth1, data, diffLen, Peaks, Ts, pats = inst.BackDiffused_constraints()
# #        depth1, data, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=False, newDelta=0.005)
#
#         Npeaks[i] = len(Peaks)
#         diffLens[i] = diffLen
#         #depths_BD.append(depth1)
#         #datas_BD.append(data)
#
#     df_Site = pd.DataFrame({'diffLens':diffLens, 'deltas':delta_arr})
#
#     df_Site.to_csv(pathResults+site+'_DiffLensVdelta_InterpBF_const.txt',sep='\t', index=False)




print('\t\t####################################')
print('\t\t### INTERPOLATION AF/BF EFFECTS ####')
print('\t\t###### INTERPOLATION AF DECON ######')
print('\t\t####################################')



delta_arr_in = np.arange(0.01,0.075,0.0005)


for j in range(len(sites)):

        # Load data
    site = sites[j]
    N_InInt = 33

    print(f'\n {site}')
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
    isoData = data_d18O
    def avg(a):
        return a[a > 0].mean()
    def std(a):
        return a[a>0].std()

    try:
        pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'
        data = pd.read_csv(pathResults + site + '_ALT_FullCore_Pshift_'+str(int(shift_in))+'_lSecs_'+str(lSecs_in)+'.csv')
        print('ALT file exists. Loading ALT data.')


        lDCT = np.asarray(data['lDCT']);lNDCT = np.asarray(data['lNDCT']);lFFT = np.asarray(data['lFFT']);
        vals_use = data['depth']

        lks = np.c_[lDCT,lNDCT,lFFT]
        lks_LT = lks[(vals_use>=dTamb)&(vals_use<=dLaki)]

        l_LT = avg(lks_LT)
        lStd_LT = std(lks_LT)

        lMean = data['lMean']
        lStd = data['lStd']
        vals_use = data['depth']

        # Otherwise compute ALTs
    except:
        print('ALT file does NOT exist. Computing ALT for core.')

        depth_ALT = np.asarray(isoData['depth'])
        d18O_ALT = np.asarray(isoData['d18O'])

            # Create annual layer thickness instance
        inst_ALT = AnnualLayerThick(depth_ALT, d18O_ALT, lSecs_in)
            # Compute ALT for entire core.
        fksMax, ls, lMean, lStd, vals_use = inst_ALT.ALT_fullCore_seq(shift=shift_in, printItes=False)
        lks_LT = ls[(vals_use>=dTamb)&(vals_use<=dLaki)]

        l_LT = avg(lks_LT)
        lStd_LT = std(lks_LT)
            # Compute an estimate for ALT at LT depth
        #l_LT = np.mean(lMean[(vals_use > self.depthMin) & (vals_use < self.depthMax)])

    ALT_LT = l_LT



    delta_arr = delta_arr_in
    diffLens = []
    #    depths = []
    #    datas = []
    #    peakss = []

    inst = BackDiffuse(site, data_d18O_LT, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)

    for i in range(len(delta_arr)):

        print(f'\n\t\tRun {i} of {len(delta_arr)}')
        depth1, data, diffLen, Peaks, Ts, pats = inst.BackDiffused_constraints(interpAfterDecon=True, newDelta=delta_arr[i])
        depth1, data, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True, newDelta=delta_arr[i])
    #        depths.append(depth1)
    #        datas.append(data)
        diffLens.append(diffLen)

    df_Site = pd.DataFrame({'diffLens':diffLens, 'deltas':delta_arr})

    df_Site.to_csv(pathResults+site+'_DiffLensVdelta_InterpAF_const.txt',sep='\t', index=False)

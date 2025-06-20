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

from src.BackDiffuse_LT import BackDiffuse
from src.GetCoreData_fct import GetCoreData
from src.Interpolation_Class import Interpolation


'''
        Function to test diff len vs. interp. delta after and before signal analysis and
        deconvolution.


        *********** TODO************

        - Make sure that final diff len corresponds to 32 peaks!
        - If diff len < 0, then don't continue!!!
'''


def getDiffLen_V_Npeaks(site_in, diffLen_start, diffLen_end, yrsInSec=33, max_Npeaks = 500, interpType='CubicSpline'):
    site = site_in

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']

    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']


    DataAll = GetCoreData(site)

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
    data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
    data_dens = DataAll[4]; data_dens_LT = DataAll[5]
    data_diff = DataAll[6]; data_diff_LT = DataAll[7]


    depth_LT = data_d18O_LT['depth']
    d18O_LT = data_d18O_LT['d18O']

    interval = np.array([min(depth_LT), max(depth_LT)])
    interpTypeAll = interpType

    inst = Interpolation(depth_LT, pd.Series(d18O_LT), interval, interpTypeAll)
    depth_LT_int1, d18O_LT_int1, Delta = inst()


    dataAll = pd.DataFrame({'depth':depth_LT_int1,'d18O':d18O_LT_int1}, index=None)

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, ysInSec=500, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)

    diffLen = inst.spectralEstimate()
    difflenEstHL = inst.diffLenEstimateHL()
    _,_,_,_, arr_DiffLens, arr_Npeaks, _,_ = inst.backDiffused(theoDiffLen=False,print_Npeaks=False, diffLenStart_In=diffLen_start, diffLenEnd_In=diffLen_end)

    idx = np.where(np.asarray(arr_Npeaks) <= yrsInSec)[0]
    Npeaks_upTo32 = np.asarray(arr_Npeaks)[idx]
    DiffLens_upTo32 = np.asarray(arr_DiffLens)[idx]


    df_Site = pd.DataFrame({'diffLen':arr_DiffLens, 'Npeaks':arr_Npeaks})
    df_Site_upTo32 = pd.DataFrame({'diffLen_upTo32':DiffLens_upTo32,'Npeaks_upTo32':Npeaks_upTo32})

    df_Site.to_csv('../Data/'+site+'_DiffLensVpeaks.txt',sep='\t', index=False)
    df_Site_upTo32.to_csv('../Data/'+site+'_DiffLensVpeaks_upTo32.txt',sep='\t', index=False)

    return


def getInterpBFdata(site_in, delta_arr_in, yrsInSec=33, interpType='CubicSpline'):
    site = site_in

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']

    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']


    DataAll = GetCoreData(site)

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
    data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
    data_dens = DataAll[4]; data_dens_LT = DataAll[5]
    data_diff = DataAll[6]; data_diff_LT = DataAll[7]


    depth_LT = data_d18O_LT['depth']
    d18O_LT = data_d18O_LT['d18O']

    delta_arr = delta_arr_in

    interval = np.array([min(depth_LT), max(depth_LT)])
    interpTypeAll = interpType

    diffLens = np.zeros(len(delta_arr))
    Npeaks = np.zeros(len(delta_arr))
    depths_BD = []
    datas_BD = []
#    depth_ints = []
#    data_ints = []

    for i in range(len(delta_arr)):
        print(f'\nRun {i}')
        print(f'Delta: {delta_arr[i]:.3f}\n')
        inst = Interpolation(depth_LT, pd.Series(d18O_LT), interval, interpTypeAll, DeltaInput=True, samplingSize=delta_arr[i])
        depth_LT_int1, d18O_LT_int1, Delta = inst()

#        depth_ints.append(depth_LT_int1)
#        data_ints.append(d18O_LT_int1)

        dataAll = pd.DataFrame({'depth':depth_LT_int1,'d18O':d18O_LT_int1}, index=None)

        inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, yrsInSec, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=2)
        diffLen = inst.spectralEstimate()
        difflenEstHL = inst.diffLenEstimateHL()
        depth1, data, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=False, newDelta=0.005)

        Npeaks[i] = len(peaks)
        diffLens[i] = diffLen
        depths_BD.append(depth1)
        datas_BD.append(data)

    df_Site = pd.DataFrame({'diffLens':diffLens, 'deltas':delta_arr})

    df_Site.to_csv('../Data/'+site+'_DiffLensVdelta_InterpBF.txt',sep='\t', index=False)

    return




def getInterpAFdata(site_in, delta_arr_in, yrsInSec=33):
    site = site_in

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']

    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']


    DataAll = GetCoreData(site)

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
    data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
    data_dens = DataAll[4]; data_dens_LT = DataAll[5]
    data_diff = DataAll[6]; data_diff_LT = DataAll[7]


    depth_LT = data_d18O_LT['depth']
    d18O_LT = data_d18O_LT['d18O']

    delta_arr = delta_arr_in
    diffLens = []
#    depths = []
#    datas = []
#    peakss = []


    inst = BackDiffuse(site, data_d18O_LT, CoresSpecs, dTamb, dLaki, yrsInSec, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)

    for i in range(len(delta_arr)):

        print(f'\n\t\tRun {i} of {len(delta_arr)}')
        depth1, data, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True, newDelta=delta_arr[i])
#        depths.append(depth1)
#        datas.append(data)
        diffLens.append(diffLen)

    df_Site = pd.DataFrame({'diffLens':diffLens, 'deltas':delta_arr})

    df_Site.to_csv('../Data/'+site+'_DiffLensVdelta_InterpAF.txt',sep='\t', index=False)
    return

'''
    Generate data for diff. len. V. N peaks.
'''
# sites = ['Crete', 'SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
# for i in range(len(sites)):
#     print('\nSite: ' + sites[i] + '\n')
#     getDiffLen_V_Npeaks(site_in=sites[i], diffLen_start=0.005, diffLen_end=0.15)


'''
    Generation of data with interpolation AFTER deconvolution (SAVE AND USE)
'''
# CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
#
# coreNames = CoresSpecs['CoreName']
#
# sites = ['SiteA', 'SiteB', 'SiteE', 'SiteG']
# delta_arr = np.arange(0.01,0.075,0.0005)
#
# for site in sites:
#     print('\n\n\n###############')
#     print('#### ' + site + ' ####')
#     print('###############')
#     getInterpAFdata(site_in=site, delta_arr_in=delta_arr)





'''
    Interpolation BEFORE deconvolution (DO NOT USE)
'''
# CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
#
# coreNames = CoresSpecs['CoreName']
#
# sites = ['Crete', 'SiteA','SiteB', 'SiteE', 'SiteG']
# deltaMins = [0.02,0.022,0.01,0.02,0.02]
# deltaMaxs = [0.13,0.12,0.14,0.12,0.11]
# yrs = [32,32,33,32,32]
#
# i = 0
# for i in range(0,len(sites)):
#     print('\n\n\n###############')
#     print('#### ' + sites[i] + ' ####')
#     print('###############')
#     delta_arr = np.linspace(deltaMins[i],deltaMaxs[i],100)
#     getInterpBFdata(site_in=sites[i], delta_arr_in=delta_arr,yrsInSec=yrs[i])
#     i += 1

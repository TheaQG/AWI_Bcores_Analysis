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
        Function to test diff len vs. interp. delta AFTER signal analysis and
        deconvolution.


        *********** TODO************

        - Make sure that final diff len corresponds to 32 peaks!
        - If diff len < 0, then don't continue!!!
'''






def getInterpAFdata(site_in, delta_arr_in):
    site = site_in

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
    depths = []
    datas = []
    peakss = []


    inst = BackDiffuse(site, data_d18O_LT, CoresSpecs, dTamb, dLaki, 32, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)

    for i in range(len(delta_arr)):

        print(f'\n\t\tRun {i} of {len(delta_arr)}')
        depth1, data, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True, newDelta=delta_arr[i])
        depths.append(depth1)
        datas.append(data)
        diffLens.append(diffLen)

    df_Site = pd.DataFrame({'diffLens':diffLens, 'deltas':delta_arr})

    df_Site.to_csv('../Data/'+site+'_DiffLensVdelta_InterpAF.txt',sep='\t', index=False)
    return




CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

coreNames = CoresSpecs['CoreName']

sites = ['Crete', 'SiteA', 'SiteB', 'SiteE', 'SiteG']
delta_arr = np.arange(0.004,0.075,0.0005)

for site in sites:
    print('\n\n\n###############')
    print('#### ' + site + ' ####')
    print('###############')
    getInterpAFdata(site_in=site, delta_arr_in=delta_arr)

import numpy as np
import os
import pandas as pd
import openpyxl
import scipy as sp

from GetCoreData_fct import GetCoreData

import sys
sys.path.append('../')

from BackDiffuse_LT import BackDiffuse

from Interpolation_Class import Interpolation

from HL_AnalyticThea_class import HL_Thea
from DiffusionProfiles_calculations import DiffusionLength

from sigmaSolver import sigma_Solver


def TempEst_analytical(site, N_InInt):
        # Read and define data of specific interest
    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']


    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']
    accum0 = CoreSpecs['Accum0']
    Temp0 = CoreSpecs['T0']

    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
    #data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
    data_dens = DataAll[4]; data_dens_LT = DataAll[5]
    data_diff = DataAll[6]; data_diff_LT = DataAll[7]


    depth_LT = data_d18O_LT['depth']
    d18O_LT = data_d18O_LT['d18O']


        # Compute diffusion length estimate interval to result in N peaks
    dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    inst = BackDiffuse(site, data_d18O_LT, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)

    depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    idxN = np.where(np.asarray(arr_Npeaks) == N_InInt)[0]
    idxNCut = idxN[:-1]

    diffLensN = np.asarray(arr_DiffLens)[idxNCut]



        # Determine temperature interval estimate (analytical solution to diffusion equation)
    accum = accum0

    sigmaSolver_inst = sigma_Solver()

    T_intEst = np.zeros(len(diffLensN))

    for i in range(len(diffLensN)):
        T_est = sigmaSolver_inst.solveTemp(sigma_data = diffLensN[i], accum = accum)# *(804.3/917.)
        T_intEst[i] = T_est

    return T_intEst, diffLensN

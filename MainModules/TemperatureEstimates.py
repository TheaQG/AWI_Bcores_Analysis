import numpy as np
import os
import pandas as pd
import openpyxl
import scipy as sp

from src.GetCoreData_fct import GetCoreData

import sys
sys.path.append('../')

from src.BackDiffuse_LT import BackDiffuse

from src.Interpolation_Class import Interpolation

from src.HL_AnalyticThea_class import HL_Thea
from src.DiffusionProfiles_calculations import DiffusionLength

from src.sigmaSolver import sigma_Solver


def TempEst_analytical(site, N_InInt, Accum_in = 0, T_in = 100):
        # Read and define data of specific interest
    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']


    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']

    if T_in != 100:
        Temp0 = T_in
    else:
        Temp0 = CoreSpecs['T0']

    if Accum_in != 0:
        accum0 = Accum_in
    else:
        accum0 = CoreSpecs['Accum0']


    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]
    #data_ECM = DataAll[2]; data_ECM_LT = DataAll[3]
    data_dens = DataAll[4]; data_dens_LT = DataAll[5]
    data_diff = DataAll[6]; data_diff_LT = DataAll[7]

    depth_LT = data_d18O_LT['depth']
    d18O_LT = data_d18O_LT['d18O']

    try:
        dens_LT = data_dens_LT['HLmodelOpti']
        densDepth_LT = data_dens_LT['depth']
    except:
        dens_LT = data_dens_LT['HLmodel']
        densDepth_LT = data_dens_LT['depth']

    dens_LT_ave = np.mean(dens_LT)*1000

    if (dens_LT_ave < 804.3):
        rhoMean = dens_LT_ave
    elif (dens_LT_ave >= 804.3):
        rhoMean = dens_LT_ave#804.3


        # Compute diffusion length estimate interval to result in N peaks
    dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    inst = BackDiffuse(site, data_d18O_LT, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT)

    depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    idxN = np.where(np.asarray(arr_Npeaks) == N_InInt)[0]
    idxNCut = idxN[:-1]

    diffLensN = np.asarray(arr_DiffLens)[idxNCut]

        # Compute sigma due to discrete sampling (taking an average of the sampling sizes in the interval)
    from src.DiffusionProfiles_calculations import sampling_sigma
    dz_ave = np.mean(np.diff(depth_LT))
    dz_std = np.std(np.diff(depth_LT))
    diffLen_sample = sampling_sigma(dz_ave)

    print(f'Average sampling size in interval: {dz_ave:.4f} +/- {dz_std:.4f}')
    diffLensN_firn = np.sqrt(diffLensN**2 - diffLen_sample**2)

        # Determine temperature interval estimate (analytical solution to diffusion equation)
    accum = accum0

    sigmaSolver_inst = sigma_Solver()

        # Compute for both firn and total diffusion length estimate
    T_intEst = np.zeros(len(diffLensN))
    T_firn_intEst = np.zeros(len(diffLensN))

    for i in range(len(diffLensN)):
        T_est = sigmaSolver_inst.solveTemp(sigma_data = diffLensN[i], accum = accum, rho_CO=rhoMean)
        T_firn_est = sigmaSolver_inst.solveTemp(sigma_data = diffLensN_firn[i], accum = accum, rho_CO=rhoMean)

        T_intEst[i] = T_est
        T_firn_intEst[i] = T_firn_est

    return T_intEst, diffLensN, T_firn_intEst, diffLensN_firn

def TempEst_analytical_arr(diffLens_in = np.array([0.08]), Accum_in = 0.3, rhoMeans_in = np.array([804.3])):
    diffLens = diffLens_in
    sigmaSolver_inst = sigma_Solver()
    T_intEst = np.zeros(len(diffLens))

    for i in range(len(diffLens)):
        T_est = sigmaSolver_inst.solveTemp(sigma_data = diffLens[i], accum = Accum_in, rho_CO = rhoMeans_in[i])
        T_intEst[i] = T_est
    return T_intEst

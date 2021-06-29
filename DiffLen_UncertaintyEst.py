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
from scipy import signal
from GetCoreData_fct import GetCoreData

from Decon import SpectralDecon
from BackDiffuse_LT import BackDiffuse

from Interpolation_Class import Interpolation

from HL_AnalyticThea_class import HL_Thea
from DiffusionProfiles_calculations import DiffusionLength

from sigmaSolver import sigma_Solver




def Calc_diffLen_RandLT(site_in, N_InInt, CoresSpecs):
    site = site_in
        # Get Laki and Tambora positions along with other core specs.
    coreNames = CoresSpecs['CoreName']
    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb_in = CoreSpecs['dTamb']
    dLaki_in = CoreSpecs['dLaki']
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

    depth_ECM = np.asarray(data_ECM['depth']); ECM = np.asarray(data_ECM['ECM'])
    depth_ECM_LT = np.asarray(data_ECM_LT['depth']); ECM_LT = np.asarray(data_ECM_LT['ECM'])


        # Define Laki and Tambora location and eruption width

    lenTamb = CoreSpecs['lenTamb']/100
    lenLaki = CoreSpecs['lenLaki']/100

    startTamb = dTamb_in - lenTamb/2; endTamb = dTamb_in + lenTamb/2
    startLaki = dLaki_in - lenLaki/2; endLaki = dLaki_in + lenLaki/2


        # Define starting ECM and d18O data
    depth_ECM_Laki = depth_ECM[(depth_ECM > startLaki) & (depth_ECM < endLaki)]
    ECM_Laki = ECM[(depth_ECM > startLaki) & (depth_ECM < endLaki)]

    depth_ECM_Tamb = depth_ECM[(depth_ECM > startTamb) & (depth_ECM < endTamb)]
    ECM_Tamb = ECM[(depth_ECM > startTamb) & (depth_ECM < endTamb)]


    depth_Laki = depth[(depth > startLaki) & (depth < endLaki)]
    d18O_Laki = d18O[(depth > startLaki) & (depth < endLaki)]

    depth_Tamb = depth[(depth > startTamb) & (depth < endTamb)]
    d18O_Tamb = d18O[(depth > startTamb) & (depth < endTamb)]


        ###############################################################
        ### Make a block here to compute distribution to draw from! ###
        ###############################################################

        # Define new 'random' (uniform) estimate of Laki and Tambora positions, within peak width

    maxTamb = lenTamb/2
    maxLaki = lenLaki/2

    randTamb = np.random.uniform(-maxTamb, maxTamb)
    randLaki = np.random.uniform(-maxLaki, maxLaki)

    dTamb = dTamb_in + randTamb
    dLaki = dLaki_in + randLaki


        # Define new sets of data in new LT_rand interval

    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]
    data_ECM = DataAll[2]
    data_dens = DataAll[4]
    data_diff = DataAll[6]


    depth = data_d18O['depth']
    d18O = data_d18O['d18O']

    depth_LT = depth[(depth >= dTamb) & (depth <= dLaki)]
    d18O_LT = d18O[(depth >= dTamb) & (depth <= dLaki)]



    depthDiff = data_diff['Depth']
    diff = data_diff['sigma_o18']

    depthDiff_LT = depthDiff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    diff_LT = diff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]



    depthDens = data_dens['depth']
    dens = data_dens['HLmodel']

    depthDens_LT = depthDens[(depthDens >= dTamb) & (depthDens <= dLaki)]
    dens_LT = dens[(depthDens >= dTamb) & (depthDens <= dLaki)]







    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]
    data_ECM = DataAll[2]
    data_dens = DataAll[4]
    data_diff = DataAll[6]


    depth = data_d18O['depth']
    d18O = data_d18O['d18O']

    depth_LT = depth[(depth >= dTamb) & (depth <= dLaki)]
    d18O_LT = d18O[(depth >= dTamb) & (depth <= dLaki)]



    depthDiff = data_diff['Depth']
    diff = data_diff['sigma_o18']

    depthDiff_LT = depthDiff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    diff_LT = diff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]



    depthDens = data_dens['depth']
    dens = data_dens['HLmodel']

    depthDens_LT = depthDens[(depthDens >= dTamb) & (depthDens <= dLaki)]
    dens_LT = dens[(depthDens >= dTamb) & (depthDens <= dLaki)]


    data_dens_LT = pd.DataFrame({'depth': depthDens_LT, 'HLmodel': dens_LT})
    data_diff_LT = pd.DataFrame({'Depth': depthDiff_LT, 'sigma_o18': diff_LT})


        # Compute diffusion length estimate:
    dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    return diffLen










def Calc_diffLen_Gauss(site_in, N_InInt, CoresSpecs):
    site = site_in
        # Get Laki and Tambora positions along with other core specs.
    coreNames = CoresSpecs['CoreName']
    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb_in = CoreSpecs['dTambCor']
    dLaki_in = CoreSpecs['dLakiCor']
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

    depth_ECM = np.asarray(data_ECM['depth']); ECM = np.asarray(data_ECM['ECM'])
    depth_ECM_LT = np.asarray(data_ECM_LT['depth']); ECM_LT = np.asarray(data_ECM_LT['ECM'])


        # Define Laki and Tambora location and eruption width

    lenTamb = CoreSpecs['lenTambCor']/100
    lenLaki = CoreSpecs['lenLakiCor']/100

    startTamb = dTamb_in - lenTamb/2; endTamb = dTamb_in + lenTamb/2
    startLaki = dLaki_in - lenLaki/2; endLaki = dLaki_in + lenLaki/2


        # Define starting ECM and d18O data
    depth_ECM_Laki = depth_ECM[(depth_ECM > startLaki) & (depth_ECM < endLaki)]
    ECM_Laki = ECM[(depth_ECM > startLaki) & (depth_ECM < endLaki)]

    depth_ECM_Tamb = depth_ECM[(depth_ECM > startTamb) & (depth_ECM < endTamb)]
    ECM_Tamb = ECM[(depth_ECM > startTamb) & (depth_ECM < endTamb)]


    depth_Laki = depth[(depth > startLaki) & (depth < endLaki)]
    d18O_Laki = d18O[(depth > startLaki) & (depth < endLaki)]

    depth_Tamb = depth[(depth > startTamb) & (depth < endTamb)]
    d18O_Tamb = d18O[(depth > startTamb) & (depth < endTamb)]


        # Define new 'random' (uniform) estimate of Laki and Tambora positions, within peak width

    maxTamb = lenTamb/2
    maxLaki = lenLaki/2

    randTamb = np.random.normal(dTamb_in, lenTamb/4)
    randLaki = np.random.normal(dLaki_in, lenLaki/5)

    dTamb = randTamb
    dLaki = randLaki


        # Define new sets of data in new LT_rand interval

    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]
    data_ECM = DataAll[2]
    data_dens = DataAll[4]
    data_diff = DataAll[6]


    depth = data_d18O['depth']
    d18O = data_d18O['d18O']

    depth_LT = depth[(depth >= dTamb) & (depth <= dLaki)]
    d18O_LT = d18O[(depth >= dTamb) & (depth <= dLaki)]



    depthDiff = data_diff['Depth']
    diff = data_diff['sigma_o18']

    depthDiff_LT = depthDiff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    diff_LT = diff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]



    depthDens = data_dens['depth']
    dens = data_dens['HLmodel']

    depthDens_LT = depthDens[(depthDens >= dTamb) & (depthDens <= dLaki)]
    dens_LT = dens[(depthDens >= dTamb) & (depthDens <= dLaki)]







    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]
    data_ECM = DataAll[2]
    data_dens = DataAll[4]
    data_diff = DataAll[6]


    depth = data_d18O['depth']
    d18O = data_d18O['d18O']

    depth_LT = depth[(depth >= dTamb) & (depth <= dLaki)]
    d18O_LT = d18O[(depth >= dTamb) & (depth <= dLaki)]



    depthDiff = data_diff['Depth']
    diff = data_diff['sigma_o18']

    depthDiff_LT = depthDiff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    diff_LT = diff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]



    depthDens = data_dens['depth']
    dens = data_dens['HLmodel']

    depthDens_LT = depthDens[(depthDens >= dTamb) & (depthDens <= dLaki)]
    dens_LT = dens[(depthDens >= dTamb) & (depthDens <= dLaki)]


    data_dens_LT = pd.DataFrame({'depth': depthDens_LT, 'HLmodel': dens_LT})
    data_diff_LT = pd.DataFrame({'Depth': depthDiff_LT, 'sigma_o18': diff_LT})


        # Compute diffusion length estimate:
    dataAll = pd.DataFrame({'depth':depth,'d18O':d18O}, index=None)

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    depthOpt, dataOpt, diffLen, Peaks, Ts, pats = inst.BackDiffused_constraints(LayerThickness=0, N_summers=0, N_winters=0, Amplitude=0, N=2000, print_Npeaks=True, theoDiffLen=True, diffLenStart_In=0, diffLenEnd_In=0.1, interpAfterDecon=True, newDelta=0, interpBFDecon=True)

    #dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    #inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    #depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    return diffLen, dTamb, dLaki



sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']

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

    np.savetxt(site+'diffLens_GaussDistwDepths.csv', np.array([diffLens,dTambs,dLakis]))

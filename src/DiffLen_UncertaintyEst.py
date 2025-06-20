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
import time
from src.GetCoreData_fct import GetCoreData
from src.SignalAttenuation import Attenuation, AnnualLayerThick

from src.Decon import SpectralDecon
from src.BackDiffuse_LT import BackDiffuse

from src.Interpolation_Class import Interpolation

from src.HL_AnalyticThea_class import HL_Thea
from src.DiffusionProfiles_calculations import DiffusionLength

from src.sigmaSolver import sigma_Solver




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










def Calc_diffLen_Gauss(site_in, N_InInt, CoresSpecs, section = 'LT', mu1 = 0, mu2 = 0, sig1 = 1, sig2 = 1):
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

    if section == 'LT':
        randTamb = np.random.normal(dTamb_in, lenTamb/4)
        randLaki = np.random.normal(dLaki_in, lenLaki/5)
    else:
        randTamb = np.random.normal(mu1, sig1)
        randLaki = np.random.normal(mu2,sig2)

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


    data_dens_LT = pd.DataFrame({'depth': depthDens_LT, 'HLmodel': dens_LT})
    data_diff_LT = pd.DataFrame({'Depth': depthDiff_LT, 'sigma_o18': diff_LT})


        # Compute diffusion length estimate:
    dataAll = pd.DataFrame({'depth':depth,'d18O':d18O}, index=None)

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    depthOpt, dataOpt, diffLen, Peaks, Ts, pats = inst.BackDiffused_constraints()

    #dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    #inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    #depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    return diffLen, dTamb, dLaki


def Calc_diffLen_Gauss_MonthVar(site_in, N_InInt, CoresSpecs, lsecs = 7, shift_in = 1.5, Nmonths = 2, transType_in='DCT', constraints = True, timeIt=False):
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

    lenLT = dLaki_in - dTamb_in
    isoData = data_d18O
    def avg(a):
        return a[a > 0].mean()
    def std(a):
        return a[a>0].std()

    try:
        pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'
        data = pd.read_csv(pathResults + site + '_ALT_FullCore_Pshift_'+str(int(1.5))+'_lSecs_'+str(7)+'.csv')
        print('ALT file exists. Loading ALT data.')


        lDCT = np.asarray(data['lDCT']);lNDCT = np.asarray(data['lNDCT']);lFFT = np.asarray(data['lFFT']);
        vals_use = data['depth']

        lks = np.c_[lDCT,lNDCT,lFFT]
        lks_LT = lks[(vals_use>=dTamb_in)&(vals_use<=dLaki_in)]

        l_LT = avg(lks_LT)
        lStd_LT = std(lks_LT)

        # Otherwise compute ALTs
    except:
        print('ALT file does NOT exist. Computing ALT for core.')

        depth_ALT = np.asarray(isoData['depth'])
        d18O_ALT = np.asarray(isoData['d18O'])

            # Create annual layer thickness instance
        inst_ALT = AnnualLayerThick(depth_ALT, d18O_ALT, lsecs)
            # Compute ALT for entire core.
        fksMax, ls, lMean, lStd, vals_use = inst_ALT.ALT_fullCore_seq(shift=shift_in, printItes=False)
        lks_LT = ls[(vals_use>=dTamb_in)&(vals_use<=dLaki_in)]

        l_LT = avg(lks_LT)
        lStd_LT = std(lks_LT)
            # Compute an estimate for ALT at LT depth
        #l_LT = np.mean(lMean[(vals_use > self.depthMin) & (vals_use < self.depthMax)])

    MLT_LT = (l_LT/(12/Nmonths)) / 2






    randTamb = np.random.normal(dTamb_in, MLT_LT)
    randLaki = np.random.normal(dLaki_in, MLT_LT)

    dTamb = randTamb
    dLaki = randLaki


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

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30, transType=transType_in)

    if constraints:
        if timeIt:
            t0 = time.time()

        depthOpt, dataOpt, diffLen, Peaks, Ts, pats = inst.BackDiffused_constraints()

        if timeIt:
            t1 = time.time()

            t_tot = t1-t0
    else:
        if timeIt:
            t0 = time.time()

        depthOpt, dataOpt, diffLen, Peaks, arr_diffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(print_Npeaks=False)

        if timeIt:
            t1 = time.time()

            t_tot = t1-t0

    #dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    #inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    #depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    if timeIt:
        return t_tot, diffLen, dTamb, dLaki
    else:
        return diffLen, dTamb, dLaki





def Calc_diffLen_Gauss_1const(site_in, N_InInt, CoresSpecs, eruption='Laki'):
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

    lenLT = dLaki_in - dTamb_in


    if eruption == 'Laki':
        lenLaki = CoreSpecs['lenLakiCor']/100
        startLaki = dLaki_in - lenLaki/2; endLaki = dLaki_in + lenLaki/2
        maxLaki = lenLaki/2
        randLaki = np.random.normal(dLaki_in, lenLaki/5)
        dLaki = randLaki
        dTamb = dLaki - lenLT
    elif eruption == 'Tambora':
        lenTamb = CoreSpecs['lenTambCor']/100
        startTamb = dTamb_in - lenTamb/2; endTamb = dTamb_in + lenTamb/2
        maxTamb = lenTamb/2
        randTamb = np.random.normal(dTamb_in, lenTamb/4)
        dTamb = randTamb
        dLaki = dTamb + lenLT

        # Define Laki and Tambora location and eruption width

        # Define new 'random' (uniform) estimate of Laki and Tambora positions, within peak width

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


    data_dens_LT = pd.DataFrame({'depth': depthDens_LT, 'HLmodel': dens_LT})
    data_diff_LT = pd.DataFrame({'Depth': depthDiff_LT, 'sigma_o18': diff_LT})


        # Compute diffusion length estimate:
    dataAll = pd.DataFrame({'depth':depth,'d18O':d18O}, index=None)

    inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    depthOpt, dataOpt, diffLen, Peaks, Ts, pats = inst.BackDiffused_constraints()

    #dataAll = pd.DataFrame({'depth':depth_LT,'d18O':d18O_LT}, index=None)

    #inst = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30)

    #depthOpt, dataOpt, diffLen, peaks, arr_DiffLens, arr_Npeaks, arr_depth, arr_data = inst.backDiffused(theoDiffLen=True,print_Npeaks=False, diffLenStart_In=0.005, diffLenEnd_In=0.15, interpAfterDecon=True)

    return diffLen, dTamb, dLaki













# sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
#
# for i in range(len(sites)):
#     site = sites[i]
#     print('\n##########'+site+'##########\n')
#     N_InInt = 33
#
#     CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
#
#     N = 500
#     diffLens = np.zeros(N)
#     dTambs = np.zeros(N)
#     dLakis = np.zeros(N)
#
#     for i in range(N):
#         print(i)
#         diffLens[i], dTambs[i], dLakis[i] = Calc_diffLen_Gauss(site, 33, CoresSpecs)
#
#     np.savetxt(site+'diffLens_GaussDistwDepths.csv', np.array([diffLens,dTambs,dLakis]))

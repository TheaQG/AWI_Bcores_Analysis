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
    - Spectral transform effects (FFT/DCT/NDCT)
        - Sigma estimate
        - Speed
"""


def Calc_diffLen_spectralTransform(site_in, N_InInt, CoresSpecs, Nmonths = 1):
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
        inst_ALT = AnnualLayerThick(depth_ALT, d18O_ALT, lSecs)
            # Compute ALT for entire core.
        fksMax, ls, lMean, lStd, vals_use = inst_ALT.ALT_fullCore_seq(shift=shift_in, printItes=False)
        lks_LT = ls[(vals_use>=self.depthMin)&(vals_use<=self.depthMax)]

        l_LT = avg(lks_LT)
        lStd_LT = std(lks_LT)
            # Compute an estimate for ALT at LT depth
        #l_LT = np.mean(lMean[(vals_use > self.depthMin) & (vals_use < self.depthMax)])

    MLT_LT = l_LT/(12/Nmonths)






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

    trans = ['DCT', 'NDCT', 'FFT']

    t0 = time.time()

    instDCT = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30, transType=trans[0])
    depthDCT, dataDCT, diffLenDCT, PeaksDCT, TsDCT, patsDCT = instDCT.BackDiffused_constraints()
    t1 = time.time()

    totalDCT = t1-t0


    t0 = time.time()
    instNDCT = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30, transType=trans[1])
    depthNDCT, dataNDCT, diffLenNDCT, PeaksNDCT, TsNDCT, patsNDCT = instNDCT.BackDiffused_constraints()
    t1 = time.time()

    totalNDCT = t1-t0


    t0 = time.time()
    instFFT = BackDiffuse(site, dataAll, CoresSpecs, dTamb, dLaki, N_InInt, diffLenData=data_diff_LT[['Depth','sigma_o18']], densData=data_dens_LT, Dist=30, transType=trans[2])
    depthFFT, dataFFT, diffLenFFT, PeaksFFT, TsFFT, patsFFT = instFFT.BackDiffused_constraints()
    t1 = time.time()

    totalFFT = t1-t0

    return dTamb, dLaki, diffLenDCT, diffLenNDCT, diffLenFFT, totalDCT, totalNDCT, totalFFT




pathResults = '/home/thea/MesterTesen/Analysis/ResultsGeneration/ResultsData/'

sites = ['SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
print('\t\t####################################')
print('\t\t### SPECTRAL TRANSFORMS EFFECTS ####')
print('\t\t###### 1 MONTH VARIATION ONLY ######')
print('\t\t####################################')


for j in range(len(sites)):
    site = sites[j]
    print('\n##########'+site+'##########\n')
    N_InInt = 33

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    N = 200
    nMonths_in = 1
    dTambs = np.zeros(N)
    dLakis = np.zeros(N)
    diffLenDCTs = np.zeros(N)
    diffLenNDCTs = np.zeros(N)
    diffLenFFTs = np.zeros(N)
    totalDCTs = np.zeros(N)
    totalNDCTs = np.zeros(N)
    totalFFTs = np.zeros(N)


    for i in range(N):
        print(i)
        dTambs[i], dLakis[i], diffLenDCTs[i], diffLenNDCTs[i], diffLenFFTs[i], totalDCTs[i], totalNDCTs[i], totalFFTs[i] = Calc_diffLen_spectralTransform(site, N_InInt, CoresSpecs, Nmonths=nMonths_in)

    np.savetxt(pathResults + site+'_diffLens_SpecTransEffect_wTiming_varyLandT.csv', np.array([dTambs,dLakis, diffLenDCTs, diffLenNDCTs, diffLenFFTs, totalDCTs, totalNDCTs, totalFFTs]))

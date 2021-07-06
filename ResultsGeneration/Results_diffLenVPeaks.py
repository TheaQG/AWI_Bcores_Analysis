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


"""
    - N peaks v. sigma
        - Pattern/No pattern
"""
sites = ['SiteA']#, 'SiteB', 'SiteD', 'SiteE', 'SiteG']
diffLens = np.linspace(0.01,0.15,200)


for i in range(len(sites)):

        # Load data
    site = sites[i]
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
        inst_ALT = AnnualLayerThick(depth_ALT, d18O_ALT, lSecs)
            # Compute ALT for entire core.
        fksMax, ls, lMean, lStd, vals_use = inst_ALT.ALT_fullCore_seq(shift=shift_in, printItes=False)
        lks_LT = ls[(vals_use>=self.depthMin)&(vals_use<=self.depthMax)]

        l_LT = avg(lks_LT)
        lStd_LT = std(lks_LT)
            # Compute an estimate for ALT at LT depth
        #l_LT = np.mean(lMean[(vals_use > self.depthMin) & (vals_use < self.depthMax)])

    ALT_LT = l_LT


    inst = BackDiffuse(site, data_d18O, CoresSpecs, dTamb, dLaki, 33, diffLenData=data_diff_LT, densData=data_dens_LT)

    patterns = np.zeros(len(diffLens))
    N_Ts = np.zeros(len(diffLens))
    N_Ps = np.zeros(len(diffLens))


    for i in range(len(diffLens)):
        #print('\n\n###############')
        if i%10 == 0:
            print(str(i) + f'/{len(diffLens)}')
            print(f'Sigma input: {diffLens[i]*100:.2f} [cm]')

        newDepth, newData, pattern, Ps, Ts, diffLen = inst.BackDiffuse_manuel_constrained(sigma=diffLens[i], ALT_LT_in = ALT_LT)
        patterns[i] = pattern
        N_Ts[i] = len(Ts)
        N_Ps[i] = len(Ps)


    np.savetxt('ResultsData/'+site+'diffLensVNpeaks_constrained.csv', np.array([diffLens,N_Ps,N_Ts,patterns]))

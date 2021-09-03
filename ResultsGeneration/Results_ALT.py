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





'''
    ESTIMATING ANNUAL LAYER THICKNESS
'''
"""
    - Amplitude attenuation
"""
"""
    - Annual Layer Thicknes
"""
sites = ['Crete','SiteA', 'SiteB', 'SiteD', 'SiteE', 'SiteG']
shiftIn = 4
lSecsIn = 7
for i in range(len(sites)):
    site = sites[i]

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')

    coreNames = CoresSpecs['CoreName']

    print(site)
    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]
    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']

    DataAll = GetCoreData(site, 'Alphabet')

    data_d18O = DataAll[0]; data_d18O_LT = DataAll[1]

    depth = np.asarray(data_d18O['depth'])
    d18O = np.asarray(data_d18O['d18O'])

    ALT_inst = AnnualLayerThick(depth, d18O, lSecsIn)
    fks, ls, lMean, lStd, secs = ALT_inst.ALT_fullCore_seq(shift=shiftIn, printItes=False)

    l_LT = np.mean(lMean[(secs > dTamb) & (secs < dLaki)])

    allData = np.c_[fks,ls,lMean,lStd,secs]

    np.savetxt('ResultsData/'+site+'_ALT_FullCore_Pshift_'+str(int(shiftIn))+'_lSecs_'+str(lSecsIn)+'.csv', allData, delimiter=",", header="fDCT,fNDCT,fFFT,lDCT,lNDCT,lFFT,lMean,lStd,depth", comments='')

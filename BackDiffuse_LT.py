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

from Decon import SpectralDecon


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{textcomp}',
    r'\usepackage{wasysym}']
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'STIXGeneral'

class BackDiffuse():

    def __init__(self, coreName, d18OData, coreSpecs, depthMin, depthMax, interpAll = False, diffLenData = None, densData = None):
        self.coreName = coreName
        self.d18OData = d18OData
        self.coreSpecs = coreSpecs
        self.depthMin = depthMin
        self.depthMax = depthMax
        self.interpAll = interpAll

        self.densData = densData
        self.diffLenData = diffLenData

        return




    def __call__():
        return


    def interpCores(self, pad = 1):
        isoData = self.d18OData
        d_in = isoData['depth']
        x_in = isoData['d18O']


        if self.interpAll:
            valMin = d_in.min()
            valmax = d_in.max()
        else:
            valMin = self.depthMin + pad
            valMax = self.depthMax + pad


        d = d_in[(d_in >= valMin) & (d_in <= valMax)]
        x = x_in[(d_in >= valMin) & (d_in <= valMax)]

        diff = np.diff(d)
        Delta = round(min(diff), 3)

        d_min = Delta * np.ceil(d.values[0]/Delta)
        d_max = Delta * np.floor(d.values[-1]/Delta)

        n = int(1 + (d_max - d_min)/Delta)

        j_arr = np.linspace(0,n,n)
        dhat = d_min + (j_arr - 1)*Delta

        f = interpolate.CubicSpline(d,x)

        xhat = f(dhat)

        return dhat, xhat, Delta




    def diffProfile():
        if not diffLenData==None:
            diffDepth = diffLenData['Depth']
            diffData = diffLenData['sigma_o18']
            return diffDepth, diffData
        else:
            print('Compute diff len profile first!')

            return [], []

    def densProfile():
        if not densData==None:
            densDepth = densData['depth']
            densHL = densData['HLmodel']

            return densDepth, densHL

        else:
            print('Compute density profile first"')

            return [], []

    def spectralEstimate():
        dInt, d18OInt, Delta = self.interpCores()

        decon_inst = SpectralDecon(dInt, d18OInt, 2000)
        a,b = decon_inst.
        return

    def diffLenEstimate():
        return

    def backDiffused():
        return

    def interpDiffData():
        return

    def countPeaks():
        return

coreName = 'SiteG'

d18OData = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+coreName+'_det.txt',' ')
densities = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/'+coreName+'DepthDens_w_Models.txt','\t')
diffLens = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/'+coreName+'_DepthDiff.txt','\t')
specsCores = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt',',')


inst = BackDiffuseLT(coreName, d18OData, specsCores,60.50,69.40, diffLenData=diffLens[['Depth','sigma_o18']], densData=densities)
test, test2, test3 = inst.interpCores()
print(test2)

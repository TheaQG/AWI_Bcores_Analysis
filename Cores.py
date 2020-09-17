import numpy as np
import pandas as pd
from pandas import ExcelWriter
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
'''
    Class to handle North Greenland Transverse B cores from AWI.
    Needs data for at least density, d18O isotope and volcanic eruption depth locations
    and takes ECM, DEP, chemical data if available.

    Returns dataframe of data given between the two volcanic eruptions Laki and Tambora.

'''
class Cores():

    def __init__(self, name, df_dens, df_d18O, df_ECM, df_DEP, volcWE):
        self.name = name
        self.df_dens = df_dens
        self.df_d18O = df_d18O
        self.df_ECM = df_ECM
        self.df_DEP = df_DEP
        self.volcWE = volcWE
        return

    def FindVolcErup(self):
        volcWE_use = []

        for i in range(len(self.volcWE)):
            if isinstance(self.volcWE[i][0],float):
                volcWE_use.append(self.volcWE[i][0])

        volcWE_use = np.asarray(volcWE_use)[~np.isnan(volcWE_use)]
        return volcWE_use

    def volcIceDepth(self):
        depthIce = np.asarray(self.df_dens['iceDepth'])
        density = np.asarray(self.df_dens['density'])
        depthWE = np.asarray(depthIce*density/1000)
        volcWE_use = self.FindVolcErup()

        idx = []
        for i in range(len(volcWE_use)):
            idx.append(int(np.where(np.logical_and(depthWE>=volcWE_use[i]-0.1, depthWE<=volcWE_use[i]+0.1))[0][0]))

        volc_depthIce = depthIce[idx]
        return volc_depthIce

    def plotCore(self):

        variables = [self.df_d18O, self.df_DEP, self.df_ECM]
        varNames = ['d18O', 'cond', 'ECM']
        yLabels = ['d18O', 'DEP', 'ECM, Conductivity']

        volcDepthIce = self.volcIceDepth()
        ns = len(variables)
        delItems = []

        for i in range(len(variables)):
            if not isinstance(variables[i], pd.DataFrame):
                ns -= 1
                delItems.append(i)

        for i in range(len(delItems)):
            del variables[delItems[i]]; del varNames[delItems[i]]; del yLabels[delItems[i]]#yLabels.remove(yLabels[variables.index(a)])


        figCore, axCore = plt.subplots(ns, figsize=(18,10))
        for i in range(ns):
            axCore[i].plot(variables[i]['depth'], variables[i][varNames[i]])
            axCore[i].set(xlabel='Depth [m]', ylabel=yLabels[i], xlim=(variables[0]['depth'].iloc[0], variables[0]['depth'].iloc[-1]))
            for j in range(len(volcDepthIce)):
                axCore[i].axvline(x = volcDepthIce[j],color='k')

        figCore.tight_layout()
        figCore.savefig('Core'+self.name+'.png', dpi=600)

        return

    def getData_LakiToTambora(self):
        volcDepthIce = self.volcIceDepth()
        loc_Tambora = volcDepthIce[1]
        loc_Laki = volcDepthIce[2]

        variables = [self.df_d18O, self.df_DEP, self.df_ECM]
        varNames = ['d18O', 'cond', 'ECM']
        yLabels = ['d18O', 'DEP', 'ECM, Conductivity']

        ns = len(variables)
        delItems = []

        for i in range(len(variables)):
            if not isinstance(variables[i], pd.DataFrame):
                ns -= 1
                delItems.append(i)

        for i in range(len(delItems)):
            del variables[delItems[i]]; del varNames[delItems[i]]; del yLabels[delItems[i]]

        LT_idx = []; LT_idxPlot = []
        dfs_LT = []; dfs_LTplot = []
        fig2Core, ax2Core = plt.subplots(ns, figsize=(18,10))

        for i in range(ns):
            LT1 = np.nonzero(variables[i]['depth'].gt(loc_Tambora)); LT1plot = np.nonzero(variables[i]['depth'].gt(loc_Tambora-1))
            LT2 = np.nonzero(variables[i]['depth'].lt(loc_Laki)); LT2plot = np.nonzero(variables[i]['depth'].lt(loc_Laki+1))
            LT_idx.append(np.intersect1d(LT1,LT2)); LT_idxPlot.append(np.intersect1d(LT1plot,LT2plot))
            dfs_LT.append(variables[i].iloc[list(LT_idx[i])]); dfs_LTplot.append(variables[i].iloc[list(LT_idxPlot[i])])

            ax2Core[i].plot(dfs_LTplot[i]['depth'], dfs_LTplot[i][varNames[i]])
            ax2Core[i].set(xlabel='Depth [m]', ylabel=yLabels[i])
            ax2Core[i].axvline(x=loc_Tambora, color='k')
            ax2Core[i].axvline(x=loc_Laki, color='k')
#        print(variables[2].iloc[list(LT_idx)])
        fig2Core.tight_layout()
        fig2Core.savefig('Core_LT_'+self.name+'.png', dpi=600)
        return dfs_LT


#coreName = 'B18'
#core_B16 = Cores(name=coreName, df_dens=pd.read_excel('DepthDensity_Bcores_lowRes.xlsx', sheet_name=coreName, index=False),
#                 df_d18O = pd.read_excel('Depth_d18O__Bcores.xlsx', sheet_name=coreName, index=False),
#                 df_ECM = pd.read_excel('DepthECM__B16_B18_B21.xlsx', sheet_name=coreName, index=False),
#                 df_DEP = pd.read_excel('DepthDEP__Bcores.xlsx', sheet_name=coreName, index=False),
#                 df_chem = pd.read_excel('DepthChem__B16_B18_B21_B29.xlsx', sheet_name=coreName, index=False),
#                 volcWE = np.asarray(pd.read_excel('VolcanicEruptions__WE_Depth.xlsx', 'Sheet1', usecols=[coreName])))
# B16_dfs_LT = core_B16.plotCore()

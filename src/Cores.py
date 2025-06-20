import numpy as np
import pandas as pd
from pandas import ExcelWriter
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class Cores():

    """
        Class to handle North Greenland Transverse B cores from AWI.
        Needs data for at least density, d18O isotope and volcanic eruption depth locations
        and takes ECM, DEP, chemical data if available.

        Returns dataframe of measured data(depth density, d18O, ECM/DEP and volcanic eruptions)
        given between the two volcanic eruptions Laki and Tambora.

    """

    """
        Methods available:

            FindVolcErup(self):
                    Searches through the passed array of volcanic eruptions
                    (measured in W.E.) and finds all non NaN values and passes
                    floats to output.

            volcIceDepth(self):
                    Calculates the location of passed volcanic eruptions (W.E.)
                    in ice depth.

            plotCore(self, saveFig=False, plotFig=True)):
                    Plots (and saves) figure of entire core data, depth vs d18O data,
                    ECM/DEP/both and est. locations of eruptions.

            getData_LakiToTambora(self, saveFig=False, plotFig=True):
                Plots (and saves) figure of core data, depth vs d18O data,
                ECM/DEP/both, between (and a little further) the eruptions Laki and Tambora,
                along with est. locations of eruptions.

            SampleResolution(self, dataSlice):
                Gives the resolution of the given dataslice: all sample sizes, the unique sample sizes,
                the maximum and the minimum sample size.
    """
    def __init__(self, name, df_dens, df_d18O, df_ECM, df_DEP, volcWE):
        self.name = name
        self.df_dens = df_dens
        self.df_d18O = df_d18O
        self.df_ECM = df_ECM
        self.df_DEP = df_DEP
        self.volcWE = volcWE
        return

    def FindVolcErup(self):
        """
            Searches through the passed array of volcanic eruptions(measured in W.E.) and finds
            all non NaN values and passes floats to output.

                Arguments:
                ----------
                    None

                Returns:
                --------
                    volcWE_use      [array of floats] Containing all available estimated
                                                      volcanic eruptions for the given core.
        """
        volcWE_use = []

        for i in range(len(self.volcWE)):
            if isinstance(self.volcWE[i][0],float):
                volcWE_use.append(self.volcWE[i][0])

        volcWE_use = np.asarray(volcWE_use)[~np.isnan(volcWE_use)]
        return volcWE_use

    def volcIceDepth(self):
        """
            Calculates the location of passed volcanic eruptions (W.E.) in ice depth.

                Arguments:
                ---------
                    None

                Returns:
                --------
                    volc_depthIce   [array of floats] Containing ice depth location
                                    of all available estimated volcanic eruptions for
                                    given core.
        """
        depthIce = np.asarray(self.df_dens['iceDepth'])
        density = np.asarray(self.df_dens['density'])
        depthWE = np.asarray(self.df_dens['weDepth'])
        volcWE_use = self.FindVolcErup()

        idx = []
        for i in range(len(volcWE_use)):
            idx.append(int(np.where(np.logical_and(depthWE>=volcWE_use[i]-0.1, depthWE<=volcWE_use[i]+0.1))[0][0]))

        volc_depthIce = depthIce[idx]
        return volc_depthIce

    def plotCore(self, saveFig=False, plotFig=True):

        """
            Plots (and saves) figure of entire core data, depth vs d18O data, ECM/DEP/both and
            est. locations of eruptions.

                Arguments:
                ---------
                    saveFig:        [bool] Default: False. Save figure? Only if plotFig == True
                    plotFig:        [bool] Default: True. Plot figure?

                Returns:
                --------
                    None
        """
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

        if plotFig:
            import matplotlib as mpl
            mpl.rcParams['mathtext.fontset'] = 'stix'
            mpl.rcParams['font.size'] = 24
            mpl.rcParams['font.family'] = 'STIXGeneral'

            figCore, axCore = plt.subplots(ns, figsize=(12,10))
            axCore[0].set_title(self.name + ' full dataset')
            for i in range(ns):
                axCore[i].plot(variables[i]['depth'], variables[i][varNames[i]])
                axCore[i].set(xlabel='Depth [m]', ylabel=yLabels[i], xlim=(variables[0]['depth'].iloc[0], variables[0]['depth'].iloc[-1]))
                for j in range(len(volcDepthIce)):
                    axCore[i].axvline(x = volcDepthIce[j],color='k')

            figCore.tight_layout()
            if saveFig:
                figCore.savefig('Figures/Core'+self.name+'.eps')

        return

    def getData_LakiToTambora(self, saveFig=False, plotFig=True):
        """
            Plots (and saves) figure of core data, depth vs d18O data, ECM/DEP/both, between
            (and a little further) the eruptions Laki and Tambora, along with est. locations of eruptions.

            Arguments:
            ---------
                saveFig:        [bool] Default: False. Save figure? Only if plotFig == True
                plotFig:        [bool] Default: True. Plot figure?

            Returns:
            --------
                dfs_LT:         [list] List containing data, depth, d18O data, ECM/DEP/both, between
                                Laki and Tambora.
        """
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
        if plotFig:
            import matplotlib as mpl
            mpl.rcParams['mathtext.fontset'] = 'stix'
            mpl.rcParams['font.size'] = 28
            mpl.rcParams['font.family'] = 'STIXGeneral'
            fig2Core, ax2Core = plt.subplots(ns, figsize=(12,9),sharex=True)
            ax2Core[0].set_title(self.name, fontsize=40)


        for i in range(ns):
            LT1 = np.nonzero(variables[i]['depth'].gt(loc_Tambora)); LT1plot = np.nonzero(variables[i]['depth'].gt(loc_Tambora-1))
            LT2 = np.nonzero(variables[i]['depth'].lt(loc_Laki)); LT2plot = np.nonzero(variables[i]['depth'].lt(loc_Laki+1))
            LT_idx.append(np.intersect1d(LT1,LT2)); LT_idxPlot.append(np.intersect1d(LT1plot,LT2plot))
            dfs_LT.append(variables[i].iloc[list(LT_idx[i])]); dfs_LTplot.append(variables[i].iloc[list(LT_idxPlot[i])])
            if plotFig:
                ax2Core[i].plot(dfs_LTplot[i]['depth'], dfs_LTplot[i][varNames[i]])
                ax2Core[i].set(ylabel=yLabels[i])
                ax2Core[i].axvline(x=loc_Tambora, color='k', alpha=0.2, ls='--')
                ax2Core[i].axvline(x=loc_Laki, color='k', alpha=0.2, ls='--')

        if plotFig:
            ax2Core[-1].set(xlabel='Depth [m]')
            fig2Core.tight_layout()
            if saveFig:
                fig2Core.savefig('Figures/Core_LT_'+self.name+'.eps')
        return dfs_LT

    def SampleResolution(self, dataSlice):
        """
            Gives the resolution of the given dataslice: all sample sizes, the unique sample
            sizes, the maximum and the minimum sample size.

            Arguments:
            ----------
                dataslice:      [array of two floats] Containiig two depth values, describing the
                                desired depth array to examine resolution for.

            Returns:
            --------
                diff:           [list of 4 arrays] Contains all sample sizes, the unique
                                sample sizes, the maximum and the minimum sample size
        """
        depth = self.df_d18O['depth']
        depthSlice = depth[(depth >= dataSlice[0]) & (depth <= dataSlice[1])]

        diffDepth = np.round(np.diff(depthSlice), decimals=4)
        diffUnique = np.unique(diffDepth)
        diffMax = diffDepth.max()
        diffMin = diffDepth.min()

        diff = [diffDepth, diffUnique, diffMax, diffMin]
        return diff

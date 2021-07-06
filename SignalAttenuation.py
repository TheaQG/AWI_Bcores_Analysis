import numpy as np
import matplotlib.pyplot as plt
from MEM_class import MEM
from Decon import SpectralDecon
from scipy import signal
from scipy import interpolate

class Attenuation():
    '''
        Signal attenuation class to estimate attenuation of isotopic signal based
        on power spectral densities by DCT, NDCT, FFT or MEM methods.
        Can be used to estimate annual layer thickness of section in ice.
        Finds resonance peaks in frequency spectrum and sorts them. If wMin != 0,
        the maximal spectral peak found is found only at frequencies > wMin.


        Methods available:
            __init__(self, t, y, wMin, M = 30, N = 4000, PSD_Type = 'DCT'):
                    Initializes class intance given two required arguments (the
                    depth series) and four optional.

            __call__(self):
                    Calls the instance and determines the maximal resonance peak
                    and its magnitude.

            DCT_PSD(self):
                    Computes the PSD from DCT and returns frequency and magnitude (w,P).
                    (Depth series is interpolated if non-uniformly sampled.)

            NDCT_PSD(self):
                    Computes the PSD from NDCT and returns frequency and magnitude (w,P).

            FFT_PSD(self):
                    Computes the PSD from FFT and returns frequency and magnitude (w,P).
                    (Depth series is interpolated if non-uniformly sampled.)

            MEM_PSD(self):
                    Computes the PSD from MEM and returns frequency and magnitude (w,P).
                    (Depth series is interpolated if non-uniformly sampled.)

            PeakProps(self):
                    DCT/NDCT/FFT:
                        Computes the maximum spectral peak above the passed wMin and returns
                        the magnitude and frequency, fk, Pk.
                    mem:
                        Analytically computes location and Pk of all significant
                        peaks in spectrum

            findPeaks(self, wMin):
                    Finds peaks in spectrum with w >= wMin and their properties.
                    Sorts them peaks and properties from highest to lowest and saves
                    in matrix of data (dataHtoL) and list of keys (keys).

            interpData(self, DeltaInput=False, DeltaIn=0.):
                    Interpolates the passed depth series (t,y) with the smallest sampling
                    distance found in the original data set as new sampling size.
    '''
    def __init__(self, t, y, wMin=0, M = 30, N = 4000, PSD_Type = 'DCT'):
        '''
            Initializes class instance. When initializing: perform PSD computation.

            Arguments:
            ----------
                t:          [array of floats] Depth data of depth sereis
                y:          [array of floats] Isotopic value data of depth series.
                wMin:       [float] Minimum frequency to search for peaks above
                M:          [int] Number of poles in MEM analysis
                N:          [int] Number of frequency points in MEM analysis.
                PSD_Type:   [str] Type of spectral transform to use. Must be in ['DCT', 'NDCT', 'FFT', 'MEM'].

            Returns:
            --------
                None
        '''


        self.t = t
        self.y = y
        self.wMin = wMin
        self.M = M
        self.N = N
        self.PSD_Type = PSD_Type
        self.dtMean = np.mean(np.diff(t))
        self.w = 0.
        self.P = 0.


        if self.PSD_Type == 'DCT':
            self.DCT_PSD()
        elif self.PSD_Type == 'NDCT':
            self.NDCT_PSD()
        elif self.PSD_Type == 'FFT':
            self.FFT_PSD()
        elif self.PSD_Type == 'MEM':
            self.MEM_PSD()

        return

    def __call__(self):
        fk, Pk = self.PeakProps()

        return fk, Pk

    def DCT_PSD(self):
        '''
            Power spectral density estimation through discrete cosine transform.

            Arguments:
            ----------
                None

            Returns:
            --------
                None
        '''
        tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

        specInst = SpectralDecon(tHat, yHat, self.N, transType='DCT')

        w, P =specInst.dct_psd()

        self.w = w
        self.P = P
        return

    def NDCT_PSD(self):
        '''
            Power spectral density estimation through nonuniform discrete cosine transform.

            Arguments:
            ----------
                None

            Returns:
            --------
                None
        '''
        specInst = SpectralDecon(self.t, self.y, self.N, transType='NDCT')

        w, P =specInst.dct_psd()

        self.w = w
        self.P = P
        return

    def FFT_PSD(self):
        '''
            Power spectral density estimation through fast Fourier transform.

            Arguments:
            ----------
                None

            Returns:
            --------
                None
        '''
        tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

        specInst = SpectralDecon(tHat, yHat, self.N, transType='FFT')

        w, P =specInst.dct_psd()

        self.w = w
        self.P = P
        return

    def MEM_PSD(self):
        '''
            Power spectral density estimation through Burg's Maximum Entropy Method.

            Arguments:
            ----------
                None

            Returns:
            --------
                None
        '''
        tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

        specInst = MEM(tHat, yHat, self.M)

        w, P = specInst(tHat, yHat, self.M, self.N, view=False)

        self.w = w
        self.P = P
        return


    def PeakProps(self):
        '''
            Find spectral resonance peaks above certain frequency limit, wMin.

            Arguments:
            ----------
                None

            Returns:
            --------
                data_all:       [matrix of floats] Containing all peaks properties
                                sorted from highest peak to lowest peak.
                keys:           [lst of str] Containing the name of the peak properties.
        '''
        # print(self.wMin)
        w = self.w
        P = self.P

        if self.PSD_Type == 'MEM':
            tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

            specInst = MEM(tHat, yHat, self.M)

            roots = specInst._find_roots()
            fk0 = specInst._peak_frequencies(roots)
            Pk0 = specInst._peak_power(roots)

            fk = fk0[fk0>0]
            Pk = Pk0[fk0>0]
            keys = ['fk', 'Pk']
            data_all = np.asarray([fk, Pk])
        else:

            #wMax, PMax, idFull, wFull, PFull, area = self.findMaxPeak()
            #fk = wMax
            #Pk = area
            data_all, keys = self.findPeaks(wMin=self.wMin)

        return data_all, keys

    def findPeaks(self, wMin):
        '''
            Finding resonance spectral peaks and all their properties above frequency
            limit wMin. Sorts peaks from highest to lowest.

            Arguments:
            ----------
                None

            Returns:
            --------
                None
        '''
        w = self.w
        P = self.P

        Pnew = np.copy(P)
        Pnew[w<wMin] = 0

        peaks, props = signal.find_peaks(Pnew, height=0.025, width=0.001)

        if bool(props):
            data_all = np.zeros((len(props.keys())+1,len(peaks)))
            data_all[0,:] = peaks
            keys = ['peaks_idx']
            for key, i in zip(props.keys(), range(1,len(peaks))):
                keys.append(key)
                data_all[i,:] = props[key]

        dataLtoH = data_all[:, np.argsort(data_all[1,:])]
        dataHtoL = np.flip(dataLtoH, axis=1)

        return dataHtoL, keys

    def interpData(self, DeltaInput=False, DeltaIn=0.):
        '''
            Interpolation of (t,y) depth series.
            Uses cubic spline interpolation with a new sample size equal to the
            minimal sample size in the original data set.

            Arguments:
            ----------
                DeltaInput:     [bool] Default = False. Is inputted new sampling size? True/False.
                DeltaIn:        [float] Default = 0. If new sampling size is inputted,
                                        what is its value?

            Returns:
            --------
                None
        '''

        d = self.t
        x = self.y

        if DeltaInput:
            Delta = DeltaIn
        else:
            diff = np.diff(d)
            Delta = round(min(diff), 3)

        d_min = Delta * np.ceil(d[0]/Delta)
        d_max = Delta * np.floor(d[-1]/Delta)

        n = int(1 + (d_max - d_min)/Delta)

        j_arr = np.linspace(0,n,n)
        dhat0 = d_min + (j_arr - 1)*Delta

        f = interpolate.CubicSpline(d,x)

        xhat0 = f(dhat0)

        dhat = dhat0[(dhat0 >= min(d)) & (dhat0 <= max(d))]
        xhat = xhat0[(dhat0 >= min(d)) & (dhat0 <= max(d))]

        return dhat, xhat, Delta




class AnnualLayerThick():

    def __init__(self, t, y, lenSecs):
        self.t = t
        self.y = y
        self.lenSecs = lenSecs

        return

    def __call__(self):
        return


    def ALT_section(self, tSec, ySec, wMinDCT=0., wMinNDCT=0., wMinFFT=0.):
        AttInstDCT = Attenuation(tSec, ySec, wMinDCT, PSD_Type = 'DCT')
        dataDCT, keysDCT = AttInstDCT()

        AttInstNDCT = Attenuation(tSec, ySec, wMinNDCT, PSD_Type = 'NDCT')
        dataNDCT, keysNDCT = AttInstNDCT()

        AttInstFFT = Attenuation(tSec, ySec, wMinFFT, PSD_Type = 'FFT')
        dataFFT, keysFFT = AttInstFFT()

#        AttInst4 = Attenuation(depthTop, d18OTop, wMin, PSD_Type = 'MEM')
#        dataMEM, keysMEM = AttInst4()



        wDCT = AttInstDCT.w
        PDCT = AttInstDCT.P

        wNDCT = AttInstNDCT.w
        PNDCT = AttInstNDCT.P

        wFFT = AttInstFFT.w
        PFFT = AttInstFFT.P

        #wMEM = AttInst4.w
        #PMEM = AttInst4.P



        fkDCT = wDCT[dataDCT[0].astype(int)]
        hsDCT = dataDCT[1]


        fkNDCT = wNDCT[dataNDCT[0].astype(int)]
        hsNDCT = dataNDCT[1]


        fkFFT = wFFT[dataFFT[0].astype(int)]
        hsFFT = dataFFT[1]
        if len(fkDCT) == 0:
            fkDCTMax = -1
        else:
            fkDCTMax = fkDCT[0]
        if len(fkNDCT) == 0:
            fkNDCTMax = -1
        else:
            fkNDCTMax = fkNDCT[0]
        if len(fkFFT) == 0:
            fkFFTMax = -1
        else:
            fkFFTMax = fkFFT[0]

        fksMax = np.asarray([fkDCTMax, fkNDCTMax, fkFFTMax])

        return fksMax, fkDCT, fkNDCT, fkFFT

    def ALT_fullCore(self):
        t = self.t
        y = self.y
        lSecs = self.lenSecs

        tMin = np.floor(min(t))
        tMax = np.ceil(max(t))

        vals = np.arange(tMin,tMax,lSecs)

        wMinDCT_in = 0
        wMinNDCT_in = 0
        wMinFFT_in = 0

        fksMax_all = [np.array([0,0,0])]


        for i in range(len(vals) - 1):
            cutOff_High = vals[i]
            cutOff_Low = vals[i+1]

            tSec = t[(t >= cutOff_High) & (t <= cutOff_Low)]
            ySec = y[(t >= cutOff_High) & (t <= cutOff_Low)]

            fksMax, _, _, _ = self.ALT_section(tSec, ySec, wMinDCT_in, wMinNDCT_in, wMinFFT_in)

            fksMax_all.append(fksMax)
            fksMax_allArr = np.asarray(fksMax_all)

            if fksMax[0] < 0:
                wMinDCT_in =  fksMax_allArr[fksMax_allArr[:,0] > 0 , 0][-1] - 0.2 #fksMax_all[-2][0]
            else:
                wMinDCT_in = fksMax[0] - 0.2
            if fksMax[1] < 0:
                wMinNDCT_in = fksMax_allArr[fksMax_allArr[:,1] > 0 , 1][-1] - 0.2 #fksMax_all[-2][1]
            else:
                wMinNDCT_in = fksMax[1] - 0.1
            if fksMax[2] < 0:

                wMinFFT_in = fksMax_allArr[fksMax_allArr[:,2] > 0 , 2][-1] - 0.2 #fksMax_all[-2][2]
            else:
                 wMinFFT_in = fksMax[2] - 0.1

        fksMax_all = fksMax_all[1:]
        ls_all = 1/np.asarray(fksMax_all)

        def avg(a):
            return a[a > 0].mean()
        def std(a):
            return a[a>0].std()

        lMean = np.apply_along_axis(avg, 1, ls_all)
        lStd = np.apply_along_axis(std, 1, ls_all)

#        lMean = np.mean(ls_all, axis = 1)
#        lStd = np.std(ls_all, axis = 1)

        return np.asarray(fksMax_all), ls_all, lMean, lStd, vals


    def ALT_fullCore_seq(self, shift = 0.1, printItes = False):
        t = self.t
        y = self.y
        lsecs = self.lenSecs

        Tmax = max(t)

        N = int(np.floor((Tmax - t[0] - lsecs)/shift))
        tmins = np.arange(0, N + 1)*shift + t[0]
        tmaxs = tmins+lsecs


        wMinDCT_in = 0
        wMinNDCT_in = 0
        wMinFFT_in = 0
        tNew = []
        fksMax_all = [np.array([0,0,0])]

        print(f'Entire core: {t[0]}-{t[-1]} [m]\n')
        for i in range(len(tmins)):#idxs:

            cutOff_High = tmins[i]#t[i]
            cutOff_Low = tmaxs[i]#t[t <= (cutOff_High + lsecs)][-1]
            if printItes and i%100 == 0:
                print(f'ITERATION # {i}')
                print(f'Depth: {cutOff_High:.2f}-{cutOff_Low:.2f} [m]')

            tSec = t[(t >= cutOff_High) & (t <= cutOff_Low)]
            ySec = y[(t >= cutOff_High) & (t <= cutOff_Low)]
            tNew.append(tSec[0] + (tSec[-1] - tSec[0])/2)
            fksMax, _, _, _ = self.ALT_section(tSec, ySec, wMinDCT_in, wMinNDCT_in, wMinFFT_in)

            fksMax_all.append(fksMax)
            fksMax_allArr = np.asarray(fksMax_all)

            if fksMax[0] < 0:
                wMinDCT_in =  fksMax_allArr[fksMax_allArr[:,0] > 0 , 0][-1] - 0.2 #fksMax_all[-2][0]
            else:
                wMinDCT_in = fksMax[0] - 0.2
            if fksMax[1] < 0:
                wMinNDCT_in = fksMax_allArr[fksMax_allArr[:,1] > 0 , 1][-1] - 0.2 #fksMax_all[-2][1]
            else:
                wMinNDCT_in = fksMax[1] - 0.1
            if fksMax[2] < 0:

                wMinFFT_in = fksMax_allArr[fksMax_allArr[:,2] > 0 , 2][-1] - 0.2 #fksMax_all[-2][2]
            else:
                 wMinFFT_in = fksMax[2] - 0.1

        fksMax_all = fksMax_all[1:]
        ls_all = 1/np.asarray(fksMax_all)

        def avg(a):
            return a[a > 0].mean()
        def std(a):
            return a[a>0].std()

        lMean = np.apply_along_axis(avg, 1, ls_all)
        lStd = np.apply_along_axis(std, 1, ls_all)

#        lMean = np.mean(ls_all, axis = 1)
#        lStd = np.std(ls_all, axis = 1)


        return np.asarray(fksMax_all), ls_all, lMean, lStd, tNew




    # def closestTwo(self, lst, K):
    #     lstIn = np.copy(lst)
    #     lstUse = np.copy(lst)
    #
    #     Close1_ = lstUse[min(range(len(lstUse)), key = lambda i: abs(lstUse[i]-K))]
    #     id1 = np.where(lstIn == Close1_)[0][0]
    #
    #
    #     if Close1_ > K:
    #         idDel = np.where(lstIn > K)
    #         newLst = np.delete(lstUse, idDel)
    #     elif Close1_ < K:
    #         idDel = np.where(lstIn < K)
    #         newLst = np.delete(lstUse, idDel)
    #     #newLst = np.delete(lstUse, id1)
    #
    #     Close2_ = newLst[min(range(len(newLst)), key = lambda i: abs(newLst[i] - K))]
    #     id2 = np.where(lstIn == Close2_)[0][0]
    #
    #     if id1 < id2:
    #         Close2 = Close1_
    #         Close1 = Close2_
    #     elif id1 > id2:
    #         Close1 = Close1_
    #         Close2 = Close2_
    #
    #     vals = [Close1, Close2]
    #     ids = np.sort([id1, id2])
    #     return vals, ids
    #
    #
    # def findMaxPeak(self):
    #     w = self.w
    #     P = self.P
    #
    #     idPeaks, propsPeaks = signal.find_peaks(P)
    #     idTroughs, propsTroughs = signal.find_peaks(-P)
    #
    #     wPeaks = w[idPeaks]
    #     PPeaks = P[idPeaks]
    #
    #     wTroughs = w[idTroughs]
    #     PTroughs = P[idTroughs]
    #
    #     idMax = np.where(PPeaks == max(PPeaks))[0][0]
    #     wMax = wPeaks[idMax]
    #     PMax = PPeaks[idMax]
    #
    #
    #     wClose, idClose = self.closestTwo(wTroughs, wMax)
    #
    #     idFull =  np.where((w >= wTroughs[idClose[0]]) & (w <= wTroughs[idClose[1]]))
    #     wFull = w[idFull]
    #     PFull = P[idFull]
    #
    #     dw = np.mean(np.diff(wFull))
    #     area = sum(PFull)*dw
    #
    #     return wMax, PMax, idFull, wFull, PFull, area

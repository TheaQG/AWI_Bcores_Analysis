import numpy as np
import matplotlib.pyplot as plt
from MEM_class import MEM
from Decon import SpectralDecon
from scipy import signal
from scipy import interpolate

class Attenuation():

    def __init__(self, t, y, M = 30, N = 4000, PSD_Type = 'DCT'):

        self.t = t
        self.y = y
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
        tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

        specInst = SpectralDecon(tHat, yHat, self.N, transType='DCT')

        w, P =specInst.dct_psd()

        self.w = w
        self.P = P
        return

    def NDCT_PSD(self):
        specInst = SpectralDecon(self.t, self.y, self.N, transType='NDCT')

        w, P =specInst.dct_psd()

        self.w = w
        self.P = P
        return

    def FFT_PSD(self):
        tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

        specInst = SpectralDecon(tHat, yHat, self.N, transType='FFT')

        w, P =specInst.dct_psd()

        self.w = w
        self.P = P
        return

    def MEM_PSD(self):
        tHat, yHat, dt = self.interpData(DeltaInput = True, DeltaIn = self.dtMean)

        specInst = MEM(tHat, yHat, self.M)

        w, P = specInst(tHat, yHat, self.M, self.N, view=False)

        self.w = w
        self.P = P
        return


    def PeakProps(self):

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
            keys = ['fk', 'peak_heights','Pk']
            data_all = np.asarray([fk, Pk])
        else:

            #wMax, PMax, idFull, wFull, PFull, area = self.findMaxPeak()
            #fk = wMax
            #Pk = area
            data_all, keys = self.findPeaks()

        return data_all, keys


    def closestTwo(self, lst, K):
        lstIn = np.copy(lst)
        lstUse = np.copy(lst)

        Close1_ = lstUse[min(range(len(lstUse)), key = lambda i: abs(lstUse[i]-K))]
        id1 = np.where(lstIn == Close1_)[0][0]


        if Close1_ > K:
            idDel = np.where(lstIn > K)
            newLst = np.delete(lstUse, idDel)
        elif Close1_ < K:
            idDel = np.where(lstIn < K)
            newLst = np.delete(lstUse, idDel)
        #newLst = np.delete(lstUse, id1)

        Close2_ = newLst[min(range(len(newLst)), key = lambda i: abs(newLst[i] - K))]
        id2 = np.where(lstIn == Close2_)[0][0]

        if id1 < id2:
            Close2 = Close1_
            Close1 = Close2_
        elif id1 > id2:
            Close1 = Close1_
            Close2 = Close2_

        vals = [Close1, Close2]
        ids = np.sort([id1, id2])
        return vals, ids

    def interpData(self, DeltaInput=False, DeltaIn=0.):

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

    def findMaxPeak(self):
        w = self.w
        P = self.P

        idPeaks, propsPeaks = signal.find_peaks(P)
        idTroughs, propsTroughs = signal.find_peaks(-P)

        wPeaks = w[idPeaks]
        PPeaks = P[idPeaks]

        wTroughs = w[idTroughs]
        PTroughs = P[idTroughs]

        idMax = np.where(PPeaks == max(PPeaks))[0][0]
        wMax = wPeaks[idMax]
        PMax = PPeaks[idMax]


        wClose, idClose = self.closestTwo(wTroughs, wMax)

        idFull =  np.where((w >= wTroughs[idClose[0]]) & (w <= wTroughs[idClose[1]]))
        wFull = w[idFull]
        PFull = P[idFull]

        dw = np.mean(np.diff(wFull))
        area = sum(PFull)*dw

        return wMax, PMax, idFull, wFull, PFull, area


    def findPeaks(self):
        w = self.w
        P = self.P

        peaks, props = signal.find_peaks(P, height=0.025, width=0.001)

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

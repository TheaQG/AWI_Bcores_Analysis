import numpy as np
import scipy as sp
import pandas as pd
import math


class transforms():

    def __init__(self,t_data = np.array([]), x_data = np.array([]), w_data = np.array([]), \
                X_data = np.array([]), N = 200):

        self.t_data = t_data
        self.x_data = x_data
        self.w_data = w_data
        self.X_data = X_data
        self.N = N
        return

    def NDCT(self):

        def Dmat_ij(w_i, t_j, N):
            D_ij = np.cos(np.pi * w_i * (t_j + 1/(2 * N)))
            return D_ij

        def Dmat(w, t):
            D = np.cos(np.outer(np.pi * w, t + 1/(2*self.N)))
            return D

        x = self.x_data
        t = self.t_data

        N = self.N

        dt = (t[-1] - t[0]) / len(t)

        if np.size(self.w_data) == 0:

            w = np.fft.fftfreq(2*N, dt)[:(2*N)//2]
        else:
            w = self.w_data

        D = Dmat(w, t)

        X = np.dot(D,x)

        X[0] = X[0] * np.sqrt(1/(4*N))
        X[1:] = X[1:] * np.sqrt(1/(2*N))
        return w, X

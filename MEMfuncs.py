# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 16:27:02 2019

@author: thea
"""
import numpy as np
import matplotlib.pyplot as plt

#file = open("/home/thea/Documents/PUK_blok1_2019/Data_Files/EGRIP/egrip22pbag_401_420.txt","r")

def load_files(file_in, O17 = False):
    runID = []
    bag = []
    sample = []
    depth = []
    d18 = []
    d17 = []
    dD = []

    with file_in as infile:
        header = infile.readline()
        counter = 0
        for line in infile:
            columns = line.strip().split()
            counter += 1
            runID.append(columns[0])
            bag.append(int(columns[1]))
            sample.append(int(columns[2]))
            depth.append(float(columns[3]))

            if O17:
                d17.append(float(columns[4]))
                d18.append(float(columns[5]))
                dD.append(float(columns[6]))
            else:
                d18.append(float(columns[4]))
                dD.append(float(columns[5]))
    if O17:
        return depth, d17, d18, dD
    else:
        return depth, d18, dD


# data = load_files(file)
# t_data1 = data[0]
# y1 = data[1]
# dD = data[2]
# #%%


def coef(y_in, t_in, M):
    y = y_in
    y -= np.average(y)
    N = np.size(y)

    P0 = np.sum(y**2)/N

    m = 0

    b1 = np.copy(y[0:-1])
    b2 = np.copy(y[1:N])

    nom = np.sum(b1[0:N-1]*b2[0:N-1])
    denom = np.sum(b1[0:N-1]**2 + b2[0:N-1]**2)

    a = np.array(())
    ref = np.array(())

    a = np.append(a,2*nom/denom)
    ref = np.append(ref,2*nom/denom)

    P = np.array(())

    P = np.append(P,P0*(1 - a[m]**2))

    for m in range(1,M):
        aa = np.array(())
        for k in range(0,m):
            aa = np.append(aa,a[k])

        for k in range(0,N - (m+1)):
            b1[k] = b1[k] - aa[m-1]*b2[k]
            b2[k] = b2[k+1] - aa[m-1]*b1[k+1]

        nom = np.sum(b1[0:N - (m+1)] * b2[0:N - (m+1)])
        denom = np.sum(b1[0:N - (m+1)]**2 + b2[0:N - (m+1)]**2)

        a = np.append(a, 2*nom/denom)
        ref = np.append(ref, 2*nom/denom)
        P = np.append(P,P[m-1]*(1 - a[m]**2))

        for kk in range(0,m):
            a[kk] = aa[kk] - a[m]*aa[m-kk-1]

    return (P, a, ref)

def power(y_in, t_in, N, M):

    P = coef(y_in, t_in, M)[0]
    a = coef(y_in, t_in, M)[1]
    t_data = t_in

    Pm = P[-1]
    dt = t_data[1] - t_data[0]

    i = 0. + 1.j

    nq = 1/(2*dt)

    f = np.arange(N/2)*2*nq/N

    nomF = Pm*dt

    sumF = 0
    for k in range(1,len(a)+1):
        sumF += a[k-1]*np.exp(-2*np.pi*i*f*k*dt)

    denomF =np.abs(1-sumF)**2

    Pf = nomF/denomF

    return (f, Pf)

#%%

def MEM_func(y_in, t_in, M):
    y = y_in
    t_data = t_in
    N = len(y_in)

    Power = power(y, t_data, N, M)
    coefs = coef(y,t_data,M)

    return Power, coefs


#MEM_func(y1, t_data1)

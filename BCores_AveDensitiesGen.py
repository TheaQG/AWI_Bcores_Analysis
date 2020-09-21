
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import ExcelWriter
import openpyxl

def calc_AveDens(path, filename, N_dataPoints):
    data = np.genfromtxt(path + filename, delimiter='\t', skip_header=16, names=True, dtype={'names': ('Age [ka BP]', 'Age [a AD]', 'Depth ice/snow [m]', 'Depth w.e. top [m]', 'Acc rate ice [kg/m**2/a]', 'δ18O H2O [‰ SMOW]'),
          'formats': (np.float, np.float, np.float, np.float, np.float, np.float)},)
    names = data.dtype.names

    N_Ave = np.int(len(data[names[2]])/N_dataPoints)
    print(N_Ave)

    depthIce = data[names[2]]; deltaDepthIce = np.diff(depthIce)
    depthIceMat = np.pad(depthIce.astype(float), (0, N_Ave - depthIce.size%N_Ave), mode='constant', constant_values=np.NaN).reshape(-1, N_Ave)
    depthIceNew = depthIceMat[:,0] + (depthIceMat[:,-1] - depthIceMat[:,0])/2
    depthIceNew[np.isnan(depthIceNew)] = max(depthIceNew) + np.mean(deltaDepthIce)

    depthWE = data[names[3]]; deltaDepthWE = np.diff(depthWE)
    depthWEMat = np.pad(depthWE.astype(float), (0, N_Ave - depthWE.size%N_Ave), mode='constant', constant_values=np.NaN).reshape(-1, N_Ave)
    depthWENew = depthWEMat[:,0] + (depthWEMat[:,-1] - depthWEMat[:,0])/2
    depthWENew[np.isnan(depthWENew)] = max(depthWENew) + np.mean(deltaDepthWE)

    dens = deltaDepthWE / deltaDepthIce

    densAve = np.nanmean(np.pad(dens.astype(float), (0, N_Ave - dens.size%N_Ave), mode='constant', constant_values=np.NaN).reshape(-1, N_Ave), axis=1)
    if N_Ave > 1:
        densSTD = np.nanstd(np.pad(dens.astype(float), (0, N_Ave - dens.size%N_Ave), mode='constant', constant_values=np.NaN).reshape(-1, N_Ave), axis=1)
    else:
        densSTD = np.zeros_like(densAve)

    if len(densAve) != len(depthIceNew) and len(densAve) != len(depthWENew):
        depthIceNew = depthIceNew[1:]
        depthWENew = depthWENew[1:]
    return densAve*1000, densSTD*1000, depthIceNew, depthWENew, N_Ave



'''
    Use above function to calculate densities for all B-cores.
'''

dens = []
path = "../Data/datasets/B_cores_AWI/Densities/"
filenames_all = ["B29_2_acc_rate_d18O.tab", "B17_2_acc_rate_d18O.tab", "B18_2_acc_rate_d18O.tab", "B19_2_acc_rate_d18O.tab", "B20_2_acc_rate_d18O.tab", "B21_2_acc_rate_d18O.tab", "B22_2_acc_rate_d18O.tab",
            "B23_2_acc_rate_d18O.tab", "B26_2_acc_rate_d18O.tab", "B27_2_acc_rate_d18O.tab", "B28_3_acc_rate_d18O.tab", "B29_2_acc_rate_d18O.tab", "B30_2_acc_rate_d18O.tab"]
#filenames = ["B16_2_acc_rate_d18O.tab", "B17_2_acc_rate_d18O.tab", "B18_2_acc_rate_d18O.tab", "B21_2_acc_rate_d18O.tab"]
filenames_chosen = filenames_all
i = 0
N_points = [100,50,70,100,70,100,70,80,70,70,50,80,80]
for file in filenames_chosen:
    dens.append(calc_AveDens(path, file, N_points[i]))
    i += 1



'''
    Create and save files with depth-density profile data
'''

coreNames_all = ['B16', 'B17','B18', 'B19', 'B20', 'B21', 'B22', 'B23', 'B26', 'B27', 'B28', 'B29', 'B30']
coreNames_chosen = coreNames_all

f_save = 'DepthDensity_Bcores_lowResAve.xlsx'
writer = pd.ExcelWriter(f_save, engine='xlsxwriter')
writer.save()
df = pd.DataFrame({'iceDepth': dens[0][2], 'density': dens[0][0], 'STD': dens[0][1], 'weDepth': dens[0][3], 'N_Ave': dens[0][4]*np.ones_like(dens[0][0])})
writer = pd.ExcelWriter(f_save, engine='openpyxl')
df.to_excel(writer, sheet_name=coreNames_chosen[0], index=False)
writer.save()

for i in range(1, len(filenames_chosen)):
    df = pd.DataFrame({'iceDepth': dens[i][2], 'density': dens[i][0], 'STD': dens[i][1], 'weDepth': dens[i][3], 'N_Ave': dens[i][4]*np.ones_like(dens[i][0])})
    writer = pd.ExcelWriter(f_save, engine='openpyxl', mode='a')
    df.to_excel(writer, sheet_name=coreNames_chosen[i], index=False)
    writer.save()

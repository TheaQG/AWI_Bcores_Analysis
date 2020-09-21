import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import ExcelWriter
import openpyxl

'''
    Define function to compute density from ice and w.e. depths
'''
def calc_dens(path, filename):
    data = np.genfromtxt(path + filename, delimiter='\t', skip_header=16, names=True, dtype={'names': ('Age [ka BP]', 'Age [a AD]', 'Depth ice/snow [m]', 'Depth w.e. top [m]', 'Acc rate ice [kg/m**2/a]', 'δ18O H2O [‰ SMOW]'),
          'formats': (np.float, np.float, np.float, np.float, np.float, np.float)},)
    names = data.dtype.names

    depthIce = data[names[2]]; deltaDepthIce = np.diff(depthIce)
    depthWE = data[names[3]]; deltaDepthWE = np.diff(depthWE)
    dens = deltaDepthWE / deltaDepthIce
    #print(depthIce, deltaDepthIce)
    return dens*1000, depthIce[:-1], depthWE[:-1]


'''
    Use above function to calculate densities for all B-cores.
'''

dens = []
path = "../Data/datasets/B_cores_AWI/Densities/"
filenames_all = ["B16_2_acc_rate_d18O.tab", "B17_2_acc_rate_d18O.tab", "B18_2_acc_rate_d18O.tab", "B19_2_acc_rate_d18O.tab", "B20_2_acc_rate_d18O.tab", "B21_2_acc_rate_d18O.tab", "B22_2_acc_rate_d18O.tab",
            "B23_2_acc_rate_d18O.tab", "B26_2_acc_rate_d18O.tab", "B27_2_acc_rate_d18O.tab", "B28_3_acc_rate_d18O.tab", "B29_2_acc_rate_d18O.tab", "B30_2_acc_rate_d18O.tab"]
#filenames = ["B16_2_acc_rate_d18O.tab", "B17_2_acc_rate_d18O.tab", "B18_2_acc_rate_d18O.tab", "B21_2_acc_rate_d18O.tab"]
filenames_chosen = filenames_all

for file in filenames_chosen:
    dens.append(calc_dens(path, file))



'''
    Create and save files with depth-density profile data
'''

coreNames_all = ['B16', 'B17','B18', 'B19', 'B20', 'B21', 'B22', 'B23', 'B26', 'B27', 'B28', 'B29', 'B30']
coreNames_chosen = coreNames_all

f_save = 'DepthDensity_Bcores_lowRes.xlsx'
writer = pd.ExcelWriter(f_save, engine='xlsxwriter')
writer.save()

df = pd.DataFrame({'iceDepth': dens[0][1], 'density': dens[0][0], 'weDepth': dens[0][2]})
writer = pd.ExcelWriter(f_save, engine='openpyxl')
df.to_excel(writer, sheet_name=coreNames_chosen[0], index=False)
writer.save()

for i in range(1, len(filenames_chosen)):
    df = pd.DataFrame({'iceDepth': dens[i][1], 'density': dens[i][0], 'weDepth': dens[i][2]})
    writer = pd.ExcelWriter(f_save, engine='openpyxl', mode='a')
    df.to_excel(writer, sheet_name=coreNames_chosen[i], index=False)
    writer.save()

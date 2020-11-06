import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import ExcelWriter
import openpyxl

'''
    Loading d18O data and generating cleaned up files.
'''

path = "../../Data/datasets/B_cores_AWI/"
filename = "NGT_d18O_highres.xlsx"

NGT_d18O_highres_data = pd.read_excel(path + filename, sheet_name=None)

headers = list(NGT_d18O_highres_data.keys())

df_B16_d18O = NGT_d18O_highres_data[headers[0]]
df_B17_d18O = NGT_d18O_highres_data[headers[1]]
df_B18_d18O = NGT_d18O_highres_data[headers[2]]
df_B19_d18O = NGT_d18O_highres_data[headers[3]]
df_B20_d18O = NGT_d18O_highres_data[headers[4]]
df_B21_d18O = NGT_d18O_highres_data[headers[5]]
df_B22_d18O = NGT_d18O_highres_data[headers[6]]
df_B23_d18O = NGT_d18O_highres_data[headers[7]]
df_B26_d18O = NGT_d18O_highres_data[headers[8]]
df_B27_d18O = NGT_d18O_highres_data[headers[9]]
df_B28_d18O = NGT_d18O_highres_data[headers[10]]
df_B29_d18O = NGT_d18O_highres_data[headers[11]]
df_B30_d18O = NGT_d18O_highres_data[headers[12]]

coreNames_d18O = ['B16', 'B17', 'B18', 'B19', 'B20', 'B21', 'B22', 'B23', 'B26', 'B27', 'B28', 'B29', 'B30']

path_save = "../../Data/datasets/B_cores_AWI/AWI_Bcores__Cleanded_xlsx/"

f_d18O_save = path_save + 'Depth_d18O.xlsx'
writer = pd.ExcelWriter(f_d18O_save, engine='xlsxwriter')
writer.save()

df16_d18O = pd.DataFrame({'depth': df_B16_d18O[df_B16_d18O.keys()[0]], 'd18O': df_B16_d18O[df_B16_d18O.keys()[2]]})
writer = pd.ExcelWriter(f_d18O_save, engine='openpyxl')
df16_d18O.to_excel(writer, sheet_name=coreNames_d18O[0], index=False)
writer.save()

df17_d18O = pd.DataFrame({'depth': df_B17_d18O[df_B17_d18O.keys()[0]], 'd18O': df_B17_d18O[df_B17_d18O.keys()[2]]})
df18_d18O = pd.DataFrame({'depth': df_B18_d18O[df_B18_d18O.keys()[1]], 'd18O': df_B18_d18O[df_B18_d18O.keys()[3]]})
df19_d18O = pd.DataFrame({'depth': df_B19_d18O[df_B19_d18O.keys()[0]], 'd18O': df_B19_d18O[df_B19_d18O.keys()[2]]})
df20_d18O = pd.DataFrame({'depth': df_B20_d18O[df_B20_d18O.keys()[1]], 'd18O': df_B20_d18O[df_B20_d18O.keys()[3]]})
df21_d18O = pd.DataFrame({'depth': df_B21_d18O[df_B21_d18O.keys()[0]], 'd18O': df_B21_d18O[df_B21_d18O.keys()[2]]})
df22_d18O = pd.DataFrame({'depth': df_B22_d18O[df_B22_d18O.keys()[0]], 'd18O': df_B22_d18O[df_B22_d18O.keys()[2]]})
df23_d18O = pd.DataFrame({'depth': df_B23_d18O[df_B23_d18O.keys()[0]], 'd18O': df_B23_d18O[df_B23_d18O.keys()[2]]})
df26_d18O = pd.DataFrame({'depth': df_B26_d18O[df_B26_d18O.keys()[0]], 'd18O': df_B26_d18O[df_B26_d18O.keys()[2]]})
df27_d18O = pd.DataFrame({'depth': df_B27_d18O[df_B27_d18O.keys()[0]], 'd18O': df_B27_d18O[df_B27_d18O.keys()[1]]})
df28_d18O = pd.DataFrame({'depth': df_B28_d18O[df_B28_d18O.keys()[0]], 'd18O': df_B28_d18O[df_B28_d18O.keys()[2]]})
df29_d18O = pd.DataFrame({'depth': df_B29_d18O[df_B29_d18O.keys()[0]], 'd18O': df_B29_d18O[df_B29_d18O.keys()[2]]})
df30_d18O = pd.DataFrame({'depth': df_B30_d18O[df_B30_d18O.keys()[0]], 'd18O': df_B30_d18O[df_B30_d18O.keys()[2]]})

with pd.ExcelWriter(f_d18O_save, engine='openpyxl', mode='a') as writer:
    df17_d18O.to_excel(writer, sheet_name=coreNames_d18O[1], index=False)
    df18_d18O.to_excel(writer, sheet_name=coreNames_d18O[2], index=False)
    df19_d18O.to_excel(writer, sheet_name=coreNames_d18O[3], index=False)
    df20_d18O.to_excel(writer, sheet_name=coreNames_d18O[4], index=False)
    df21_d18O.to_excel(writer, sheet_name=coreNames_d18O[5], index=False)
    df22_d18O.to_excel(writer, sheet_name=coreNames_d18O[6], index=False)
    df23_d18O.to_excel(writer, sheet_name=coreNames_d18O[7], index=False)
    df26_d18O.to_excel(writer, sheet_name=coreNames_d18O[8], index=False)
    df27_d18O.to_excel(writer, sheet_name=coreNames_d18O[9], index=False)
    df28_d18O.to_excel(writer, sheet_name=coreNames_d18O[10], index=False)
    df29_d18O.to_excel(writer, sheet_name=coreNames_d18O[11], index=False)
    df30_d18O.to_excel(writer, sheet_name=coreNames_d18O[12], index=False)

#####################################################################################################################################33

'''
    Loading ECM, DEP and chem data and generating cleaned up files.
'''

path2 = "../../Data/datasets/B_cores_AWI/"
filename2 = "chemistry_conductivity_NGT93_95_all.xlsx"

B_Cores_chem_cond__data = pd.read_excel(path2 + filename2, sheet_name=None)

headers2 = list(B_Cores_chem_cond__data.keys())

df_B16_chem = B_Cores_chem_cond__data[headers2[0]]
df_B16_ECM = B_Cores_chem_cond__data[headers2[1]]
df_B17_DEP = B_Cores_chem_cond__data[headers2[2]]
df_B18_chem = B_Cores_chem_cond__data[headers2[3]]
df_B18_DEP = B_Cores_chem_cond__data[headers2[4]]
df_B18_ECM = B_Cores_chem_cond__data[headers2[5]]
df_B19_DEP = B_Cores_chem_cond__data[headers2[6]]
df_B20_DEP = B_Cores_chem_cond__data[headers2[7]]
df_B21_ECM = B_Cores_chem_cond__data[headers2[8]]
df_B21_chem = B_Cores_chem_cond__data[headers2[9]]
df_B22_DEP = B_Cores_chem_cond__data[headers2[10]]
df_B23_DEP = B_Cores_chem_cond__data[headers2[11]]
df_B26_DEP = B_Cores_chem_cond__data[headers2[12]]
df_B27_DEP = B_Cores_chem_cond__data[headers2[13]]
df_B28_DEP = B_Cores_chem_cond__data[headers2[14]]
df_B29_chem = B_Cores_chem_cond__data[headers2[15]]


'''
    Generating files with depth and DEP data
'''

coreNames_DEP = ['B17', 'B18', 'B19', 'B20', 'B22', 'B23', 'B26', 'B27', 'B28']

f_DEP_save = path_save + 'DepthDEP.xlsx'
writer = pd.ExcelWriter(f_DEP_save, engine='xlsxwriter')
writer.save()

df17_DEP = pd.DataFrame({'depth': df_B17_DEP[df_B17_DEP.keys()[0]], 'g_250k': df_B17_DEP[df_B17_DEP.keys()[1]], 'cp_250k': df_B17_DEP[df_B17_DEP.keys()[2]]})
writer = pd.ExcelWriter(f_DEP_save, engine='openpyxl')
df17_DEP.to_excel(writer, sheet_name=coreNames_DEP[0], index=False)
writer.save()

df18_DEP = pd.DataFrame({'depth': df_B18_DEP[df_B18_DEP.keys()[4]], 'g_250k': df_B18_DEP[df_B18_DEP.keys()[5]], 'cp_250k': df_B18_DEP[df_B18_DEP.keys()[6]]})
df19_DEP = pd.DataFrame({'depth': df_B19_DEP[df_B19_DEP.keys()[0]], 'g_250k': df_B19_DEP[df_B19_DEP.keys()[1]], 'cp_250k': df_B19_DEP[df_B19_DEP.keys()[2]]})
df20_DEP = pd.DataFrame({'depth': df_B20_DEP[df_B20_DEP.keys()[0]], 'dens': df_B20_DEP[df_B20_DEP.keys()[1]], 'cond': df_B20_DEP[df_B20_DEP.keys()[2]]})
df22_DEP = pd.DataFrame({'depth': df_B22_DEP[df_B22_DEP.keys()[0]], 'dens': df_B22_DEP[df_B22_DEP.keys()[1]], 'cond': df_B22_DEP[df_B22_DEP.keys()[2]]})
df23_DEP = pd.DataFrame({'depth': df_B23_DEP[df_B23_DEP.keys()[0]], 'dens': df_B23_DEP[df_B23_DEP.keys()[1]], 'cond': df_B23_DEP[df_B23_DEP.keys()[2]]})
df26_DEP = pd.DataFrame({'depth': df_B26_DEP[df_B26_DEP.keys()[1]], 'sig_250k': df_B26_DEP[df_B26_DEP.keys()[2]]})
df27_DEP = pd.DataFrame({'depth': df_B27_DEP[df_B27_DEP.keys()[0]], 'sig_250k': df_B27_DEP[df_B27_DEP.keys()[1]]})
df28_DEP = pd.DataFrame({'depth': df_B28_DEP[df_B28_DEP.keys()[0]], 'sig_250k': df_B28_DEP[df_B28_DEP.keys()[1]]})


with pd.ExcelWriter(f_DEP_save, engine='openpyxl', mode='a') as writer:
    df18_DEP.to_excel(writer, sheet_name=coreNames_DEP[1], index=False)
    df19_DEP.to_excel(writer, sheet_name=coreNames_DEP[2], index=False)
    df20_DEP.to_excel(writer, sheet_name=coreNames_DEP[3], index=False)
    df22_DEP.to_excel(writer, sheet_name=coreNames_DEP[4], index=False)
    df23_DEP.to_excel(writer, sheet_name=coreNames_DEP[5], index=False)
    df26_DEP.to_excel(writer, sheet_name=coreNames_DEP[6], index=False)
    df27_DEP.to_excel(writer, sheet_name=coreNames_DEP[7], index=False)
    df28_DEP.to_excel(writer, sheet_name=coreNames_DEP[8], index=False)


#Creating file with only actual useful DEP data
coreNames_DEP2 = ['B18', 'B19', 'B20', 'B22', 'B23']

f_DEP_save2 = path_save + 'DepthDEP__Clean.xlsx'
writer = pd.ExcelWriter(f_DEP_save2, engine='xlsxwriter')
writer.save()

df18_DEP = pd.DataFrame({'depth': df_B18_DEP[df_B18_DEP.keys()[4]], 'cond': df_B18_DEP[df_B18_DEP.keys()[6]]})
writer = pd.ExcelWriter(f_DEP_save, engine='openpyxl')
df18_DEP.to_excel(writer, sheet_name=coreNames_DEP[0], index=False)
writer.save()

df19_DEP = pd.DataFrame({'depth': df_B19_DEP[df_B19_DEP.keys()[0]], 'cond': df_B19_DEP[df_B19_DEP.keys()[2]]})
df20_DEP = pd.DataFrame({'depth': df_B20_DEP[df_B20_DEP.keys()[0]], 'cond': df_B20_DEP[df_B20_DEP.keys()[2]]})
df22_DEP = pd.DataFrame({'depth': df_B22_DEP[df_B22_DEP.keys()[0]], 'cond': df_B22_DEP[df_B22_DEP.keys()[2]]})
df23_DEP = pd.DataFrame({'depth': df_B23_DEP[df_B23_DEP.keys()[0]], 'cond': df_B23_DEP[df_B23_DEP.keys()[2]]})


with pd.ExcelWriter(f_DEP_save2, engine='openpyxl', mode='a') as writer:
    df19_DEP.to_excel(writer, sheet_name=coreNames_DEP2[1], index=False)
    df20_DEP.to_excel(writer, sheet_name=coreNames_DEP2[2], index=False)
    df22_DEP.to_excel(writer, sheet_name=coreNames_DEP2[3], index=False)
    df23_DEP.to_excel(writer, sheet_name=coreNames_DEP2[4], index=False)

'''
    Generating files with depth and chemical measurements
'''

coreNames_chem = ['B16', 'B18', 'B21', 'B29']


f_chem_save = path_save + 'DepthChem.xlsx'
writer = pd.ExcelWriter(f_chem_save, engine='xlsxwriter')
writer.save()

df16_chem = pd.DataFrame({'depth': df_B16_chem[df_B16_chem.keys()[0]], 'MSA': df_B16_chem[df_B16_chem.keys()[3]], 'Cl': df_B16_chem[df_B16_chem.keys()[4]],
                        'Clx': df_B16_chem[df_B16_chem.keys()[5]], 'Br': df_B16_chem[df_B16_chem.keys()[6]], 'NO3': df_B16_chem[df_B16_chem.keys()[7]],
                        'SO4': df_B16_chem[df_B16_chem.keys()[8]], 'SO4nss': df_B16_chem[df_B16_chem.keys()[9]], 'Na': df_B16_chem[df_B16_chem.keys()[10]],
                        'NH4': df_B16_chem[df_B16_chem.keys()[11]], 'K': df_B16_chem[df_B16_chem.keys()[12]], 'Mg': df_B16_chem[df_B16_chem.keys()[13]],
                        'Ca': df_B16_chem[df_B16_chem.keys()[14]], 'conductivity': df_B16_chem[df_B16_chem.keys()[18]], 'pH': df_B16_chem[df_B16_chem.keys()[17]]})
writer = pd.ExcelWriter(f_chem_save, engine='openpyxl')
df16_chem.to_excel(writer, sheet_name=coreNames_chem[0], index=False)
writer.save()

df18_chem = pd.DataFrame({'depth': df_B18_chem[df_B18_chem.keys()[0]], 'MSA': df_B18_chem[df_B18_chem.keys()[4]], 'Cl': df_B18_chem[df_B18_chem.keys()[5]],
                        'Clx': df_B18_chem[df_B18_chem.keys()[6]], 'Br': df_B18_chem[df_B18_chem.keys()[7]], 'NO3': df_B18_chem[df_B18_chem.keys()[8]],
                        'SO4': df_B18_chem[df_B18_chem.keys()[9]], 'SO4nss': df_B18_chem[df_B18_chem.keys()[10]], 'Na': df_B18_chem[df_B18_chem.keys()[11]],
                        'NH4': df_B18_chem[df_B18_chem.keys()[12]], 'K': df_B18_chem[df_B18_chem.keys()[13]], 'Mg': df_B18_chem[df_B18_chem.keys()[14]],
                        'Ca': df_B18_chem[df_B18_chem.keys()[15]], 'F': df_B18_chem[df_B18_chem.keys()[16]]})
df21_chem = pd.DataFrame({'depth': df_B21_chem[df_B21_chem.keys()[0]], 'MSA': df_B21_chem[df_B21_chem.keys()[4]], 'Cl': df_B21_chem[df_B21_chem.keys()[5]],
                        'Clx': df_B21_chem[df_B21_chem.keys()[6]], 'Br': df_B21_chem[df_B21_chem.keys()[7]], 'NO3': df_B21_chem[df_B21_chem.keys()[8]],
                        'SO4': df_B21_chem[df_B21_chem.keys()[9]], 'SO4nss': df_B21_chem[df_B21_chem.keys()[10]], 'Na': df_B21_chem[df_B21_chem.keys()[11]],
                        'NH4': df_B21_chem[df_B21_chem.keys()[12]], 'K': df_B21_chem[df_B21_chem.keys()[13]], 'Mg': df_B21_chem[df_B21_chem.keys()[14]],
                        'Ca': df_B21_chem[df_B21_chem.keys()[15]], 'F': df_B21_chem[df_B21_chem.keys()[16]], 'pH': df_B21_chem[df_B21_chem.keys()[17]]})
df29_chem = pd.DataFrame({'depth': df_B29_chem[df_B29_chem.keys()[0]], 'Ca': df_B29_chem[df_B29_chem.keys()[1]], 'NH4': df_B29_chem[df_B29_chem.keys()[3]],
                        'CH2O': df_B29_chem[df_B29_chem.keys()[5]], 'H2O2': df_B29_chem[df_B29_chem.keys()[7]], 'density': df_B29_chem[df_B29_chem.keys()[11]],
                        'accumulation': df_B29_chem[df_B29_chem.keys()[1]]})

with pd.ExcelWriter(f_chem_save, engine='openpyxl', mode='a') as writer:
    df18_chem.to_excel(writer, sheet_name=coreNames_chem[1], index=False)
    df21_chem.to_excel(writer, sheet_name=coreNames_chem[2], index=False)
    df29_chem.to_excel(writer, sheet_name=coreNames_chem[3], index=False)



coreNames_ECM = ['B16', 'B18', 'B21']

f_ECM_save = path_save + 'DepthECM.xlsx'
writer = pd.ExcelWriter(f_ECM_save, engine='xlsxwriter')
writer.save()

df16 = pd.DataFrame({'depth': df_B16_ECM[df_B16_ECM.keys()[0]], 'ECM': df_B16_ECM[df_B16_ECM.keys()[-1]]})
writer = pd.ExcelWriter(f_ECM_save, engine='openpyxl')
df16.to_excel(writer, sheet_name=coreNames_ECM[0], index=False)
writer.save()

df18 = pd.DataFrame({'depth': df_B18_ECM[df_B18_ECM.keys()[1]], 'ECM': df_B18_ECM[df_B18_ECM.keys()[0]]})
df21 = pd.DataFrame({'depth': df_B21_ECM[df_B21_ECM.keys()[5]], 'ECM': df_B21_ECM[df_B21_ECM.keys()[-1]]})

with pd.ExcelWriter(f_ECM_save, engine='openpyxl', mode='a') as writer:
    df18.to_excel(writer, sheet_name=coreNames_ECM[1], index=False)
    df21.to_excel(writer, sheet_name=coreNames_ECM[2], index=False)

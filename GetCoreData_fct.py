import numpy as np
import os
import pandas as pd
from pandas import ExcelWriter



def GetCoreData(site_in):
    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
    coreNames = CoresSpecs['CoreName']
    site = site_in

    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]

    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']
    Accum0 = CoreSpecs['Accum0']
    T0 = CoreSpecs['T0']
    dens0 = CoreSpecs['dens0']
    z0 = CoreSpecs['z0']


    try:
        site_d18O =  pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+ site + '_det.txt', ',')
    except:
        print('No d18O file found, setting empty df instead')
        site_d18O = pd.DataFrame(index=np.arange(0,1000),columns=['depth','d18O'])
        site_d18O = site_d18O.fillna(0)

    try:
        site_ECM = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetECM/'+ site + '_ECM.txt', ',')
    except:
        print('No ECM file found, setting empty df instead')
        site_ECM = pd.DataFrame(index=np.arange(0,1000),columns=['depth','ECM'])
        site_ECM = site_ECM.fillna(0)

    try:
            site_Dens = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/'+ site + 'DepthDens_w_Models.txt', '\t')
    except:
        print('No density file found, setting empty df instead')
        site_Dens = pd.DataFrame(index=np.arange(0,1000),columns=['depth', 'dens', 'rhoMeas', 'HLmodel', 'HLmodelOpti'])
        site_Dens = site_Dens.fillna(0)

    try:
        site_Diff = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/'+ site + '_DepthDiff.txt', '\t')
    except:
        print('No diffusion file found, setting empty df instead')
        site_Diff = pd.DataFrame(index=np.arange(0,1000),columns=['Depth','Density', 'sigma_D', 'sigma_o18', 'sigma_o17'])
        site_Diff = site_Diff.fillna(0)

        # Define d18O details data. Define btw. Laki and Tamb.
    depth = site_d18O['depth'][1:]
    d18O = site_d18O['d18O'][1:]
    depth_LT = site_d18O['depth'][(site_d18O['depth'] >= dTamb) & (site_d18O['depth'] <= dLaki)]
    d18O_LT = site_d18O['d18O'][(site_d18O['depth'] >= dTamb) & (site_d18O['depth'] <= dLaki)]

    data_d18O = pd.DataFrame(data = {'depth': depth, 'd18O': d18O})
    data_d18O_LT = pd.DataFrame(data = {'depth': depth_LT, 'd18O': d18O_LT})

        # Define ECM data. Define btw. Laki and Tamb.
    depthECM = site_ECM['depth']
    ECM = site_ECM['ECM']
    depthECM_LT = depthECM[(depthECM >= dTamb) & (depthECM <= dLaki)]
    ECM_LT = ECM[(depthECM >= dTamb) & (depthECM <= dLaki)]

    data_ECM = pd.DataFrame(data = {'depth': depthECM, 'ECM': ECM})
    data_ECM_LT = pd.DataFrame(data = {'depth': depthECM_LT, 'ECM': ECM_LT})

        # Define density measurements: raw, model and fudged model. Define btw. Laki and Tamb.
    if len(site_Dens.columns) == 2:
        depthRho = site_Dens['depth']
        HLmodel = site_Dens['HLmodel']
        depthRho_LT = depthRho[(depthRho >= dTamb) & (depthRho <= dLaki)]
        HLmodel_LT = HLmodel[(depthRho >= dTamb) & (depthRho <= dLaki)]

        data_dens = pd.DataFrame(data = {'depth': depthRho, 'HLmodel': HLmodel})
        data_dens_LT = pd.DataFrame(data = {'depth': depthRho_LT, 'HLmodel': HLmodel_LT})
    elif len(site_Dens.columns) == 4:
        depthRho = site_Dens['depth']
        rhoMeas = site_Dens['rhoMeas']
        HLmodel = site_Dens['HLmodel']
        HLmodelOpti = site_Dens['HLmodelOpti']
        depthRho_LT = depthRho[(depthRho >= dTamb) & (depthRho <= dLaki)]
        rhoMeas_LT = rhoMeas[(depthRho >= dTamb) & (depthRho <= dLaki)]
        HLmodel_LT = HLmodel[(depthRho >= dTamb) & (depthRho <= dLaki)]
        HLmodelOpti_LT =HLmodelOpti[(depthRho >= dTamb) & (depthRho <= dLaki)]

        data_dens = pd.DataFrame(data = {'depth': depthRho, 'rhoMeas': rhoMeas, 'HLmodel': HLmodel, 'HLmodelOpti': HLmodelOpti})
        data_dens_LT = pd.DataFrame(data = {'depth': depthRho_LT, 'rhoMeas': rhoMeas_LT, 'HLmodel': HLmodel_LT,
                                         'HLmodelOpti': HLmodelOpti_LT})

        # Define diffusion length measurements: o18, o17 and D. Define btw. Laki and Tamb.
    depthDiff = site_Diff['Depth']
    sigma_D = site_Diff['sigma_D']
    sigma_o18 = site_Diff['sigma_o18']
    sigma_o17 = site_Diff['sigma_o17']
    depthDiff_LT = depthDiff[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    sigma_o18_LT = sigma_o18[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    sigma_o17_LT = sigma_o17[(depthDiff >= dTamb) & (depthDiff <= dLaki)]
    sigma_D_LT = sigma_D[(depthDiff >= dTamb) & (depthDiff <= dLaki)]

    data_diff = pd.DataFrame(data = {'Depth': depthDiff, 'sigma_D': sigma_D, 'sigma_o18': sigma_o18,
                                     'sigma_o17': sigma_o17})
    data_diff_LT = pd.DataFrame(data = {'Depth': depthDiff_LT, 'sigma_D': sigma_D_LT,
                                     'sigma_o18': sigma_o18_LT, 'sigma_o17': sigma_o17_LT})

#    df = pd.DataFrame(data={'d18O': data_d18O, 'd18O_LT': data_d18O_LT, 'ECM': data_ECM, 'ECM_LT': data_ECM_LT,
#                            'dens': data_dens, 'dens_LT': data_dens_LT, 'diff': data_diff, 'diff_LT': data_diff_LT})
    return  data_d18O, data_d18O_LT, data_ECM, data_ECM_LT, data_dens, data_dens_LT, data_diff, data_diff_LT

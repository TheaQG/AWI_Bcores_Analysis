import numpy as np
import os
import pandas as pd
from pandas import ExcelWriter
from HL_AnalyticThea_class import HL_Thea
from DiffusionProfiles_calculations import DiffusionLength

'''
        *********** TODO************

        - (26/01/21) Fix function to return corespecs of specific core as well as all data.
        - (12/02/21) Add function for density profile creation
        - (12/02/21) Add function for diffusion profile creation
'''

def GetCoreData(site_in, type_in='Alphabet'):
    '''
        Function to locate data files of given core and construct dataframes consisting
        of all data of interest: isotopes, ECM, denisty and diffusion length.
        If no file with data exists, an empty dataframe is produced instead.
        Consideres both files with measured density data and files with only modelled
        density data.

        Arguments:
        ----------
            site_in:            [str] Name of core drilling site.

        returns:
        --------
            data_d18O:          [pd.DataFrame]  Measured isotope data of full core, 'depth', 'd18O'
            data_d18O_LT:       [pd.DataFrame]  Measured isotope data of LT section, 'depth', 'd18O'
            data_ECM:           [pd.DataFrame]  Measured ECM data of full core, 'depth', 'ECM'
            data_ECM_LT:        [pd.DataFrame]  Measured ECM data of LT section, 'depth', 'ECM'
            data_dens:          [pd.DataFrame]  Measured density data of full core, 'depth',
                                                ('rhoMeas'), 'HLmodel', ('HLmodelOpti').
            data_dens_LT:       [pd.DataFrame]  Measured density data of LT section, 'depth',
                                                ('rhoMeas'), 'HLmodel', ('HLmodelOpti').
            data_diff:          [pd.DataFrame]  Measured diffusion data of full core, 'Depth', 'sigma_D',
                                                'sigma_o18', 'sigma_o17'
            data_diff_LT:       [pd.DataFrame]  Measured diffusion data of LT section, 'Depth', 'sigma_D',
                                                'sigma_o18', 'sigma_o17'

    '''

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

        # Load isotope data. If no file to be found, create empty data frame.
    try:
        if type_in == 'Alphabet':
            site_d18O =  pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+ site + '_det.txt', ',')
        elif type_in == 'AWI_Bcores':
            site_d18O =  pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/B_cores_AWI/Bcoresd18O/Depth_d18O__'+ site + '.txt', '\t')
    except:
        print('No d18O file found, setting empty df instead')
        site_d18O = pd.DataFrame(index=np.arange(0,1000),columns=['depth','d18O'])
        site_d18O = site_d18O.fillna(0)


        # Load ECM data. If no file to be found, create empty data frame.
    try:
        if type_in == 'Alphabet':
            site_ECM = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetECM/'+ site + '_ECM.txt', ',')
        elif type_in == 'AWI_Bcores':
            site_ECM =  pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/B_cores_AWI/BcoresECMDEP/'+ site + '_ECM.txt', '\t')
    except:
        print('No ECM file found, setting empty df instead')
        site_ECM = pd.DataFrame(index=np.arange(0,1000),columns=['depth','ECM'])
        site_ECM = site_ECM.fillna(0)


        # Load density data. If no file to be found, create empty data frame.
    try:
        if type_in == 'Alphabet':
            site_Dens = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/'+ site + 'DepthDens_w_Models.txt', '\t')
        elif type_in == 'AWI_Bcores':
            site_Dens =  pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/B_cores_AWI/BcoresDens/'+ site + 'DepthDens_w_Models.txt', '\t')

    except:
        print('No density file found, setting empty df instead')
        site_Dens = pd.DataFrame(index=np.arange(0,1000),columns=['depth', 'dens', 'rhoMeas', 'HLmodel', 'HLmodelOpti'])
        site_Dens = site_Dens.fillna(0)


        # Load diffusion length data. If no file to be found, create empty data frame.
    try:
        if type_in == 'Alphabet':
            site_Diff = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/'+ site + '_DepthDiff.txt', '\t')
        elif type_in == 'AWI_Bcores':
            site_Diff =  pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/B_cores_AWI/BcoresDiff/'+ site + '_DepthDiff.txt', '\t')
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

        # Define density measurements: raw, model and fudged model. Check the length
        # of density - use only modelled, if only two columns. Define btw. Laki and Tamb.
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

    return  data_d18O, data_d18O_LT, data_ECM, data_ECM_LT, data_dens, data_dens_LT, data_diff, data_diff_LT











def GetDensProfile(site_in, path_densMeas, delim_densMeas, path_isoMeas, delim_isoMeas, path_outFile, delim_outFile, area_in='Alphabet', zMeas_str = 'depth', densMeas_str = 'density'):
    import HL_AnalyticThea_class

    from HL_AnalyticThea_class import HL_Thea


    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
    coreNames = CoresSpecs['CoreName']
    site = site_in

    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]

    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']
    bdot0 = CoreSpecs['Accum0']
    Temp0 = CoreSpecs['T0']+273.15
    dens0 = CoreSpecs['dens0']
    z0 = CoreSpecs['z0']


    try:
        densMeas = pd.read_csv(path_densMeas,delim_densMeas)
        densMeas_in = True
    except:
        print('No density measurements. Creating purely analytical profile.')
        densMeas_in = False
        dens0 = 350

    if densMeas_in:
        z_vec = densMeas[zMeas_str]
        rho_vec = densMeas[densMeas_str].astype('float64')

        hl_inst = HL_Thea(z_meas = z_vec, rho_meas = rho_vec,\
                             Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = False)
        hl_model = hl_inst.model(z_vec)

        hl_instOpti = HL_Thea(z_meas = z_vec, rho_meas = rho_vec,\
                             Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = True)

        hl_modelOpti = hl_instOpti.model(z_vec)

        data = pd.DataFrame({'depth':z_vec,'rhoMeas':rho_vec,'HLmodel': hl_model['rhoHL'],'HLmodelOpti': hl_modelOpti['rhoHL']})
        data.to_csv(path_outFile,index=None,sep=delim_outFile)

    else:
        z_vec = np.asarray(pd.read_csv(path_isoMeas,delim_isoMeas)[zMeas_str])

        hl_inst = HL_Thea(Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = False)
        hl_model = hl_inst.model(z_vec)

        data = pd.DataFrame({'depth':z_vec,'HLmodel': hl_model['rhoHL']})
        data.to_csv(path_outFile,index=None,sep=delim_outFile)

    return

'''
    Density profile creation:
'''
# from GetCoreData_fct import GetDensProfile
# print('Creating density profile')
# core = 'SiteH'
#
# pathDens = '/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/' + core + '_DepthDens.txt'
# pathOut = '/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/' + core + 'DepthDens_w_Models.txt'
# pathIso = '/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+core+'_det.txt'
#
# delimIso = ','
# delimOut = '\t'
# delimDens = '\t'
#
# GetDensProfile(site_in = core, path_densMeas=pathDens,delim_densMeas=delimDens, path_isoMeas=pathIso, delim_isoMeas=delimIso, path_outFile=pathOut, delim_outFile=delimOut, area_in='Alphabet', zMeas_str = 'depth', densMeas_str = 'density')



def GetDiffProfile(site_in, path_outFile, delim_outFile, path_densMeas,delim_densMeas, path_isoMeas, delim_isoMeas, densMeas_str='density', zMeas_str='depth', dz_in=0.55):
    import HL_AnalyticThea_class
    from HL_AnalyticThea_class import HL_Thea
    from DiffusionProfiles_calculations import DiffusionLength

    CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
    coreNames = CoresSpecs['CoreName']
    site = site_in

    core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
    CoreSpecs = CoresSpecs.iloc[core_idx]

    dTamb = CoreSpecs['dTamb']
    dLaki = CoreSpecs['dLaki']
    bdot0 = CoreSpecs['Accum0']
    Temp0 = CoreSpecs['T0']+273.15
    dens0 = CoreSpecs['dens0']
    z0 = CoreSpecs['z0']


    try:
        densMeas = pd.read_csv(path_densMeas,delim_densMeas)
        densMeas_in = True
    except:
        print('No density measurements. Creating purely analytical profile.')
        densMeas_in = False
        dens0 = 350

    if densMeas_in:
        z_vec = densMeas[zMeas_str]
        rho_vec = densMeas[densMeas_str].astype('float64')

        hl_inst = HL_Thea(z_meas = z_vec, rho_meas = rho_vec,\
                             Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = False)
        hl_model = hl_inst.model(z_vec)

        hl_instOpti = HL_Thea(z_meas = z_vec, rho_meas = rho_vec,\
                             Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = True)

        hl_modelOpti = hl_instOpti.model(z_vec)
        f0_fin = hl_modelOpti['f0_fin']; f1_fin = hl_modelOpti['f1_fin']

        sigma_arr = diffProfileCalc(P = 0.75, temp = Temp0, accum = bdot0, rho_surf = dens0, f0 = f0_fin, f1 = f1_fin,\
                           dz = dz_in, z_final = max(z_vec), fileout = path_outFile)

#        sigma_inst = SigmaToolbox()
#        sigma_arr = sigma_inst.experiment2(P = 0.75, temp = Temp0, accum = bdot0, rho_o = dens0, \
#                        fo = f0_fin, f1 = f1_fin, dz = dz_in, z_final = max(z_vec), fileout =path_outFile)

    else:
        z_vec = np.asarray(pd.read_csv(path_isoMeas,delim_isoMeas)[zMeas_str])


        hl_inst = HL_Thea(Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = False)
        hl_model = hl_inst.model(z_vec)
        f0_fin = hl_model['f0_fin']; f1_fin = hl_model['f1_fin']

        sigma_arr = diffProfileCalc(P = 0.75, temp = Temp0, accum = bdot0, rho_surf = dens0, f0 = f0_fin, f1 = f1_fin,\
                           dz = dz_in, z_final = max(z_vec), fileout = path_outFile)
#        sigma_inst = SigmaToolbox()
#        sigma_arr = sigma_inst.experiment2(P = 0.75, temp = Temp0, accum = bdot0, rho_o = dens0, \
#                        fo = f0_fin, f1 = f1_fin, dz = dz_in, z_final = max(z_vec), fileout =path_outFile)


    return

def diffProfileCalc(P = 0.75, temp = 244.15, accum = 0.025, rho_surf = 350.0, f0 = 1, f1 = 1,\
                   dz = 1, z_final = 100, fileout = False):
    sigma_o17_num = []
    sigma_o18_num = []
    sigma_D_num = []

    herron_model = HL_Thea(Temp_0 = temp, Acc_0 = accum, rho_0 = rho_surf, f0_init = f0, f1_init = f1).model(np.arange(0, z_final, dz))
    rhos = 1000*herron_model['rhoHL']
    depths = herron_model['z']

    sigmaInst = DiffusionLength(P = P, rho_surf = rho_surf, f0 = f0, f1 = f1)

#    for rho in rhos:
#        sigma_all = sigmaInst.semi_analytical_HL(rho = rho, T = temp, accum = accum)
#        sigma_o17_num = np.append(sigma_o17_num, sigma_all[2])
#        sigma_o18_num = np.append(sigma_o18_num, sigma_all[1])
#        sigma_D_num = np.append(sigma_D_num, sigma_all[0])

    sigma_D_analyt, sigma_o18_analyt, sigma_o17_analyt = sigmaInst.analytical_HL(rhos, T = temp, accum = accum)

    if not fileout == False:
        f = open(fileout, "w")
        f.write("Depth\tDensity\tsigma_D\tsigma_o18\tsigma_o17\n")
        data_out = np.transpose(np.vstack((depths, rhos, sigma_D_analyt, sigma_o18_analyt, sigma_o17_analyt)))
        np.savetxt(f, data_out, delimiter = "\t", fmt = ("%0.2f", "%0.3f", "%0.6f", "%0.6f", "%0.6f"))
        f.close()

    return




'''
    Diffusion profile creation:
'''

#from GetCoreData_fct import GetDiffProfile

#core = 'SiteB'
# print('Creating diffusion length profile')
# PathDens = '/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDens/' + core + '_DepthDens.txt'
# PathOut = '/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/AlphabetDiff/' + core + '_DepthDiff.txt'
# PathIso = '/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/Alphabet_cores/Alphabetd18O/'+core+'_det.txt'
#
# DelimIso = ','
# DelimOut = '\t'
# DelimDens = '\t'
#
# GetDiffProfile(site_in=core, path_outFile=PathOut, delim_outFile=DelimOut, path_densMeas=PathDens, delim_densMeas=DelimDens, path_isoMeas=PathIso, delim_isoMeas=DelimIso)




'''
    OLD MODULE! DO NOT USE. LINKS TO VASILEIOS GKINIS sigma.py
'''
# def GetDiffProfile(site_in, path_outFile, delim_outFile, path_densMeas,delim_densMeas, path_isoMeas, delim_isoMeas, densMeas_str='density', zMeas_str='depth', dz_in=0.55):
#     import HL_AnalyticThea_class
#     from HL_AnalyticThea_class import HL_Thea
#
#     import sigma
#     from sigma import SigmaToolbox
#
#     CoresSpecs = pd.read_csv('/home/thea/Documents/KUFysik/MesterTesen/Data/CoreSpecs.txt', ',')
#     coreNames = CoresSpecs['CoreName']
#     site = site_in
#
#     core_idx = coreNames[CoresSpecs['CoreName'] == site].index[0]
#     CoreSpecs = CoresSpecs.iloc[core_idx]
#
#     dTamb = CoreSpecs['dTamb']
#     dLaki = CoreSpecs['dLaki']
#     bdot0 = CoreSpecs['Accum0']
#     Temp0 = CoreSpecs['T0']+273.15
#     dens0 = CoreSpecs['dens0']
#     z0 = CoreSpecs['z0']
#
#
#     try:
#         densMeas = pd.read_csv(path_densMeas,delim_densMeas)
#         densMeas_in = True
#     except:
#         print('No density measurements. Creating purely analytical profile.')
#         densMeas_in = False
#         dens0 = 350
#
#     if densMeas_in:
#         z_vec = densMeas[zMeas_str]
#         rho_vec = densMeas[densMeas_str].astype('float64')
#
#         hl_inst = HL_Thea(z_meas = z_vec, rho_meas = rho_vec,\
#                              Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = False)
#         hl_model = hl_inst.model(z_vec)
#
#         hl_instOpti = HL_Thea(z_meas = z_vec, rho_meas = rho_vec,\
#                              Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = True)
#
#         hl_modelOpti = hl_instOpti.model(z_vec)
#         f0_fin = hl_modelOpti['f0_fin']; f1_fin = hl_modelOpti['f1_fin']
#
#
#         sigma_inst = SigmaToolbox()
#         sigma_arr = sigma_inst.experiment2(P = 0.75, temp = Temp0, accum = bdot0, rho_o = dens0, \
#                         fo = f0_fin, f1 = f1_fin, dz = dz_in, z_final = max(z_vec), fileout =path_outFile)
#
#     else:
#         z_vec = np.asarray(pd.read_csv(path_isoMeas,delim_isoMeas)[zMeas_str])
#
#
#         hl_inst = HL_Thea(Acc_0 = bdot0, Temp_0 = Temp0, rho_0 = dens0, opti = False)
#         hl_model = hl_inst.model(z_vec)
#         f0_fin = hl_model['f0_fin']; f1_fin = hl_model['f1_fin']
#
#         sigma_inst = SigmaToolbox()
#         sigma_arr = sigma_inst.experiment2(P = 0.75, temp = Temp0, accum = bdot0, rho_o = dens0, \
#                         fo = f0_fin, f1 = f1_fin, dz = dz_in, z_final = max(z_vec), fileout =path_outFile)
#
#
#     return

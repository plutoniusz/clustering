import os
import nitorch.tools.qmri.io as qio
from nitorch.tools.qmri.param import ParameterMap
from nitorch.tools.qmri.relax import GREEQOptions
from nitorch.io import savef
from nitorch.tools.qmri.relax._mpm._greeq import greeq
from pathlib import Path
from datetime import datetime

# which echoes to estimate
model = 'chi'
noise_model = 'chi'
contrasts = ["pd", "t1", "mt"] # order of contrasts
penalty = dict(pd=10, r1=10, mt=5, r2s=5)
# dataset choice
datasets = ['cl'] #'mpm'/ 'cl'/ 'mc'
cl = bool(datasets==["cl"])
if cl:
    datasets = ['112111', '128221', '130519', '170192', '176117', '208010', '211787', '214685', '232237', '308597', '324038', '330406', '346878']
    datasets = datasets[11:12]
    print(datasets)

datime = str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

# paths to the datasets
for dataset in datasets:
    # for naming the result folders by the estimation date
    
    print(f"date: {datime}, dataset: {dataset}")

    # paths to dataset
    cwd = os.getcwd()
    if dataset == 'mc':
        pth_mris = [os.path.join(cwd, '4John_Klara/mtw_mfc_3dflash_v1k_180deg_RR_0038'),
                    os.path.join(cwd, '4John_Klara/pdw_mfc_3dflash_v1k_RR_0036'),
                    os.path.join(cwd, '4John_Klara/t1w_mfc_3dflash_v1k_RR_0034')]
    elif dataset == 'mpm':
        pth_mris = [os.path.join(cwd, 'MPM/pdw_mfc_3dflash_v1i_R4_0009'),
                    os.path.join(cwd, 'MPM/t1w_mfc_3dflash_v1i_R4_0015'),
                    os.path.join(cwd, 'MPM/mtw_mfc_3dflash_v1i_R4_0012')]
        for filename in os.listdir(str(os.path.join(cwd, "MPM/pdw_mfc_3dflash_v1i_R4_0009/Results/Supplementary"))):
            if filename.endswith("B1map.nii"):
                b1p_map = str(os.path.join(cwd, "MPM/pdw_mfc_3dflash_v1i_R4_0009/Results/Supplementary", filename))
            elif filename.endswith("B1ref.nii"):
                b1p_ref = str(os.path.join(cwd, "MPM/pdw_mfc_3dflash_v1i_R4_0009/Results/Supplementary", filename))
        b1m_map = [0]*len(contrasts)
        for contrast_number,contrast in enumerate(contrasts):
            b1m_map[contrast_number] = str(os.path.join(cwd, "MPM/pdw_mfc_3dflash_v1i_R4_0009/Results/Supplementary/sensMap_HC_over_BC_division_" + contrast.upper() + ".nii"))
    elif cl:
        #gen_path = "/data/underworld/01_DATA/00_ANALYSIS/02_KBAS/03_data/source/mri/"
        gen_path = "/data/underworld/kbas/03_data/source_dif/mri/"
        path_dataset = os.path.join(gen_path, dataset)
        path_scans = os.listdir(path_dataset)
        path_scans = os.path.join(path_dataset, path_scans[0])
        pth_mris = [os.path.join(path_scans, "13"), #pdw
                    os.path.join(path_scans, "10"), #t1w
                    os.path.join(path_scans, "16")] #mtw
        #path_dataset = os.path.join("/data/underworld/01_DATA/00_ANALYSIS/02_KBAS/03_data/raw/", dataset)
        path_dataset = os.path.join("/data/underworld/kbas/03_data/raw_dif", dataset)
        path_pom = os.listdir(path_dataset)
        path_dataset = os.path.join(path_dataset, path_pom[0])
        path_suppl = os.path.join(path_dataset, "anat/Results/Supplementary")
        for filename in os.listdir(str(path_suppl)):
            if filename.endswith("B1map.nii"):
                b1p_map = str(os.path.join(path_suppl, filename))
            elif filename.endswith("B1ref.nii"):
                b1p_ref = str(os.path.join(path_suppl, filename))
        b1m_map = [0]*len(contrasts)
        for contrast_number, contrast in enumerate(contrasts):
            b1m_map[contrast_number] = os.path.join(path_dataset, "anat/Results/Supplementary/sensMap_HC_over_BC_division_" + contrast.upper() + ".nii")

    # paths to mris
    mris = [0]*len(contrasts)
    for contrast_number, contrast in enumerate(contrasts):
        mris_pom = []
        for filename in os.listdir(str(pth_mris[contrast_number])):
            if filename.endswith(".nii") & filename.startswith("sMP"):
                mris_pom.append(str(os.path.join(pth_mris[contrast_number], filename)))
            else:
                continue
        mris[contrast_number] = mris_pom

    print("precomputind field maps")
    if dataset == 'mc':
        transmit = None
        receive = [None]*3
    else:
        transmit = qio.PrecomputedFieldMap(b1p_map, magnitude=b1p_ref)
        receive = [qio.PrecomputedFieldMap(b1m_map[0],  unit='a.u'), #pd
                    qio.PrecomputedFieldMap(b1m_map[1], unit='a.u'), #t1
                    qio.PrecomputedFieldMap(b1m_map[2], unit='a.u')] #mt

    # folder for saving the results of this estimation
    save_folder = '/' + datime + '/' + dataset
    folder_make = Path(cwd + save_folder).mkdir(parents=True, exist_ok=True)

    # logfile
    log_txt = f"date: {datime}\n dataset: {dataset}\n penalty: {penalty}\n models: {model}\n save_folder {save_folder}\n\n"
    with open(cwd + save_folder + '/log_est.txt', 'a') as f:
        f.write(log_txt)


    # preparing path lists of observations

    mris_le = [0]*len(contrasts)
    for contrast_number, contrast in enumerate(contrasts):
        mris_le[contrast_number] = qio.GradientEchoMulti.from_fnames(mris[contrast_number], mt=bool(contrast=="mt"))

    # options for parameter estimation
    mris_le[2].mt = True
    opt = GREEQOptions()
    opt.recon.space = 1
    opt.backend.device = 'cuda'
    opt.preproc.register = False #is this correct
    opt.optim.nb_levels = 1
    opt.verbose = 1
    opt.likelihood = model
    opt.noisemodel = model
    opt.penalty.factor = penalty
    print("start greeq")
    pd, r1, r2s, mt = greeq(mris_le, transmit = transmit, receive = receive, opt=opt)

    # saving parameter maps
    pm_path = Path(cwd + save_folder + '/parameter_map').mkdir(parents=True, exist_ok=True)
    savef(pd.volume, os.path.join(cwd + save_folder + '/parameter_map' +'/pd'+'.nii'), affine= pd.affine)
    savef(r1.volume, os.path.join(cwd + save_folder + '/parameter_map' +'/r1'+'.nii'), affine= r1.affine)
    savef(r2s.volume, os.path.join(cwd + save_folder + '/parameter_map' +'/r2s'+'.nii'), affine= r2s.affine)
    savef(mt.volume, os.path.join(cwd + save_folder + '/parameter_map' +'/mt'+'.nii'), affine= mt.affine)

    # svaing noise estimates
    noise_txt = f"model: {model}, noise: {pd.noise}, dof: {pd.dof}\n"
    with open(cwd + save_folder + '/noise_' + dataset + '_' + datime + '.txt', 'a') as f:
        f.write(noise_txt)

        
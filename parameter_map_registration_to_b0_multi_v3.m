%-----------------------------------------------------------------------
% Job saved on 22-May-2023 15:11:43 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

subjects = [112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 324038 330406 346878];

for i = 1:numel(subjects)

    subject = num2str(subjects(i));

    spm('defaults','fmri');
    spm_jobman('initcfg');

    % b0 images
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/data/underworld/kbas/03_data/processed_dif'};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^y_1_000.._sub-' subject '_ses-.*_desc-average-sbref_b0_mb\.nii$'];
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {['/data/underworld/kbas/03_data/source_dif/mri/' subject]};
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sMP.*-0010-.*01\.nii';
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';

    matlabbatch{3}.spm.util.defs.comp{1}.def(1) =cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{3}.spm.util.defs.comp{2}.id.space(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sMP.*-0010-.*01\.nii)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

    matlabbatch{3}.spm.util.defs.out{1}.savedef.ofname = ['r_deformation_' subject '_qmri'];
    matlabbatch{3}.spm.util.defs.out{1}.savedef.savedir.saveusr = {['/data/underworld/kbas/03_data/processed_dif/' subject]};

    % resliced deformations
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/data/underworld/home/kbas/03_data/processed_dif/' subject};
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['y_r_deformation_' subject '_qmri.nii$'];
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';

    % read mt map from a directory
    matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^mt.nii$';
    matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    % apply deformation and reslice to diffusion resolution
    matlabbatch{6}.spm.util.defs.comp{1}.def(1) = cfg_dep('File Selector (Batch Mode): Selected Files ([y_r_deformation subject .nii])', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{6}.spm.util.defs.comp{2}.id.space(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{6}.spm.util.defs.out{1}.push.fnames(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^(mt)\.nii$)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{6}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{6}.spm.util.defs.out{1}.push.savedir.saveusr(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    % template
    matlabbatch{6}.spm.util.defs.out{1}.push.fov.file = {'/data/underworld/kbas/03_data/processed_dif/mu_mb.nii'};
    matlabbatch{6}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{6}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{6}.spm.util.defs.out{1}.push.prefix = 'average_';

    % read pd map from a directory
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^pd.nii$';
    matlabbatch{7}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    % apply deformation and reslice to diffusion resolution
    matlabbatch{8}.spm.util.defs.comp{1}.def(1) = cfg_dep('File Selector (Batch Mode): Selected Files ([y_r_deformation subject .nii])', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.util.defs.comp{2}.id.space(1) =  cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.util.defs.out{1}.push.fnames(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^(pd)\.nii$)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{8}.spm.util.defs.out{1}.push.savedir.saveusr(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    % template
    matlabbatch{8}.spm.util.defs.out{1}.push.fov.file = {'/data/underworld/kbas/03_data/processed_dif/mu_mb.nii'};
    matlabbatch{8}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{8}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{8}.spm.util.defs.out{1}.push.prefix = 'average_';

    % read r1 map from a directory
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^r1.nii$';
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    % apply deformation and reslice to diffusion resolution
    matlabbatch{10}.spm.util.defs.comp{1}.def(1) = cfg_dep('File Selector (Batch Mode): Selected Files ([y_r_deformation subject .nii])', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{10}.spm.util.defs.comp{2}.id.space(1) =  cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{10}.spm.util.defs.out{1}.push.fnames(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^t1.nii$)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{10}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{10}.spm.util.defs.out{1}.push.savedir.saveusr(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    % template
    matlabbatch{10}.spm.util.defs.out{1}.push.fov.file = {'/data/underworld/kbas/03_data/processed_dif/mu_mb.nii'};
    matlabbatch{10}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{10}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{10}.spm.util.defs.out{1}.push.prefix = 'average_';

    % read r2 map from a directory
    matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^r2s.nii$';
    matlabbatch{11}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    % apply deformation and reslice to diffusion resolution
    matlabbatch{12}.spm.util.defs.comp{1}.def(1) = cfg_dep('File Selector (Batch Mode): Selected Files ([y_r_deformation subject .nii])', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{12}.spm.util.defs.comp{2}.id.space(1) =  cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{12}.spm.util.defs.out{1}.push.fnames(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^r2s.nii$)', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{12}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{12}.spm.util.defs.out{1}.push.savedir.saveusr(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
    % template
    matlabbatch{12}.spm.util.defs.out{1}.push.fov.file = {'/data/underworld/kbas/03_data/processed_dif/mu_mb.nii'};
    matlabbatch{12}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{12}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{12}.spm.util.defs.out{1}.push.prefix = 'average_';

    spm_jobman('run',matlabbatch);
end

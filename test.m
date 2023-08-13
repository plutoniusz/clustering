subjects = [112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 324038 330406 346878];
for i = 1:numel(subjects)

    subject = num2str(subjects(i));

    spm('defaults','fmri');
    spm_jobman('initcfg');

    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {['/data/underworld/home/kbas/03_data/derivatives_dif/' subject]};
    %matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^sub-' subject '_ses-.*_space-orig_desc-dwi-skullstripped_b0\.nii$'];
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^sub-' subject '_ses-.*_desc-average-sbref_b0\.nii$'];
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
%     matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/data/underworld/home/kbas/03_data/source_b0'};
%     matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^y_1_000.._sub-' subject '_ses-.*_space-orig_desc-dwi-skullstripped_b0_mb\.nii$']; y_template112111.nii
    
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/data/underworld/home/kbas/03_data/processed_dif'};
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^y_1_000.._sub-' subject '_ses-.*_desc-average-sbref_b0_mb\.nii$']; %['^y_1_000.._sub-' subject '_ses-.*_space-orig_desc-dwi-skullstripped_b0_mb2\.nii$']; %['^y_template' subject '.nii$'];y_1_00004_sub-170192_ses-20190111_desc-average-sbref_b0_mb
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
    
    matlabbatch{3}.spm.util.defs.comp{1}.def(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    %matlabbatch{3}.spm.util.defs.out{1}.push.fnames(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub- subject _ses-.*_space-orig_desc-dwi-skullstripped_b0\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{3}.spm.util.defs.out{1}.push.fnames(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{3}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{3}.spm.util.defs.out{1}.push.savedir.saveusr = {'/data/underworld/home/kbas/03_data/processed_dif'};
    %matlabbatch{3}.spm.util.defs.out{1}.push.fov.file = {'/data/underworld/home/kbas/spm12/toolbox/mb/mu_X.nii'};
    matlabbatch{3}.spm.util.defs.out{1}.push.fov.file = {'/data/underworld/kbas/03_data/processed_dif/mu_mb.nii'};
    %matlabbatch{3}.spm.util.defs.out{1}.push.fov.bbvox.bb = [NaN NaN NaN
    %                                                         NaN NaN NaN];
    %matlabbatch{3}.spm.util.defs.out{1}.push.fov.bbvox.vox = [NaN NaN NaN];
    matlabbatch{3}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{3}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{3}.spm.util.defs.out{1}.push.prefix = 'warped_';


    spm_jobman('run',matlabbatch);
end

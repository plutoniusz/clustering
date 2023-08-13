subjects = [112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 324038 330406 346878];
for subject=subjects

    subject = num2str(subject);

    spm('defaults','fmri');
    spm_jobman('initcfg');


    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {['/data/underworld/kbas/03_data/derivatives_dif/' subject]};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^sub-' subject '_ses-\d+_desc-average-sbref_b0\.nii$'];
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';

    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/data/underworld/kbas/03_data/processed_dif'};
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = ['^y_1_000.._sub-' subject '_ses-.*_desc-average-sbref_b0_mb\.nii$'];
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

    matlabbatch{3}.spm.util.defs.comp{1}.def(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^y_1_000.._sub- subject _ses-.*_desc-average-sbref_b0_mb\.nii$)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

    matlabbatch{3}.spm.util.defs.comp{2}.id.space(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub- subject _ses-\d+_desc-average-sbref_b0\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

    matlabbatch{3}.spm.util.defs.out{1}.savedef.ofname = ['r_deformation_' subject '_dif'];
    matlabbatch{3}.spm.util.defs.out{1}.savedef.savedir.saveusr = {['/data/underworld/kbas/03_data/processed_dif/' subject]};

    spm_jobman('run',matlabbatch);
end
spm('defaults','fmri');
spm_jobman('initcfg');

% read mt map from a directory
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/112111']};
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^average_r1.nii$';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';

matlabbatch{2}.spm.tools.mb.fil.images(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^r1.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

spm_jobman('run',matlabbatch);

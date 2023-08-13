%-----------------------------------------------------------------------
% Job saved on 22-Mar-2023 13:34:21 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% generates deformation files between subjects
% based on the b0 files
% uses existing template mu_X.nii

spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/data/underworld/home/kbas/03_data/derivatives_dif'};
%matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sub-(112111|128221|130519|170192|176117|208010|210022|211787|214685|232237|308597|324038|330406|346878)_ses-\d+_space-orig_desc-dwi-skullstripped_b0\.nii$';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sub-(112111|128221|130519|170192|176117|208010|211787|214685|232237|308597|324038|330406|346878)_ses-\d+_desc-average-sbref_b0\.nii$';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
% matlabbatch{2}.spm.tools.mb.run.mu.exist = {'/home/kbas/Documents/spm12/toolbox/mb/data/mu_X.nii'};

matlabbatch{2}.spm.tools.mb.run.mu.create.K = 5;
matlabbatch{2}.spm.tools.mb.run.mu.create.vx = 1.71;
matlabbatch{2}.spm.tools.mb.run.mu.create.mu_settings = [1e-05 0.5 0];
matlabbatch{2}.spm.tools.mb.run.mu.create.issym = 0;

matlabbatch{2}.spm.tools.mb.run.aff = 'SE(3)';
matlabbatch{2}.spm.tools.mb.run.v_settings = [0.0001 0 0.4 0.1 0.4];
matlabbatch{2}.spm.tools.mb.run.del_settings = Inf;
matlabbatch{2}.spm.tools.mb.run.onam = 'mb';
matlabbatch{2}.spm.tools.mb.run.odir = {'/data/underworld/home/kbas/03_data/processed_dif'};
matlabbatch{2}.spm.tools.mb.run.cat = {{}};
matlabbatch{2}.spm.tools.mb.run.gmm.chan.images(1) = cfg_dep( 'File Selector (Batch Mode): Selected Files (^sub-(112111|128221|130519|170192|176117|208010|211787|214685|232237|308597|324038|330406|346878)_ses-\d+_desc-average-sbref_b0\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.tools.mb.run.gmm.chan.inu.inu_reg = 10000;
matlabbatch{2}.spm.tools.mb.run.gmm.chan.inu.inu_co = 40;
matlabbatch{2}.spm.tools.mb.run.gmm.chan.modality = 1;
matlabbatch{2}.spm.tools.mb.run.gmm.labels.false = [];
matlabbatch{2}.spm.tools.mb.run.gmm.pr.file = {};
%matlabbatch{2}.spm.tools.mb.run.gmm.pr.hyperpriors = [];
matlabbatch{2}.spm.tools.mb.run.gmm.pr.hyperpriors = {
                                                      'b0_priors'
                                                      {
                                                      0.01
                                                      0.01
                                                      }
                                                      };

matlabbatch{2}.spm.tools.mb.run.gmm.tol_gmm = 0.0005;
matlabbatch{2}.spm.tools.mb.run.gmm.nit_gmm_miss = 32;
matlabbatch{2}.spm.tools.mb.run.gmm.nit_gmm = 8;
matlabbatch{2}.spm.tools.mb.run.gmm.nit_appear = 8;
matlabbatch{2}.spm.tools.mb.run.accel = 0.8;
matlabbatch{2}.spm.tools.mb.run.min_dim = 8;
matlabbatch{2}.spm.tools.mb.run.tol = 0.001;
matlabbatch{2}.spm.tools.mb.run.sampdens = 2;
matlabbatch{2}.spm.tools.mb.run.save = true;
matlabbatch{2}.spm.tools.mb.run.nworker = 0;

spm_jobman('run',matlabbatch);

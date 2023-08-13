%-----------------------------------------------------------------------
% Job saved on 19-May-2023 19:30:35 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

subjects = [112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 324038 330406 346878];

for i = 1:numel(subjects)

    subject = num2str(subjects(i));

    spm('defaults','fmri');
    spm_jobman('initcfg');

    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {['/data/underworld/home/kbas/03_data/source_dif/mri/' subject]};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sMP.*-(0010|0013|0016)-.*\.nii$';
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_file_split.name = 'mpm_data';
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_file_split.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sMP.*-(0010|0013|0016)-.*\.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_file_split.index = {
                                                                         1
                                                                         [2 3 4 5 6 7 8]
                                                                         9
                                                                         [10 11 12 13 14 15 16]
                                                                         17
                                                                         [18 19 20 21 22]
                                                                         };
    matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {['/data/underworld/home/kbas/03_data/derivatives_dif/' subject]};
    matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sub-\d+_ses-\d+_desc-average-sbref_b0\.nii$';
    matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {['/data/underworld/home/kbas/03_data/raw_dif/' subject]};
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sensMap.*_(MT|PD|T1).nii$';
    matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';
    matlabbatch{5}.cfg_basicio.file_dir.file_ops.cfg_file_split.name = 'field_map';
    matlabbatch{5}.cfg_basicio.file_dir.file_ops.cfg_file_split.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sensMap.*_(MT|PD|T1).nii$)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{5}.cfg_basicio.file_dir.file_ops.cfg_file_split.index = {
                                                                        1
                                                                        2
                                                                        3
                                                                        };
    matlabbatch{6}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub-\d+_ses-\d+_desc-average-sbref_b0\.nii$)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{6}.spm.spatial.coreg.estimate.source(1) = cfg_dep('File Set Split: mpm_data (1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{1}));
    matlabbatch{6}.spm.spatial.coreg.estimate.other(1) = cfg_dep('File Set Split: mpm_data (2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{2}));
    matlabbatch{6}.spm.spatial.coreg.estimate.other(2) = cfg_dep('File Set Split: field_map (3)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{3}));
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    matlabbatch{7}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub-\d+_ses-\d+_desc-average-sbref_b0\.nii$)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{7}.spm.spatial.coreg.estimate.source(1) = cfg_dep('File Set Split: mpm_data (3)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{3}));
    matlabbatch{7}.spm.spatial.coreg.estimate.other(1) = cfg_dep('File Set Split: mpm_data (4)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{4}));
    matlabbatch{7}.spm.spatial.coreg.estimate.other(2) = cfg_dep('File Set Split: field_map (2)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{2}));
    matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    matlabbatch{8}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub-\d+_ses-\d+_desc-average-sbref_b0\.nii$)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.spatial.coreg.estimate.source(1) = cfg_dep('File Set Split: mpm_data (5)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{5}));
    matlabbatch{8}.spm.spatial.coreg.estimate.other(1) = cfg_dep('File Set Split: mpm_data (6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{6}));
    matlabbatch{8}.spm.spatial.coreg.estimate.other(2) = cfg_dep('File Set Split: field_map (1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{1}));
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    % label the T1w image
    matlabbatch{9}.spm.tools.mb.fil.images(1) = cfg_dep('File Set Split: mpm_data (1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('{}',{1}));
    % reslice the labels using label option
    matlabbatch{10}.spm.spatial.coreg.write.ref(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sub-\d+_ses-\d+_desc-average-sbref_b0\.nii$)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{10}.spm.spatial.coreg.write.source(1) = cfg_dep('Image Labelling: Labelled brains', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','labels', '()',{':'}));
    matlabbatch{10}.spm.spatial.coreg.write.roptions.interp = -1;
    matlabbatch{10}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{10}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{10}.spm.spatial.coreg.write.roptions.prefix = 'r';
    % calculate thalamus mask
    matlabbatch{11}.spm.util.imcalc.input(1) = cfg_dep('Coregister: Reslice: Resliced Images', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
    matlabbatch{11}.spm.util.imcalc.output = 'thalamus_mask_diff';
    matlabbatch{11}.spm.util.imcalc.outdir = {['/data/underworld/home/kbas/03_data/processed_dif/' subject]};
    matlabbatch{11}.spm.util.imcalc.expression = '(i1==59) + (i1==60)';
    matlabbatch{11}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{11}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{11}.spm.util.imcalc.options.mask = 0;
    matlabbatch{11}.spm.util.imcalc.options.interp = 1;
    matlabbatch{11}.spm.util.imcalc.options.dtype = 4;

    % increase the mask area
    matlabbatch{12}.spm.spatial.smooth.data = cfg_dep('Image Calculator: ImCalc Computed Image: thalamus_mask_diff', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));;
    matlabbatch{12}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{12}.spm.spatial.smooth.dtype = 0;
    matlabbatch{12}.spm.spatial.smooth.im = 0;
    matlabbatch{12}.spm.spatial.smooth.prefix = 's';

    
    matlabbatch{13}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{13}.spm.util.imcalc.output = 'thalamus_mask_diff_inc';
    matlabbatch{13}.spm.util.imcalc.outdir = {['/data/underworld/home/kbas/03_data/processed_dif/' subject]};
    matlabbatch{13}.spm.util.imcalc.expression = 'i1>0.05';
    matlabbatch{13}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{13}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{13}.spm.util.imcalc.options.mask = 0;
    matlabbatch{13}.spm.util.imcalc.options.interp = 1;
    matlabbatch{13}.spm.util.imcalc.options.dtype = 4;

    % zip the thalamus mask
    matlabbatch{14}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(1) = cfg_dep('Image Calculator: ImCalc Computed Image: thalamus_mask_diff_inc', substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{14}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = {''};
    matlabbatch{14}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep = true;

    % calculate amygdala mask
    matlabbatch{15}.spm.util.imcalc.input(1) = cfg_dep('Coregister: Reslice: Resliced Images', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
    matlabbatch{15}.spm.util.imcalc.output = 'amygdala_mask_diff';
    matlabbatch{15}.spm.util.imcalc.outdir = {['/data/underworld/home/kbas/03_data/processed_dif/' subject]};
    matlabbatch{15}.spm.util.imcalc.expression = '(i1==31) + (i1==32)';
    matlabbatch{15}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{15}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{15}.spm.util.imcalc.options.mask = 0;
    matlabbatch{15}.spm.util.imcalc.options.interp = 1;
    matlabbatch{15}.spm.util.imcalc.options.dtype = 4;

    % increase the mask area
    matlabbatch{16}.spm.spatial.smooth.data = cfg_dep('Image Calculator: ImCalc Computed Image: amygdala_mask_diff', substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));;
    matlabbatch{16}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{16}.spm.spatial.smooth.dtype = 0;
    matlabbatch{16}.spm.spatial.smooth.im = 0;
    matlabbatch{16}.spm.spatial.smooth.prefix = 's';

    matlabbatch{17}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{17}.spm.util.imcalc.output = 'amygdala_mask_diff_inc';
    matlabbatch{17}.spm.util.imcalc.outdir = {['/data/underworld/home/kbas/03_data/processed_dif/' subject]};
    matlabbatch{17}.spm.util.imcalc.expression = 'i1>0.05';
    matlabbatch{17}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{17}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{17}.spm.util.imcalc.options.mask = 0;
    matlabbatch{17}.spm.util.imcalc.options.interp = 1;
    matlabbatch{17}.spm.util.imcalc.options.dtype = 4;
    

    % zip amygdala mask
    matlabbatch{18}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(1) = cfg_dep('Image Calculator: ImCalc Computed Image: amygdala_mask_diff_inc', substruct('.','val', '{}',{17}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{18}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = {''};
    matlabbatch{18}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep = true;
    spm_jobman('run',matlabbatch);

    
    spm_jobman('run',matlabbatch);
end
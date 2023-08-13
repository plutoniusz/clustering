subjects = [112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 346878];
maps = ["mt", "pd", "r1", "r2s"];

for i = 1:numel(subjects)
    subject = num2str(subjects(i));

    for j = 1:numel(maps);
        map = convertStringsToChars(maps(j));
        text = strcat('^average_', map, '.nii$');
        thalamus_mask = convertStringsToChars(['thalamus_' map '_qmri_' subject]);
        amygdala_mask = convertStringsToChars(['amygdala_' map '_qmri_' subject]);

        spm('defaults','fmri');
        spm_jobman('initcfg');

        % read map from a directoryS
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/112111']};
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = 'label_average_r1.nii$';
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';

        % read map from a directory
        matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.dir(1) = {['/data/underworld/home/kbas/clustering/2023-07-26_17-42-58/' subject]};
        matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.filter = text;
        matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPListRec';

        % calculate thalamus mask
        matlabbatch{3}.spm.util.imcalc.input(1) = cfg_dep('File Selector (Batch Mode): Selected Files (label_average_r1.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{3}.spm.util.imcalc.input(2) = cfg_dep('File Selector (Batch Mode): Selected Files (^average_ map \.nii$)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{3}.spm.util.imcalc.output = thalamus_mask;
        matlabbatch{3}.spm.util.imcalc.outdir = {['/data/underworld/home/kbas/03_data/processed_dif/' subject]};
        matlabbatch{3}.spm.util.imcalc.expression = '((i1==59) + (i1==60)).*i2';
        matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{3}.spm.util.imcalc.options.mask = 0;
        matlabbatch{3}.spm.util.imcalc.options.interp = 1;
        matlabbatch{3}.spm.util.imcalc.options.dtype = 4;

        % calculate amygdala mask
        matlabbatch{4}.spm.util.imcalc.input(1) = cfg_dep('File Selector (Batch Mode): Selected Files (label_average_r1.nii$)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{4}.spm.util.imcalc.input(2) = cfg_dep('File Selector (Batch Mode): Selected Files (^average_ map \.nii$)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{4}.spm.util.imcalc.output =amygdala_mask;
        matlabbatch{4}.spm.util.imcalc.outdir = {['/data/underworld/home/kbas/03_data/processed_dif/' subject]};
        matlabbatch{4}.spm.util.imcalc.expression = '((i1==31) + (i1==32)).*i2';
        matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{4}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{4}.spm.util.imcalc.options.mask = 0;
        matlabbatch{4}.spm.util.imcalc.options.interp = 1;
        matlabbatch{4}.spm.util.imcalc.options.dtype = 4;

        spm_jobman('run',matlabbatch);

    end
end
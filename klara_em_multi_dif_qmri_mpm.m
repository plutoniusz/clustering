% X - 2d intensity histogram, dimenstions MxN,
%     where N is the number of voxels in seed region
%     and M is the number of target voxels
% P - matrix of means of multinominal distributions MxK
% g - mixing proportions of the clusters
% Z - latent variable matrix such that p(z_n|g) = \prod_{k=1}^K g_k^{z_{kn}}
% K - number of clusters
% N - number of seed voxels
% r - responsibilities kxn

%%
if true
    subjects = [112111 128221 130519 170192 176117 208010 211787 214685 232237 308597 324038 330406 346878];
    gen_path =  "/data/underworld/kbas/03_data";
    sub_path = [];
    for subject=subjects
        directories = dir(gen_path + '/derivatives_mpm/' + subject);
        for i=1:numel(directories)
            if contains(directories(i).name, '20')
                sub_path = [sub_path; string(directories(i).name)];
                break
            end
        end
    end

    maps_gen = '/data/underworld/kbas/clustering/2024-02-05_10-17-29';


    gen_list = [];
    path_list = [];
    def_list = [];
    b0_list = [];
    mt_list = [];
    pd_list = [];
    r1_list = [];
    r2_list = [];

    for i=1:numel(subjects)

        num = num2str(sprintf('%05d', i));

        % diffusion folders
        path = strcat(gen_path, '/derivatives_mpm/', string(subjects(i)), '/', sub_path(i), '/dwi/fsl-probtrackx-1');
        gen_list = [gen_list; path];
        % connectivity matrix files
        path = strcat(path, '/fdt_matrix2.dot');
        path_list = [path_list; path];
        % deformations path
        path = strcat(gen_path, '/processed_mpm/', 'y_1_', num, '_sub-', string(subjects(i)), '_ses-', sub_path(i), '_desc-average-sbref_b0_mb_mpm.nii');
        % path = strcat(gen_path, '/processed_mpm/', string(subjects(i)), '/y_r_deformation_', string(subjects(i)) ,'_dif', '.nii');
        def_list = [def_list; path];
        % b0 files path
        path = strcat(gen_path, '/derivatives_mpm/', string(subjects(i)), '/', sub_path(i), '/dwi/qmap-preproc-b0/sub-', string(subjects(i)), '_ses-', sub_path(i), '_desc-average-sbref_b0.nii');
        b0_list = [b0_list; path];
%         % parameter maps path
%         path = strcat(maps_gen, '/', string(subjects(i)), '/average_mt.nii');
%         mt_list = [mt_list; path];
%         path = strcat(maps_gen, '/', string(subjects(i)), '/average_pd.nii');
%         pd_list = [pd_list; path];
%         path = strcat(maps_gen, '/', string(subjects(i)), '/average_r1.nii');
%         r1_list = [r1_list; path];
%         path = strcat(maps_gen, '/', string(subjects(i)), '/average_r2s.nii');
%         r2_list = [r2_list; path];
        % parameter maps path
        %/data/underworld/kbas/03_data/processed_mpm/112111/amygdala_mt_mpm_112111.nii
        path = strcat(gen_path, '/processed_mpm/', string(subjects(i)), '/amygdala_mt_mpm_', string(subjects(i)), '.nii');
        mt_list = [mt_list; path];
        path = strcat(gen_path, '/processed_mpm/', string(subjects(i)), '/amygdala_pd_mpm_', string(subjects(i)), '.nii');
        pd_list = [pd_list; path];
        path = strcat(gen_path, '/processed_mpm/', string(subjects(i)), '/amygdala_r1_mpm_', string(subjects(i)), '.nii');
        r1_list = [r1_list; path];
        path = strcat(gen_path, '/processed_mpm/', string(subjects(i)), '/amygdala_r2s_mpm_', string(subjects(i)), '.nii');
        r2_list = [r2_list; path];





    end

    %%
    x_list = cell(numel(subjects),1);
    xq_list = cell(numel(subjects),1);
    xq_whole_list = cell(numel(subjects),1);
    whole_list = cell(numel(subjects),1);
    phi_list = cell(numel(subjects),1);
    phi_list_2 = cell(numel(subjects),1);
    addpath([getenv('FSLDIR') '/etc/matlab']);
    
    %figure;
    ind_all = [];
    ind_wb_all = [];
    ind_all_q = [];

    %maska = niftiread([gen_list{i} '/fdt_paths.nii.gz']);
    [mask,~,scales] = read_avw([gen_list{1} '/fdt_paths.nii.gz']);
    [mask_t,~,scales_t] = read_avw('/data/underworld/kbas/03_data/processed_mpm/mpm_warped_112111sMP02874-0010-00001-000224-01.nii');

    %save_avw(mask, [gen_list{i} '/' num2str(i) '_mask.nii'] ,'i',scales);
    %save_avw(mask_t, [gen_list{i} '/' num2str(i) '_mask_t.nii'] ,'i',scales_t);

    for i=1:numel(subjects)
        disp(i)
        % connectivity matrix
        x = spconvert(load(path_list(i)))';
        % converting deformation to sparse format
        coord = load([gen_list{i} '/coords_for_fdt_matrix2'])+1; % correcting for matlab indexing
        ind   = sub2ind(size(mask),coord(:,1),coord(:,2),coord(:,3));
        coord_wb = load([gen_list{i} '/tract_space_coords_for_fdt_matrix2'])+1;
        ind_wb = sub2ind(size(mask), coord_wb(:,1), coord_wb(:,2), coord_wb(:,3));

        pom = kb_invdef2sparse(num2str(def_list(i)), num2str(b0_list(i)), '/data/underworld/kbas/03_data/processed_mpm/softmax_mb_mpm.nii');
        pom = load(pom.sparse{1});
        phi_list{i} = pom.Phi;
        % phi_list{i} = pom.Phi(:,ind_wb);
        % phi_list_2{i} = pom.Phi(:,ind);
        clear pom;

        [a1,a2] = size(phi_list{i});

        % %maska = niftiread([gen_list{i} '/fdt_paths.nii.gz']);
        % [mask,~,scales] = read_avw([gen_list{i} '/fdt_paths.nii.gz']);
        % [mask_t,~,scales_t] = read_avw('/data/underworld/kbas/03_data/processed_mpm/mpm_warped_112111sMP02874-0010-00001-000224-01.nii');

        % %save_avw(mask, [gen_list{i} '/' num2str(i) '_mask.nii'] ,'i',scales);
        % %save_avw(mask_t, [gen_list{i} '/' num2str(i) '_mask_t.nii'] ,'i',scales_t);
        % mask = 0*mask;
        coord = load([gen_list{i} '/coords_for_fdt_matrix2'])+1; % correcting for matlab indexing
        ind   = sub2ind(size(mask),coord(:,1),coord(:,2),coord(:,3));
        coord_wb = load([gen_list{i} '/tract_space_coords_for_fdt_matrix2'])+1;
        ind_wb = sub2ind(size(mask), coord_wb(:,1), coord_wb(:,2), coord_wb(:,3));

        whole_list{i} = sparse(a2, a2);
        whole_list{i}(ind_wb, ind) = x;
        clear x;
        %mask = mask*0;
        %mask(ind_wb) = 1;
        %save_avw(mask, [gen_list{i} '/' num2str(i) '_whole_list.nii'] ,'i',scales);
        %spy(whole_list{i});
        %text = string(['whole_list_' i '.fig'])
        %saveas(gcf,text);

        %whole_list{i} = phi_list{i}*whole_list{i}(ind_wb, ind)*phi_list_2{i}';
        whole_list{i} = phi_list{i}*whole_list{i}*phi_list{i}';

        %phi = full(sum(phi_list{i},1));
        %phi= reshape(phi, [140,136,90]);
        %text = join( ["phi_", string(subjects(i)), ".nii"])
        %niftiwrite(phi, text);

        %phi = full(sum(phi_list{i},2));
        %size(phi);
        %phi= reshape(phi, [227,422,463]);
        %text = join( ["phi_b_", string(subjects(i)), ".nii"])
        %niftiwrite(phi, text);

        %sum_a = full(sum(whole_list{i},1));
        %sum_b = full(sum(whole_list{i},2));
        %sum_a = reshape(sum_a, [227,422,463]);
        %sum_b = reshape(sum_b, [227,422,463]);
        %text = join( ["sum_a_", string(subjects(i)), ".nii"])
        %niftiwrite(sum_a, text);
        %text = join(["sum_b_", string(subjects(i)), ".nii"])
        %niftiwrite(sum_b, text);

        xq_whole_list{i} = permute(cat(4,niftiread(pd_list(i)), niftiread(mt_list(i)), niftiread(r1_list(i)), niftiread(r2_list(i))), [4,1,2,3]);


        %[a,b] = find(whole_list{i});
        % mask_t = mask_t*0;
        % mask_t(a) = 1;
        % save_avw(mask_t, [gen_list{i} '/' num2str(i) '_whole_list_warped.nii'] ,'i',scales_t);


        %ind_all = [ind_all; b];
        %ind_wb_all = [ind_wb_all; a];
        path = convertStringsToChars([fullfile(gen_list{i}, 'whole_list.mat')]);
        x = whole_list{i}; %(logical(whole_list{i}));
        save(path, 'x', '-v7.3');
        whole_list{i} = fullfile(gen_list{i}, 'whole_list.mat');

        %aq = find(xq_whole_list{i}(1,:,:,:));
        % %disp(aq)
        %ind_all_q = [ind_all_q; aq];
        % %save("aq.mat","aq", '-v7.3');
        % %disp("aq")

        path = convertStringsToChars([fullfile(gen_list{i}, 'xq_whole_list.mat')]);
        xq = xq_whole_list{i};
        save(path, 'xq', '-v7.3');
        xq_whole_list{i} = fullfile(gen_list{i}, 'xq_whole_list.mat');


        % % mask files
        % mask = 0*mask;
        % mask(ind) = 1;
        % save_avw(mask, [gen_list{i} '/' num2str(i) '_ind.nii'] ,'i',scales);
         
        % mask = 0*mask;
        % mask(ind_wb) = 2;
        % save_avw(mask, [gen_list{i} '/' num2str(i) '_ind_wb.nii'] ,'i',scales);
         
        % %mask = zeros(227,422,463);
        % mask_t = 0*mask_t;
        % mask_t(b) = 3;
        % save_avw(mask_t, [gen_list{i} '/' num2str(i) '_ind_warped.nii'] ,'i',scales_t);

        % mask_t = 0*mask_t;
        % mask_t(a) = 4;
        % save_avw(mask_t, [gen_list{i} '/' num2str(i) '_ind_wb_warped.nii'] ,'i',scales_t);       

        % mask_t = 0*mask_t;
        % mask_t(:) = xq_whole_list{i}(1,:);
        % save_avw(mask_t, [gen_list{i} '/' num2str(i) '_pd.nii'] ,'i',scales_t);

        clear phi_list{i} x xq

    end

%     %%
%     ind_all = unique(ind_all);
%     ind_wb_all = unique (ind_wb_all);
%     %figure;
%     for i=1:numel(subjects)
%         x_list{i} = whole_list{i}(ind_wb_all, ind_all);
%         xq_list{i} = zeros(4,numel(ind_all));
%         for j=1:4
%             xq_list{i}(j,:) = xq_whole_list{i}(j, ind_all);
% 
%             mask = 0*mask;
%             mask(:) = xq_whole_list{i}(j,:);
%             save_avw(mask, [gen_list{i} '/' num2str(i) '_map_' num2str(j) '.nii'] ,'i',scales);
%         end
%         %spy(x_list{i});
%         %hold on;
%     end
%     %hold off;
    %%

    aq = find(xq_whole_list{1}(1,:,:,:));
    [a,b] = find(whole_list{1});
    ind_all = b; [ind_all; b];
    ind_wb_all = a; [ind_wb_all; a];
    %ind_all = unique(ind_all);
    %ind_wb_all = unique (ind_wb_all);
    ind_all_q = aq;
    %ind_all_q = unique(ind_all_q);
    %save("ind_all_q.mat","ind_all_q", '-v7.3');
    %disp("size ind_all_q")
    %disp(size(ind_all_q))

    %figure;
    for i=1:numel(subjects)
        x_list{i} = whole_list{i}(ind_wb_all, ind_all_q);
        %disp(size(x_list{i}));
        xq_list{i} = zeros(4,numel(ind_all_q));
        for j=1:4
            xq_list{i}(j,:) = xq_whole_list{i}(j,ind_all_q);
            %disp(size(xq_list{i}));
            %xq_list{i}(j,:) = xq_whole_list{i}(j, ind_all);

            mask_t = 0*mask_t;
            mask_t(:) = xq_whole_list{i}(j,:);
            save_avw(mask_t, [gen_list{i} '/' num2str(i) '_map_' num2str(j) '.nii'] ,'i',scales_t);
        end
        %spy(x_list{i});
        %hold on;
    end
    %hold off;
    %%

   
else
    K  = 10;
    N  = 1000;
    M  = 20;
    g0 = softmax(randn(K,1)*0.5);
    P0 = softmax(randn(M,K),1);
    R0 = multinom_random(repmat(g0,[1 N]),1);
    X  = multinom_random(P0*R0,1000);
end

% %%
% for n=1:13
%     figure;
%     h1 = histogram(xq_list{n}(1,:))
% end

%%
K = 6 ;

%save("x_list.mat","x_list", '-v7.3');
%save("xq_list.mat","xq_list", '-v7.3');

[ELBO,Alpha,Beta, list_R, g]  = mixcode1joint_working(x_list,xq_list,K);
%[ELBO,Alpha,Beta, list_R]  = mixcode1qmri_common(x_list,xq_list,K);
%[ELBO,Alpha,Beta,alpha0,beta0,K,list_R,lG]  = mixcode2joint(x_list,xq_list,K);

save("list_R_common.mat","list_R", '-v7.3');
% save("ind_all.mat")
% save("K.mat")
% save("mask.mat")
% save("gen_list.mat")
% save("scales")
% save("g")

%%
% separate file for each cluster for each subject
for s=1:13
    for i=1:K
        %disp(size(mask))
        %mask = zeros(227,422,463);
        %mask = zeros(140,143,103);
        mask_t = mask_t*0;
        %disp(sum(list_R{s}(i,:)));
        mask_t(ind_all_q) = full(list_R{s}(i, :));
        save_avw(mask_t, [gen_list{s} '/clusters_joint_working_' num2str(s) '_' num2str(i)] ,'i',scales_t);
        %niftiinfo([gen_list{s} '/clusters_template_new_' num2str(s) '_' num2str(i)]).raw = niftiinfo(b0_list{i}).raw;        
        %niftiwrite(mask, [gen_list{s} '/clusters_template_new_test' num2str(s) '_' num2str(i)])%,  niftiinfo(b0_list{i}))nif
    end
end

%%
% for each subject single file with most probable cluster
for s=1:13
    %mask = zeros(140,143,103);
    %mask = zeros(227,422,463);
    [~, idx] = max(full(list_R{s}), [], 1);
    mask_t = mask_t*0;
    mask_t(ind_all_q) = idx;
    text = [gen_list{s} '/clusters_joint_working_' num2str(s)];
    save_avw(mask_t, text ,'i',scales_t);
    gunzip([text '.nii.gz']);
end

%%
% for each cluster a tissue prior
%G = exp(g);
exp(sum(g,1))
for i=1:K
        %mask = zeros(140,143,103);
        %mask = zeros(227,422,463);
        mask_t = mask_t*0;
        G = exp(g);
        mask_t(ind_all_q) = G(i, :);
        save_avw(mask_t, [gen_list{1} '/clusters_joint_working_G_' num2str(i)] ,'i',scales_t);
end


%%
% separate file for each cluster for each subject
%// for s=1:13
%//     for i=1:56
%//         mask = zeros(140,143,103);
%//         mask(ind_all) = full(list_R{s}(i, :));
%//         save_avw(mask, [gen_list{s} '/clusters_joint_optimised_' num2str(s) '_' num2str(i)] ,'i',scales);
%//         %niftiinfo([gen_list{s} '/clusters_template_new_' num2str(s) '_' num2str(i)]).raw = niftiinfo(b0_list{i}).raw;        
%//         %niftiwrite(mask, [gen_list{s} '/clusters_template_new_test' num2str(s) '_' num2str(i)])%,  niftiinfo(b0_list{i}))nif
%//     end
%// end

%%
% for each subject single file with most probable cluster
% for s=1:13
%     gunzip([gen_list{s} '/clusters_joint_common_' num2str(s) '.nii.gz']);
% end








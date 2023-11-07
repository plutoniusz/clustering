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
        directories = dir(gen_path + '/derivatives_dif/' + subject);
        for i=1:numel(directories)
            if contains(directories(i).name, '20')
                sub_path = [sub_path; string(directories(i).name)];
                break
            end
        end
    end

    maps_gen = '/data/underworld/kbas/clustering/2023-07-26_17-42-58';


    gen_list = [];
    path_list = [];
    def_list = [];
    b0_list = [];
    mt_list = [];
    pd_list = [];
    r1_list = [];
    r2_list = [];

    for i=1:numel(subjects)

        % diffusion folders
        path = strcat(gen_path, '/derivatives_dif/', string(subjects(i)), '/', sub_path(i), '/dwi/fsl-probtrackx-1');
        gen_list = [gen_list; path];
        % connectivity matrix files
        path = strcat(path, '/fdt_matrix2.dot');
        path_list = [path_list; path];
        % deformations path
        path = strcat(gen_path, '/processed_dif/', string(subjects(i)), '/y_r_deformation_', string(subjects(i)) ,'_dif', '.nii');
        def_list = [def_list; path];
        % b0 files path
        path = strcat(gen_path, '/derivatives_dif/', string(subjects(i)), '/', sub_path(i), '/dwi/qmap-preproc-b0/sub-', string(subjects(i)), '_ses-', sub_path(i), '_desc-average-sbref_b0.nii');
        b0_list = [b0_list; path];
        % parameter maps path
        path = strcat(maps_gen, '/', string(subject), '/average_mt.nii');
        mt_list = [mt_list; path];
        path = strcat(maps_gen, '/', string(subject), '/average_pd.nii');
        pd_list = [pd_list; path];
        path = strcat(maps_gen, '/', string(subject), '/average_r1.nii');
        r1_list = [r1_list; path];
        path = strcat(maps_gen, '/', string(subject), '/average_r2s.nii');
        r2_list = [r2_list; path];

    end

    %%
    x_list = cell(numel(subjects),1);
    xq_list = cell(numel(subjects),1);
    xq_whole_list = cell(numel(subjects),1);
    whole_list = cell(numel(subjects),1);
    phi_list = cell(numel(subjects),1);
    addpath([getenv('FSLDIR') '/etc/matlab']);
    
    %figure;
    ind_all = [];
    ind_wb_all = [];

    for i=1:numel(subjects)

        x = spconvert(load(path_list(i)))';

        pom = kb_def2sparse(num2str(def_list(i)),num2str(b0_list(i)),'/data/underworld/kbas/03_data/processed_dif/softmax_mb.nii');
        pom = load(pom.sparse{1});
        phi_list{i} = pom.Phi;

        [a1,~] = size(phi_list{i});

        maska = niftiread([gen_list{i} '/fdt_paths.nii.gz']);
        [mask,~,scales] = read_avw([gen_list{i} '/fdt_paths.nii.gz']);

%         save_avw(mask, [gen_list{i} '/' num2str(i) '_mask.nii'] ,'i',scales);
        mask = 0*mask;
        coord = load([gen_list{i} '/coords_for_fdt_matrix2'])+1; % correcting for matlab indexing
        ind   = sub2ind(size(mask),coord(:,1),coord(:,2),coord(:,3));
        coord_wb = load([gen_list{i} '/tract_space_coords_for_fdt_matrix2'])+1;
        ind_wb = sub2ind(size(mask), coord_wb(:,1), coord_wb(:,2), coord_wb(:,3));

        whole_list{i} = sparse(a1, a1);
        whole_list{i}(ind_wb, ind) = x;
        whole_list{i} = phi_list{i}'*whole_list{i}*phi_list{i};

        [a,b] = find(whole_list{i});

        ind_all = [ind_all; b];
        ind_wb_all = [ind_wb_all; a];
        
        xq_whole_list{i} = permute(cat(4,niftiread(pd_list(i)), niftiread(mt_list(i)), niftiread(r1_list(i)), niftiread(r2_list(i))), [4,1,2,3]);

        % mask files
%         mask = 0*mask;
%         mask(ind) = 1;
%         save_avw(mask, [gen_list{i} '/' num2str(i) '_ind.nii'] ,'i',scales);
%         
%         mask = 0*mask;
%         mask(ind_wb) = 2;
%         save_avw(mask, [gen_list{i} '/' num2str(i) '_ind_wb.nii'] ,'i',scales);
%         
%         mask = zeros(140,143,103);
%         mask = 0*mask;
%         mask(b) = 3;
%         save_avw(mask, [gen_list{i} '/' num2str(i) '_ind_warped.nii'] ,'i',scales);
% 
%         mask = 0*mask;
%         mask(a) = 4;
%         save_avw(mask, [gen_list{i} '/' num2str(i) '_ind_wb_warped.nii'] ,'i',scales);       

%         mask = 0*mask;
%         mask(:) = xq_whole_list{i}(1,:);
%         save_avw(mask, [gen_list{i} '/' num2str(i) '_pd.nii'] ,'i',scales);

    end

    %%
    ind_all = unique(ind_all);
    ind_wb_all = unique (ind_wb_all);
    %figure;
    for i=1:numel(subjects)
        x_list{i} = whole_list{i}(ind_wb_all, ind_all);
        xq_list{i} = zeros(4,numel(ind_all));
        for j=1:4
            xq_list{i}(j,:) = xq_whole_list{i}(j, ind_all);

%             mask = 0*mask;
%             mask(:) = xq_whole_list{i}(j,:);
%             save_avw(mask, [gen_list{i} '/' num2str(i) '_map_' num2str(j) '.nii'] ,'i',scales);
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

%%
K = 20;
ll_list = [];
%disp(size(x_list))
%disp(size(xq_list))
%disp(K)
[ELBO,Alpha,Beta, list_R]  = mixcode1joint(x_list,xq_list,K);

%%
%figure;
%plot(ll_list);
%savefig('log_likelihood_1');

for s=1:13
    for i=1:K
        mask = mask*0;
        mask(ind_all) = full(list_R{s}(i, :));
        save_avw(mask, [gen_list{s} '/clusters_joint_' num2str(s) '_' num2str(i)] ,'i',scales);
        %niftiinfo([gen_list{s} '/clusters_template_new_' num2str(s) '_' num2str(i)]).raw = niftiinfo(b0_list{i}).raw;        
        %niftiwrite(mask, [gen_list{s} '/clusters_template_new_test' num2str(s) '_' num2str(i)])%,  niftiinfo(b0_list{i}))nif
    end
end







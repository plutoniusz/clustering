% X - 2d intensity histogram, dimenstions MxN,
%     where N is the number of voxels in seed region
%     and M is the number of target voxels
% P - matrix of means of multinominal distributions MxK
% g - mixing proportions of the clusters
% Z - latent variable matrix such that p(z_n|g) = \prod_{k=1}^K g_k^{z_{kn}}
% K - number of clusters
% N - number of seed voxels
% r - responsibilities kxn

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
    gen_list = [];
    path_list = [];
    def_list = [];
    b0_list = [];
    for i=1:numel(subjects)
        path = strcat(gen_path, '/derivatives_dif/', string(subjects(i)), '/', sub_path(i), '/dwi/fsl-probtrackx-1');
        gen_list = [gen_list; path];
        path = strcat(path, '/fdt_matrix2.dot');
        path_list = [path_list; path];
        path = strcat(gen_path, '/processed_dif/', string(subjects(i)), '/y_r_deformation_', string(subjects(i)) ,'_dif', '.nii');
        def_list = [def_list; path];
        path = strcat(gen_path, '/derivatives_dif/', string(subjects(i)), '/', sub_path(i), '/dwi/qmap-preproc-b0/sub-', string(subjects(i)), '_ses-', sub_path(i), '_space-orig_desc-dwi-skullstripped_b0.nii');
        b0_list = [b0_list; path];        
    end

    x_list = cell(14,1);
    whole_list = cell(14,1);
    phi_list = cell(14,1);
    addpath([getenv('FSLDIR') '/etc/matlab']);
    
    %figure;
    ind_all = [];
    ind_wb_all = [];
    for i=1:13
        x = spconvert(load(path_list(i)))';
        %spy(x);
        disp(num2str(def_list(i)))
        [phi_list{i},~,~] = spm_def2sparse(num2str(def_list(i)), num2str(b0_list(i)));
        [a1,~] = size(phi_list{i});

        [mask,~,scales] = read_avw([gen_list{i} '/fdt_paths.nii.gz']);
        mask = 0*mask;
        coord = load([gen_list{i} '/coords_for_fdt_matrix2'])+1; % correcting for matlab indexing
        ind   = sub2ind(size(mask),coord(:,1),coord(:,2),coord(:,3));
        coord_wb = load([gen_list{i} '/tract_space_coords_for_fdt_matrix2'])+1;
        ind_wb = sub2ind(size(mask), coord_wb(:,1), coord_wb(:,2), coord_wb(:,3));

        whole_list{i} = sparse(a1, a1);
        whole_list{i}(ind_wb, ind) = x;
        whole_list{i} = phi_list{i}'*whole_list{i}*phi_list{i};

        %spy(whole_list{i});
        %hold on;
        
        [a,b] = find(whole_list{i});
        ind_all = [ind_all; b];
        ind_wb_all = [ind_wb_all; a];

    end
    %hold off
    %disp(size(ind_all));
    ind_all = unique(ind_all);
    %disp(size(ind_all))
    %disp(size(ind_wb_all))
    ind_wb_all = unique (ind_wb_all);
    %disp(size(ind_wb_all))
    %figure;
    for i=1:13
        x_list{i} = whole_list{i}(ind_wb_all, ind_all);
        %spy(x_list{i});
        %hold on;
    end
    %hold off;
   

else
    K  = 10;
    N  = 1000;
    M  = 20;
    g0 = softmax(randn(K,1)*0.5);
    P0 = softmax(randn(M,K),1);
    R0 = multinom_random(repmat(g0,[1 N]),1);
    X  = multinom_random(P0*R0,1000);
end

K = 20;
[M,N] = size(x_list{i});
P     = exp(randn(M,K)*0.01);
P     = log(P./sum(P,1));
g     = log(ones(K,1)/K);

ll_list = [];
[P,g,R,ll,ll_list] = em(x_list,P,g,1000,ll_list);
figure;
plot(ll_list);
savefig('log_likelihood_1');


for s=1:13
    for i=1:K
        mask(ind_all) = full(R{s}(i, :));
        save_avw(mask, [gen_list{s} '/' num2str(s) '_' num2str(i)] ,'i',scales);
        %niftiinfo([gen_list{s} '/clusters_template_new_' num2str(s) '_' num2str(i)]).raw = niftiinfo(b0_list{i}).raw;        
        %niftiwrite(mask, [gen_list{s} '/clusters_template_new_test' num2str(s) '_' num2str(i)])%,  niftiinfo(b0_list{i}))nif
    end
end




function Z = multinom_random(P,Nsamp)
    if nargin<2
        N = 1;
    end
    [M,N] = size(P);
    U     = cumsum(P,1); %cumulative column sum  % Upper value
    L     = [zeros(1,N); U(1:(end-1),:)]; % Lower values
    Z     = zeros(M,N);
    for i=1:Nsamp
        t = rand(1,N);
        s = sum((t>=L & t<U).*(1:M)',1);
        Z = Z + full(sparse(s,1:N,1,M,N));
    end
end


function [P,g,R,ll,ll_list] = em(x_list,P,g,nit, ll_list)
    ll = -Inf;
    for iter= 1:nit
        disp(iter)
        r = 0;
        xr = 0;
        for i=1:13
            R{i} = e_step(x_list{i},P,g);
            xr = xr + x_list{i}*R{i}';
            r = r + R{i};
        end
        [P,g] = m_step(xr,r);
        ll_o  = ll;
        ll_list = [ll_list, ll_o];
        ll    = loglikelihood(x_list,P,g,R);
        if abs(ll-ll_o) < abs(ll*1e-9); break; end
    end
%     for i =1:13
%         ll = ll + factorial_stuff(x_list{i});
%     end
end


function ll = loglikelihood(x_list,P,g,R)
    ll = 0;
    for i=1:13
        ll = ll + sum(LSE((P'*x_list{i} + g)'*R{i},1),2);

    end
end


function R = e_step(xs,P,g)
    %R = (softmax(log(P)'*xs + log(g),1));
    R = sparse(softmax(P'*xs + g,1));
end


function [P, g] = m_step(xr,r)
    alpha0 = 1e-3; % Behaves like having a Dirichlet prior
    alpha = xr + alpha0;
    P = psi(alpha) - psi(sum(alpha,1));
    %P = P./sum(P,1); %dlaczego tak a nie *p(old)
    %g = r + 1e-10; % To prevent numerical problems
    beta0 = 1e-3;
    beta = r + beta0;
    %g = g/sum(g);
    g = psi(beta)-psi(sum(beta,1));
end


function S = softmax(A,dim)
    % exp(A)./sum(exp(A),dim)
    if nargin<2
        dim = 1;
    end
    t = max(A,[],dim);
    T = exp(A - t);
    S = T./sum(T,dim);
end

function L = LSE(A,dim)
    % log(sum(exp(A)),dim)
    % https://en.wikipedia.org/wiki/LogSumExp
    if nargin<2
        dim = 1;
    end
    t = max(A,[],dim);
    L = log(sum(exp(A - t))) + t;
end


function ll_const = factorial_stuff(X)
    ll_const = sum(logfactorial(sum(X,1)) - sum(logfactorial(X),1));
end


function L = logfactorial(X)
    % log(X!)
    L = gammaln(sparse(X)+1);
end



function [ELBO,Alpha,Beta, list_R, lG] = fit_model(X,Xq,K,alpha0,beta0)
    if nargin<3, K      = 10;    end
    if nargin<4, alpha0 =  0.1; end
    if nargin<5, beta0  =  0.1; end

    if not(numel(X)==numel(Xq)), disp('missing data for one or more subjects'); end

    % Initialise model
    N      = numel(X);
    [M,J]  = size(X{1});
    alpha0 = alpha0*ones(M,1);
    beta0  =  beta0*ones(K,1);

    % Different initialisations will give different results
    Alpha  = exp(randn(M,K)*0.001);
    Beta   = exp(randn(K,J)*0.001);

    % qmi priors
    nu0 = ((100-3-eps)*rand(1)+3+eps);
    L0 = tril(randn(4,4));
    W0  = L0*L0'*0.001;
    c0  = (rand(1)+eps)*0.001;
    mu0 = randn(4,1)*0.001;

    % initialise update equations and sufficient statistics
    R = rand(K,J);
    R = R./sum(R,2);
    [R_sum, x_mean, S_mean] = suffsta(R, Xq{1}, K, J);
    [c,nu,W,mu] = update_quations(c0, nu0, W0, mu0, R_sum, x_mean, S_mean, K);
    nu = repmat(nu,[1,1,N]);
    c = repmat(c,[1,1,N]);
    W = repmat(W,[1,1,1,N]);
    mu = repmat(mu,[1,1,N]);
   

    elbo0 = 0;
    % for n=1:N
    %     elbo0 = elbo0 + sum(C1(X{n})); % memory issues
    % end
    elbo0 = elbo0 - N*K*logB(W0,nu0);

    ELBO   = -Inf;
    for it=1:1000
        ELBO_prev = ELBO;
        lP  = psi(Alpha) - psi(sum(Alpha, 1)); % E[ln P]
        lG  = psi(Beta)  - psi(sum(Beta,  1)); % E[ln G]
        
        A   = 0;
        B   = 0;

        lse = 0;
        elbo_pom = 0;

        XML = []; %equation 64
        lL = []; %equation 65
        R = [];
        
        for n=1:N
            
            for j=1:J
                for k=1:K
                    XML{n}(k,j) = 4*c(:,k,n)^(-1) + nu(:,k,n)*(Xq{n}(:,j)- mu(:,k,n))'*W(:,:,k,n)*(Xq{n}(:,j)-mu(:,k,n));
                    lL{n}(k)  = Elogdet(W(:,:,k,n),nu(:,k,n));
                end
            end

            [R{n},lse_n] = softmax(lP'*X{n} + lG + lL{n}'./2 - XML{n}./2-2*log(2*pi)); % R = E[Z]
            A   = A + X{n}*R{n}';
            B   = B + R{n};
            lse = lse + sum(lse_n);

        end

        elbo71 = 0;
        elbo74 = 0;
        elbo77 = 0;

        R_sum = [];
        x_mean = [];
        S_mean = [];

        for n=1:N

            x_mean{n} = zeros(4,K);
            S_mean{n} = zeros(4,4,k);

            [R_sum{n}, x_mean{n}, S_mean{n}] = suffsta(R{n}, Xq{n}, K, J);

            for k=1:K
            
                %E[ln(p(xq|z,mu,L))]
                elbo71 = elbo71 + R_sum{n}(:,k)*(lL{n}(:,k)-4*c(:,k,n)^(-1)-nu(:,k,n)*trace(S_mean{n}(:,:,k)*W(:,:,k,n))-nu(:,k,n)*(x_mean{n}(:,k)-mu(:,k,n))'*W(:,:,k,n)*(x_mean{n}(:,k)-mu(:,k,n))-4*log(2*pi))/2;

                %E[ln(p(mu,L))]
                elbo74 = elbo74 + (4*log(c0/(2*pi)) + lL{n}(:,k) - 4*c0/c(:,k,n) - c0*nu(:,k,n)*(mu(:,k,n)-mu0)'*W(:,:,k,n)*(mu(:,k,n)-mu0) + (nu0-5)*lL{n}(:,k)-nu(:,k,n)*trace(W0^(-1)*W(:,:,k,n)))/2;

                %E[ln(q(mu,L))]
                elbo77 = elbo77 + lL{n}(:,k)/2 + 2*log(c(:,k,n)/(2*pi)) - 2 - H_Wishart(W(:,:,k,n), nu(:,k,n));
            end
            
            elbo_pom = elbo_pom +sum(R{N}.*(lL{n}'./2 - XML{n}./2-2*log(2*pi)), "all");

        end

        ELBO = elbo0 + sum(lse) + ...
        sum(sum((alpha0 - Alpha).*lP, 1) - C(Alpha) + C(alpha0), 2) +...
        sum(sum((beta0  - Beta ).*lG, 1) - C(Beta) + C(beta0), 2) +...
        elbo71+elbo74-elbo77-elbo_pom;
%%

        Alpha = A + alpha0;
        Beta  = B +  beta0;
        for n=1:N                
            c(:,:,n) = c0 + R_sum{n};
            nu(:,:,n) = nu0 + R_sum{n};

            for k=1:K              
                W(:,:,k,n) = (W0^(-1) + R_sum{n}(:,k)*S_mean{n}(:,:,k) + c0*R_sum{n}(:,k)*(x_mean{n}(:,k)-mu0)*(x_mean{n}(:,k)-mu0)'/(c0+R_sum{n}(:,k)))^(-1);
                mu(:,k,n) = (c0*mu0+x_mean{n}(:,k)*R_sum{n}(:,k))/c(:,k,n);
            end            
        end

        epsilon = eps(ELBO*0+numel(lP));
       
        if true
            % Display stuff for debugging purposes
            if ELBO<ELBO_prev, str = '<='; else str=''; end
            fprintf('%d %g%s\n', it, ELBO,str);
            %subplot(2,2,1); imagesc(Alpha./sum(Alpha)); drawnow
            %subplot(2,2,2); imagesc(Beta./sum(Beta)); drawnow
        end
        if ELBO - ELBO_prev < abs(ELBO*epsilon), break; end
    end
end

function [R,lse] = softmax(A, dim)
    if nargin<2, dim=1; end
    mx  = max(A,[],dim);
    R   = exp(A-mx);
    se  = sum(R,dim);
    R   = R./se;
    lse = log(se) + mx;
end

function lp = C(alpha)
    lp = gammaln(sum(alpha)) - sum(gammaln(alpha));
end

function lp = C1(alpha)
    lp = gammaln(sum(alpha,1)+1) - sum(gammaln(alpha+1),1);
end

function [c,nu,W,mu] = update_quations(c0, nu0, W0, mu0, R_sum, x_mean, S_mean, K)
    c(:,:) = c0 + R_sum;
    nu(:,:) = nu0 + R_sum;
        for k=1:K              
            W(:,:,k) = (W0^(-1) + R_sum(:,k)*S_mean(:,:,k) + c0*R_sum(:,k)*(x_mean(:,k)-mu0)*(x_mean(:,k)-mu0)'/(c0+R_sum(:,k)))^(-1);
            mu(:,k) = (c0*mu0+x_mean(:,k)*R_sum(:,k))/c(:,k);
        end
end

function [c,nu,W,mu] = update_quations_n(c0, nu0, W0, mu0, R_sum, x_mean, S_mean, K, n)
    c(:,:,n) = c0 + R_sum{n};
    nu(:,:,n) = nu0 + R_sum{n}
        for k=1:K              
            W(:,:,k,n) = (W0^(-1) + R_sum{n}(:,k)*S_mean{n}(:,:,k) + c0*R_sum{n}(:,k)*(x_mean{n}(:,k)-mu0)*(x_mean{n}(:,k)-mu0)'/(c0+R_sum{n}(:,k)))^(-1);
            mu(:,k,n) = (c0*mu0+x_mean{n}(:,k)*R_sum{n}(:,k))/c(:,k,n);
        end            
end

function [R_sum, x_mean, S_mean] = suffsta(R, Xq, K, J)
    x_mean = zeros(4,K);
    S_mean = zeros(4,4,K);
    R_sum = sum(R',1)+eps;
    for k=1:K
        x_mean(:,k) = Xq*R(k,:)'/(R_sum(:,k));
        for j=1:J
            S_mean(:,:,k) = S_mean(:,:,k) + R(k,j)*(Xq(:,j)-x_mean(:,k))*(Xq(:,j)-x_mean(:,k))'/(R_sum(:,k));
        end
    end
end

function h = H_Wishart(W,nu)
    D = size(W,1);
    h = -logB(W,nu) - (nu-D-1)*Elogdet(W,nu)/2 + nu*D/2;
end

function lb = logB(W,nu)
    D  = size(W,1);
    s  = 0;
    for i=1:D
        s = s + gammaln((nu+1-i)/2);
    end
    lb = -nu*logdet(W)/2 - (nu*D*log(2)/2 + D*(D-1)*log(pi)/4 + s);
end

function ld = Elogdet(W,nu)
    D  = size(W,1);
    ld = D*log(2) + logdet(W);
    for i=1:D
        ld = ld + psi((nu+1-i)/2);
    end
end

function ld = logdet(W)
    ld = sum(log(diag(chol(W))))*2;
end

function [ELBO,Alpha,Beta, list_R] = fit_model(X,Xq,K,alpha0,beta0)
    if nargin<3, K      = 10;    end
    if nargin<4, alpha0 =  0.1; end
    if nargin<5, beta0  =  0.1; end

    if not(numel(X)==numel(Xq)), disp('missing data for one or more subjects'); end

    % Initialise model
    N      = numel(X);
    N = 13;
    [M,J]  = size(X{1});
    % alpha0 = alpha0*ones(M,1,'like',X{1});
    % beta0  =  beta0*ones(K,1,'like',X{1});
    alpha0 = alpha0*ones(M,1);%,'like',X{1});
    beta0  =  beta0*ones(K,1);%,'like',X{1});

    % Different initialisations will give different results
    % Alpha  = exp(randn(M,K,'like',X{1})*0.001);
    % Beta   = exp(randn(K,J,'like',X{1})*0.001);
    Alpha  = exp(randn(M,K)*0.001);
    Beta   = exp(randn(K,J)*0.001);

    nu0 = (100-3-eps)*rand(1)+3+eps;
    L0 = tril(randn(4,4));
    W0  = L0*L0';
    c0  = rand(1)+eps;
    mu0 = randn(4,1)*0.001;

    %nu = max(exp(randn(1,K,N)),3.001);
    nu = (100-3-eps)*rand(1,K,N)+3+eps;
    for n=1:N
        for k=1:K
            L(:,:,k,n)=tril(randn(4,4));
            W(:,:,k,n) = L(:,:,k,n)*L(:,:,k,n)';
        end
    end
    c  = rand(K,1,N)+eps;
    mu =   randn(4,K,N)*0.001;

    elbo0 = 0;
%     for n=1:N
%         elbo0 = elbo0 + sum(C1(X{n})); % memory issues
%     end
    elbo_b = N*K*Bf(W0,nu0);

    ELBO   = -Inf;
    for it=1:1000
        ELBO_prev = ELBO;
        lP  = psi(Alpha) - psi(sum(Alpha, 1)); % E[ln P]
        lG  = psi(Beta)  - psi(sum(Beta,  1)); % E[ln G]
        

        A   = 0;
        B   = 0;
        lse = 0;
        XML = cell(N,1);
        lL = cell(N,1);


        for n=1:N    
            for j=1:J
                for k=1:K
                    det_W(k) = det(W(:,:,k,n));
                    XML{n}(k,j) = 4*c(k,:,n) + nu(:,k,n)*(Xq{n}(:,j)-mu(:,k,n))'*W(:,:,k,n)*(Xq{n}(:,j)-mu(:,k,n));
                end
            end
            sum_psi = psi((nu(:,:,n))/2)+psi((nu(:,:,n) - 1)/2)+psi((nu(:,:,n) - 2)/2)+psi((nu(:,:,n) -3 )/2);
            lL{n}  = sum_psi + log(16) + log(det_W);
            
            [R,lse_n] = softmax(lP'*X{n} + lG + lL{n}'./2 - XML{n}./2); % R = E[Z]
            %if n==1, subplot(2,2,3); imagesc(R); end
            A   = A + X{n}*R';
            B   = B + R;
            lse = lse + sum(lse_n);
            list_R{n} = R;
        end

        elbo71 = 0;
        elbo74 = 0;
        elbo77 = 0;
        for n=1:N
            %disp(n)
            R_sum{n} = sum(list_R{n},2)+eps;
            x_mean{n} = sum(Xq{n}*list_R{n}', 1)./R_sum{n}';
            for k=1:K
                %disp(k)
                S_mean{n}(k) = sum(list_R{n}(k,:)*(Xq{n}-x_mean{n}(k))'*(Xq{n}-x_mean{n}(k)),2);
                
                %E[ln(p(xq|z,mu,L))]
                elbo71 = elbo71 + R_sum{n}(k)*(lL{n}(:,k)-4*c(k,:,n)^(-1)-nu(:,k,n)*trace(S_mean{n}(k)*W(:,:,k,n))-nu(:,k,n)*(x_mean{n}(:,k)-mu(:,k,n))'*W(:,:,k,n)*(x_mean{n}(:,k)-mu(:,k,n))-4*log(2*pi))/2;
                %E[ln(p(mu,L))]
                elbo74 = elbo74 + 4*log(c0/(2*pi)) + lL{n}(:,k) - 4*c0/c(k,:,n) - c0*nu(:,k,n)*(mu(:,k,n)-mu0)'*W(:,:,k,n)*(mu(:,k,n)-mu0) + (nu0-5)*lL{n}(:,k)/2-nu(:,k,n)*trace(W0^(-1)*W(:,:,k,n))/2;
                %E[ln(q(mu,L))]
                elbo77 = elbo77 + lL{n}(:,k)/2 + 2*log(c(k,:,n)/(2*pi)) - 2 - H(lL{n}(k), W(:,:,k,n), nu(:,k,n));

            end
        end
        
        ELBO  = elbo0 + sum(lse) +...
        sum(sum((alpha0 - Alpha).*lP, 1) - C(Alpha) + C(alpha0), 2) +...
        sum(sum((beta0  - Beta ).*lG, 1) - C(Beta) + C(beta0), 2) +...
        elbo71+elbo74-elbo77-elbo_b;
%%

        Alpha = A + alpha0;
        Beta  = B +  beta0;
        for n=1:N
            for k=1:K              
                W(:,:,k,n) = (W0^(-1) + R_sum{n}(k)*S_mean{n}(k) + c0*R_sum{n}(k)*(x_mean{n}(k)-mu0)*(x_mean{n}(k)-mu0)')^-1;
            end
            c(:,:,n) = c0 + R_sum{n};
            nu(:,:,n) = nu0 + R_sum{n}';
            mu(:,:,n) = (c0*mu0+x_mean{n}*R_sum{n})./c(:,:,n)';
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

function lp = Bf(alpha, beta)
    lp = det(alpha)^(-beta/2)*(2^(2*beta)*pi^(3)+gamma(beta/2)+gamma((beta-1)/2)+gamma((beta-2)/2)+ gamma((beta-3)/2))^(-1);

end

function lp = H(alpha,beta,delta)
    disp("H")
    lp = -log(Bf(beta,delta))
    lp =  lp - (delta-5)/2*alpha+2*delta 
end


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

    nu0 = ((100-3-eps)*rand(1)+3+eps);
    L0 = tril(randn(4,4));
    W0  = L0*L0'*0.001;
    c0  = (rand(1)+eps)*0.001;
    mu0 = randn(4,1)*0.001;

    nu = (100-3-eps)*rand(1,K,N)+3+eps;
    for n=1:N
        for k=1:K
            L(:,:,k,n)=tril(randn(4,4));
            W(:,:,k,n) = L(:,:,k,n)*L(:,:,k,n)';
        end
    end
    c  = rand(1,K,N)+eps;
    mu =   randn(4,K,N)*0.001;

    elbo0 = 0;
%     for n=1:N
%         elbo0 = elbo0 + sum(C1(X{n})); % memory issues
%     end
    elbo_b = N*K*logB(W0,nu0);

    ELBO   = -Inf;
    for it=1:1000
        ELBO_prev = ELBO;
        lP  = psi(Alpha) - psi(sum(Alpha, 1)); % E[ln P]
        lG  = psi(Beta)  - psi(sum(Beta,  1)); % E[ln G]
        

        A   = 0;
        B   = 0;

        lse = 0;
        pom = 0;

        XML = []; %cell(N,1);
        lL = []; %cell(N,1);
        list_R = [];
        
        for n=1:N
            
            for j=1:J
                for k=1:K
                    XML{n}(k,j) = 4*c(:,k,n);
                    XML{n}(k,j) = XML{n}(k,j) + nu(:,k,n)*(Xq{n}(:,j)- mu(:,k,n))'*W(:,:,k,n)*(Xq{n}(:,j)-mu(:,k,n));
                    lL{n}  = Elogdet(W(:,:,k,n),nu(:,:,n));
                end
            end
            %disp(sum(lP'*X{n} + lG + lL{n}'./2 - XML{n}./2-2*log(2*pi),2))


            [R,lse_n] = softmax(lP'*X{n} + lG + lL{n}'./2 - XML{n}./2-2*log(2*pi)); % R = E[Z]
            A   = A + X{n}*R';
            B   = B + R;
            lse = lse + sum(lse_n);
            list_R{n} = R;
            pom = pom + sum((lL{n}'./2 - XML{n}./2-2*log(2*pi))*R');
        end

        elbo71 = 0;
        elbo74 = 0;
        elbo77 = 0;

        for n=1:N

            R_sum{n} = sum(list_R{n}',1)+eps;
            %disp(R_sum{n})

            x_mean{n} = zeros(4,K);
            S_mean{n} = zeros(4,4,k);

            for k=1:K

                x_mean{n}(:,k) = Xq{n}*list_R{n}(k,:)'/(R_sum{n}(:,k));
               
                for j=1:J
                    S_mean{n}(:,:,k) = S_mean{n}(:,:,k) + list_R{n}(k,j)*(Xq{n}(:,j)-x_mean{n}(:,k))*(Xq{n}(:,j)-x_mean{n}(:,k))'/(R_sum{n}(:,k));
                end

                %E[ln(p(xq|z,mu,L))]
                elbo71 = elbo71 + R_sum{n}(:,k)*(lL{n}(:,k)-4*c(:,k,n)^(-1)-nu(:,k,n)*trace(S_mean{n}(:,:,k)*W(:,:,k,n))-nu(:,k,n)*(x_mean{n}(:,k)-mu(:,k,n))'*W(:,:,k,n)*(x_mean{n}(:,k)-mu(:,k,n))-4*log(2*pi))/2;

                %E[ln(p(mu,L))]
                elbo74 = elbo74 + (4*log(c0/(2*pi)) + lL{n}(:,k) - 4*c0/c(:,k,n) - c0*nu(:,k,n)*(mu(:,k,n)-mu0)'*W(:,:,k,n)*(mu(:,k,n)-mu0) + (nu0-5)*lL{n}(:,k)-nu(:,k,n)*trace(W0^(-1)*W(:,:,k,n)))/2;

                %E[ln(q(mu,L))]
                elbo77 = elbo77 + lL{n}(:,k)/2 + 2*log(c(:,k,n)/(2*pi)) - 2 - H_Wishart(W(:,:,k,n), nu(:,k,n));

            end
        end

        ELBO  = elbo0 + sum(lse) + sum(pom) + ...
        sum(sum((alpha0 - Alpha).*lP, 1) - C(Alpha) + C(alpha0), 2) +...
        sum(sum((beta0  - Beta ).*lG, 1) - C(Beta) + C(beta0), 2) +...
        elbo71+elbo74-elbo77-elbo_b;

%         disp(elbo0 + sum(lse))
%         disp(sum(sum((alpha0 - Alpha).*lP, 1) - C(Alpha) + C(alpha0), 2))
%         disp(elbo71)
%         disp(elbo74)
%         disp(elbo77)
%         disp(elbo_b)

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

function lp = Bf(alpha, beta)
    lp = det(alpha)^(-beta/2);
    lp = lp*(2^(2*beta)*pi^(3)*gamma(beta/2)*gamma((beta-1)/2)*gamma((beta-2)/2)*gamma((beta-3)/2))^(-1);
end

function lp = H(alpha,beta,delta)
    lp = -log(Bf(beta,delta));
    lp =  lp - (delta-5)/2*alpha+2*delta;
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

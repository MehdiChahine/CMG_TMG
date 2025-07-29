function X = CMG_TMG_STC(Q,mu,R,c,Tc,Tbin)
    
    % Get size variables
    T = Tbin + Tc;
    [M,N] = size(R);

    % Mean shift
    c = c -R*mu;

    % Compute A and b
    A = R';
    b = c;

    % Intialize scale parameter s
    s = 1/N;
    
    % Initialize variables
    W = zeros(M,1);
    X = zeros(N,T);
    STOP = 0; ACC = 0; t = 1; 
    ind = [];
    NPID = 50;ACCPID = 0;

       
    
    %% Gibbs Sampler Initialization
    X_ = (1/N)*ones(N,1)-1/N^2;
    X(:,1) = (X_ - mu);
    U = A'*X(:,1) - b;
    if any(U <= 0)
        warning('Wrong initialization!!')
    end

    for m = 1:M
        W(m) = gigrnd_ONE(1/2,1/s,abs(U(m)^2)/s);        
    end
    
    tic;
    %% Gibbs Sampler Loop
    while ~STOP
        t = t + 1;
        for m = 1:M
%             W(m,t) = gigrnd_ONE(1/2,1/s,abs(U(m)^2)/s);
            W(m) = 1/inverseGaussian(1/abs(U(m)), 1/s);
        end
    
        nu = W(:) + b;
    	Omega = 1./(W(:));
    
        F = chol(Q + A*(Omega.*A')/s);
    
    	X(:,t) = F\(randn(N,1) + F'\A*(Omega.*nu)/s);

        U = A'*X(:,t) - b;
        if t > Tbin
            if ~any(U <= 0) 
                ind = [ind;t];
                ACC = ACC + 1;
                if ACC > Tc
                    STOP = 1;
                end
            end
        else
            if ~any(U <= 0)
                ACCPID = ACCPID + 1;
            end
            
            if mod(t,NPID) == 0
                rho_1 = zeros(N,1);
                IDX = t-NPID+1:t;
                for nrho = 1:N
                    rho_1(nrho) = acf_lag1(X(nrho,IDX)');
                end
                ESSR_1 = min((1-rho_1)./(1+rho_1));
                ACCRATE = ACCPID/NPID;

                e_pid = ACCRATE - ESSR_1;

                Kp = (Tbin-t)^2/Tbin^2 * 1/N;

                s = max(1e-8,s + Kp * e_pid);
                ACCPID = 0;
            end
        end
        
        if ~mod(t,5000)
            fprintf('Iter N : %d || %d \n',t,ACC);
        end
    end
    % Acceptation
    X = X(:,ind);
    
    % Transormation of X
    X = X + mu;




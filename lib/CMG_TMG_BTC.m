function X = CMG_TMG_BTC(Q,mu,R,c,d,Tc,Tbin)
%CMG_TMG_STC generates samples distributed according to a Truncated
%Multivariate Gaussian (TMG) restricted by linear inqualities (c <= Rx <= d) as defined in eq (1) 
%of the paper "A Continuous Gaussian Mixture Approach to Sample
%Multivariate Gaussians constrained by linear inequalities"


%% The inputs are as follows:
% --- Q    NxN precision matrix,
% --- mu   Nx1 mean vector,
% --- R    MxN constraints matrix,
% --- c    Mx1 constraints vector,
% --- Tbin Number of burn-in samples
% --- Tc   Number of required samples

%% the output is:
% --- X an NxTc matrix of samples distributed according to the target TMG


%% Main Code
    % Get size variables
    T = Tbin + Tc;
    [M,N] = size(R);
    
    % Mean shift, get c tilde : 
    c = c -R*mu;
    d = d - R*mu;


    % Compute A and b
    Delta = diag(1./(d-c));

    A = (Delta*R)';
    b = Delta*c;

    % Intialize scale parameter s
    s = 1/N;

    % Initialize variables
    W = zeros(M,1);
    X = zeros(N,T);
    STOP = 0; ACC = 0; t = 1; 
    ind = [];
    NPID = 50;ACCPID = 0;


   
    %% Gibbs sampler
    % Initialization
    X_ = (1/N)*ones(N,1)-1/N^2;
    X(:,1) = (X_ - mu);
    U = A'*X(:,1) - b;
    if any(abs(U-1/2) >= 1/2)
        warning('Wrong initialization!!')
    end

    for m = 1:M
        W(m) = p_W(U(m),sqrt(s));
    end
    
    %% Loop
    while ~STOP
        t = t + 1;
%         tStart = tic;
        for m = 1:M
            W(m) = p_W(U(m),sqrt(s));
        end
    
        nu = W + b;
        Omega = (1./(W.*(1 - W)));
        
        F = chol(Q + A*(Omega.*A')/s);
    	X(:,t) = F\(randn(N,1) + F'\A*(Omega.*nu)/s);
	
        U = A'*X(:,t) - b;
        if t > Tbin
            if ~any(abs(U-1/2) >= 1/2) 
                ind = [ind;t];
                ACC = ACC + 1;
                if ACC > Tc
                    STOP = 1;
                end
            end
        else
            if ~any(abs(U-1/2) >= 1/2) 
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

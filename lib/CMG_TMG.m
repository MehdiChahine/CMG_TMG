function X = CMG_TMG(Q,mu,R0,c0,R1,c1,d1,Tc,Tbin)
%CMG_TMG_STC generates samples distributed according to a Truncated
%Multivariate Gaussian (TMG) restricted by linear inqualities (c <= Rx) as defined in eq (1)
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
    [M0,N] = size(R0);

    [M1,N] = size(R1);

    % Mean shift
    c0 = c0 -R0*mu; 
    c1 = c1 -R1*mu; 
    d1 = d1 -R1*mu; 
    
    Delta = diag(1./(d1-c1));
    
    DR1 = (Delta*R1);
    Dc1 = Delta*c1;

    % Intialize scale parameter s
    s = 1/N;
    
    % Initialize variables
    W = zeros(M0,1);
    V = zeros(M1,1);
    
    X = zeros(N,T);
    STOP = 0; ACC = 0; t = 1; 
    ind = [];
    NPID = 50;ACCPID = 0;

    %% Gibbs Sampler Initialization
    % Initialization
    X_ = zeros(N,1);
    X(:,1) = (X_ - mu);
    Uw = R0*X(:,t) - c0;
    Uv = DR1*X(:,t) - Dc1;
    U = [Uw;Uv];
    if any(U <= 0)
        warning('Wrong initialization!!')
    end
    
    for m = 1:M0
%         W(m) = gigrnd_ONE(1/2,1/s^2,abs(Uw(m)^2)/s^2);
        W(m) = 1/inverseGaussian(1/abs(Uw(m)), 1/s);
    end
    
    for m = 1:M1
        V(m,t) = p_W(Uv(m),s);
    end
    %% Gibbs Sampler Loop
    while ~STOP
        t = t + 1;
        for m = 1:M0
%             W(m,t) = gigrnd_ONE(1/2,1/s,abs(U(m)^2)/s);
            W(m) = 1/inverseGaussian(1/abs(Uw(m)), 1/s);
        end
    
        nuW = W + c0;
    	OmegaW = 1./(W);

        for m = 1:M1
            V(m) = p_W(Uv(m),sqrt(s));
        end

        nuV = V + Dc1;
        OmegaV = (1./(V.*(1 - V)));


        F = chol(Q + R0'*(OmegaW.*R0/s) + DR1'*(OmegaV.*DR1)/s);
        X(:,t) = F\randn(N,1) + F\(F'\(R0'*(OmegaW.*nuW)/s + DR1'*(OmegaV.*nuV)/s));
    
        Uw = R0*X(:,t) - c0;
        Uv = DR1*X(:,t) - Dc1;
        
        if t > Tbin
            if ~any(abs(Uv-1/2) >= 1/2) && ~any(Uw <= 0)     
                ind = [ind;t];
                ACC = ACC + 1;
                if ACC > Tc
                    STOP = 1;
                end
            end
        else
            if ~any(abs(Uv-1/2) >= 1/2) && ~any(Uw <= 0)     
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




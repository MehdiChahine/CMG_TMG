function [Sample,ACC] = GLSWRAND_DEV(x_0,s)

% Case 1 : x_0 in [0,1]
if 0 <= x_0 && x_0 <= 1
    %Special cases x_ = {0,1}
    if x_0 == 1 || x_0 == 0
        EX = (-1)^(1-x_0);
        x_0 = 1;
        %Recurent computations
        x_1 = 2*x_0 - 1;
        % Mode
        M = asin(x_1);

        % compute a and b (-s and t in the paper by Devroye)
        a = -M + 2*atan((1-sqrt(2*s^2-4*x_0^2+4*x_0))/(sqrt(2)*s+x_1));
         % Algorithm
        %nu = 1;% 1/s^2 * ((x_1)*sec(b+M) - tan(b+M))^2; % -psi(-t)
        %eta = -2*(sin(b+M)-(x_1))*((x_1)*sin(b+M)-1)/(s^2*cos(b+M)^3); % -dpsi(t)

        theta = 1;%1/s^2 * ((x_1)*sec(a+M) - tan(a+M))^2; % -psi(-s)
        xi = 2*(sin(a+M)-(x_1))*((x_1)*sin(a+M)-1)/(s^2*cos(a+M)^3); % -dpsi(t)

        p = 1/xi;
        %r = 1/eta;

        %tpr = b - r*nu;
        spr = -a - p*theta;
        q = pi/2 + spr;

        ACC = 0;
        while 1
            ACC = ACC + 1;
            U = rand;
            V = rand;
            if U < q/(p+q)
                X = -spr + q*V;
            else
                X = -spr + p*log(V);
            end


            if -spr <= X && X <= pi/2
                XI = 1;
            else
                XI = exp(-theta + xi*(X-a));
            end
            %XI = (-spr <= X && X <= tpr ) + ( X > tpr ) * exp(-nu - eta*(X-b)) + (X < -spr) * exp(-theta + xi*(X-a));
            if rand*XI <= exp(-1/(2*s^2) * ((x_1)*sec(X+M) - tan(X+M))^2) && (abs(X+M) <= pi/2)
                Sample = EX * (X+M);
                Sample = 1/2*(sin(Sample)+1);
                break;
            end
        end
    else
        %Recurent computations
        x_1 = 2*x_0 - 1;
        % Mode
        M = asin(x_1);

        % compute a and b (-s and t in the paper by Devroye)
        a = -M + 2*atan((1-sqrt(2*s^2-4*x_0^2+4*x_0))/(sqrt(2)*s+x_1));
        b = -M + 2*atan((sqrt(2*s^2-4*x_0^2+4*x_0)-1)/(sqrt(2)*s-x_1));

        % Algorithm
        nu = 1;% 1/s^2 * ((x_1)*sec(b+M) - tan(b+M))^2; % -psi(-t)
        eta = -(sin(b+M)-(x_1))*((x_1)*sin(b+M)-1)/(s^2*cos(b+M)^3); % -dpsi(t)

        theta = 1;%1/s^2 * ((x_1)*sec(a+M) - tan(a+M))^2; % -psi(-s)
        xi = (sin(a+M)-(x_1))*((x_1)*sin(a+M)-1)/(s^2*cos(a+M)^3); % -dpsi(t)

        p = 1/xi;
        r = 1/eta;

        tpr = b - r*nu;
        spr = -a - p*theta;
        q = tpr + spr;

        ACC = 0;
        while 1
            ACC = ACC + 1;
            U = rand;
            V = rand;
            if U < q/(p+q+r)
                X = -spr + q*V;
            elseif U <  (q+r)/(p+q+r)
                X = tpr - r*log(V);
            else
                X = -spr + p*log(V);
            end


            if -spr <= X && X <= tpr
                XI = 1;
            elseif tpr < X 
                XI = exp(-nu - eta*(X-b));
            elseif  X < -spr
                XI = exp(-theta + xi*(X-a));
            end
            %XI = (-spr <= X && X <= tpr ) + ( X > tpr ) * exp(-nu - eta*(X-b)) + (X < -spr) * exp(-theta + xi*(X-a));
            if rand*XI <= exp(-1/(2*s^2) * ((x_1)*sec(X+M) - tan(X+M))^2) && (abs(X+M) <= pi/2)
                Sample = X+M;
                Sample = 1/2*(sin(Sample)+1);
                break;
            end
        end
    end
elseif x_0 < 0% Case 2 : x_0 < 0
    %Recurent computations
    x_1 = 2*x_0 - 1;
    
    % Mode
    M = asin(1/x_1);
    
    % compute a and b (-s and t in the paper by Devroye)
    z = sqrt(2*s^2+((x_1)*sec(M) - tan(M))^2);
    a = -M + 2*atan((1+sqrt(z^2-4*x_0^2+4*x_0))/(-z+x_1));   
    b = -M + 2*atan((sqrt(z^2-4*x_0^2+4*x_0)-1)/(z-x_1));
    
    % Algorithm
    nu = 1;% 1/s^2 * ((x_1)*sec(b+M) - tan(b+M))^2; % -psi(-t)
    eta = -(sin(b+M)-(x_1))*((x_1)*sin(b+M)-1)/(s^2*cos(b+M)^3); % -dpsi(t)
    
    theta = 1;%1/s^2 * ((x_1)*sec(a+M) - tan(a+M))^2; % -psi(-s)
    xi = (sin(a+M)-(x_1))*((x_1)*sin(a+M)-1)/(s^2*cos(a+M)^3); % -dpsi(t)
    
    p = 1/xi;
    r = 1/eta;
    
    tpr = b - r*nu;
    spr = -a - p*theta;
    q = tpr + spr;
    
    ACC = 0;
    while 1
        ACC = ACC + 1;
        U = rand;
        V = rand;
        if U < q/(p+q+r)
            X = -spr + q*V;
        elseif U <  (q+r)/(p+q+r)
            X = tpr - r*log(V);
        else
            X = -spr + p*log(V);
        end
        
        
        if -spr <= X && X <= tpr
            XI = 1;
        elseif tpr < X 
            XI = exp(-nu - eta*(X-b));
        elseif  X < -spr
            XI = exp(-theta + xi*(X-a));
        end
        %XI = (-spr <= X && X <= tpr ) + ( X > tpr ) * exp(-nu - eta*(X-b)) + (X < -spr) * exp(-theta + xi*(X-a));
        if rand*XI <= exp(-1/(2*s^2) * ((x_1)*sec(X+M) - tan(X+M))^2 + 1/(2*s^2) * ((x_1)*sec(M) - tan(M))^2) && (abs(X+M) <= pi/2)
            Sample = X+M;
            Sample = 1/2*(sin(Sample)+1);
            break;
        end
    end
else %Case 3: x_0 > 1
    %Recurent computations
    x_1 = 2*x_0 - 1;
    
    % Mode
    M = asin(1/x_1);
    
    % compute a and b (-s and t in the paper by Devroye)
    z = -sqrt(2*s^2+((x_1)*sec(M) - tan(M))^2);
    a = -M + 2*atan((1+sqrt(z^2-4*x_0^2+4*x_0))/(-z+x_1));   
    b = -M + 2*atan((sqrt(z^2-4*x_0^2+4*x_0)-1)/(z-x_1));
    
    % Algorithm
    nu = 1;% 1/s^2 * ((x_1)*sec(b+M) - tan(b+M))^2; % -psi(-t)
    eta = -(sin(b+M)-(x_1))*((x_1)*sin(b+M)-1)/(s^2*cos(b+M)^3); % -dpsi(t)
    
    theta = 1;%1/s^2 * ((x_1)*sec(a+M) - tan(a+M))^2; % -psi(-s)
    xi = (sin(a+M)-(x_1))*((x_1)*sin(a+M)-1)/(s^2*cos(a+M)^3); % -dpsi(t)
    
    p = 1/xi;
    r = 1/eta;
    
    tpr = b - r*nu;
    spr = -a - p*theta;
    q = tpr + spr;
    
    ACC = 0;
    while 1
        ACC = ACC + 1;
        U = rand;
        V = rand;
        if U < q/(p+q+r)
            X = -spr + q*V;
        elseif U <  (q+r)/(p+q+r)
            X = tpr - r*log(V);
        else
            X = -spr + p*log(V);
        end
        
        
        if -spr <= X && X <= tpr
            XI = 1;
        elseif tpr < X 
            XI = exp(-nu - eta*(X-b));
        elseif  X < -spr
            XI = exp(-theta + xi*(X-a));
        end
        %XI = (-spr <= X && X <= tpr ) + ( X > tpr ) * exp(-nu - eta*(X-b)) + (X < -spr) * exp(-theta + xi*(X-a));
        if rand*XI <= exp(-1/(2*s^2) * ((x_1)*sec(X+M) - tan(X+M))^2 + 1/(2*s^2) * ((x_1)*sec(M) - tan(M))^2) && (abs(X+M) <= pi/2)
            Sample = X+M;
            Sample = 1/2*(sin(Sample)+1);
            break;
        end
    end
end
end

% Case 1 :
% psi = -1/(2*s^2) * ((x_1)*sec(x+M) - tan(x+M))^2;
% dpqi = (sin(x+M)-(x_1))*((x_1)*sin(x+M)-1)/(s^2*cos(x+M)^3);

% Case 2 :
% psi = -1/(2*s^2) * ((x_1)*sec(x+M) - tan(x+M))^2 + 1/s^2 * ((x_1)*sec(M) - tan(M))^2;
% dpqi = (sin(x+M)-(x_1))*((x_1)*sin(x+M)-1)/(s^2*cos(x+M)^3);

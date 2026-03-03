function D_end = runSEIRD(theta, T)
% runSEIRD  SEIRD simulator returning cumulative deaths at time T.
%  theta length 6  → basic model
%  theta length 9  → intervention (β reduced by k after t_int)

if nargin<2, T = 160; end
beta  = theta(1); sigma = theta(2);
gamma = theta(3); mu    = theta(4);
N     = theta(5); I0    = theta(6);

% Possible intervention
hasInt = numel(theta)==9;
if hasInt
    beta2  = theta(7); t_int = theta(8); k = theta(9);  %#ok<NASGU>
end

% Initial state
S0 = N - 2*I0;  E0 = I0;  R0 = 0; D0 = 0;
y0 = [S0;E0;I0;R0;D0];

rhs = @(t,y) seirdRHS(t,y,theta);  % nested RHS function
[tout, yout] = ode45(rhs, [0 T], y0, odeset('RelTol',1e-8,'AbsTol',1e-10));
D_end = yout(end,5);  % deaths at T

    function dy = seirdRHS(t,y,theta)
        b = beta;  % default
        if hasInt && t>=theta(8)  % t_int
            b = theta(9)*beta;  % k*beta
        end
        S=y(1);E=y(2);I=y(3);R=y(4);D=y(5);  %#ok<NASGU>
        dS = -b*I/N*S;
        dE =  b*I/N*S - sigma*E;
        dI =  sigma*E - (gamma+mu)*I;
        dR =  gamma*I;
        dD =  mu*I;
        dy = [dS;dE;dI;dR;dD];
    end
end

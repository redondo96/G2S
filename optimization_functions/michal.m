function [y] = michal(xx, m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MICHALEWICZ FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% xx = [x1, x2]
% m = constant (optional), with default value 10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 1)
    m = 10;
end

d = length(xx);
sum = 0;

for ii = 1:d
	xi = xx(ii);
	new = sin(xi) * (sin(ii*xi^2/pi))^(2*m);
	sum  = sum + new;
end

y = -sum;


% GRADIENT
% 
% evalPts = xx;    
% 
% N = size(evalPts,1);
% gradMag = zeros(N,1);
% 
% iIdx = (1:size(evalPts,2)).';
% two_m = 2*m;
% for k = 1:N
%     x = evalPts(k,:).';
%     a = sin(x);
%     b = sin(iIdx .* x.^2 / pi);
%     b2m    = b.^two_m;
%     dbdx   = cos(iIdx .* x.^2 / pi) .* (2*iIdx .* x / pi);
%     dti_dx = cos(x) .* b2m + a .* two_m .* b.^(two_m-1) .* dbdx;
%     g      = -dti_dx;
%     gradMag(k) = norm(g);
% end


% d = numel(xx);
% g = zeros(d,1);
% 
% for ii = 1:d
%     xi   = xx(ii);
% 
%     Si   = sin(ii * xi^2 / pi);  % S_i
%     Ci   = cos(ii * xi^2 / pi);  % cos(arg) for derivative
% 
%     term1 =  cos(xi) * Si^(2*m);  %  cos(xi) · S_i^{2m}
%     term2 =  sin(xi) * (2*m) * Si^(2*m-1) ...
%              * Ci * (2*ii*xi/pi);  % sin(xi) * (2m) * S^{2m-1} * cos(arg) * d(arg)/dx
% 
%     g = -( term1 + term2 );  % negative sign from f definition
% end


% [N, d] = size(xx);
% two_m  = 2*m;
% 
% % Pre-allocate gradient matrix: G(k,i)
% G = zeros(N, d);
% 
% for i = 1:d
%     xi     = xx(:, i);                       % N×1
%     a      = sin(xi);                        % sin(x_i)
%     arg    = i .* xi.^2 ./ pi;               % i*xi^2/π
%     b      = sin(arg);                       % S_i
%     b2m    = b.^two_m;                       % S_i^{2m}
%     dbdx   = cos(arg) .* (2*i .* xi ./ pi);  % dS_i/dx_i
% 
%     % Gradient in x_i for all points (N×1)
%     G(:, i) = -( cos(xi).*b2m ...
%                + a .* two_m .* b.^(two_m-1) .* dbdx );
% end
% 
% % Gradient magnitude (N×1)
% gradMag = sqrt( sum(G.^2, 2) );

end

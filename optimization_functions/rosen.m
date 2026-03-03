function [y] = rosen(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROSENBROCK FUNCTION
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
% INPUT:
%
% xx = [x1, x2, ..., xd]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = length(xx);
sum = 0;
for ii = 1:(d-1)
	xi = xx(ii);
	xnext = xx(ii+1);
	new = 100*(xnext-xi^2)^2 + (xi-1)^2;
	sum = sum + new;
end

y = sum;



% % GRADIENT
% 
% evalPts = xx;    
% 
% N = size(evalPts,1);
% gradMag = zeros(N,1);
% 
% for k = 1:N
%     x = evalPts(k,:).';
%     d = numel(x);
%     g = zeros(d,1);
% 
%     g(1) = -400*x(1)*(x(2)-x(1)^2) + 2*(x(1)-1);
%     g(d) =  200*(x(d) - x(d-1)^2);
%     for j = 2:d-1
%         g(j) = 200*(x(j)-x(j-1)^2) ...
%               -400*x(j)*(x(j+1)-x(j)^2) ...
%               +2*(x(j)-1);
%     end
%     gradMag(k) = norm(g);
% end

end

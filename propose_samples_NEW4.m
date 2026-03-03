function [Xnew, Xnew_val] = propose_samples_NEW4(samples_tri, ...
    values,new_budget,discret_spl)

% Proposal for 2-D / 3-D subspaces
% Inputs:
%   samples_tri  (N×d) : existing coordinates (d = 2 or 3)
%   values       (N×1) : associated values  (higher -> more interesting)
%   new_budget         : how many new points to return
%   discret_spl        : discretized values of the subspace dimensions
%
% Outputs:
%   Xnew      q×(d+1) : newly proposed coordinates
%   Xnew_val


% Mixing coefficient 0<T≤1 (Temperature)
Tmix = 1;
% Kernel width as fraction of box diagonal
sigmaFac = 0.02;

% Proportion of box to drop
edgeCut = 0.025;  % Default: outer 2.5% per dimension

[N,d] = size(samples_tri);
assert(d==2 || d==3, 'Sampler designed for 2- or 3-D spaces.');

% Lower and upper bounds as 1×d row-vectors
lb = cellfun(@(v) v(1) , discret_spl)';   % first element of each grid
ub = cellfun(@(v) v(end), discret_spl)';  % last  element of each grid


if size(values,2)>=2
    % GRADIENT norm
    values = vecnorm(values, 2, 2);  % vectorized norm computation
end


%% 1. Normalize coordinates to [0,1] box
Xn = (samples_tri - lb) ./ (ub - lb);
boxDiag = sqrt(d);  % diag length of [0,1]^d
sigma = sigmaFac * boxDiag;  % isotropic Gaussian width


%% 2. Rescale 'values' to positive weights, then compress with k-means
w  = values - min(values);  % shift to non-negative
if all(w==0), w = ones(size(w)); end
w  = w / sum(w);  % normalize

% Filter out edge points
maskCore = all( Xn > edgeCut & Xn < (1-edgeCut), 2 );  % N×1 logical

% If the filter removes every point, return to all points
if ~any(maskCore)
    maskCore = true(size(Xn,1),1);
end

Xcore = Xn(maskCore,:);      % coordinates without border points
wcore = w(maskCore);         % matching weights
wcore = wcore / sum(wcore);  % renormalize to sum=1

% Weighted subsampling to get M centres
nCenters = chooseNumCentres(wcore, new_budget);

M = min(nCenters, N);  % cannot exceed N
idxC = randsample(size(Xcore,1), M, true, wcore);  % with replacement
centres = Xcore(idxC,:);  % M×d (in unit cube)
wC = wcore(idxC);  % weights
wC = wC / sum(wC);  % renormalize


%% 3. Draw 'new_budget' points from the mixture (Tmix⋅KDE + (1-Tmix)⋅Uniform)
U  = rand(new_budget,1);

Xnew = zeros(new_budget,d);
for i = 1:new_budget
    if U(i) < Tmix
        % Draw from KDE component
        k     = randsample(numel(wC),1,true,wC);
        mu    = centres(k,:);
        Xcand = mu + sigma*randn(1,d);
        Xcand = max(min(Xcand,1),0);  % keep inside [0,1]
    else
        % Draw uniform exploration sample
        Xcand = rand(1,d);
    end
    Xnew(i,:) = lb + Xcand.*(ub-lb);
end


%% 4. Interpolate values
switch d
    case 2
        F = scatteredInterpolant(samples_tri(:,1), ...
                                 samples_tri(:,2), ...
                                 values, 'natural', 'nearest');
        
        vhat = F(Xnew(:,1), Xnew(:,2));
    
    case 3
        F = scatteredInterpolant(samples_tri(:,1), ...
                                 samples_tri(:,2), ...
                                 samples_tri(:,3), ...
                                 values, 'natural', 'nearest');

        vhat = F(Xnew(:,1), Xnew(:,2), Xnew(:,3));

    otherwise
        error('Sampler designed for 2-D or 3-D subspaces.');
end

Xnew_val = [Xnew, vhat(:)];


end

function [samples_pos_tri, weighted_means_pos] = ...
         discretiseCoords(samples_tri, weighted_means, discret_spl, dim_spl)

%   Snap each row of samples_tri to an integer grid defined by dim_spl
%   and average weighted_means when multiple rows land on the same bin.
%
% Inputs:
%   samples_tri      N×D  continuous coordinates
%   weighted_means   N×1  values attached to every row of samples_tri
%   discret_spl      1×D  cell; discret_spl{d}(1) = min_d, end = max_d
%   dim_spl          1×D  number of discrete indices per dimension
%
% Outputs:
%   samples_pos_tri      M×D  integer indices (1...dim_spl(d))
%   weighted_means_pos   M×1  averaged values for each unique row


[N, D] = size(samples_tri);
assert(numel(discret_spl) == D, 'discret_spl must have D cells');
assert(numel(dim_spl) == D, 'dim_spl must have length D');
assert(size(weighted_means,1) == N, 'weighted_means must match rows');


%% 1. Map every coordinate to its integer bin index
idxMat = zeros(N, D);

for d = 1:D
    lo = discret_spl{d}(1);    % lower bound
    hi = discret_spl{d}(end);  % upper bound
    M = dim_spl(d);            % number of bins
    step = (hi-lo) / (M-1);    % uniform spacing

    % Nearest index in 1...M
    idx = round((samples_tri(:,d) - lo) / step) + 1;
    idx = min(max(idx,1), M);  % clamp
    idxMat(:,d) = idx;
end


%% 2. Collapse duplicates & average values
[samples_pos_tri, ~, ic] = unique(idxMat, 'rows');   % ic: group id (N×1)
weighted_means_pos = accumarray(ic, weighted_means, [], @mean);

end

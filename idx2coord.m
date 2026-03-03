function samples_tri = ...
         idx2coord(samples_pos_tri, discret_spl, dim_spl)

% Convert the integer grid indices in samples_pos_tri back to their
% coordinates. The row order is preserved, so the output can stay aligned
% with weighted_means_pos (or any other row-wise data).
%
% Inputs:
%   samples_pos_tri   M×D   integer indices
%   discret_spl       1×D   cell; discret_spl{d}(1)=lb , end=ub
%   dim_spl           1×D   number of grid nodes per dimension
%
% Output:
%   samples_tri       M×D   coordinates in the original domain


[M, D] = size(samples_pos_tri);
assert(numel(discret_spl) ==D, 'discret_spl must have D cells');
assert(numel(dim_spl) ==D, 'dim_spl length must be D');

samples_tri = zeros(M, D);

for d = 1:D
    lo = discret_spl{d}(1);     % lower bound of dim d
    hi = discret_spl{d}(end);   % upper bound
    Mbin = dim_spl(d);          % number of grid points (1...Mbin)
    step = (hi-lo) / (Mbin-1);  % uniform spacing

    idx = samples_pos_tri(:,d);  % integer indices (1...Mbin)
    samples_tri(:,d) = lo + (idx - 1) * step;  % map back to coordinate
end

end

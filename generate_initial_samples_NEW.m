function samples_comp_orig = generate_initial_samples_NEW(splits, dim, ...
    initial_budget, discret, opt_function)

% Inputs :  splits – 1×nSub cell, splits{s} = column indices of subspace s
%           dim
%           initial_budget
%           discret – 1×D bounds
% 
% Output:   samples_comp_orig – initial_budget × D matrix


% Multiplier for Sobol candidate pool
overshoot = 8;
% Initial radius as a fraction of box diagonal
minDist0  = 0.15;

lb = cellfun(@(v) v(1) , discret)';  % first element of each grid
ub = cellfun(@(v) v(end), discret)';

D    = numel(dim);     % global dimensionality
nSub = numel(splits);  % number of subspaces


% Enumerate all subspace corners (each is a full-D point)

mid = (lb + ub) / 2;  % filler for other coords (mid-value)

cornerList = [];
for s = 1:nSub
    J = splits{s};                     % indices of subspace s
    ds = numel(J);                     % 2 or 3
    combos = dec2bin(0:2^ds-1) - '0';  % 0→lower 1→upper
    for c = 1:size(combos,1)
        p = mid;                         % start with mid-values
        choose = combos(c,:)==0;
        p(J(choose)) = lb(J(choose));    % set lower bounds
        p(J(~choose)) = ub(J(~choose));  % set upper bounds
        cornerList = [cornerList; p];    %#ok<AGROW>
    end
end
cornerList       = unique(cornerList, 'rows');  % just in case of overlap
Ncorner          = size(cornerList,1);

if Ncorner > initial_budget
    error('initial_budget (%d) smaller than #required corners (%d).', ...
           initial_budget, Ncorner);
end

% Sobol candidate pool
remain  = initial_budget - Ncorner;
sob     = scramble(sobolset(D,'Skip',1000,'Leap',100), ...
    'MatousekAffineOwen');
Npool   = overshoot * remain;
pool    = net(sob, Npool);
pool    = bsxfun(@plus, lb, bsxfun(@times, pool, (ub-lb)));

% Maximin radius
boxDiag = norm(ub-lb);
radius  = minDist0 * boxDiag;

% Hash of already-taken points (string key)
row2key = @(q) sprintf('%.12g_', q);
taken   = containers.Map('KeyType','char','ValueType','logical');
for k = 1:Ncorner
    taken(row2key(cornerList(k,:))) = true;
end

accepted = cornerList;  % start with corners
keep = Ncorner;
idx = 1;

% Adaptively decreasing radius
while keep < initial_budget && idx <= Npool
    x = pool(idx,:);   idx = idx+1;

    % Distance test in every subspace
    ok = true;
    for s = 1:nSub
        J    = splits{s};
        if size(accepted,1)>0
            [~,dist] = knnsearch(accepted(:,J), x(J), 'K',1);
            if dist < radius
                ok = false; break
            end
        end
    end

    if ok
        keep = keep + 1;
        accepted(keep,:) = x;
        taken(row2key(x)) = true;
    end

    % If the pool is exhausted early, shrink radius and extend pool
    if idx > Npool && keep < initial_budget
        radius = radius * 0.7;
        extra  = net(sob, overshoot*remain);
        extra  = bsxfun(@plus, lb, bsxfun(@times, extra, (ub-lb)));
        pool   = [pool; extra];  %#ok<AGROW>
        Npool  = size(pool,1);
    end
end
Xinit = accepted(1:keep,:);
if keep < initial_budget
    warning('Only %d/%d samples produced; increase ''overshoot''.', ...
             keep, initial_budget);
end


% Evaluate function for all samples (vectorized)
s = arrayfun(@(idx) opt_function(Xinit(idx, :)), 1:keep).';


samples_comp_orig = [Xinit, s];

end

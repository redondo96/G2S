function M = chooseNumCentres(w, new_budget)

% Adaptive number of KDE centres
%
%  Inputs:
%    w          : N×1 positive weights, already normalised (sum=1)
%    new_budget : batch size we intend to propose
%
%  Output:
%    M : number of KDE centres


% Coverage: cumulative weight threshold
covT  = 0.9;
% Target q/M ratio
alpha = 2;

wSorted = sort(w,'descend');
kCov = find(cumsum(wSorted) >= covT, 1,'first');

if isempty(kCov), kCov = numel(w); end  % extremely flat weights

% combine both criteria
M = min( ceil( new_budget / alpha ), kCov );
end

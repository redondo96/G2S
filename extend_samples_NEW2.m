function [sample_comp, taken] = extend_samples_NEW2(dim,sample,splits, ...
    taken,discret,i,store,st)

% Extend a point in one subspace) to the full D–dimensional space.
%
% Inputs:
%   dim             full global dimensionality D
%   sample          1×d_sub     coordinates in the current subspace
%   splits          1×nSub cell, splits{s} = indices of subspace s
%   taken           samples already taken
%   discret         1×D   cell, discret{d}(1)=lb, end=ub
%   i               subspace id
%   store           1×nSub cell, each cell holds .samples{st}  (M×(d+1))
%   st              strategy index
%
% Outputs:
%   sample_comp     1×D   full-dimensional point (new and unique)
%   taken           samples already taken (updated)


subIdx = i;

if length(dim) <= 3
    % If there are no subspaces, the sample should not be extended
    sample_comp = sample;
else

    row2key = @(p) sprintf('%.12g_', p);

    % Precompute full dimension D
    nSub = numel(splits);
    % D = sum(cellfun(@numel, splits));
    
    lb  = cellfun(@(v) v(1),  discret)';   % 1×D
    ub  = cellfun(@(v) v(end), discret)';  % 1×D
    mid = (lb + ub)/2;

    % Try until we find a NEW full-dimensional point
    maxTrials = 100;
    trial = 0;
    
    while true
        trial = trial + 1;
        if trial > maxTrials
            error('extendSample:CouldNotFindUnique', ...
                  'Failed to build a unique sample after %d trials.', ...
                  maxTrials);
        end
    
        xfull = mid;  % start with mids everywhere
        xfull(splits{subIdx}) = sample;  % insert the known chunk
    
        % Fill every other subspace
        for j = 1:nSub
            if j == subIdx, continue; end
    
            Mj = store{j}.samples{st};  % Mj×(d_j+1) or []

            if isempty(Mj)
                % No information. Keep the mid-point already in xfull
                continue
            else
                vals = abs(Mj(:,end));  % weights (use abs in case of neg.)
                if all(vals==0), vals = ones(size(vals)); end

                % vals = vals / sum(vals);
    
                idxRow = randsample(numel(vals),1,true,vals);  % weighted pick
                chunk  = Mj(idxRow,1:end-1);  % 1×d_j
                xfull(splits{j}) = chunk;  % insert chunk
            end
        end
    
        % Uniqueness check
        key = row2key(xfull);
        if ~isKey(taken,key)
            sample_comp     = xfull;
            return
        end
        % Otherwise retry with a different weighted pick
    end
    
end

end


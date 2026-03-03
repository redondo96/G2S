function [samples_pos_tri_run, ...
    samples_OUTPUT, ...
    taken, ...
    samples_pos_comp_TRACK, ...
    samples_tri_run, ...
    samples_tri_sts, samples_comp_sts, ...
    samples_comp_eval, ...
    samples_pos_err_subesp, err, ...
    samples_pos_sts, ...
    strategy1, strategy1_vals, ...
    strategy2, strategy2_vals, ...
    strategy3, strategy3_vals, ...
    samples_tri_wm, samples_pos_err_wm, ...
    store] = execute_strategy(st, ...
                    ndim_spl, dim_spl, discret_spl, ...
                    new_budget, ...
                    samples_pos_comp_TRACK, ...
                    samples_pos_tri_run, ...
                    i, ...
                    spl, ...
                    dims_left, dim, ...
                    discret, ...
                    opt_function, ...
                    samples_OUTPUT, ...
                    taken, ...
                    samples_tri_run, ...
                    samples_tri_sts, samples_comp_sts, ...
                    s, ...
                    splits, ...
                    models_sts, ...
                    error_function, ...
                    samples_pos_err_subesp, ...
                    samples_pos_sts, ...
                    strategy1, strategy1_vals, ...
                    strategy2, strategy2_vals, ...
                    strategy3, strategy3_vals, ...
                    samples_tri_wm, samples_pos_err_wm, ...
                    store)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

samples_tri = samples_OUTPUT(:,spl);
weighted_means = samples_OUTPUT(:,end);

% 1) Unique rows. ic  tells to which row each original row belongs
[samples_tri, ~, ic] = unique(samples_tri, 'rows');
% 2) Average the values
weighted_means = accumarray(ic, weighted_means, [], @mean);

switch st
    case 1
        % STRATEGY 1: Gradient Estimation

        T = delaunayTriangulation(samples_tri);

        grad_v = estimate_gradient(T,[samples_tri(:,1:ndim_spl) weighted_means]);  % samples_tri(:,end)
        values = grad_v;

        [samples_tri_prop, samples_tri_prop_val] = propose_samples_NEW4(samples_tri, ...
            values,new_budget,discret_spl);

        % strategy1{i,end+1} = [samples_tri(:,1:ndim_spl) vecnorm(grad_v,2,2)];
        % strategy1_vals{i,end+1} = samples_pos_tri_prop;

    case 2
        [samples_pos_tri, weighted_means_pos] = ...
            discretiseCoords(samples_tri, weighted_means, discret_spl, dim_spl);

        % STRATEGY 2: Criticality
        [criticality] = main_sampling_base(dim_spl,samples_pos_tri,weighted_means_pos);  % base for strategies 2 and 3
        % TO USE THE WHOLE TENSOR
        % [~,criticality] = main_sampling_base(dim,samples_pos_comp,samples_comp);
        values = criticality;

        samples_tri = idx2coord(samples_pos_tri, discret_spl, dim_spl);

        [samples_tri_prop, samples_tri_prop_val] = propose_samples_NEW4(samples_tri, ...
            values,new_budget,discret_spl);

        % strategy2{i,end+1} = [samples_tri(:,1:ndim_spl) criticality];
        % strategy2_vals{i,end+1} = samples_pos_tri_prop;

    case 3
        [samples_pos_tri, weighted_means_pos] = ...
            discretiseCoords(samples_tri, weighted_means, discret_spl, dim_spl);
        
        % STRATEGY 3: Gradient of Criticality
        [criticality] = main_sampling_base(dim_spl,samples_pos_tri,weighted_means_pos);  % base for strategies 2 and 3
        
        % (CHECKK)
        % if all(criticality == criticality(1))
        %     grad_cv = ones(numel(criticality), ndim_spl);
        % else

        samples_tri = idx2coord(samples_pos_tri, discret_spl, dim_spl);

        T = delaunayTriangulation(samples_tri);

        grad_cv = estimate_gradient(T, ...
            [samples_tri(:,1:ndim_spl) weighted_means_pos],criticality);
        values = grad_cv;

        [samples_tri_prop, samples_tri_prop_val] = propose_samples_NEW4(samples_tri, ...
            values,new_budget,discret_spl);

        % strategy3{i,end+1} = [samples_tri(:,1:ndim_spl) vecnorm(grad_cv,2,2)];
        % strategy3_vals{i,end+1} = samples_pos_tri_prop;
end

samples_pos_sts{st,i} = [samples_pos_sts{st,i}; samples_tri_prop];  %%%%%%%

% Save information for next RUNS
samples_pos_tri_run{i} = unique([samples_pos_tri_run{i}; samples_tri_prop],'rows','stable');  %%%%%%%

samples_comp_eval = zeros(new_budget, length(spl)+length(dims_left)+1);

for sam=1:size(samples_tri_prop,1)

    sample = samples_tri_prop(sam,:);
    [sample_comp, taken] = extend_samples_NEW2(dim,sample,splits, ...
        taken,discret,i,store,st);

    % Evaluate function
    sample_comp_eval = [sample_comp, opt_function(sample_comp)];
    samples_comp_eval(sam,:) = sample_comp_eval;

    % Each time the function is evaluated, the samples are accumulated
    % in the output variable
    samples_OUTPUT = cat(1,samples_OUTPUT,sample_comp_eval);

end

% We store the samples proposed by the strategy
samples_comp_sts{s,i} = cat(1,samples_comp_sts{s,i},samples_comp_eval);  %%%%%%%

% ERROR of the samples
errs = error_function(models_sts{s},samples_comp_eval);
% And store them together with the positions
samples_err = [samples_comp_eval(:,1:end-1), errs];


% i) Most recent
store{i}.samples{st} = samples_tri_prop_val;

% ii) All
% store{i}.samples{st} = [store{i}.samples{st}; samples_tri_prop_val];

% iii) Mix
% S = store{i}.samples{st};
% N  = size(samples_tri_prop_val, 1);
% k  = min(N, size(S,1));  % Cannot take more than available
% 
% [~, order] = sortrows(S, size(S,2), 'descend');
% topRows    = S(order(1:k), :);
% store{i}.samples{st} = [topRows; samples_tri_prop_val];


samples_pos_err_subesp{i} = samples_err;

samples_pos_err_wm{i,end+1} = samples_err;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ERROR

% Total ERROR
err = mean(samples_err(:,end-1));  % difference with estimates

% Search minimum value
% err = min(error_function(cat(1,samples_tri_orig{i},samples_tri_sts{s,i})));

% Tensor estimation
% err = error_function(cat(1,samples_comp_orig,samples_comp_sts{s,i}), ...
%     samples_comp_test,samples_pos_test,dim);

end


function models_sts = update_models(n_strategies,models_sts, ...
    samples_comp_orig,samples_comp_extra,samples_comp_sts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

parfor s=1:n_strategies
    
    samples_comp_model = cat(1,samples_comp_orig,samples_comp_extra, ...
        vertcat(samples_comp_sts{s,:}));
    
    % Create new model
    X = samples_comp_model(:,1:end-1);
    y = samples_comp_model(:,end);

    % % Test with other kernels if necessary:
    % % 'exponential', 'matern32' or 'matern52'
    % gprMdl = fitrgp(X,y,'KernelFunction', 'matern32', ...
    %     'FitMethod','fic','ActiveSetSize',25);  % 'PredictMethod','fic'
    % % cgprMdl = compact(gprMdl);
    
    Mdl = fitrensemble(X,y,'NumLearningCycles',50);
    models_sts{s} = Mdl;
    
end

end


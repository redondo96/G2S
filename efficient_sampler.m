function elapsed = efficient_sampler(ex,func_name,dims,total_budget)

% continuous_complex_2D, epidemic_3D, rosenbrock_3D

% func_name = "well";
% dims = 24;
% total_budget = 50000;

st = dbstack;
f_name = st.name;

tic

% Fix the seed for reproducibility
rng(42+ex, 'twister');  % 51

%% 1. INPUT PARAMETERS

new_budget = 64;  % budget for next iterations - 4
% n_runs = 12;  % 200

fixed_strategy = 0;  % 0 (no), 1, 2, 3

% input parameter: err_min_value | err_tensor_estimation
err_func_name = 'err_samples_prediction_stdev';

opt_function = str2func(func_name);
name = strcat(func_name,"_",int2str(dims));


%% 2. INTERNAL PARAMETERS

n_strategies = 3;
if fixed_strategy ~= 0
    n_strategies = 1;
end

analysis = logical(1-0);  % USE false  1-1

epsilon = 0.1;  % for epsilon-greedy strategy


%% 3. INITIALIZATION

error_function = str2func(err_func_name);

% paramspace not used
[~,discret,steps,dim,n_dims,splits] = initialize_function(name);

initial_budget = min(2000*numel(splits), floor(0.1*total_budget));  % 20*|50*

% Structural information about subspaces
% n_dims_left_run, dims_before_run & dims_after_run not used
[ndims_run,dims_run,steps_run,discrets_run,~,~,~,dims_left_run] = ...
    initialize_splits(splits,dim,steps,discret,n_dims);
clear steps

% We will execute each strategy at least once.
% After that we concatenate the new results
err_means = cellfun(@(x) zeros(n_strategies,n_strategies), splits, ...
    'UniformOutput', false);

errs_sts = cell(n_strategies,length(splits));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checkings
samples_pos_sts = cell(n_strategies,length(splits));
strategy1 = cell(numel(splits),1);
strategy1_vals = cell(numel(splits),1);
strategy2 = cell(numel(splits),1);
strategy2_vals = cell(numel(splits),1);
strategy3 = cell(numel(splits),1);
strategy3_vals = cell(numel(splits),1);
samples_tri_wm = cell(numel(splits),1);
samples_pos_err_wm = cell(numel(splits),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUN variables are used to share information of the different subspaces
% across the different runs
samples_pos_tri_run = cell(size(splits));
samples_tri_run = cell(size(splits));
% triangs_run = cell(size(splits));

% Variable that stores the samples proposed by each strategy in each subspace
samples_tri_sts = cell(n_strategies,length(splits));
% Variable that stores the samples proposed by each strategy
samples_comp_sts = cell(n_strategies,length(splits));

% To share information
% samples_pos_err_subesp_prev = cell(size(splits));
store = cell(1, numel(splits));

for s = 1:numel(splits)
    d      = numel(splits{s});  % 2 or 3 for this subspace
    empty  = zeros(0, d+1);     % 0 rows, (d coords + val) cols
    
    % One cell per strategy, every cell starts with an empty matrix
    store{s}.samples = repmat({empty}, 1, n_strategies);

    % (optional) Keep meta-data
    % store{s}.dims    = d;           % dimensionality
    % store{s}.idx     = splits{s};   % columns of the global point
end

% To add extra close samples
samples_comp_extra = [];

% For evaluating fit
% samples_comp_eval = [];

% For analysis and representation purposes (plot, etc.)
if analysis
    samples_timeline = cell(numel(splits),1);
    for i = 1:numel(splits)
        samples_timeline{i} = zeros(n_strategies,1);
    end

    % times = zeros(1,n_runs);

    % Save other variables if necessary
end


%% 4. INITIAL SAMPLES

samples_comp_orig = generate_initial_samples_NEW(splits, dim, ...
    initial_budget, discret, opt_function);

% % To avoid duplicates when proposing samples
samples_pos_comp_TRACK = [];

% Output variable
samples_OUTPUT = samples_comp_orig;

% Build a fast hash of already-used points
row2key = @(p) sprintf('%.12g_', p);
taken   = containers.Map('KeyType','char','ValueType','logical');
for k = 1:size(samples_OUTPUT,1)
    taken(row2key(samples_OUTPUT(k,:))) = true;
end

if size(samples_OUTPUT,1) < total_budget
    
    % Create the model with the initial samples
    X = samples_comp_orig(:,1:end-1);
    y = samples_comp_orig(:,end);
    
    Mdl = fitrensemble(X,y,'NumLearningCycles',50);
    models_sts = repmat({Mdl}, 1, n_strategies);  % save the models
    clear X y Mdl
    
    
    %% 5. FIRST RUNS (EACH STRATEGY ONCE)
    
    samples_pos_tri_copy = samples_pos_tri_run;
    samples_tri_copy = samples_tri_run;
    
    samples_pos_err_copy = cell(size(splits));
    
    st = fixed_strategy;
    
    for s=1:n_strategies
        
        % % Time performance
        % if analysis
        %     tic
        % end
    
        if fixed_strategy == 0
            st = s;
        end
        
        % To share information
        samples_pos_err_subesp = cell(size(splits));
    
        for i=1:length(splits)
    
            % Structural information about the subspace
            % (does not change over time)
            % steps_spl not used
            [spl,ndim_spl,dim_spl,~,discret_spl,dims_left] = ...
                get_split_information(splits,ndims_run,dims_run, ...
                steps_run,discrets_run,dims_left_run,i);
    
            % Information about the samples we have so far
            % samples_pos_tri = samples_pos_tri_run{i};
            % samples_tri = samples_tri_run{i};
            % T = triangs_run{i};
            
            % Execute Strategy
            [samples_pos_tri_copy, ...
             samples_OUTPUT, ...
             taken, ...
             samples_pos_comp_TRACK, ...
             samples_tri_copy, ...
             samples_tri_sts, samples_comp_sts, ...
             ~, ...
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
                samples_pos_tri_copy, ...
                i, ...
                spl, ...
                dims_left, dim, ...
                discret, ...
                opt_function, ...
                samples_OUTPUT, ...
                taken, ...
                samples_tri_copy, ...
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
                store);
    
            % Save the information for the end of the initial strategy loop
            samples_pos_err_copy{i} = [samples_pos_err_copy{i}; ...
                samples_pos_err_subesp{i}];
    
            % We store the error in the strategy and subspace compartment
            % in order to calculate the mean for MAB
            errs_sts{s,i} = err;
    
            % We store the error of the first rounds
            err_means{i}(1:n_strategies,s) = err;
    
        end  % [SPLIT/SUBSPACE LOOP]
        
        % % Time performance
        % if analysis
        %     tc = toc;
        %     times(s) = tc;
        %     % fprintf('%i - %f\n', s, tc);
        % end
    
    end  % [STRATEGY LOOP]
    
    % Update the strategy model with new samples
    models_sts = update_models(n_strategies,models_sts, ...
        samples_comp_orig,samples_comp_extra,samples_comp_sts);
    
    % We keep the updated samples
    samples_pos_tri_run = samples_pos_tri_copy;
    samples_tri_run = samples_tri_copy;
    
    % samples_pos_err_subesp = samples_pos_err_copy;  % In reality, they are no longer positions (RENAME)
    
    clear samples_pos_tri_copy samples_tri_copy samples_pos_err_copy
    
    
    % st_ex = zeros(length(splits), n_runs);
    % sam_ex = zeros(length(splits), n_runs);
    % time_ex = zeros(length(splits), n_runs);
    % cls = zeros(1, n_runs);
    % cmb = zeros(1, n_runs);
    % sam_mdl = zeros(n_strategies, n_runs);
    % time_mdl = zeros(n_strategies, n_runs);
    

    %% 6. NEXT RUNS
    
    r = n_strategies+1;
    % fixed_strategy = 3;
    
    % for r=n_strategies+1:6
    while size(samples_OUTPUT,1) < total_budget
    
        % fprintf('Run %i\n', r)
        
        % % Time performance
        % if analysis
        %     tic
        % end
    
        % Duplicate errors
        for j=1:length(splits)
            err_means{j} = [err_means{j}; err_means{j}(end,:)];
        end
    
        % To share information
        % samples_pos_err_subesp_prev = samples_pos_err_subesp;
        samples_pos_err_subesp = cell(size(splits));
    
        % For evaluating fit
        % samples_comp_eval = [];
    
        for i=1:length(splits)
    
            % Structural information about the subspace
            % (does not change over time)
            % steps_spl not used
            [spl,ndim_spl,dim_spl,~,discret_spl,dims_left] = ...
                get_split_information(splits,ndims_run,dims_run, ...
                steps_run,discrets_run,dims_left_run,i);

    
            % STRATEGY SELECTION
            if fixed_strategy == 0
                argmax = epsilon_greedy(epsilon,err_means{i}(r,:), ...
                    n_strategies);
                
                % If the strategy is not fixed, the index in which
                % to store the results will be the one indicated by the
                % selected strategy (argmax)
                s = argmax;
            else
                argmax = fixed_strategy;
                % If the strategy is fixed, the index in which
                % to store the results will be 1 because there is only
                % one strategy
                s = 1;
            end
    
            % fprintf('Split %i. Strategy %i\n', i, argmax)
            
            % tic
            % Execute Strategy
            [samples_pos_tri_run, ...
             samples_OUTPUT, ...
             taken, ...
             samples_pos_comp_TRACK, ...
             samples_tri_run, ...
             samples_tri_sts, samples_comp_sts, ...
             ~, ...
             samples_pos_err_subesp, err, ...
             samples_pos_sts, ...
             strategy1, strategy1_vals, ...
             strategy2, strategy2_vals, ...
             strategy3, strategy3_vals, ...
             samples_tri_wm, samples_pos_err_wm, ...
             store] = execute_strategy(argmax, ...
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
                store);
            
            % tc_ex = toc;
            % fprintf('%i samples. Execution time: %fs\n', ...
            %     size(samples_tri,1), tc_ex)
            % st_ex(i,r) = argmax;
            % sam_ex(i,r) = size(samples_tri,1);
            % time_ex(i,r) = tc_ex;
            
            % We store the error in the strategy and subspace compartment
            % in order to calculate the mean for MAB
            errs_sts{s,i} = cat(1,errs_sts{s,i},err);
    
            % The error is updated with the information that is obtained.
            % That is, the new error is not concatenated, but the average
            % of the errors so far is concatenated
            err_means{i}(r,s) = mean(errs_sts{s,i});
    
            % Information for analysis
            if analysis
                % strategy 'argmax' chosen
                samples_timeline{i} = [samples_timeline{i}; argmax];
            end
    
        end  % [SPLIT/SUBSPACE LOOP]
        
        % Update the strategy model with new samples
        models_sts = update_models(n_strategies,models_sts, ...
            samples_comp_orig,samples_comp_extra,samples_comp_sts);
    
        r = r+1;
        
        % % Time performance
        % if analysis
        %     tc = toc;
        %     times(r) = tc;
        %     fprintf('%i - %f\n', r, tc);
        % end
        
        fprintf('%i samples generated so far\n', size(samples_OUTPUT,1))
        % disp(' ')
    
    end  % [RUN LOOP]

end  % [IF initial_samples LOOP]

elapsed = toc;

% strategies = [samples_timeline{1}, samples_timeline{1}];
% fn = sprintf("experiments/ex%d_strategies.xlsx",ex);
% writematrix(strategies, fn);

filename = sprintf("experiments/%d/ex%d_%s_%s_%d_%d.mat",ex,ex, ...
    f_name,func_name,dims,total_budget);
% _150_50_101_2
% _100_10_101_5_02_max_max
dataset = samples_OUTPUT;
save(filename, "dataset")

end

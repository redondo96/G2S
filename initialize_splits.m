function [ndims_run, dims_run, steps_run, discrets_run, ...
    ndims_left_run, dims_before_run, dims_after_run, dims_left_run] = ...
    initialize_splits(splits, dim,steps, discret, ndims)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

num_splits = length(splits);
dims_run = cell(num_splits,1);
ndims_run = cell(num_splits,1);
steps_run = cell(num_splits,1);
discrets_run = cell(num_splits,1);

ndims_left_run = cell(num_splits,1);
dims_before_run = cell(num_splits,1);
dims_after_run = cell(num_splits,1);
dims_left_run = cell(num_splits,1);

% corners_run = cell(num_splits,1);
% n_corners_run = cell(num_splits,1);

for i=1:length(splits)

    spl = splits{i}; ndims_spl = length(spl);
    dim_spl = dim(spl); steps_spl = steps(spl);
    discret_spl = cell(ndims_spl,1);
    for s=1:ndims_spl
        discret_spl{s} = discret{spl(s)};
    end
    
    ndims_left = ndims-ndims_spl;
    
    dims_before = [];
    dims_after = [];
    for ii=1:length(splits)
        if ii<i
            dims_before = cat(2,dims_before,splits{ii});
        elseif ii>i
            dims_after = cat(2,dims_after,splits{ii});
        end
    end
    
    dims_before = sort(dims_before);
    dims_after = sort(dims_after);
    dims_left = cat(2,dims_before,dims_after);

    % % CORNERS
    % % if ndims_spl==1
    % %     corners = [1;dim_spl];
    % % end
    % if ndims_spl==2
    % corners = [1 1;
    %            1 dim_spl(2);
    %            dim_spl(1) 1;
    %            dim_spl(1) dim_spl(2)];
    % 
    % elseif ndims_spl==3
    % corners = [1 1 1;
    %            1 dim_spl(2) 1;
    %            1 1 dim_spl(3);
    %            1 dim_spl(2) dim_spl(3);
    %            dim_spl(1) 1 1;
    %            dim_spl(1) dim_spl(2) 1;
    %            dim_spl(1) 1 dim_spl(3);
    %            dim_spl(1) dim_spl(2) dim_spl(3)];
    % end
    % n_corners = size(corners,1);
    
    % Outputs
    
    ndims_run{i} = ndims_spl;
    dims_run{i} = dim_spl;
    steps_run{i} = steps_spl;
    discrets_run{i} = discret_spl;
    
    ndims_left_run{i} = ndims_left;
    dims_before_run{i} = dims_before;
    dims_after_run{i} = dims_after;
    dims_left_run{i} = dims_left;

    % corners_run{i} = corners;
    % n_corners_run{i} = n_corners;
end

end


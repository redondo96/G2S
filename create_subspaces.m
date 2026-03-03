function split = create_subspaces(ndims,manual_splits)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (0<ndims) && (ndims<=3)

    split = {1:ndims};

elseif ndims > 3 && nargin == 1  % No manual splits

    [nsub3,nsub2] = number_subsp(ndims);

    split = cell(nsub3+nsub2,1);

    for i=1:nsub3
        split{i} = [3*i-2 3*i-1 3*i];
    end
    
    if nsub2==1
        split{nsub3+nsub2} = [3*nsub3+1 3*nsub3+2];
    elseif nsub2==2
        split{nsub3+nsub2-1} = [3*nsub3+1 3*nsub3+2];
        split{nsub3+nsub2} = [3*nsub3+3 3*nsub3+4];
    end

elseif ndims > 3 && nargin > 1  % Manual splits

    % We check that the manual splits do not share
    % dimensions with each other
    elems = [];
    for i=1:length(manual_splits)
        elems = cat(2,elems,manual_splits{i});
    end
    elems_u = unique(elems);
    if length(elems) ~= length(elems_u)
        error('Subspaces cannot share dimensions. Check manual dimensions')
    end

    
    dims_taken = [];
    for i=1:length(manual_splits)
        dims_taken = cat(2,dims_taken,manual_splits{i});
    end
    dims_taken = sort(dims_taken);

    ndims_left = ndims - length(dims_taken);

    [nsub3,nsub2] = number_subsp(ndims_left);

    autom_splits = cell(nsub3+nsub2,1);

    candidates = setdiff(1:ndims,dims_taken);

    for i=1:nsub3
        autom_splits{i} = [candidates(3*i-2) ...
                           candidates(3*i-1) ...
                           candidates(3*i)];
    end
    if nsub2==1
        autom_splits{nsub3+nsub2} = [candidates(end-1) candidates(end)];
    elseif nsub2==2
        autom_splits{nsub3+nsub2-1} = [candidates(end-3) candidates(end-2)];
        autom_splits{nsub3+nsub2} = [candidates(end-1) candidates(end)];
    end

    split = [manual_splits'; autom_splits];

else
    error('The function cannot have a negative number of dimensions')
end

end


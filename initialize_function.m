function varargout = initialize_function(name)

nOutputs = nargout;
varargout = cell(1,nOutputs);

% rosen, rosen_2, rosen_10, rosen_6
% ackley_6, rastrigin_6, levy_6, schwefel_6, powell_3


%% TEST

if name == "continuous_complex_2D" && nOutputs == 6
    % Load Space
    load Space_Slanted_Border_Function.mat Space X Y

    % Dimensions of the space
    dim = size(Space);
    ndims = size(dim,2);

    discret = cell(2,1);
    discret{1} = X(1,:);
    discret{2} = Y(:,1)';
    
    % Step sizes of the discrete space
    stepsizeX = zeros(dim(1)-1,1);
    for i=1:dim(1)-1
        stepsizeX(i) = X(1,i+1)-X(1,i);
    end
    stepsizeX = mean(stepsizeX);
    
    stepsizeY = zeros(dim(2)-1,1);
    for i=1:dim(2)-1
        stepsizeY(i) = Y(i+1,1)-Y(i,1);
    end
    stepsizeY = mean(stepsizeY);

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    
    % OUTPUT
    paramspace = cell(1,3);
    paramspace{1} = Space;
    paramspace{2} = X; paramspace{3} = Y;
    
    varargout{1} = paramspace;
    % varargout{2} = Space;
    % varargout{3} = X; varargout{4} = Y;
    varargout{2} = discret;

    % varargout{3} = [X(1,1),stepsizeX,X(1,end);
    %                 Y(1,1),stepsizeY,Y(end,1)];
    varargout{3} = [stepsizeX; stepsizeY];
    
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;
    
    % support = [X(1,:); Y(:,1)'];

    % Plot function
    figure
    surf(X,Y,Space,'EdgeColor','none')


elseif name == "epidemic_3D" && nOutputs == 6
    % Load Space
    load ParaSpaceG.mat ParaSpaceG
    start = 0.025; step = 0.025; finish = 1;
    [X,Y,Z] = ndgrid(start:step:finish);
    
    % Dimensions of the space
    dim = size(ParaSpaceG);
    ndims = size(dim,2);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    
    % OUTPUT
    paramspace = cell(1,4);
    paramspace{1} = ParaSpaceG;
    paramspace{2} = X; paramspace{3} = Y; paramspace{4} = Z;

    varargout{1} = paramspace;
    % varargout{2} = ParaSpaceG;
    % varargout{3} = X; varargout{4} = Y; varargout{5} = Z;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);  % repmat([start,step,finish],ndims,1)
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;
    
    % support = [start:step:finish; start:step:finish; start:step:finish];


elseif name == "rosenbrock_3D" && nOutputs == 6
    start = -2; step = 0.1; finish = 2;
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    ndims = 3;

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


%% COMPETITORS

elseif startsWith(name,'rosen_') && nOutputs == 6
    start = -2.048; step = 0.0816; finish = 2.048;  %% 0.05
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    if ndims == 10
        splits = cell(1,4);
        splits{1} = [1 2 3];
        splits{2} = [4 6];
        splits{3} = [5 7];
        splits{4} = [8 9 10];
        % splits{5} = [11 12];
    elseif ndims == 18
        splits = cell(1,6);
        splits{1} = [1 2 3];
        splits{2} = [4 5 6];
        splits{3} = [7 8 9];
        splits{4} = [10 11 12];
        splits{5} = [13 14 15];
        splits{6} = [16 17 18];
    end
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'normal_') && nOutputs == 6

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});
    
    if ndims == 2
        start = -3; step = 1e-3; finish = 3;
    else
        start = -3; step = 1e-2; finish = 3;  % 0.075 | 0.25
    end
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % if ndims == 6
    %     splits = cell(1,3);
    %     splits{1} = [1 2];
    %     splits{2} = [3 4];
    %     splits{3} = [5 6];
    if ndims == 10
        splits = cell(1,4);
        splits{1} = [1 2 3];
        splits{2} = [4 6];
        splits{3} = [5 7];
        splits{4} = [8 9 10];
        % splits{5} = [11 12];
    elseif ndims == 18
        splits = cell(1,6);
        splits{1} = [1 2 3];
        splits{2} = [4 5 6];
        splits{3} = [7 8 9];
        splits{4} = [10 11 12];
        splits{5} = [13 14 15];
        splits{6} = [16 17 18];
    end
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif (startsWith(name,'peak_') || ...
        startsWith(name,'well_') || ...
        startsWith(name,'crossVall_')) && nOutputs == 6

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});
    
    if ndims == 2
        start = -2; step = 1e-3; finish = 2;
    else
        start = -2; step = 0.08; finish = 2;  % 0.075 | 1e-2 | 0.04
    end
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims

    % splits = cell(1,3);
    % splits{1} = [1 2];
    % splits{2} = [3 4];
    % splits{3} = [5 6];

    % if ndims == 6
        % splits = cell(1,2);
        % % splits{1} = [1 2 3]; splits{2} = [4 5 6];
        % % splits{1} = [1 2 5]; splits{2} = [3 4 6];
        % % splits{1} = [1 3 4]; splits{2} = [2 5 6];
        % % splits{1} = [1 4 5]; splits{2} = [2 3 6];
        % splits{1} = [1 5 6]; splits{2} = [2 3 4];
    % end
    
    % if ndims == 10
    %     splits = cell(1,4);
    %     splits{1} = [1 2 3];
    %     splits{2} = [4 6];
    %     splits{3} = [5 7];
    %     splits{4} = [8 9 10];
    %     % splits{5} = [11 12];
    % elseif ndims == 18
    %     splits = cell(1,6);
    %     splits{1} = [1 2 3];
    %     splits{2} = [4 5 6];
    %     splits{3} = [7 8 9];
    %     splits{4} = [10 11 12];
    %     splits{5} = [13 14 15];
    %     splits{6} = [16 17 18];
    % end
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'ackley_') && nOutputs == 6
    start = -32.768; step = 0.65; finish = 32.768;  % 0.1
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % splits = cell(1,4);
    % splits{1} = [1 2 3];
    % splits{2} = [4 6];
    % splits{3} = [5 7];
    % splits{4} = [8 9 10];
    % % splits{5} = [11 12];
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'griewank_') && nOutputs == 6
    start = -600; step = 24; finish = 600;  % 6

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % if ndims == 6
        % splits = cell(1,2);
        % splits{1} = [1 2 3]; splits{2} = [4 5 6];
        % % splits{1} = [1 2 5]; splits{2} = [3 4 6];
        % % splits{1} = [1 3 4]; splits{2} = [2 5 6];
        % % splits{1} = [1 4 5]; splits{2} = [2 3 6];
        % % splits{1} = [1 5 6]; splits{2} = [2 3 4];
    % end
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'langer_') && nOutputs == 6
    start = 0; step = 0.2; finish = 10;

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;

elseif startsWith(name,'rothyp_') && nOutputs == 6
    start = -65.536; step = 2.6; finish = 65.536;  % 2

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'zakharov_') && nOutputs == 6
    start = -5; step = 0.3; finish = 10;  % 0.15

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'michal_') && nOutputs == 6
    start = 0; step = pi/50; finish = pi;

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims

    % splits = cell(1,3);
    % splits{1} = [1 2];
    % splits{2} = [3 4];
    % splits{3} = [5 6];

    % if ndims == 6
        % splits = cell(1,2);
        % splits{1} = [1 2 3]; splits{2} = [4 5 6];
        % % splits{1} = [1 2 5]; splits{2} = [3 4 6];
        % % splits{1} = [1 3 4]; splits{2} = [2 5 6];
        % % splits{1} = [1 4 5]; splits{2} = [2 3 6];
        % % splits{1} = [1 5 6]; splits{2} = [2 3 4];
    % end

    % if ndims == 12
    %     splits = cell(1,6);
    %     splits{1} = [1 2];
    %     splits{2} = [3 4];
    %     splits{3} = [5 6];
    %     splits{4} = [7 8];
    %     splits{5} = [9 10];
    %     splits{6} = [11 12];
    % end

    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'stybtang_') && nOutputs == 6
    start = -5; step = 0.2; finish = 5;

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

%% OTHERS

elseif startsWith(name,'rastrigin_') && nOutputs == 6
    start = -5; step = 0.25; finish = 5;
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % splits = cell(1,4);
    % splits{1} = [1 2 3];
    % splits{2} = [4 6];
    % splits{3} = [5 7];
    % splits{4} = [8 9 10];
    % % splits{5} = [11 12];
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'levy_') && nOutputs == 6
    start = -10; step = 0.5; finish = 10;
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % splits = cell(1,4);
    % splits{1} = [1 2 3];
    % splits{2} = [4 6];
    % splits{3} = [5 7];
    % splits{4} = [8 9 10];
    % % splits{5} = [11 12];
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'schwefel_') && nOutputs == 6
    start = -500; step = 10; finish = 500;
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % splits = cell(1,4);
    % splits{1} = [1 2 3];
    % splits{2} = [4 6];
    % splits{3} = [5 7];
    % splits{4} = [8 9 10];
    % % splits{5} = [11 12];
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


elseif startsWith(name,'powell_') && nOutputs == 6
    start = -4; step = 0.1; finish = 5;
    % [X1,X2,X3,X4,X5,X6] = ndgrid(start:step:finish);

    % Dimensions of the space
    s = split(name,'_');
    ndims = str2double(s{2});

    dim = size(start:step:finish,2)+zeros(1,ndims);

    ssf = repmat([start,step,finish],ndims,1);
    discret = cell(ndims,1);
    for i = 1:ndims
        discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
    end
    
    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    % splits = create_subspaces(ndims,{[1 2 5],[7,13]});  % custom dims
    
    % splits = cell(1,4);
    % splits{1} = [1 2 3];
    % splits{2} = [4 6];
    % splits{3} = [5 7];
    % splits{4} = [8 9 10];
    % % splits{5} = [11 12];
    
    % OUTPUT
    % For now paramspace = [discret, ndims]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = ndims;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = repmat(step,ndims,1);
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;


% elseif name == "sin_1D" && nOutputs == 6
% 
%     start = 0; step = 0.02; finish = 1.6;
% 
%     % Dimensions of the space
%     ndims = 1;
% 
%     dim = size(start:step:finish,2)+zeros(1,ndims);
% 
%     ssf = repmat([start,step,finish],ndims,1);  % start step finish
%     discret = cell(ndims,1);
%     for i = 1:ndims
%         discret{i} = ssf(i,1):ssf(i,2):ssf(i,3);
%     end
% 
%     % SPLITS
%     % We prioritize spaces of dimension 3. This could be changed
%     splits = create_subspaces(ndims);
% 
%     % OUTPUT
%     % For now paramspace = [discret, ndims]
%     paramspace = cell(1,2);
%     paramspace{1} = discret;
%     paramspace{2} = ndims;
% 
%     varargout{1} = paramspace;
%     varargout{2} = discret;
%     varargout{3} = repmat(step,ndims,1);
%     varargout{4} = dim; varargout{5} = ndims;
%     varargout{6} = splits;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************************************************************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SURVEYS

elseif name == "survey_2D" && nOutputs == 7

    %% READ SURVEY DATA

    % zipFileName = fullfile('2023-eleciones-autonomicas', ...
    %                        '2023-05-06-filtered-Madrid-200.zip');
    % zipFileContent = zip_getContent(zipFileName);
    % char({zipFileContent.file_name}) % displays the file names
    % 
    % 
    % 
    % filename = fullfile('2023-eleciones-autonomicas', ...
    %                     '2023-05-06-filtered-Madrid-200', ...
    %                     '2023-05-06-Madrid-200.csv');
    % 
    % opts = detectImportOptions(filename);
    % preview(filename,opts)
    % 
    % dbstop if error
    % survey = readtable(filename);  % 'PreserveVariableNames',true
    % 
    % survey = readcell(filename);
    
    filename = 'survey-Madrid.csv';
    % survey = readtable(filename);  % ERROR - ASK
    survey_cell = readcell(filename);
    survey = cell2table(survey_cell(2:end,:));
    % Variable names
    v_names = strings(1,size(survey_cell,2));
    for v=1:size(survey_cell,2)
        v_names(v) = survey_cell{1,v};
    end
    survey.Properties.VariableNames = v_names;

    % survey-Valencia.csv -> survey (table)
    % GENDER | YEAR
    % PP | medico

    % Dimensions of the space
    ndims = 2;

    discret = cell(ndims,1);

    discret{1} = ["female" "male"];  % Gender
    discret{2} = [1950 1960 1970 1980 1990 2000];  % Year

    parties = ["PP.1" "MM.1" "PS.1" "Vox.1" "UP.1"];
    other_voting = ["blanco.1" "no.vota.1" "no.decidido.1"];
    control = [];  % "medico"

    target_vars = cell(3,1);
    target_vars{1} = parties;
    target_vars{2} = other_voting;
    target_vars{3} = control;
    
    dim = [];
    for i=1:ndims
        dim = cat(2,dim,length(discret{i}));
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    
    % OUTPUT
    % For now paramspace = [discret, survey (table)]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = survey;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = ones(ndims,1);  % step = 1
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;
    varargout{7} = target_vars;


elseif name == "survey_3D" && nOutputs == 7

    % READ SURVEY DATA
    
    filename = 'survey-Valencia-3D.csv';
    % survey = readtable(filename);  % ERROR - ASK
    survey_cell = readcell(filename);
    survey = cell2table(survey_cell(2:end,:));
    % Variable names
    v_names = strings(1,size(survey_cell,2));
    for v=1:size(survey_cell,2)
        v_names(v) = survey_cell{1,v};
    end
    survey.Properties.VariableNames = v_names;

    % survey-Valencia-3D.csv -> survey (table)
    % GENDER | INCOME | YEAR
    % PP | PSOE | UP | medico

    % Dimensions of the space
    ndims = 3;

    discret = cell(ndims,1);

    discret{1} = ["female" "male"];  % Gender
    discret{2} = ["lower" "middle" "high" "prefer_not_to_say"];  % Income
    discret{3} = [1950 1960 1970 1980 1990 2000];  % Year

    target_v = ["PP.2" "PS.2" "UP.2"];
    control_v = "medico";

    n_vars = size([target_v control_v],2);
    
    dim = [];
    for i=1:ndims
        dim = cat(2,dim,length(discret{i}));
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    
    % OUTPUT
    % For now paramspace = [discret, survey (table)]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = survey;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = ones(ndims,1);  % step = 1
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;
    varargout{7} = n_vars;


elseif name == "survey_3DEdu" && nOutputs == 7

    % READ SURVEY DATA
    
    filename = 'survey-Valencia-3DEdu.csv';
    % survey = readtable(filename);  % ERROR - ASK
    survey_cell = readcell(filename);
    survey = cell2table(survey_cell(2:end,:));
    % Variable names
    v_names = strings(1,size(survey_cell,2));
    for v=1:size(survey_cell,2)
        v_names(v) = survey_cell{1,v};
    end
    survey.Properties.VariableNames = v_names;

    % survey-Valencia-3DEdu.csv -> survey (table)
    % GENDER | EDUCATION | YEAR
    % PP | PSOE | UP | medico

    % Dimensions of the space
    ndims = 3;

    discret = cell(ndims,1);

    discret{1} = ["female" "male"];  % Gender
    discret{2} = ["elementary_school" "middle_school" "high_school" ...
        "vocational_technical_college" "university" ...
        "postgraduate"];  % Education
    discret{3} = [1950 1960 1970 1980 1990 2000];  % Year

    target_v = ["PP.2" "PS.2" "UP.2"];
    control_v = "medico";

    n_vars = size([target_v control_v],2);
    
    dim = [];
    for i=1:ndims
        dim = cat(2,dim,length(discret{i}));
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    splits = create_subspaces(ndims);
    
    % OUTPUT
    % For now paramspace = [discret, survey (table)]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = survey;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = ones(ndims,1);  % step = 1
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;
    varargout{7} = n_vars;


elseif name == "survey_5D" && nOutputs == 7

    % READ SURVEY DATA
    
    filename = 'survey-Valencia-5D.csv';
    % survey = readtable(filename);  % ERROR - ASK
    survey_cell = readcell(filename);
    survey = cell2table(survey_cell(2:end,:));
    % Variable names
    v_names = strings(1,size(survey_cell,2));
    for v=1:size(survey_cell,2)
        v_names(v) = survey_cell{1,v};
    end
    survey.Properties.VariableNames = v_names;

    % survey-Valencia-5D.csv -> survey (table)
    % GENDER | EDUCATION | INCOME | EMPLOYMENT | YEAR
    % PP | PSOE | UP | medico

    % Dimensions of the space
    ndims = 5;

    discret = cell(ndims,1);

    discret{1} = ["female" "male"];  % Gender
    discret{2} = ["elementary_school" "middle_school" "high_school" ...
        "vocational_technical_college" "university" ...
        "postgraduate"];  % Education
    discret{3} = ["lower" "middle" "high" "prefer_not_to_say"];  % Income
    discret{4} = ["employed_for_wages" "self_employed" "homemaker" ...
        "unemployed_looking" "unemployed_not_looking" "student" ...
        "retired" "unable_to_work" "other"];  % Employment
    discret{5} = [1950 1960 1970 1980 1990 2000];  % Year

    target_v = ["PP.2" "PS.2" "UP.2"];
    control_v = "medico";

    n_vars = size([target_v control_v],2);
    
    dim = [];
    for i=1:ndims
        dim = cat(2,dim,length(discret{i}));
    end

    % SPLITS
    % We prioritize spaces of dimension 3. This could be changed
    % splits = create_subspaces(ndims);

    splits = cell(1,2);
    splits{1} = [1 5];
    splits{2} = [2 3 4];
    
    % OUTPUT
    % For now paramspace = [discret, survey (table)]
    paramspace = cell(1,2);
    paramspace{1} = discret;
    paramspace{2} = survey;

    varargout{1} = paramspace;
    varargout{2} = discret;
    varargout{3} = ones(ndims,1);  % step = 1
    varargout{4} = dim; varargout{5} = ndims;
    varargout{6} = splits;
    varargout{7} = n_vars;


else
    error(['Function non implemented or ' ...
        'number of arguments do not match space size'])
end

end

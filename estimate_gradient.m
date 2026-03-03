function grad_v = estimate_gradient(varargin)
       % [areas_vols,grad,grad_v,grad_v_norm]
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Replacing `estimate_grad_fun` and `estimate_grad_crit`

% nargin = 2;  % 3

T = varargin{1};
samples = varargin{2};
if nargin == 3
    criticality = varargin{3};
end

ndims = size(samples,2)-1;
n_triang = size(T,1);  % number of triangles

% x = samples(:,1);
% y = samples(:,2);
% z = samples(:,3);
% 
% TR = triangulation(T.ConnectivityList,x,y,z);
% F = featureEdges(TR,pi/6)';
% figure
% plot3(x(F),y(F),z(F),'k');

% Precompute value sources based on nargin
if nargin == 2
    % Gradient of the function values
    get_value = @(vertex_idx) samples(vertex_idx, ndims+1);  % get from samples
elseif nargin == 3
    % Gradient of the criticalities
    get_value = @(vertex_idx) criticality(vertex_idx);  % get from criticality
end

if ndims == 2

    % samples = [x y value]
    vertices = T.ConnectivityList;

    % Extract coordinates of vertices (i, j, k) for all triangles
    % It is possible to change samples to T.Points
    i_coords = samples(vertices(:,1), 1:ndims);  % [n_triang x ndims]
    j_coords = samples(vertices(:,2), 1:ndims);  % [n_triang x ndims]
    k_coords = samples(vertices(:,3), 1:ndims);  % [n_triang x ndims]

    % Extract function/criticality values for all vertices
    i_f = get_value(vertices(:,1));  % [n_triang x 1]
    j_f = get_value(vertices(:,2));  % [n_triang x 1]
    k_f = get_value(vertices(:,3));  % [n_triang x 1]

    % Compute edge vectors (i.e., triangle sides)
    ik = i_coords - k_coords;  % [n_triang x ndims]
    ji = j_coords - i_coords;  % [n_triang x ndims]

    % Compute triangle areas
    areas = 0.5 * ...
        abs((j_coords(:,1) - i_coords(:,1)) .* (k_coords(:,2) - i_coords(:,2)) - ...
            (k_coords(:,1) - i_coords(:,1)) .* (j_coords(:,2) - i_coords(:,2)));  % [n_triang x 1]

    % Alternative way
    % area_t = 0.5 * abs(i(1)*j(2) + j(1)*k(2) + k(1)*i(2) - j(1)*i(2) - k(1)*j(2) - i(1)*k(2));
    
    % Rotation matrix for 90 degrees
    rotation_matrix = [0 -1; 1 0];  % 90-degree rotation in 2D

    % Rotate vectors using the rotation matrix
    ikr = ik * rotation_matrix';  % Rotate ik (results in [n_triang x 2])
    jir = ji * rotation_matrix';  % Rotate ji (results in [n_triang x 2])

    % ang = (pi/180) * 90;  % 90 degrees
    % % NORM IS THE SAME AFTER ROTATION
    % % Rotate vector 1
    % [th,r] = cart2pol(ik(1),ik(2));
    % [xr1,yr1] = pol2cart(th+ang, r);
    % % Rotate vector 2
    % [th,r] = cart2pol(ji(1),ji(2));
    % [xr2,yr2] = pol2cart(th+ang, r);
    % 
    % ikr = [xr1 yr1];
    % jir = [xr2 yr2];

    grad = (j_f - i_f) .* (ikr ./ (2 * areas)) + ...
           (k_f - i_f) .* (jir ./ (2 * areas));

    % disp(grad)

elseif ndims == 3
    
    % volums = zeros(n_triang,1);
    grad = zeros(n_triang,ndims);

    % samples = [x y z value]
    vertices = T.ConnectivityList;

    % Extract coordinates of vertices (i, j, k) for all triangles
    % It is possible to change samples to T.Points
    i_coords = samples(vertices(:,1), 1:ndims);  % [n_triang x ndims]
    j_coords = samples(vertices(:,2), 1:ndims);  % [n_triang x ndims]
    k_coords = samples(vertices(:,3), 1:ndims);  % [n_triang x ndims]
    h_coords = samples(vertices(:,4), 1:ndims);  % [n_triang x ndims]

    % Extract function/criticality values for all vertices
    i_f = get_value(vertices(:,1));  % [n_triang x 1]
    j_f = get_value(vertices(:,2));  % [n_triang x 1]
    k_f = get_value(vertices(:,3));  % [n_triang x 1]
    h_f = get_value(vertices(:,4));  % [n_triang x 1]
    
    % Compute edge vectors (i.e., triangle sides)
    ji = j_coords - i_coords;  % [n_triang x ndims]
    ki = k_coords - i_coords;  % [n_triang x ndims]
    hi = h_coords - i_coords;  % [n_triang x ndims]

    ik = i_coords - k_coords;
    hk = h_coords - k_coords;
    ih = i_coords - h_coords;
    jh = j_coords - h_coords;

    ji_f = j_f - i_f;
    ki_f = k_f - i_f;
    hi_f = h_f - i_f;

    % Number of triangles
    parfor t = 1:n_triang
        
        volum_t = (1/6) * abs(det([ji(t,:); ki(t,:); hi(t,:)]));

        if volum_t ~= 0
            % volums(t) = volum_t;
            
            grad(t,:) = (ji_f(t)) * (cross(ik(t,:),hk(t,:)) / (2*volum_t)) + ...
                        (ki_f(t)) * (cross(ih(t,:),jh(t,:)) / (2*volum_t)) + ...
                        (hi_f(t)) * (cross(ki(t,:),ji(t,:)) / (2*volum_t));

            % If the volume is 0, we leave the gradient in the tetrahedron
            % with value 0
            % because the solid angle will be 0 (or practically 0)
            % anyway and we can ignore its value

            % When the 4 points are on the same plane,
            % the volume could be 0 (det=0)

            % % Example:
            % P1 = [1.8000, 0.9500, 1.1000];
            % P2 = [1.6500, 0.9500, 0.9000];
            % P3 = [1.8500, 0.9000, 0.7500];
            % P4 = [2.0000, 0.9000, 0.9500];
            % 
            % X = [1.8000, 1.6500, 1.8500, 2.0000];
            % Y = [0.9500, 0.9500, 0.9000, 0.9000];
            % Z = [1.1000, 0.9000, 0.7500, 0.9500];
            % 
            % % Points
            % scatter3(X,Y,Z,'.','LineWidth',100)
            % xlim([1.5, 2]);
            % ylim([0.8, 1]);
            % zlim([0.6, 1.2]);
            % 
            % % Vectors
            % PA = [P1;P2];
            % PB = [P1;P3];
            % PC = [P1;P4];
            % 
            % plot3(PA(:,1),PA(:,2),PA(:,3),'r','LineWidth',2)
            % hold on
            % plot3(PB(:,1),PB(:,2),PB(:,3),'r','LineWidth',2)
            % hold on
            % plot3(PC(:,1),PC(:,2),PC(:,3),'r','LineWidth',2)
            % grid on
            % xlim([1.5, 2]);
            % ylim([0.8, 1]);
            % zlim([0.6, 1.2]);
        end
        
    end  % [parfor]
    
    % disp(grad)
end

% % TO RETURN
% if ndims == 2
%     areas_vols = areas;
% elseif ndims == 3
%     areas_vols = volums;
% end

% PLOT
% if ndims == 2
%     for p = 1:size(T.Points,1)
%         figure
%         triplot(T)
%         hold on  
%         triplot(T(V{p,:},:),samples(:,1),samples(:,2),'Color','r')
%         hold on
%         plot(T.Points(p,1),T.Points(p,2),'k.','MarkerSize',12)
%         hold off
%     end
% elseif ndims == 3
%     for p = 1:size(T.Points,1)
%         figure
%         tetramesh(T,'FaceAlpha',0.3)
%         hold on  
%         tetramesh(T(V{p,:},:),samples(:,1),samples(:,2),samples(:,3),'Color','r')
%         hold on
%         plot(T.Points(p,1),T.Points(p,2),T.Points(p,3),'k.','MarkerSize',12)
%         hold off
%     end
% end


%% AVERAGE GRADIENT ON STAR (AGS)

% Compute star (neighbourhood)
V = vertexAttachments(T);
% V{p,:}

% Extract points once for reuse
points = T.Points;

grad_v = zeros(size(points,1),ndims);
% grad_v_norm = zeros(size(points,1),1);

% % Create a reduced version of T.ConnectivityList for each vertex
% vertex_to_triangles = cell(size(points,1), 1);
% for p = 1:size(points,1)
%     vertex_to_triangles{p} = T.ConnectivityList(V{p,:},:);  % Extract relevant triangles
% end

if ndims == 2

    % Number of vertex
    for p = 1:size(points,1)
        % Get triangle indices connected to vertex `p`
        triangle_ids = V{p,:};
        % Triangles
        triangl = T(triangle_ids,:);
        
        % Initialize storage for angles and gradient sums
        num_triangles = size(triangl,1);
        angles = zeros(num_triangles,1);
        sum_grad = zeros(1,ndims);
    
        % Loop over connected triangles to compute neighbours angles
        for i = 1:num_triangles
            tri_id = triangle_ids(i);
            point_ids = triangl(i,:);
    
            point_ids(point_ids==p) = [];
            P1 = points(p,:);
            P2 = points(point_ids(1),:);
            P3 = points(point_ids(2),:);
            
            % In radians
            angle = atan2(2*areas(tri_id), dot(P2-P1,P3-P1));
            % Convert to degrees
            % angle = angle*180/pi;
            
            % Accumulate weighted gradient
            angles(i) = angle;
            sum_grad = sum_grad + angle*grad(tri_id,:);
        end
        
        g = sum_grad/sum(angles);  % The sum of the angles cannot be 0
        % disp(g)
    
        grad_v(p,:) = g;
        % Euclidean norm
        % grad_v_norm(p) = norm(g);
    end

elseif ndims == 3

    % Number of vertex
    for p = 1:size(points,1)
        % Get triangle indices connected to vertex `p`
        triangle_ids = V{p,:};
        % Triangles
        triangl = T(triangle_ids,:);
        
        % Initialize storage for angles and gradient sums
        num_triangles = size(triangl,1);
        angles = zeros(num_triangles,1);
        sum_grad = zeros(1,ndims);
    
        % Loop over connected triangles to compute neighbours angles
        for i = 1:num_triangles
            tri_id = triangle_ids(i);
            point_ids = triangl(i,:);
    
            point_ids(point_ids==p) = [];
            P1 = points(p,:);
            P2 = points(point_ids(1),:);
            P3 = points(point_ids(2),:);
            P4 = points(point_ids(3),:);
        
            % In steradians
            angle = abs(Solid_Angle_Triangle(P1,P2,P3,P4));  % without sign
            
            % angle = solidangle(P2-P1,P3-P1,P4-P1);

            % When the 3 edges forming the angle are on the same plane,
            % the function 'solidangle' could return complex values
            
            % % Example:
            % P1 = [1.2000, 1.7000, 1.4000];
            % P2 = [1.5000, 1.8000, 1.1000];
            % P3 = [1.3000, 1.6000, 1.4000];
            % P4 = [1.4000, 1.9000, 1.1000];
            % 
            % X = [1.2000, 1.5000, 1.3000, 1.4000];
            % Y = [1.7000, 1.8000, 1.6000, 1.9000];
            % Z = [1.4000, 1.1000, 1.4000, 1.1000];
            % 
            % % Points
            % scatter3(X,Y,Z,'.','LineWidth',100)
            % xlim([1, 1.6]);
            % ylim([1.4, 2]);
            % zlim([1, 1.6]);
            % 
            % % Vectors
            % PA = [P1;P2];
            % PB = [P1;P3];
            % PC = [P1;P4];
            % 
            % plot3(PA(:,1),PA(:,2),PA(:,3),'r','LineWidth',2)
            % hold on
            % plot3(PB(:,1),PB(:,2),PB(:,3),'r','LineWidth',2)
            % hold on
            % plot3(PC(:,1),PC(:,2),PC(:,3),'r','LineWidth',2)
            % grid on
            % xlim([1, 1.6]);
            % ylim([1.4, 2]);
            % zlim([1, 1.6]);
            
            % Accumulate weighted gradient
            angles(i) = angle;

            sum_grad = sum_grad + angle*grad(tri_id,:);
        end
        
        % To check (ndims == 3):
        % - The sum of the angles of a vertex point must sum to (4*pi)/8 = 1.5708
        % - The sum of an angle on one side of the cube must sum to (4*pi)/2 = 6.2832
        % - The sum of an interior point must sum to 4*pi = 12.5664
    
        % if ndims == 3
        %     if p<=8
        %         fprintf('corner: %f\n',sum(angles));
        %     elseif sum(angles)<12
        %         fprintf('point [%f, %f, %f]: %f\n',points(p,:),sum(angles));
        %     else
        %         fprintf('%f\n',sum(angles));
        %     end
        % end
        % The same could be done for 2 dimension


        % Note that a center point has 4π sr.
        % Therefore, a corner of the cube will have π/2 sr.
        
        sa = sum(angles);
        g = sum_grad/sa;  % The sum of the angles cannot be 0
        % disp(g)
        
        % Condition to de-emphasize parameter space edges.
        % If we are at an edge, we reduce the gradient value by 2  %%%%%%%% CHECK
        if sa < 4*pi
            grad_v(p,:) = g/2;
        else
            grad_v(p,:) = g;
        end
        % Euclidean norm
        % grad_v_norm(p) = norm(g);
    
    end

end

% CHECK:
% faceNormal
% nearestNeighbor
% vertexNormal

end


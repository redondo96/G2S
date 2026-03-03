function [criticality] = main_sampling_base(dim,samples_pos,samples)  % fit,
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if all(samples==samples(end))

    % Same information to all points
    criticality = ones(size(samples))./numel(samples);

else
    
    % Tolerance (for HOSVD calculation)
    eps = 0.0002;  % 0.0002
    
    % Create sparse tensor
    sparseT = sptensor(samples_pos,samples(:,end),dim);
    
    tensor_orig = full(sparseT);  % Dense tensor
    % T_doub = double(tensor_orig);
    
    % HOSVD estimation
    % 0.25*siz
    % 'ranks',[8, 8]
    decom = hosvd(tensor_orig,2*sqrt(eps),'verbosity',0);  % 'verbosity',1
    
    % CALCULATE FIT
    tensor_recon = full(decom);

    % fit: in reality it is loss
    % fit = norm(tensor_recon-tensor_orig)/norm(tensor_orig);
    
    n_samples = size(samples_pos,1);
    losses = zeros(n_samples,1);
    % Degree of fit
    for i = 1:n_samples
        losses(i) = abs(tensor_orig(samples_pos(i,:))-tensor_recon(samples_pos(i,:)));
    end
    
    % Criticality
    sum_losses = sum(losses,'all');
    if sum_losses ~= 0
        criticality = losses ./ sum_losses;
    else
        % Same probability for all of them
        criticality = ones(n_samples,1) * (1/n_samples);
    end
    
    % if all(losses == 0)
    %     % If the tensor has been perfectly recomposed
    %     % We give random values to 'criticality'
    %     rand_v = rand(1,n_samples)';
    %     % Normalize values to sum to 1
    %     criticality = rand_v / sum(rand_v);
    %     warning('Randomly generated run samples')
    % else
    %     criticality = losses./sum_losses;
    % end

end

% PLOT CRITICALITIES
% X, Y
% % indexes = samples_pos;
% % xx = zeros(1,n_samples);
% % yy = zeros(1,n_samples);
% % for i = 1:n_samples
% %     xx(i) = X(indexes(i,1),indexes(i,2));
% %     yy(i) = Y(indexes(i,1),indexes(i,2));
% % end
% 
% % x = [0.1750 0.4210 0.8350 0.5080 0.6310 0.3400 0.6250 0.7390 0.1120 0.6760 0.8590 0.6730 0.3130 0.7060 0.1780 0.4420 0.7960 0.3250 0.3190 0.1990 0.1120 0.6190 0.2620 0.3010 0.4750 0.1390 0.5380 0.2110 0.5500 0.6370];
% % y = [0.5180 0.7520 0.2000 0.4310 0.3110 0.2690 0.3410 0.4640 0.5030 0.6110 0.5210 0.7250 0.3560 0.8720 0.2210 0.7130 0.5180 0.6290 0.3050 0.3500 0.8120 0.9410 0.4400 0.7310 0.8720 0.8870 0.2630 0.2270 0.3290 0.8720];
% % criticality = [0.0000 0.0006 0 0 0 0 0 0 0.0001 0 0 0 0 0.0797 0 0.0051 0.0000 0.0027 0 0.6491 0.0000 0.0001 0.2592 0.0000 0.0000 0.0000 0 0 0 0.0034];
% 
% % scatter3(x,y,criticality,[],criticality,'filled')  % '.r','markersize',10
% % axis([0.1 1 0.2 1 0 1])
% 
% scatter3(samples(:,1),samples(:,2),criticality,'filled')
% % colorbar
% % colormap turbo

end


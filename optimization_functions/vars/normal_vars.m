function y = normal_vars(vars)
% Evaluate the 6D multivariate normal PDF at point x, given mean mu and cov Sigma.
%
% x     : dx1 vector
% mu    : dx1 mean vector
% Sigma : dxd covariance matrix

    d = width(vars);

    x = zeros(1,d);
    for i=1:d
        x(i) = vars.("x"+i);
    end

    mu = zeros(d,1);
    Sigma = eye(d);
    
    detSigma = det(Sigma);
    invSigma = inv(Sigma);
    
    % Normalizing constant
    normConst = 1 / sqrt((2*pi)^d * detSigma);
    
    diff = (x(:) - mu);
    exponent = -0.5 * (diff' * invSigma * diff);
    
    y = -(normConst * exp(exponent));
end

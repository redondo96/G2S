function pdf_val = normal(x)
% Evaluate the 6D multivariate normal PDF at point x, given mean mu and cov Sigma.
%
% x     : dx1 vector
% mu    : dx1 mean vector
% Sigma : dxd covariance matrix
    
    x = x(:);
    d = length(x);

    mu = zeros(d,1);
    Sigma = eye(d);
    
    detSigma = det(Sigma);
    invSigma = inv(Sigma);
    
    % Normalizing constant
    normConst = 1 / sqrt((2*pi)^d * detSigma);
    
    diff = (x - mu);
    exponent = -0.5 * (diff' * invSigma * diff);
    
    pdf_val = -(normConst * exp(exponent));
end

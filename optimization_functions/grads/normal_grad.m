function grad = normal_grad(x)

% Computes:
%   grad = gradient(fval wrt x), which is a d x 1 vector
%
% Inputs:
%   x    : d x 1 column vector
%   mu   : d x 1 column vector (default zero if omitted)
%   Sigma: d x d covariance (default identity if omitted)
    
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
    
    % pdf_val = -(normConst * exp(exponent));

    % Gradient of fval w.r.t x
    % = c * exp(exponent) * invSigma * (x - mu).
    grad = normConst * exp(exponent) * (invSigma * diff);
end

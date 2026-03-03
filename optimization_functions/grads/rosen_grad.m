function g = rosen_grad(xx)

% -------------------------------------------------------------
% Gradient of the Rosenbrock function in d dimensions.
% Returns a column vector of size d x 1.
% -------------------------------------------------------------

    d = length(xx);
    g = zeros(d,1);
    for i = 1:d
        % Part that depends on (x_i, x_{i+1}) => if i < d
        if i < d
            g(i) = g(i) - 400 * xx(i) * (xx(i+1) - xx(i)^2) + 2 * (xx(i) - 1);
        end
        
        % Part that depends on (x_{i-1}, x_i) => if i > 1
        if i > 1
            g(i) = g(i) + 200 * (xx(i) - xx(i-1)^2);
        end
    end
end

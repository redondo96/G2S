function val = crossVall(x)
% crossVall  Sum of axis-aligned cosine valleys in R^d.
%
%   val = crossVall(x)
%
%   INPUT
%     x    : vector (d×1 or 1×d), d ≥ 1
%
%   OUTPUT
%     val  : scalar f(x)  (−d*amp ≤ f(x) ≤ 0)
%
%   The domain of interest is the hypercube [-2,2]^d,
%   but the function is defined for any real x.
    
    % Depth of each valley
    amp = 1;
    % Half-width w of each valley
    width = 1;

    x = x(:);  % column vector
    v = @(u) (abs(u)<=width) .* (-0.5*amp*(cos(pi*u/width)+1));
    val = sum(arrayfun(v, x));

end

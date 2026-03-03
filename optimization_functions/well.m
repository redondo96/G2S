function val = well(x)

% Well Function:
%   “Well with plateau” function in n dimensions.
%   - a: inner radius of the flat plateau.
%   - b: outer radius where f goes to 0 (b > a).
%   - c: value on the plateau (e.g., c < 0 if a minimum is desired).
%
%   val = c                                                  when r <= a
%   val smoothly transitions from c to 0 with a half-cosine  when a < r < b
%   val = 0                                                  when r >= b

    a = 1.0;
    b = 2.0;
    c = -1.0;  % negative

    x = x(:);
    r = norm(x,2);  % Euclidean distance
    if r >= b
        val = 0;
    elseif r <= a
        val = c;
    else
        % We calculate alpha in [0,1]
        alpha = (r - a) / (b - a);
        % Half-cosine transition
        val = c + (0 - c) * (1 - cos(pi*alpha))/2;
    end
    

    % ---------------------------------------------------------------------

    % a = 1.0;
    % b = 2.0;
    % c = -1.0;
    % 
    % x   = x(:);
    % d   = numel(x);
    % val = 0.0;
    % 
    % sumsq = 0.0;
    % for i = 1:d
    %     sumsq = sumsq + x(i)^2;
    %     r     = sqrt(sumsq);
    % 
    %     if r <= a
    %         w = c;
    %     elseif r >= b
    %         w = 0.0;
    %     else
    %         alpha = (r - a) / (b - a);
    %         w = c + (1 - cos(pi*alpha))/2 * (0 - c);
    %     end
    % 
    %     val = val + w;
    % end

    
    % GRADIENT
    % 
    % a = 1.0; 
    % b = 2.0;
    % c = -1.0;  % negative
    % 
    % x = x(:);
    % r = norm(x, 2);  % Euclidean distance
    % if r >= b
    %     val = 0;
    % elseif r <= a
    %     val = 0;
    % else
    %     alpha  = (r - a) ./ (b - a);  % in (0,1)
    %     % derivative wrt r
    %     df_dr  = (0 - c) * (pi/(2*(b - a))) .* sin(pi*alpha);
    %     val = abs(df_dr);
    % end
    

    % ---------------------------------------------------------------------
    
    % well_shifted:  Flat-bottom cosine well centred at mu in R^d.
    %
    %   val = well_shifted(x, a, b, c)
    %   val = well_shifted(x, a, b, c, mu)
    %
    %   INPUT
    %     x  : d×1 or 1×d vector (evaluation point)
    %     a  : inner radius (flat inner plateau)
    %     b  : outer radius (value is 0 for r >= b), must satisfy b > a
    %     c  : plateau value (typically c < 0 for a "well")
    %     mu : (optional) centre of the well (default mu = 0.5*ones(d,1))
    %
    %   OUTPUT
    %     val : scalar function value
    %
    %   Domain of interest: [-2,2]^d (not enforced here; clamp/check if needed).
    % 
    % a = 0.5;
    % b = 1.0;
    % c = -1.0;  % negative
    % 
    % x = x(:);
    % d = numel(x);
    % mu = 0.5 * ones(d,1);
    % 
    % r = norm(x - mu, 2);  % shifted radius
    % 
    % if r >= b
    %     val = 0;
    % elseif r <= a
    %     val = c;
    % else
    %     alpha = (r - a) / (b - a);  % in (0,1)
    %     val = c + (0 - c) * (1 - cos(pi * alpha)) / 2;
    % end

end

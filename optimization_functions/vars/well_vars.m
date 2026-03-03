function val = well_vars(vars)

% Well Function:
%   “Well with plateau” function in n dimensions.
%   - a: inner radius of the flat plateau.
%   - b: outer radius where f goes to 0 (b > a).
%   - c: value on the plateau (e.g., c < 0 if a minimum is desired).
%
%   val = c                                                  when r <= a
%   val smoothly transitions from c to 0 with a half-cosine  when a < r < b
%   val = 0                                                  when r >= b

    d = width(vars);
    x = zeros(1,d);
    for i=1:d
        x(i) = vars.("x"+i);
    end

    a = 1.0; 
    b = 2.0;
    c = -1.0;  % negative

    x = x(:);
    r = norm(x, 2);  % Euclidean distance
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

end

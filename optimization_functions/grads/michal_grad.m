function [gradMag] = michal_grad(evalPts)
%
%   INPUTS
%     evalPts  :  N×d matrix – each row is a point x = [x1 ... xd].
%
%   OUTPUTS
%     gradMag  :  N×1 vector – gradient magnitude at each point.
    
    m = 10;
    
    N = size(evalPts,1);
    gradMag = zeros(N,1);

    iIdx = (1:size(evalPts,2)).';
    two_m = 2*m;
    for k = 1:N
        x = evalPts(k,:).';
        a = sin(x);
        b = sin(iIdx .* x.^2 / pi);
        b2m    = b.^two_m;
        dbdx   = cos(iIdx .* x.^2 / pi) .* (2*iIdx .* x / pi);
        dti_dx = cos(x) .* b2m + a .* two_m .* b.^(two_m-1) .* dbdx;
        g      = -dti_dx;
        gradMag(k) = norm(g);
    end

end

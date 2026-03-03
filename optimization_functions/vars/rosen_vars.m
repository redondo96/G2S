function y = rosen_vars(vars)
    % 'vars' is an array or table of variables
    % typically in MATLAB bayesopt, each dimension is a separate variable
    % e.g., vars.x1, vars.x2, ..., vars.x18

    dim = width(vars);

    x = zeros(1,dim);
    for i=1:dim
        x(i) = vars.("x"+i);
    end

    val = 0;
    for i=1:dim-1
        val = val + 100*( x(i+1) - x(i)^2 )^2 + (x(i) - 1)^2;
    end
    y = val;
end

function y = peak_vars(vars)
% Evaluate the peak function at point x.

    dim = width(vars);

    x = zeros(1,dim);
    for i=1:dim
        x(i) = vars.("x"+i);
    end

    % Centers of each peak in R^N
    c1 = -1*ones(1,dim);
    c2 =  0*ones(1,dim);
    c3 =  1*ones(1,dim);

    % Scale parameter (sigma) for peak width
    sigma = 0.15;

    % Calculation of the square distance to each center
    dist1 = sum((x - c1).^2);
    dist2 = sum((x - c2).^2);
    dist3 = sum((x - c3).^2);

    % Different amplitudes
    A1 = 1.0;
    A2 = 0.6;  % less amplitude in the center
    A3 = 1.0;

    % Gaussian distributions are added together
    val =   A1*exp(-dist1 / (2*sigma^2)) ...
          + A2*exp(-dist2 / (2*sigma^2)) ...
          - A3*exp(-dist3 / (2*sigma^2));

    % ---------------------------------------------------------------------

    % 2. We add multidimensional deterministic noise using sine
    % functions (high frequencies -> “fine” variation).
    % 
    % Example: for each dimension i, we use a different frequency.
    % We generate a frequency vector (freqs) of length d, and
    % optionally a phase vector (phases).

    baseFreq = 4;  % base frequency in the first dimension
    freqs   = (baseFreq : baseFreq + (dim-1)) * (2*pi);  % e.g. 3,4,5,... in radians
    phases  = linspace(0, pi/2, dim);  % distributed example phases

    % Overall noise amplitude
    amp = 0.1;

    noise = 0;
    for i = 1:dim
        noise = noise + sin(freqs(i)*x(i) + phases(i));
    end
    % We standardize to avoid excessive accumulation
    noise = noise / dim;
    % We scale up using amplitude
    noise = amp * noise;

    % We add the noise to the final value
    y = val + noise;

end

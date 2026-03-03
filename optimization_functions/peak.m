function val = peak(x)
    % peak: Function with peaks in N dimensions.
    % Input:
    %   x: vector of length N where each xi ∈ [-2, 2].
    %
    % Output:
    %   val: value of the peak function.

    % Ensure that x is a vector of length N.
    x = x(:).';
    d = length(x);

    % Centers of each peak in R^N
    c1 = -0.5 * ones(1,d);
    % c2 =  0 * ones(1,d);
    c2 =  0.5 * ones(1,d);
    % c3 =  0.4 * ones(1,d);

    % Scale parameter (sigma) for peak width
    sigma = 0.4;

    % Calculation of the square distance to each center
    dist1 = sum((x - c1).^2);
    dist2 = sum((x - c2).^2);
    % dist3 = sum((x - c3).^2);
    % dist4 = sum((x - c4).^2);

    % Different amplitudes
    A1 = 3.0;
    A2 = 1.0;
    % A3 = 1.0;
    % A4 = 1.0;

    % Gaussian distributions are added together
    val =   A1*exp(-dist1 / (2*sigma^2)) ...
          - A2*exp(-dist2 / (2*sigma^2));
          % - A3*exp(-dist3 / (2*sigma^2));
          % - A4*exp(-dist4 / (2*sigma^2));
    
    % ---------------------------------------------------------------------
    
    % %% 2. We add multidimensional deterministic noise using sine
    % % functions (high frequencies -> “fine” variation).
    % 
    % % Example: for each dimension i, we use a different frequency.
    % % We generate a frequency vector (freqs) of length d, and
    % % optionally a phase vector (phases).
    % 
    % baseFreq = 4;  % base frequency in the first dimension
    % freqs = (baseFreq : baseFreq + (d-1)) * (2*pi);  % e.g. 3,4,5,... in radians
    % phases = linspace(0, pi/2, d);  % distributed example phases
    % 
    % % Overall noise amplitude
    % amp = 0.1;
    % 
    % noise = 0;
    % for i = 1:d
    %     noise = noise + sin(freqs(i)*x(i) + phases(i));
    % end
    % % We standardize to avoid excessive accumulation
    % noise = noise / d;
    % % We scale up using amplitude
    % noise = amp * noise;
    % 
    % % We add the noise to the final value
    % val = val + noise;


    % sigma = 0.40;
    % A1    = 3.0;
    % A2    = 1.0;
    % mu    = 0.5;
    % 
    % x = x(:);
    % d = numel(x);
    % 
    % s1 = 0;
    % s2 = 0;
    % val = 0;
    % 
    % for i = 1:d
    %     s1 = s1 + (x(i)+mu)^2;
    %     s2 = s2 + (x(i)-mu)^2;
    % 
    %     val = val ...
    %           + A1 * exp(-s1 /(2*sigma^2)) ...
    %           - A2 * exp(-s2 /(2*sigma^2));
    % end


    
    % GRADIENT
    % 
    % sigma = 0.4;
    % evalPts = x;
    % 
    % N = size(evalPts,1);
    % gradMag = zeros(N,1);
    % 
    % % % centres and amplitudes
    % % d   = size(evalPts,2);
    % % C   = [-0.5*ones(d,1),  0.5*ones(d,1)];
    % % A   = [3,               1];
    % % 
    % % for k = 1:N
    % %     x = evalPts(k,:).';
    % %     g = zeros(d,1);
    % %     for j = 1:2
    % %         diff  = x - C(:,j);
    % %         coeff = A(j)*exp(-sum(diff.^2)/(2*sigma^2));
    % %         g     = g - (coeff/sigma^2).*diff;
    % %     end
    % %     gradMag(k) = norm(g);
    % % end
    % 
    % d   = size(evalPts,2);
    % c1  = -0.5 * ones(1,d);  % centre of large positive peak
    % c2  =  0.5 * ones(1,d);  % centre of smaller negative peak
    % A1    = 3.0;
    % A2    = 1.0;
    % invS2 = 1 / sigma^2;
    % 
    % for k = 1:N
    %     x  = evalPts(k,:);
    % 
    %     % Squared distances to each centre
    %     dist1 = sum((x - c1).^2);
    %     dist2 = sum((x - c2).^2);
    % 
    %     % gradient = -(A1/σ²)·exp(-d1/2σ²)(x-c1) + (A2/σ²)·exp(-d2/2σ²)(x-c2)
    %     g  = -(A1*invS2)*exp(-dist1/(2*sigma^2)) * (x - c1) ...
    %          + (A2*invS2)*exp(-dist2/(2*sigma^2)) * (x - c2);
    % 
    %     gradMag(k) = norm(g);
    % end

end

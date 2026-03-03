function argmax = epsilon_greedy(epsilon,err_means_v,n_strategies)

% EPSILON-GREEDY
% We use epsilon to encourage exploration.
% If we do not want so much exploration try 'epsilon-decreasing',
% 'Initial Optimism epsilon-greedy' or 'VDBE'

rand_val = rand;

if rand_val < epsilon
    % Exploration: choose a strategy at random
    argmax = randi(n_strategies);
else
    % Exploitation: choose the best known strategy
    % We want to focus on the strategy that yields the maximum errors
    % (put more samples there)
    [~, argmax] = max(err_means_v);
end

end

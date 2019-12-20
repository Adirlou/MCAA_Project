function [pi0_change] = compute_limiting_dist(initial_distributions, transition_prob, limiting_t, statespace_size)
% compute_limiting_dist Compute the evolution of given initial
% distributions over time

% Know if it is time-homogenous or not
[~, ~, interval_time] = size(transition_prob);

% Get the nmber of initial ditributions given
[num_initials,~] = size(initial_distributions);

% Create the tensor that will contain the updates for our initial vectors
pi0_change = zeros(limiting_t, statespace_size, num_initials);

% Set initial distributions
pi0_change(1, :, :) = initial_distributions;

% If the transition is not time-homogeneous, use the various transition
% matrices
if interval_time > 1
    
    % Update each distribution according to the transition matrix (heterogenous)
    for t=2:interval_time
        for s=1:num_initials
            pi0_change(t, :, s) = pi0_change(t-1, :, s) * transition_prob(:, :, t);
        end
    end
else
    
    % Update each distribution according to the transition matrix (homogenous)
    for t=2:limiting_t
        for s=1:num_initials
            pi0_change(t, :, s) = pi0_change(t-1, :, s) * transition_prob(:, :, 1);
        end
    end

end
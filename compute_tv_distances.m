function [tv_distances] = compute_tv_distances(pi0_change)
% compute_tv_distances Given the evolution of one or multiple initial
% distributions, compute the evolution of the TV distance with the value of
% distribution for largest time t

% Get the nmber of initial ditributions given
[limiting_t, ~, num_initials] = size(pi0_change);

tv_distances = zeros(limiting_t, num_initials);

for t=1:limiting_t
    for s=1:num_initials
        tv_distances(t, s) = tv_dist(pi0_change(t, :, s), pi0_change(limiting_t, :, s));
    end
end

end
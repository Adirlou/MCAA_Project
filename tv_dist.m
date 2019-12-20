function [tv_distance] = tv_dist(p, q)
% TV_dist Computes the TV distance between distributions p and q

tv_distance = 0.5 * sum(abs(p-q));

end
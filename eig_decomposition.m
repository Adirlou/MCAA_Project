function [stationary_dist, spectral_gap] = eig_decomposition(transition_prob)
% eig_decomposition Compute the eigenvalue decomposition of the given
% transition matrix, and return the stationary distribution together with the
% spectral gap

% Compute the eig decomposition (note that we keep the left eigenvectors)
[~, D, W] = eig(transition_prob);

eigenvalues = diag(D);

% Sort the eigenvalues (in order to retrieve the eigenvector with
% eigenvalue 1 easily
[~, ind] = sort(eigenvalues, 'descend');

ordered_W = W(:, ind);

% Retrieve the corresponding eigenvector and normalize it
eigenvector = ordered_W(:, 1);
stationary_dist = eigenvector / sum(eigenvector);

% To get lambda_star, we sort eigenvalues wrt their absolute value
sorted_eigenvalues_abs = sort(abs(eigenvalues), 'descend');

% Compute the spectral gap
spectral_gap = 1 - sorted_eigenvalues_abs(2);

end
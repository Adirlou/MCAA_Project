%% Testing for Time-homogeneity

% Size of the state space
statespace_size = 5;

% Number of chains to run
n_chain = 1000000;

% Amount of time 
time = 50;

% Initial distribution (uniform)
pi0 = ones(1, statespace_size) / statespace_size;

% Run the chain and get outputs
X = chain_3(n_chain, time + 1, pi0);

% Compute the transition matrix (when testing for time-homogeneity)
transition_prob = get_transition_matrix(X, n_chain, time, statespace_size, true);

% We plot the evolution of the transition matrix as a function of t
plot(reshape(transition_prob, 25, time).')
title('Evolution of the elements $p_{ij}^{(t)}$ of $P^{(t)}$ as a function of $t$','interpreter','latex')
xlabel('$$t$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "p_t_chain_3.png")

% Here, the chain is not time-homogeneous, the transition matrix actually
% alternate between two behaviors.

%% Computing the evolution of pi0 (uniform)

% Compute the evolution of pi0 using the transition matrix computed before
pi0_change = compute_limiting_dist(pi0, transition_prob, time, statespace_size);

% Plot the evolution of each element of the limiting distribution
plot(pi0_change)
title('Evolution of components of $\pi^{(t)}$ as a function of $t$', 'interpreter', 'latex')
legend({'$\pi_1^{(t)}$','$\pi_2^{(t)}$', '$\pi_3^{(t)}$', '$\pi_4^{(t)}$','$\pi_5^{(t)}$'}, 'interpreter', 'latex')
xlabel('$$t$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "pi_t_from_uniform_chain_3.png")

% We see here that there is no limiting distribution.
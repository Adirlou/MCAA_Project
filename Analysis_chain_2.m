%% Testing for Time-homogeneity

% Size of the state space
statespace_size = 5;

% Number of chains to run
n_chain = 500000;

% Amount of time 
time = 50;

% Initial distribution (uniform)
pi0 = ones(1, statespace_size) / statespace_size;

% Run the chain and get outputs
X = chain_2(n_chain, time + 1, pi0);

% Compute the transition matrix (when testing for time-homogeneity)
transition_prob = get_transition_matrix(X, n_chain, time, statespace_size, true);

% We plot the evolution of the transition matrix as a function of t
plot(reshape(transition_prob, 25, time).')
title('Evolution of the elements $p_{ij}^{(t)}$ of $P^{(t)}$ as a function of $t$','interpreter','latex')
xlabel('$$t$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "p_t_chain_2.png")

% Since the chain is time-homogeneous, we recompute the transition matrix
% more precisely, no need for a sequence of transition matrices

%% Recomputing the transition matrix more precisely

% Since time-homogeneous, only need 1 chain
n_chain = 1;

% Amount of time 
time = 1000000;

% Run the chain and get outputs
X = chain_2(n_chain, time, pi0);

% Compute the transition matrix (knowing it is time-homogeneous)
transition_prob = get_transition_matrix(X, n_chain, time, statespace_size, false);

% Save the transition matrix
%save('P_hat_chain_2.mat','transition_prob')

%% Plot Markov Chain

% Convert to Markov Chain Matlab format
markov_chain = dtmc(transition_prob);

% Plot the underlying graph
graphplot(markov_chain, 'LabelEdges',true)

% Save the graph of the chain
%saveas(gcf, "graph_chain_2.png")

%% Computing the evolution of pi0 (uniform)

% Step until which to update pi(t)
limiting_t = 100;

% Compute the evolution of pi0 using the transition matrix computed before
pi0_change = compute_limiting_dist(pi0, transition_prob, limiting_t, statespace_size);

% Plot the evolution of each element of the limiting distribution
plot(pi0_change)
title('Evolution of components of $\pi^{(t)}$ as a function of $t$', 'interpreter', 'latex')
legend({'$\pi_1^{(t)}$','$\pi_2^{(t)}$', '$\pi_3^{(t)}$', '$\pi_4^{(t)}$','$\pi_5^{(t)}$'}, 'interpreter', 'latex')
xlabel('$$t$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "pi_t_from_uniform_chain_2.png")

% We see that there is no limiting distribution.

%% Get stationary distribution and spectral gap from eigenvalue decomposition

% Compute the spectral decomposition, but return only the stationary
% distribution found (since there is no limiting distribution here, we
% don't return the spectral gap)
[stationary_from_eig, ~] = eig_decomposition(transition_prob);

%% Plotting values of the stationary distribution in bar chart

bar(categorical({'\pi_1','\pi_2','\pi_3','\pi_4', '\pi_5'}), stationary_from_eig)
title('Bar chart representing the components of the stationary distribution $\pi$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "stationary_pi_barchart_chain_2.png")
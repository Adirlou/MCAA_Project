%% Testing for Time-homogeneity

% Size of the state space
statespace_size = 5;

% Number of chains to run
n_chain = 500000;

% Amount of time 
time = 200;

% Initial distribution (uniform)
pi0 = ones(1, statespace_size) / statespace_size;

% Run the chain and get outputs
X = chain_4(n_chain, time + 1, pi0);

% Compute the transition matrix (when testing for time-homogeneity)
transition_prob = get_transition_matrix(X, n_chain, time, statespace_size, true);

% We plot the evolution of the transition matrix as a function of t
plot(reshape(transition_prob, 25, time).')
title('Evolution of the elements $p_{ij}^{(t)}$ of $P^{(t)}$ as a function of $t$','interpreter','latex')
xlabel('$$t$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "p_t_chain_4.png")

% We see here that the chain is not time-homogeneous, and moreover we see
% that each element converges to a value. So as time passes, the Markov
% chain converges to a time-homogeneous one.

%% Plot Markov Chain

% Note here that we choose to display it because this chain converges to a
% time-homogeneous Markov chain. Thus we display the Markov chain to which
% it converges

% Convert to Markov Chain Matlab format
markov_chain = dtmc(transition_prob(:, :, time));

% Plot the underlying graph
graphplot(markov_chain, 'LabelEdges',true)

% Save the graph of the chain
%saveas(gcf, "graph_chain_4.png")

%% Computing the evolution of pi0 (uniform)

% Step until which to update pi(t)
limiting_t = 200;

% Compute the evolution of pi0 using the transition matrix computed before
pi0_change = compute_limiting_dist(pi0, transition_prob, limiting_t, statespace_size);

% Plot the evolution of each element of the limiting distribution
plot(pi0_change)
title('Evolution of components of $\pi^{(t)}$ as a function of $t$', 'interpreter', 'latex')
legend({'$\pi_1^{(t)}$','$\pi_2^{(t)}$', '$\pi_3^{(t)}$', '$\pi_4^{(t)}$','$\pi_5^{(t)}$'}, 'interpreter', 'latex')
xlabel('$$t$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "pi_t_from_uniform_chain_4.png")

% We see here that there is a limiting distribution, so we save it.
% Get the limiting ditribution
limiting_pi = pi0_change(limiting_t, :);

% Here we compute it again, but with the way that was expected from the
% statement (that is using counting), to get (average) error < 10^(-3)
% Note that since X was computed before for large Nchain, we need not produce
% samples again

limiting_pi_expected = expected_limiting(X, n_chain, limiting_t, statespace_size);

% Save the the limiting distribution
%save('pi_hat_chain_4.mat','limiting_pi_expected')

%% Plotting values of the limiting distribution in bar chart

bar(categorical({'\pi_1^{(t)}','\pi_2^{(t)}','\pi_3^{(t)}','\pi_4^{(t)}', '\pi_5^{(t)}'}), limiting_pi)
title('Bar chart representing the components of the limiting distribution $\pi^{(t)}$ as $t \to \infty$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "pi_t_barchart_chain_4.png")

%% Computing TV distance and plotting

% Compute the evolution of the TV distance
tv_distances = compute_tv_distances(pi0_change);

% Plot the evolution of the TV distance
plot(tv_distances)
title('TV distance between $\pi^{(t)}$ and $\pi$ as a function of $t$ from uniform $\pi^{(0)}$','interpreter','latex')
xlabel('$$t$$', 'interpreter', 'latex')
ylabel('$$\|\cdot\|_{TV}$$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "tv_from_uniform_chain_4.png")

%% Computing evolution of pi0 when dirac

% Set the initial distributions, here starting from each possible state
% with probability 1. So this translate into an identity matrix with our
% implementation
initial_distributions = eye(statespace_size);

% Compute the evolution of all the given initial distributions
pi0_change_dirac = compute_limiting_dist(initial_distributions, transition_prob, limiting_t, statespace_size);

% Compute the evolution of the TV distance for all the given initial
% ditributions
tv_distances_dirac = compute_tv_distances(pi0_change_dirac);

% Compute numerically the value T_epsilon
t_epsilon_candidates = zeros(1, statespace_size);

for i=1:statespace_size
    t_epsilon_candidates(i) = find(tv_distances_dirac(:, i) < 0.005, 1);
end

% Find the largest T_epsilon starting from all possible states
t_epsilon = max(t_epsilon_candidates);

% PLot the evolution of the TV distances
plot(tv_distances_dirac)
xline(t_epsilon,':','$$T_\epsilon$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
title('TV distance between $\pi_{\delta_i}^{(t)}$ (starting from state $i$) and $\pi$ as a function of $t$','interpreter','latex')
xlabel('$$t$$', 'interpreter', 'latex')
ylabel('$$\|\cdot\|_{TV}$$', 'interpreter', 'latex')
legend({'$\pi_{\delta_1}^{(t)}$','$\pi_{\delta_2}^{(t)}$', '$\pi_{\delta_3}^{(t)}$', '$\pi_{\delta_4}^{(t)}$','$\pi_{\delta_5}^{(t)}$'}, 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "tv_from_diracs_chain_4.png")

%% Get stationary distribution and spectral gap from eigenvalue decomposition

% Even though this part was not asked, we thought it would be fun to make
% the analysis for when we converge to a time-homogeneous Markov chain.

% Compute the spectral decomposition, but return only the stationary
% distribution found together with the spectral gap
time_homogeneous_transition = transition_prob(:, :, time);
[stationary_from_eig, spectral_gap] = eig_decomposition(time_homogeneous_transition);

% Compute the evolution of the theoretical upperbound as a function of t
upperbound_tv = zeros(limiting_t, 1);

for i=1:limiting_t
    upperbound_tv(i) = (1-spectral_gap)^i / (sqrt(min(stationary_from_eig)));
end

% Find theoretical T_epsilon numerically
t_epsilon_theoretical = find(upperbound_tv < 0.005, 1);

% Plot the upperbound together with the TV distances found previously
plot(upperbound_tv)
hold on
plot(tv_distances)
xline(t_epsilon,':','$$T_\epsilon$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xline(t_epsilon_theoretical,':','$$T_\epsilon^{theoretical}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
hold off
title('TV distance between $\pi^{(t)}$ and $\pi$ and upperbound','interpreter','latex')
xlabel('$$t$$', 'interpreter', 'latex')
ylabel('$$\|\cdot\|_{TV}$$', 'interpreter', 'latex')
legend({'$\max_{i}\lambda_*^t/(2\sqrt{\pi_i})$', '$\|\pi^{(t)} - \pi\|_{TV}$'}, 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "tv_upperbound_chain_4.png")

%% Plotting values of the "stationary" distribution in bar chart

% This also was not explicitely asked, but since we converge to a
% time-homogeneous chain, we still display the equivalent of the stationary
% distribution
bar(categorical({'\pi_1','\pi_2','\pi_3','\pi_4', '\pi_5'}), stationary_from_eig)
title('Bar chart representing the components of the stationary distribution $\pi$', 'interpreter', 'latex')

% Save the plot
%saveas(gcf, "stationary_pi_barchart_chain_4.png")
%% Define some parameters
% Parameters to compute the limiting distribution
N_chain = 5000;
Time = 200;

% Parameters to compute the transition matrix of the a chain
N_chain_transition = 1;
Time_transition = 1e6;

% Extra parameters
statespace_size = 5;
eps = 5e-3;
%% Compute stationnary distribution of the markov chain

% We define the distribution with we want to sample from
pi_a = [16,8,4,2,1]/31;

% We define the limiting distributions matrix
% dim 1 : initial state
% dim 2 : state space
pi_lim = zeros(statespace_size,statespace_size);

for x0=1:statespace_size
    X = MP_chain_2(N_chain,Time,pi_a,x0);
    pi_lim(x0,:) = expected_limiting(X, N_chain, Time, statespace_size);
end

%% Plotting values of the limiting distribution in bar chart
figure()
bar(categorical({'\pi_1^{(t)}','\pi_2^{(t)}','\pi_3^{(t)}','\pi_4^{(t)}', '\pi_5^{(t)}'}), pi_lim.')
legend('$x_0 = 1$','$x_0 = 2$','$x_0 = 3$','$x_0 = 4$','$x_0 = 5$',"Interpreter","latex")
saveas(gcf, "b211.png")

%% Total variation distances and mixing time (numerical approximation)

% We compute transition matrix of the metropolis hasting markov chain,
% initial state is not important, we take 1
% we just take one chain but we take a long time 100000
X = MP_chain_2(N_chain_transition,Time_transition,pi_a,1);
transition_prob = get_transition_matrix(X, N_chain_transition,Time_transition, statespace_size, false);

% We compute the evolution of the distributions on this markov chains for
% each initial states until t = limiting_t
initial_distributions = eye(statespace_size);
limiting_t = 200;
pi0_change_dirac = compute_limiting_dist(initial_distributions, transition_prob, limiting_t, statespace_size);

% We compute the tv distances evolution and the practical mixing time at eps = 5e-3
tv_distances_dirac = compute_tv_distances(pi0_change_dirac);

% We compute numericaly the mixing time at eps = 5e-3
t_epsilon = find(max(tv_distances_dirac,[],2) < eps, 1);

% Compute the evolution of the theoretical upperbound as a function of t
upperbound_tv = zeros(limiting_t, 1);

% Get spectral gap of the chain
[~, spectral_gap] = eig_decomposition(transition_prob);

for i=1:limiting_t
    upperbound_tv(i) = (1-spectral_gap)^i / (sqrt(min(pi_a)));
end

% Find theoretical T_epsilon numerically
t_epsilon_theoretical = find(upperbound_tv < 0.005, 1);

figure()
plot(tv_distances_dirac) 
hold on
plot(upperbound_tv)
xline(t_epsilon,':','$$T_{\epsilon, pr}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xline(t_epsilon_theoretical,':','$$T_{\epsilon, th}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xlabel('$$t$$', 'interpreter', 'latex')
ylabel('$$\|\cdot\|_{TV}$$', 'interpreter', 'latex')
legend({'$\pi_{\delta_1}^{(t)}$','$\pi_{\delta_2}^{(t)}$', '$\pi_{\delta_3}^{(t)}$', '$\pi_{\delta_4}^{(t)}$','$\pi_{\delta_5}^{(t)}$','$\max_i \lambda_*^t/(2\sqrt{\pi_i})$'}, 'interpreter', 'latex')
saveas(gcf,'b212.png')

% Display the mixing time
fprintf('Measured mixing time: %d\n',t_epsilon);
fprintf('Theoretical mixing time: %d\n',t_epsilon_theoretical);
fprintf('Spectral gap: %f\n',spectral_gap);
%% Compute stationnary distribution of the markov chain

% We define the distribution with we want to sample from
pi_a = [1,1,4,1,1]/8;

% We define the limiting distributions matrix
% dim 1 : initial state
% dim 2 : state space
pi_lim = zeros(statespace_size,statespace_size);

for x0=1:statespace_size
    X = MP_chain_2(N_chain,Time,pi_a,x0);
    pi_lim(x0,:) = expected_limiting(X, N_chain, Time, statespace_size);
end

%% Plotting values of the limiting distribution in bar chart
figure()
bar(categorical({'\pi_1^{(t)}','\pi_2^{(t)}','\pi_3^{(t)}','\pi_4^{(t)}', '\pi_5^{(t)}'}), pi_lim.')
legend('$x_0 = 1$','$x_0 = 2$','$x_0 = 3$','$x_0 = 4$','$x_0 = 5$',"Interpreter","latex")
saveas(gcf, "b221.png")

%% Total variation distances and mixing time (numerical approximation)

% We compute transition matrix of the metropolis hasting markov chain,
% initial state is not important, we take 1
% we just take one chain but we take a long time 100000
X = MP_chain_2(N_chain_transition,Time_transition,pi_a,1);
transition_prob = get_transition_matrix(X,N_chain_transition,Time_transition, statespace_size, false);

% We compute the evolution of the distributions on this markov chains for
% each initial states until t = limiting_t
initial_distributions = eye(statespace_size);
limiting_t = 200;
pi0_change_dirac = compute_limiting_dist(initial_distributions, transition_prob, limiting_t, statespace_size);

% We compute the tv distances evolution and the practical mixing time at eps = 5e-3
tv_distances_dirac = compute_tv_distances(pi0_change_dirac);

% We compute numericaly the mixing time at eps = 5e-3
t_epsilon = find(max(tv_distances_dirac,[],2) < eps, 1);

% Compute the evolution of the theoretical upperbound as a function of t
upperbound_tv = zeros(limiting_t, 1);

% Get spectral gap of the chain
[~, spectral_gap] = eig_decomposition(transition_prob);

for i=1:limiting_t
    upperbound_tv(i) = (1-spectral_gap)^i / (sqrt(min(pi_a)));
end

% Find theoretical T_epsilon numerically
t_epsilon_theoretical = find(upperbound_tv < 0.005, 1);

figure()
plot(tv_distances_dirac) 
hold on
plot(upperbound_tv)
xline(t_epsilon,':','$$T_{\epsilon, pr}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xline(t_epsilon_theoretical,':','$$T_{\epsilon, th}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xlabel('$$t$$', 'interpreter', 'latex')
ylabel('$$\|\cdot\|_{TV}$$', 'interpreter', 'latex')
legend({'$\pi_{\delta_1}^{(t)}$','$\pi_{\delta_2}^{(t)}$', '$\pi_{\delta_3}^{(t)}$', '$\pi_{\delta_4}^{(t)}$','$\pi_{\delta_5}^{(t)}$','$\max_i \lambda_*^t/(2\sqrt{\pi_i})$'}, 'interpreter', 'latex')
saveas(gcf,'b222.png')

% Display the mixing time
fprintf('Measured mixing time: %d\n',t_epsilon);
fprintf('Theoretical mixing time: %d\n',t_epsilon_theoretical);
fprintf('Spectral gap: %f\n',spectral_gap);
%% Compute stationnary distribution of the markov chain

% We define the distribution we want to sample from
pi_a = [4,2,1,2,4]/13;

% We define the limiting distributions matrix
% dim 1 : initial state
% dim 2 : state space
pi_lim = zeros(statespace_size,statespace_size);

for x0=1:statespace_size
    X = MP_chain_2(N_chain,Time,pi_a,x0);
    pi_lim(x0,:) = expected_limiting(X, N_chain, Time, statespace_size);
end

%% Plotting values of the limiting distribution in bar chart
figure()
bar(categorical({'\pi_1^{(t)}','\pi_2^{(t)}','\pi_3^{(t)}','\pi_4^{(t)}', '\pi_5^{(t)}'}), pi_lim.')
legend('$x_0 = 1$','$x_0 = 2$','$x_0 = 3$','$x_0 = 4$','$x_0 = 5$',"Interpreter","latex")
saveas(gcf, "b231.png")

%% Total variation distances and mixing time (numerical approximation)

% We compute transition matrix of the metropolis hasting markov chain,
% initial state is not important, we take 1
% we just take one chain but we take a long time 100000
X = MP_chain_2(N_chain_transition,Time_transition,pi_a,1);
transition_prob = get_transition_matrix(X, N_chain_transition,Time_transition, statespace_size, false);

% We compute the evolution of the distributions on this markov chains for
% each initial states until t = limiting_t
initial_distributions = eye(statespace_size);
limiting_t = 200;
pi0_change_dirac = compute_limiting_dist(initial_distributions, transition_prob, limiting_t, statespace_size);

% We compute the tv distances evolution and the practical mixing time at eps = 5e-3
tv_distances_dirac = compute_tv_distances(pi0_change_dirac);

% We compute numericaly the mixing time at eps = 5e-3
t_epsilon = find(max(tv_distances_dirac,[],2) < eps, 1);

% Compute the evolution of the theoretical upperbound as a function of t
upperbound_tv = zeros(limiting_t, 1);

% Get spectral gap of the chain
[~, spectral_gap] = eig_decomposition(transition_prob);

for i=1:limiting_t
    upperbound_tv(i) = (1-spectral_gap)^i / (sqrt(min(pi_a)));
end

% Find theoretical T_epsilon numerically
t_epsilon_theoretical = find(upperbound_tv < 0.005, 1);

figure()
plot(tv_distances_dirac) 
hold on
plot(upperbound_tv)
xline(t_epsilon,':','$$T_{\epsilon, pr}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xline(t_epsilon_theoretical,':','$$T_{\epsilon, th}$$', 'interpreter', 'latex', 'LabelOrientation', 'horizontal','LabelVerticalAlignment', 'bottom');
xlabel('$$t$$', 'interpreter', 'latex')
ylabel('$$\|\cdot\|_{TV}$$', 'interpreter', 'latex')
legend({'$\pi_{\delta_1}^{(t)}$','$\pi_{\delta_2}^{(t)}$', '$\pi_{\delta_3}^{(t)}$', '$\pi_{\delta_4}^{(t)}$','$\pi_{\delta_5}^{(t)}$','$\max_i \lambda_*^t/(2\sqrt{\pi_i})$'}, 'interpreter', 'latex')
saveas(gcf,'b232.png')

% Display the mixing time
fprintf('Measured mixing time: %d\n',t_epsilon);
fprintf('Theoretical mixing time: %d\n',t_epsilon_theoretical);
fprintf('Spectral gap: %f\n',spectral_gap);
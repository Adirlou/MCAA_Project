function X = MP_chain_1(N_chain,Time,pi_a,x0)

    % Define state space size
    statespace_size = 5;
    
    % Parameters to get a nice approximation of the base chain
    % We chose nchain * time = 500'000;
    nchain = 1;
    time = 500000;
    
    % Initial distribution to compute transition matrix of base chain
    pi0 = ones(1,statespace_size) / statespace_size;
    
    % Get transition matrix of the chain 1
    out = chain_1(nchain, time, pi0);
    psi = get_transition_matrix(out, nchain, time, statespace_size, 0);
    
    % All acceptance probabilities are initialized to 1
    a = ones(statespace_size,statespace_size);
    
    for i=1:statespace_size
        for j=1:statespace_size
            a(i,j) = min(a(i,j),pi_a(j)*psi(j,i)/(pi_a(i)*psi(i,j)));
        end
    end
    
    % Create transition matrix of metropolis hasting algorithm
    P = zeros(statespace_size,statespace_size);
    
    for i=1:statespace_size
        for j=1:statespace_size
            if j ~= i
                P(i,j) = psi(i,j)*a(i,j);
            end
        end
        
        P(i,i) = 1 - sum(P(i,:),2);
    end
    
    mc = dtmc(P);
        
    X = ones(Time,N_chain);
    
    pi0 = zeros(statespace_size,1);
    pi0(x0) = 1;
    
    for n=1:N_chain
        % simulate method returns data X on random walks of length numSteps
        % (2nd argument) through sequences of states in the discrete-time Markov 
        % chain mc from an initial state pi0.
        X(:,n) = simulate(mc,Time-1,'X0',pi0);  
    end
    
end

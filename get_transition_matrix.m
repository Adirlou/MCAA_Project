function [transition_prob] = get_transition_matrix(X, n_chain, time, statespace_size, test_heterogenous) 
    % get_transition_matrix 
    % 
    % Computes the transition matrix of a Markov
    % chain given some of its outputs. 
    
    % If the parameter "test_heterogenous" is set to true, this will
    % compute one transition matrix per timestep.
    
    % If the parameter "test_heterogenous" is set to false, this will
    % compute one transition matrix per chain. Moreover, the average error
    % of the computed transitions will be displayed.
    
    % If we are testing if the matrix is heterogenous
    if test_heterogenous
        
        % A tensor that will contain, mutliple transition matrices for
        % different time t
        transition_prob = zeros(statespace_size, statespace_size, time);
        
        % Update the transition matrices
        for t= 1:time
            for n= 1:n_chain    
                transition_prob(X(t, n), X(t+1, n), t) = transition_prob(X(t, n), X(t+1, n), t) + 1;
            end
    
            % Convert counts to probabilities (normalize the rows)
            transition_prob(:, :, t) = transition_prob(:, :, t)./sum(transition_prob(:, :, t),2);
        end
    else
        % A tensor that will contain, mutliple transition matrices, one for
        % each chain that is run
        transition_prob = zeros(statespace_size, statespace_size, n_chain);
        
        % Array to keep track of how many times we were in a specific state
        state_count = zeros(1, statespace_size);   
        
        for n= 1:n_chain
            
            % For each output we have, update the transition count and
            % state count
            for t= 1:(time-1)
                state_count(X(t,n)) = state_count(X(t,n)) + 1;
                transition_prob(X(t, n), X(t+1, n), n) = transition_prob(X(t, n), X(t+1, n), n) + 1;
            end
            
            % Need to account for an additional timestep for state count
            state_count(X(time,n)) = state_count(X(time,n)) + 1;
            
            % Convert counts to probabilities (normalize the rows)
            transition_prob(:, :, n) = transition_prob(:, :, n)./sum(transition_prob(:, :, n),2);
        end
        
        % Here, we compute the average error of the entries of the
        % transition matrix (using the central limit theorem, we compute
        % the standard deviation
        
        % Compute standard deviation of the transitions
        std_trans = arrayfun(@(x) sqrt(x*(1-x)),transition_prob);
        
        % Take the square root of the state counts
        state_count = arrayfun(@(x) sqrt(x),state_count);
        
        % Divide each row of by the sqrt of the state counts
        final = bsxfun(@rdivide,std_trans, state_count);
        
        % Compute the average error
        avg_error = mean(final, "all");
        
        % Display the error
        fprintf('Average Error: %d\n',avg_error);
    end
end
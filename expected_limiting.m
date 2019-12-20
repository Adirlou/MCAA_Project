function [pi0_change] = expected_limiting(X, n_chain, time, statespace_size) 
    % expected_limiting 
    % 
    % Computes the limiting distribution from the samples 
    
    % Create the tensor that will contain the updates for our initial vectors
    pi0_change = zeros(1, statespace_size);
    
        
    % For each output we have, update the state count
    for n= 1:n_chain
        pi0_change(X(time, n)) = pi0_change(X(time, n)) + 1;
    end
    
    % Convert counts to probabilities (normalize the rows)
    pi0_change = pi0_change./sum(pi0_change);
    
    % Here, we compute the average error of the entries of the
    % limiting distribution (using the central limit theorem, we compute
    % the standard deviation
    
    % Compute standard deviation of the estimated pi_i
    std_pi = arrayfun(@(x) sqrt(x*(1-x)),pi0_change);
    
    % Divide each of them by sqrt(n_chain)
    final = arrayfun(@(x) x/sqrt(n_chain), std_pi);
    
    % Compute the average error
    avg_error = mean(final, "all");
    
    % Display the error
    fprintf('Average Error: %d\n',avg_error);
end
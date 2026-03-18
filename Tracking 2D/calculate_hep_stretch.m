function [probability, distr] = calculate_hep_stretch(M, p, N)
    if p > 1, p = p/100; end
    q = 1 - p;

    % Column-stochastic A: v_{t+1} = A * v_t
    A = zeros(N+1);
    % From S_k, k=0..N-2  (columns are "from", rows are "to")
    for k = 1:N
        A(1,   k) = q;      % to S0 on tail
        A(k+1, k) = p;      % to S_k+1 on head
    end
    
    % Absorbing state
    A(N+1, N+1) = 1;

    start = zeros(N+1,1); start(1) = 1;
    distr = A^M * start;
    probability = distr(end);

    % sanity: columns should sum to 1
    % assert(all(abs(sum(A,1)-1) < 1e-12))
end


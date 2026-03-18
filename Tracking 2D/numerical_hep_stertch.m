function numerical_hep_stertch

% ---- Parameters ----
M = 20;         % length of each sequence
p = 0.5;         % probability of heads
numTrials = 10000;

% ---- Simulate sequences: M x numTrials logical matrix (1=head, 0=tail) ----
H = rand(M, numTrials) < p;

% ---- Compute max run-length of heads for each sequence ----
maxRun = zeros(1, numTrials);
for t = 1:numTrials
    x = H(:, t);
    d = diff([0; x; 0]);          % run starts/ends of ones
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    if isempty(starts)
        maxRun(t) = 0;
    else
        runLens = ends - starts + 1;
        maxRun(t) = max(runLens);
    end
end

% ---- Counts/probabilities for "≥ N heads in a row" for N=1..M ----
histMax = accumarray(maxRun(:) + 1, 1, [M + 1, 1]);  % include 0-bin
cumDesc = cumsum(flipud(histMax(2:end)));            % cumulative from top, skip 0
counts  = flipud(cumDesc);                           % counts(N) = # with maxRun ≥ N
probs   = counts / numTrials;

% ---- Display a quick summary ----
fprintf('Empirical P(≥ N heads in a row) over %d trials (M=%d, p=%.3f):\n', numTrials, M, p);
for N = 1:min(M,10)
    fprintf('  N=%2d: count=%5d, prob=%.4f\n', N, counts(N), probs(N));
end

% Optional: plot
figure; 
plot(1:M, probs, 'LineWidth', 1.5);
xlabel('N (run length threshold)'); ylabel('Empirical P(\geq N heads)');
title(sprintf('Runs in %d sequences (M=%d, p=%.2f)', numTrials, M, p));
grid on;
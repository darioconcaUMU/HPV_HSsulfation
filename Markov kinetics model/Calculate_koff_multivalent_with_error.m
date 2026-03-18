
prob = [ ...
    0.270600109  0.270670548  0.180494004  0.090270494  0.036117597;
    0.223220635  0.263062460  0.206677004  0.121783009  0.057407843;
    0.186919266  0.247913839  0.219207908  0.145369376  0.077122238;
    0.366850062  0.170034183  0.052540309  0.012176158  0.002257449
];

toff = [ ...
    2.50E+02  3.59E+05  6.89E+08  1.49E+12  3.42E+15;
    6.67E+01  4.05E+04  3.27E+07  2.98E+10  2.89E+13;
    8.33E+00  9.22E+02  1.35E+05  2.24E+07  3.94E+09;
    8.20E-02  1.09E-01  1.48E-01  2.05E-01  2.89E-01
];

toff_lo = [ ...
    142.857142900000  16599.0952200000  2560594.53400000  444366891         82256695859.0000;
    52.6000000000000  3801.45206100000  363576.199400000  39117679.5500000  4489306967.00000;
    6.25000000000000  77                1214.18440400000  21516.2018400000  406672.764900000;
    0.0575000000000000 0.0596000000000000 0.0618000000000000 0.0641000000000000 0.0665291600000000
];

toff_hi = [ ...
    1000               1484017.76300000  2935422615.00000   6532130000000.00  1.55048000000000e+16;
    90.9090909100000  18731.7382100000  5133753.38400000   1582857761.00000  520568000000.000;
    12.5000000000000  542.023333700000  30980.3425100000   1991816.09100000  136596731.800000;
    0.142857143000000 0.162765739000000 0.186373621000000  0.214454091000000 0.247952875000000
];

frame_rate = 1/15;
dt         = 1/frame_rate;
max_time   = 600;

t0    = 3*dt;                     % 45 s
time0 = 0:dt:max_time;
idx   = time0 >= t0;              % display window
t     = time0(idx);

nMC = 5000;                       % Monte Carlo samples (adjust)
alpha = 0.05;                     % 95% band

figure; hold on

for j = 1:size(prob,1)

    p = prob(j,:);

    % --- Deterministic envelope (guaranteed band)
    D_lo_env = zeros(size(time0));
    D_hi_env = zeros(size(time0));
    for i = 1:size(prob,2)
        D_lo_env = D_lo_env + p(i) * exp(-time0 ./ toff_lo(j,i));
        D_hi_env = D_hi_env + p(i) * exp(-time0 ./ toff_hi(j,i));
    end

    % --- Monte Carlo: choose distribution inside bounds
    % Log-uniform sampling: log(tau) ~ Uniform(log(lo), log(hi))
    logLo = log(toff_lo(j,:));
    logHi = log(toff_hi(j,:));

    % Sample taus: (nMC x nComp)
    taus = exp(logLo + rand(nMC, numel(p)) .* (logHi - logLo));

    % Compute D(t) for each sample (vectorized)
    % Dk(t) = sum_i p_i * exp(-t/tau_i)
    Dmc = zeros(nMC, numel(time0));
    for i = 1:numel(p)
        Dmc = Dmc + p(i) .* exp(-time0 ./ taus(:,i));
    end

    % --- Normalization at first available point (t0), using MEAN curve baseline
    % If you want normalization based on each sample's value at t0, replace this normFactor.
    % Here we normalize everything to the "mean-parameter" curve at t0 if available,
    % otherwise to the midpoint-parameter curve.
    if exist('toff','var') && ~isempty(toff)
        D_mean = zeros(size(time0));
        for i = 1:numel(p)
            D_mean = D_mean + p(i) * exp(-time0 ./ toff(j,i));
        end
    else
        tau_mid = sqrt(toff_lo(j,:) .* toff_hi(j,:)); % geometric midpoint
        D_mean = zeros(size(time0));
        for i = 1:numel(p)
            D_mean = D_mean + p(i) * exp(-time0 ./ tau_mid(i));
        end
    end

    [~, i0] = min(abs(time0 - t0));
    normFactor = max(D_mean(i0), eps);

    D_mean     = D_mean     / normFactor;
    D_lo_env   = D_lo_env   / normFactor;
    D_hi_env   = D_hi_env   / normFactor;
    Dmc_norm   = Dmc        / normFactor;

    % --- Percentile band from MC
    loP = prctile(Dmc_norm(:,idx), 100*(alpha/2), 1);
    hiP = prctile(Dmc_norm(:,idx), 100*(1-alpha/2), 1);

    % --- Plot: mean + MC bounds as dashed (same color)
    h = plot(t, D_mean(idx), 'LineWidth', 2);
    c = h.Color;

    plot(t, loP, '--', 'Color', c, 'LineWidth', 1.2);
    plot(t, hiP, '--', 'Color', c, 'LineWidth', 1.2);

    % Optional: also plot the guaranteed envelope as dotted
    % plot(t, D_lo_env(idx), ':', 'Color', c, 'LineWidth', 1.0);
    % plot(t, D_hi_env(idx), ':', 'Color', c, 'LineWidth', 1.0);

end

xlabel('Time (s)')
ylabel('Normalized bound fraction')
xlim([t0 max_time])
ylim([0 1.05])
grid on
box on

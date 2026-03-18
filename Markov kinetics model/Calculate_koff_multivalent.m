function [k,A,C,data] = Calculate_koff_multivalent

frame_rate = 1;%1/15;
max_time   = 450;%1800;

time = 1/frame_rate*3 : 1/frame_rate : max_time;

data(:,1) = time;

A = 1726/1e6;
mu = [1.1 1.3 1.5 1.6]*1e3.*A;

for i=1:4
    prob(i,:) = poisspdf(1:5, mu(i));
end


% prob = [0.270600109 0.270670548 0.180494004 0.090270494 0.036117597;
%         0.223220635 0.263062460 0.206677004 0.121783009 0.057407843;
%         0.186919266 0.247913839 0.219207908 0.145369376 0.077122238;
%         0.366850062 0.170034183 0.052540309 0.012176158 0.002257449];

toff = [ ...
    2.50E+02  7.18E+04  2.74E+07  1.18E+10  5.41E+12;
    6.67E+01  8.09E+03  1.30E+06  2.36E+08  4.56E+10;
    8.33E+00  1.89E+02  5.58E+03  1.86E+05  6.58E+06;
    8.20E-02  8.74E-02  9.33E-02  9.97E-02  1.07E-01
];

toff2 = [ ...
    1.43E+02  1.66E+04  2.56E+06  4.44E+08  8.23E+10;
    5.26E+01  3.80E+03  3.64E+05  3.91E+07  4.49E+09;
    6.25E+00  7.70E+01  1.21E+03  2.15E+04  4.07E+05;
    5.75E-02  5.96E-02  6.18E-02  6.41E-02  6.65E-02
];

toff3 = [ ...
    1.00E+03  1.48E+06  2.94E+09  6.53E+12  1.55E+16;
    9.09E+01  1.87E+04  5.13E+06  1.58E+09  5.21E+11;
    1.25E+01  5.42E+02  3.10E+04  1.99E+06  1.37E+08;
    1.43E-01  1.63E-01  1.86E-01  2.14E-01  2.48E-01
];

figure; hold on

col_scale = [ ...
    0   0   0;
    140 125  55;
    255 160  64;
    255   0   0] / 255;

hLine = gobjects(size(prob,1),1);   % store line handles

k  = zeros(size(prob,1),3);
A  = zeros(size(prob,1),3);
C  = zeros(size(prob,1),3);

for j = 1:size(prob,1)

    nComp = size(prob,2);

    temp_dissociation  = zeros(nComp, numel(time));
    temp_dissociation2 = zeros(nComp, numel(time));
    temp_dissociation3 = zeros(nComp, numel(time));

    for i = 1:nComp
        temp_dissociation(i,:)  = prob(j,i) * exp(-time ./ toff(j,i));
        temp_dissociation2(i,:) = prob(j,i) * exp(-time ./ toff2(j,i));
        temp_dissociation3(i,:) = prob(j,i) * exp(-time ./ toff3(j,i));
    end

    dissociation  = sum(temp_dissociation,  1);
    dissociation2 = sum(temp_dissociation2, 1);
    dissociation3 = sum(temp_dissociation3, 1);

    normFactor = max(dissociation(1), eps);
    dissociation  = dissociation; % / normFactor;
    dissociation2 = dissociation2; %/ normFactor;
    dissociation3 = dissociation3; %/ normFactor;

    [k(j,1), A(j,1), C(j,1)] = fit_single_exp_offset(time, dissociation);
    [k(j,2), A(j,2), C(j,2)]  = fit_single_exp_offset(time, dissociation2);
    [k(j,3), A(j,3), C(j,3)] = fit_single_exp_offset(time, dissociation3);

    dMin = min(dissociation2, dissociation3);
    dMax = max(dissociation2, dissociation3);

    % --- shaded area (exclude from legend)
    hFill = fill([time fliplr(time)], ...
                 [dMin fliplr(dMax)], ...
                 col_scale(j,:), ...
                 'FaceAlpha', 0.20, ...
                 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');

    uistack(hFill, 'bottom');

    % --- mean line (store handle)
    hLine(j) = plot(time, dissociation, ...
                    'Color', col_scale(j,:), ...
                    'LineWidth', 2);
    data(:,j+1) = dissociation;
end

xlabel('Time (s)')
ylabel('Normalized bound fraction')
box on

legend(hLine, {'Heparin','2-O-deS','6-O-deS','N-deS'}, 'Location', 'best','Box','off')

end



function [k, A, C] = fit_single_exp_offset(time, y)
% Fit y = A * exp(-k*t) + C

    % ensure column vectors
    time = time(:);
    y    = y(:);

    % remove NaNs / infs
    idx = isfinite(time) & isfinite(y);
    time = time(idx);
    y    = y(idx);

    % --- initial guesses (important for stability)
    C0 = min(y(end-5:end));                 % late-time plateau
    A0 = max(y(1) - C0, eps);               % amplitude
    k0 = 1 / max(time(end) - time(1), eps); % rough rate

    % --- nonlinear least squares
    model = @(p,t) p(1) * exp(-p(2)*t) + p(3);   % [A k C]

    p0 = [A0, k0, C0];
    lb = [0,   0,  -inf];
    ub = [inf, inf, inf];

    opts = optimoptions('lsqcurvefit', ...
        'Display','off');

    p = lsqcurvefit(model, p0, time, y, lb, ub, opts);

    A = p(1);
    k = p(2);
    C = p(3);
end
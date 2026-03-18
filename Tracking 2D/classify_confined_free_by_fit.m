function [labels, fit_quality] = classify_confined_free_by_fit(properties, doPlot, plotType)
% Model-comparison classifier for free/confined from MSD curves
% Optionally plot only a given class (plotType)
% Usage: [labels, fit_quality] = classify_confined_free_by_fit(properties, doPlot, plotType)
% doPlot: true/false
% plotType: 'all' (default), 'free', 'confined', 'uncertain', 'undetermined'

if nargin < 2 || isempty(doPlot), doPlot = true; end
if nargin < 3 || isempty(plotType), plotType = 'all'; end

nPerFigure = 16;
r2_thresh_free = 0.9;
r2_thresh_conf = 0.8;
% deltaR2   = 0.03;

nPerFigure = (round(sqrt(nPerFigure)))^2;

nTracks = numel(properties);
labels = cell(nTracks,1);
fit_quality = struct('R2_linear',nan(nTracks,1),'R2_exp',nan(nTracks,1));

if doPlot
    plotIndices = 1;
end

for i = 1:nTracks
    msd = properties(i).msd(:);
    lags = (1:length(msd))';
    if length(msd)<6 || any(isnan(msd))
        labels{i} = 'undetermined';
        continue
    end

    % Linear fit
    pf_lin = polyfit(lags, msd, 1);
    fit_lin = polyval(pf_lin, lags);
    SSres_lin = sum((msd-fit_lin).^2);
    SStot = sum((msd-mean(msd)).^2);
    R2_lin = 1-SSres_lin/SStot;

    % Exponential fit
    exp_model = @(b, x) b(1)*(1-exp(-x/b(2))) + b(3);
    b0 = [max(msd), round(length(msd)/3), min(msd)];
    opts = statset('nlinfit');
    opts.RobustWgtFun = 'bisquare';
    try
        b = nlinfit(lags, msd, exp_model, b0, opts);
        fit_exp = exp_model(b, lags);
        SSres_exp = sum((msd-fit_exp).^2);
        R2_exp = 1-SSres_exp/SStot;
    catch
        R2_exp = -Inf;
        fit_exp = nan(size(lags));
    end

    fit_quality.R2_linear(i) = R2_lin;
    fit_quality.R2_exp(i)    = R2_exp;

    % ---- Decision logic ----
    if R2_exp < r2_thresh_conf && R2_lin < r2_thresh_free
        labels{i} = 'uncertain';
    elseif (R2_lin < r2_thresh_free) && (R2_exp > r2_thresh_conf)
        labels{i} = 'confined';
    elseif (R2_lin > r2_thresh_free)
        labels{i} = 'free';
    else
        labels{i} = 'uncertain';
    end

    % -------- Plotting block (filtered) --------
    if doPlot
        % Plot only if class matches or if 'all'
        if strcmpi(plotType,'all') || strcmpi(plotType,labels{i})
            if mod(plotIndices-1, nPerFigure) == 0
                figure('Name','MSD Model Comparison','Color','w');
                tiledlayout(sqrt(nPerFigure),sqrt(nPerFigure), 'Padding','compact','TileSpacing','compact');
            end
            nexttile;
            plot(lags, msd, 'ko-', 'LineWidth',1.2); hold on;
            plot(lags, fit_lin, 'b--', 'LineWidth',1.2);
            plot(lags, fit_exp, 'r-', 'LineWidth',1.2);
            legend({'MSD','Linear fit','Exp fit'},'Location','northwest','FontSize',8);
            title({['Track ' num2str(i) ' (' labels{i} ')'], ...
                sprintf('R^2_{lin}=%.2f, R^2_{exp}=%.2f', R2_lin, R2_exp)},'FontSize',9);
            xlabel('Lag'); ylabel('MSD');
            set(gca,'FontSize',9);

            if mod(plotIndices, nPerFigure) == 0
                drawnow;
                disp('Press any key or click in the figure to continue...');
                waitforbuttonpress;
            end
            plotIndices = plotIndices + 1;
        end
    end
end

% Final pause if a page was partially filled
if doPlot && mod(plotIndices-1, nPerFigure) ~= 0
    drawnow;
    disp('Press any key or click in the figure to continue...');
    waitforbuttonpress;
end

% Show summary
confinedCount = sum(strcmp(labels,'confined'));
freeCount     = sum(strcmp(labels,'free'));
uncertainCount = sum(strcmp(labels,'uncertain'));
undetCount = sum(strcmp(labels,'undetermined'));
fprintf('Confined: %d, Free: %d, Uncertain: %d, Undetermined: %d\n', ...
    confinedCount, freeCount, uncertainCount, undetCount);

end

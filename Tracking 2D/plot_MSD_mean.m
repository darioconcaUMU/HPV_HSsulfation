function plot_MSD_mean(properties, method)
% plot_MSD_mean_sem(properties, method)
% Plots mean and SEM of MSDs from structure array 'properties' (field 'MSD')
% Also fits for D (either to mean MSD or average of individual fits)
%
% method = 'mean' (default): fit mean MSD
% method = 'individual': fit each MSD, then average D

if nargin<2, method = 'mean'; end

nTracks = numel(properties);
msdTracks = cell(nTracks,1);

for i=1:nTracks
    msdTracks{i} = properties(i).msd(:); % ensure column
end

maxLag = max(cellfun(@length, msdTracks));
allMSD = nan(nTracks, maxLag);

for i=1:nTracks
    len = length(msdTracks{i});
    allMSD(i, 1:len) = msdTracks{i};
end

meanMSD = nanmean(allMSD,1);
semMSD  = nanstd(allMSD,0,1) ./ sqrt(sum(~isnan(allMSD),1));
lags = (1:maxLag)'; % row vector

% ----- Plot -----
figure; hold on;

% Overlay individual MSDs as faint lines (optional)
h1 = plot(lags, allMSD', 'Color', [0.7 0.7 0.7 0.3]); % faint gray lines

% Shaded SEM
h2 = fill([lags; flipud(lags)], ...
     [meanMSD+semMSD, fliplr(meanMSD-semMSD)], ...
     [0.3 0.7 1], 'FaceAlpha',0.3, 'EdgeColor','none');

% Mean MSD line
h3 = plot(lags, meanMSD, 'b-', 'LineWidth', 2);

% -------------- Diffusion Coefficient Extraction --------------
switch lower(method)
    case 'mean'
        % Fit mean MSD
        fitLags = ~isnan(meanMSD);
        coeffs = polyfit(lags(fitLags), meanMSD(fitLags), 1); % linear fit
        D = coeffs(1)/4; % for 2D: MSD = 4Dt + c
        h4 = plot(lags, polyval(coeffs, lags), 'r--', 'LineWidth', 2);
        legend([h1(1), h2, h3, h4], ...
            {'Individual MSDs','Mean \pm SEM','Mean MSD','Fit to mean'}, ...
            'Location', 'northwest');
        txt = sprintf('D from fit to mean MSD:\nD = %.3g', D);

    case 'individual'
        Ds = nan(nTracks,1);
        fit_lag_min = 1; % start fitting from first lag
        for i=1:nTracks
            thisMSD = allMSD(i,:)';
            idx = ~isnan(thisMSD);
            if sum(idx)>1
                c = polyfit(lags(idx), thisMSD(idx), 1);
                Ds(i) = c(1)/4;
            end
        end
        D = nanmean(Ds);
        % Optional: plot mean linear fit for visual clarity
        h4 = plot(lags, 4*D*lags + nanmean(meanMSD - 4*D*lags'), 'r--', 'LineWidth', 2);
        legend([h1(1), h2, h3, h4], ...
            {'Individual MSDs','Mean \pm SEM','Mean MSD','Avg. fit'}, ...
            'Location', 'northwest');
        txt = sprintf('The mean D from the individual fit:\nD = %.3g', D);

    otherwise
        error('Unknown method. Use \"mean\" or \"individual\".');
end

xlabel('Lag');
ylabel('MSD');
title('Mean MSD with SEM and Linear Fit');

% Get axes limits for top-right placement
xl = xlim;
yl = ylim;
x_ann = xl(2)*0.9; % rightmost x
y_ann = yl(2)*0.9; % top y

text(x_ann, y_ann, txt, ...
    'BackgroundColor','w','EdgeColor','k','Margin',5, ...
    'VerticalAlignment','top', 'HorizontalAlignment','right');

box on; hold off;

end

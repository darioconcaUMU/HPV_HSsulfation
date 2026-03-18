function histogram_D_from_MSD(properties)
% histogram_D_from_MSD(properties)
% Fits each MSD, extracts D, plots histogram of D values (log scale)
% Overlays lines for: mean/median of log(D)s and log(D) from mean MSD
% Writes all three values in a box in the top right

nTracks = numel(properties);
Ds = nan(nTracks,1);
maxLag = max(arrayfun(@(p) numel(p.msd), properties));
allMSD = nan(nTracks, maxLag);

% Fit each MSD individually
for i = 1:nTracks
    msd = properties(i).msd(:);    % Ensure column
    lags = (1:length(msd))';
    allMSD(i,1:length(msd)) = msd;
    if sum(~isnan(msd)) > 1
        coeffs = polyfit(lags, msd, 1);
        Ds(i) = coeffs(1)/4; % 2D: slope = 4D
    end
end

% Remove bad fits
Ds_log = log(Ds(Ds>0 & ~isnan(Ds)));

% Fit the mean MSD
meanMSD = nanmean(allMSD,1);
meanLags = (1:length(meanMSD))';
idx = ~isnan(meanMSD);
meanFit = polyfit(meanLags(idx), meanMSD(idx)', 1);
D_meanMSD = meanFit(1)/4;

% Values for the annotation box
mean_Ds_log   = mean(Ds_log);
median_Ds_log = median(Ds_log);
log_D_meanMSD = log(D_meanMSD);

% ---- Plot histogram and vertical lines ----
figure;
histogram(Ds_log, 'FaceColor', [0.3 0.7 1], 'EdgeColor', 'k');
xlabel('Ln(D (diffusion coefficient))');
ylabel('Count');
title('Distribution of D from Individual MSD fits');
hold on

yl = ylim;
xline(mean_Ds_log,   'r--', 'LineWidth',2);
xline(median_Ds_log, 'g-.', 'LineWidth',2);
xline(log_D_meanMSD, 'b-',  'LineWidth',2);

legend({'D distribution','Mean of Ds','Median of Ds','D from mean MSD'}, 'Location','best');

% --- Annotation box in the top right ---
drawnow; % Make sure axes limits are up to date
xl = xlim; yl = ylim;
x_box = xl(2) - 0.02*(xl(2)-xl(1)); % 2% from right edge
y_box = yl(2) - 0.02*(yl(2)-yl(1)); % 2% from top

annotation_txt = sprintf([ ...
    'Mean D:   %.3g\n' ...
    'Median D: %.3g\n' ...
    'D_{mean MSD}: %.3g'], ...
    exp(mean_Ds_log), exp(median_Ds_log), D_meanMSD);

text(x_box, y_box, annotation_txt, ...
    'HorizontalAlignment','right', 'VerticalAlignment','top', ...
    'FontSize',10, ...
    'BackgroundColor','w','EdgeColor','k','Margin',6);

hold off

disp(['Mean D: ' num2str(exp(mean_Ds_log))]);
disp(['Median D: ' num2str(exp(median_Ds_log))]);
disp(['D from mean MSD: ' num2str(D_meanMSD)]);

end

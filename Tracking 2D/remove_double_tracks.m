function [tracks_kept, removed_idx, report] = remove_double_tracks(tracks, opts)
% REMOVE_DOUBLE_TRACKS  Remove tracks that bounce between two distant positions.
% tracks: cell array, each cell is [N x 3] = [time, x, y]
%
% Criteria (all must hold to remove):
%   • Separation:  d12 > sepFactor * r_intra   and  d12 > absSepMin
%   • Occupancy:   min(cluster_fraction) >= minFrac
%   • Switching:   switchFrac >= switchThresh  and  numBlocks >= minBlocks
%
% Inputs:
%   tracks : cell array of [time, x, y]
%   opts   : (optional) struct with parameters:
%       .sepFactor     (default 5)     % required separation multiple
%       .absSepMin     (default 0)     % absolute minimum separation (units of x/y)
%       .minFrac       (default 0.2)   % min fraction of points in each cluster
%       .switchThresh  (default 0.5)   % fraction of time-steps that switch clusters
%       .minBlocks     (default 4)     % min # of alternating blocks (run-lengths)
%       .minPoints     (default 8)     % min points in a track to test
%       .plotMax       (default 0)     % # removed tracks to plot (0 = none)
%
% Outputs:
%   tracks_kept : same as input but with bouncing tracks removed
%   removed_idx : logical vector of removed tracks (length = numel(tracks))
%   report      : struct with per-track diagnostics

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'sepFactor'),    opts.sepFactor = 5;    end
    if ~isfield(opts,'absSepMin'),    opts.absSepMin = 0;    end
    if ~isfield(opts,'minFrac'),      opts.minFrac = 0.10;   end
    if ~isfield(opts,'switchThresh'), opts.switchThresh = 0.10; end
    if ~isfield(opts,'minBlocks'),    opts.minBlocks = 4;    end
    if ~isfield(opts,'minPoints'),    opts.minPoints = 8;    end
    if ~isfield(opts,'plotMax'),      opts.plotMax = 10;      end

    nT = numel(tracks);
    removed_idx = false(nT,1);

    % diagnostics (optional)
    d12_all = nan(nT,1);
    r_intra_all = nan(nT,1);
    fracMin_all = nan(nT,1);
    switchFrac_all = nan(nT,1);
    numBlocks_all = nan(nT,1);

    nPlotted = 0;

    for i = 1:nT

        
        T = tracks{i};
        if isempty(T) || size(T,2) < 3 || size(T,1) < opts.minPoints
            continue
        end
        X = [T(:,2), T(:,3)];
        n = size(X,1);

        % --- 2-cluster fit
        try
            [lbl, C] = kmeans(X, 2, 'Replicates', 5, 'MaxIter', 200, 'Display','off');
        catch
            [lbl, C] = kmeans(X, 2, 'Replicates', 5, 'MaxIter', 200);
        end

        % cluster sizes and fractions
        n1 = sum(lbl==1); n2 = sum(lbl==2);
        f1 = n1/n; f2 = n2/n;
        fracMin = min(f1, f2);

        % robust within-cluster spread: median radius to centroid
        r1 = median( sqrt(sum((X(lbl==1,:) - C(1,:)).^2, 2)) );
        r2 = median( sqrt(sum((X(lbl==2,:) - C(2,:)).^2, 2)) );
        r_intra = max([r1, r2, realmin]);   % guard from zero
        d12 = norm(C(1,:) - C(2,:));

        % switching along time (labels in observation order)
        switches = sum(lbl(2:end) ~= lbl(1:end-1));
        switchFrac = switches / (n-1);

        % number of alternating blocks (run-length encoding of labels)
        numBlocks = 1 + switches;

        % store diagnostics
        d12_all(i) = d12;
        r_intra_all(i) = r_intra;
        fracMin_all(i) = fracMin;
        switchFrac_all(i) = switchFrac;
        numBlocks_all(i) = numBlocks;

        % --- criteria
        condSep    = (d12 > opts.sepFactor * r_intra) && (d12 > opts.absSepMin);
        condOcc    = (fracMin >= opts.minFrac);
        %condSwitch = (switchFrac >= opts.switchThresh) && (numBlocks >= opts.minBlocks);
        condSwitch = 1;

        if condSep && condOcc && condSwitch
            removed_idx(i) = true;

            % optional plot for removed track
            if nPlotted < opts.plotMax
                nPlotted = nPlotted + 1;
                % center for visualization
                cx = mean(X(:,1),'omitnan'); cy = mean(X(:,2),'omitnan');
                Xc = X(:,1) - cx; Yc = X(:,2) - cy;

                figure('Name',sprintf('Removed track %d (d/ri=%.2f, fracMin=%.2f, switch=%.2f)', ...
                        i, d12/r_intra, fracMin, switchFrac));
                hold on; axis equal; grid on
                % plot by cluster
                plot(Xc(lbl==1), Yc(lbl==1), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
                plot(Xc(lbl==2), Yc(lbl==2), '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.2);
                % centers
                plot(C(1,1)-cx, C(1,2)-cy, 'o', 'MarkerSize', 7, 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k');
                plot(C(2,1)-cx, C(2,2)-cy, 'o', 'MarkerSize', 7, 'MarkerFaceColor',[0.00 0.45 0.74], 'MarkerEdgeColor','k');
                % separation annotation
                title(sprintf('Removed %d | d_{12}=%.3g, r_{intra}=%.3g, d/r=%.2f, minFrac=%.2f, switches=%.0f/%d', ...
                      i, d12, r_intra, d12/r_intra, fracMin, switches, n-1));
                legend({'cluster 1','cluster 2','c1','c2'}, 'Location','best');
                xlabel('x (centered)'); ylabel('y (centered)');
                hold off
            end
        end
    end

    tracks_kept = tracks(~removed_idx);

    report = struct( ...
        'sepFactor', opts.sepFactor, ...
        'absSepMin', opts.absSepMin, ...
        'minFrac', opts.minFrac, ...
        'switchThresh', opts.switchThresh, ...
        'minBlocks', opts.minBlocks, ...
        'minPoints', opts.minPoints, ...
        'removed_idx', removed_idx, ...
        'd12', d12_all, ...
        'r_intra', r_intra_all, ...
        'fracMin', fracMin_all, ...
        'switchFrac', switchFrac_all, ...
        'numBlocks', numBlocks_all );
end

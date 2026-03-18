function [tracks_clean, info] = clean_outliers(tracks, opts)
% CLEAN_OUTLIERS
% Adaptive, radius-based outlier removal with endpoint handling and iteration.
% - Flag tracks with z = (maxR-meanR)/meanR above a global z-threshold (mu + kSigma*std).
% - For flagged tracks:
%     • If first/last frame is beyond rAccept: REMOVE that frame (delete row).
%     • If interior point beyond rAccept and (optionally) neighbors are inside: replace by neighbors' mean.
%     • If total fixes (endpoint removals + interior replacements) > maxCorrections: remove track.
% - Plots each flagged track centered at (0,0) with a circle of radius rAccept; marks corrected and removed endpoints.
%
% INPUT
%   tracks : cell array; each cell is [N x 3] = [time, x, y]
%   opts   : (optional) struct:
%       .kSigma               (default 1)
%       .centerMode           (default 'mean')   % 'mean' or 'first'
%       .maxCorrections       (default 5)        % > this -> remove track
%       .maxItersSafe         (default 50)
%       .plotMax              (default inf)
%       .requireNeighborsInside (default true)   % only fix interior if neighbors are inside
%
% OUTPUT
%   tracks_clean : cleaned tracks, removed ones pruned
%   info         : struct with metrics/flags

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'kSigma'),                 opts.kSigma = 1; end
    if ~isfield(opts,'centerMode'),             opts.centerMode = 'mean'; end
    if ~isfield(opts,'maxCorrections'),         opts.maxCorrections = 5; end
    if ~isfield(opts,'maxItersSafe'),           opts.maxItersSafe = 50; end
    if ~isfield(opts,'plotMax'),                opts.plotMax = inf; end
    if ~isfield(opts,'requireNeighborsInside'), opts.requireNeighborsInside = true; end

    nT = numel(tracks);
    tracks_clean = tracks;

    meanR = nan(nT,1); maxR = nan(nT,1); z = nan(nT,1); nPts = zeros(nT,1);
    centers = nan(nT,2);

    % ---------- Pass 1: per-track center, radii, z ----------
    for i = 1:nT
        T = tracks{i};
        if isempty(T) || size(T,2) < 3 || size(T,1) < 2, continue; end

        switch lower(opts.centerMode)
            case 'first'
                cx = T(1,2); cy = T(1,3);
            otherwise % 'mean'
                cx = mean(T(:,2),'omitnan'); cy = mean(T(:,3),'omitnan');
        end
        centers(i,:) = [cx, cy];

        r = hypot(T(:,2)-cx, T(:,3)-cy);
        mR = mean(r,'omitnan'); aR = max(r,[],'omitnan');
        nPts(i) = numel(r);
        meanR(i) = mR; maxR(i) = aR;

        if ~isfinite(mR) || mR == 0
            z(i) = 0;
        else
            z(i) = (aR - mR) / mR;
        end
    end

    good = ~isnan(z);
    muZ = mean(z(good));
    sdZ = std(z(good),0);
    if sdZ == 0 || ~isfinite(sdZ), sdZ = eps; end
    zThr = muZ + opts.kSigma * sdZ;

    outlierTrack = (z > zThr) & (nPts >= 2);   % need at least 2 to do anything

    removed_idx = false(nT,1);
    correctedIdx = cell(nT,1);
    removedEndPts = cell(nT,1);   % coordinates of removed endpoints (centered), for plotting
    rAccept = nan(nT,1);

    % ---------- Pass 2: iterate fixes on flagged tracks ----------
    nPlotted = 0;
    for i = find(outlierTrack).'
        Jend = [];
        Jcand = [];
        T0 = tracks{i};
        t = T0(:,1); x = T0(:,2); y = T0(:,3);
        cx = centers(i,1); cy = centers(i,2);

        rAcc = meanR(i) * (1 + zThr);
        rAccept(i) = rAcc;

        totalFixes = 0; removed = false; iter = 0;

        while iter < opts.maxItersSafe
            iter = iter + 1;

            n = numel(x);
            if n < 2, removed = true; break; end

            r = hypot(x - cx, y - cy);
            Jall = find(r > rAcc).';           % all outlier indices (including endpoints)

            if isempty(Jall), break; end

            % --- Handle endpoint outliers: remove the frame(s) ---
            endMask = (Jall == 1) | (Jall == n);
            if any(endMask)
                Jend = Jall(endMask);

                % store endpoint coords (centered) before removal for plotting
                removedEndPts{i} = [removedEndPts{i}; [x(Jend)-cx, y(Jend)-cy]];

                % remove rows
                t(Jend) = [];
                x(Jend) = [];
                y(Jend) = [];

                totalFixes = totalFixes + numel(Jend);
                if totalFixes > opts.maxCorrections
                    removed = true; break;
                end

                % after removal, re-loop (indices shifted)
                continue;
            end

            % --- Interior outliers: optionally require neighbors inside; then average neighbors ---
            Jcand = Jall(Jall > 1 & Jall < n);
            if isempty(Jcand), break; end

            if opts.requireNeighborsInside
                kLow = ones(numel(Jcand),1)';
                kHigh = kLow;
                okNLow = (r(Jcand-kLow) <= rAcc);
                okNHigh = (r(Jcand+kHigh) <= rAcc);
                while sum(okNLow)<numel(okNLow)
                    kLow = kLow+~okNLow';
                    okNLow = (r(Jcand-kLow) <= rAcc);
                end
                while sum(okNHigh)<numel(okNHigh)
                    kHigh = kHigh+~okNHigh';
                    okNHigh = (r(Jcand+kHigh) <= rAcc);
                end
                J = Jcand;
            else
                J = Jcand;
                kLow = 1;
                kHigh = 1;
            end

            if isempty(J)
                % no interior points that pass neighbor test
                break;
            end

            % neighbor-average correction
            x(J) = 0.5*(x(J-kLow) + x(J+kHigh));
            y(J) = 0.5*(y(J-kLow) + y(J+kHigh));
            totalFixes = totalFixes + numel(J);
            correctedIdx{i} = unique([correctedIdx{i}; J(:)]); %#ok<AGROW>

            if totalFixes > opts.maxCorrections
                removed = true; break;
            end
        end

        % ---------- Plot: centered original + circle + cleaned / removed endpoints ----------
        if nPlotted < opts.plotMax
            nPlotted = nPlotted + 1;

            X0 = T0(:,2) - cx; Y0 = T0(:,3) - cy;      % original centered
            figure('Name',sprintf('Track %d (z=%.3g, thr=%.3g, r_{acc}=%.3g)', i, z(i), zThr, rAcc));
            hold on; axis equal; grid on
            th = linspace(0,2*pi,256);
            plot(rAcc*cos(th), rAcc*sin(th), 'k--', 'LineWidth', 1.2); % threshold circle
            plot(X0, Y0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1); % original

            if ~removed
                Xc = x - cx; Yc = y - cy;
                plot(Xc, Yc, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
                % corrected interior points
                J = correctedIdx{i};
                if ~isempty(J)
                    plot(Xc(J), Yc(J), 'o', 'MarkerSize', 6, ...
                        'MarkerFaceColor','y','MarkerEdgeColor','k');
                end
            end

            % removed endpoints (from original coordinates)
            if ~isempty(removedEndPts{i})
                pts = removedEndPts{i};
                plot(pts(:,1), pts(:,2), 'x', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, 'MarkerSize', 8);
            end

            xlabel('x (centered)'); ylabel('y (centered)');
            if ~removed
                ttl = sprintf('Track %d: fixes=%d (<=%d OK)', i, totalFixes, opts.maxCorrections);
            else
                ttl = sprintf('Track %d: fixes=%d > %d — REMOVED', i, totalFixes, opts.maxCorrections);
            end
            title(ttl);
            legend_name = {'max accepted radius','original','cleaned'};
            if ~isempty(Jcand)
                temp_size = numel(legend_name)+1;
                legend_name{temp_size} = 'corrected pts';
            end
            if ~isempty(Jend)
                temp_size = numel(legend_name)+1;
                legend_name{temp_size} = 'removed endpoints';
            end 

            legend(legend_name, 'Location','best');
                   
            hold off
        end

        if removed || numel(x) < 2
            tracks_clean{i} = [];      % mark for pruning
            removed_idx(i) = true;
        else
            tracks_clean{i} = [t, x, y];
        end
    end

    % prune removed tracks
    tracks_clean = tracks_clean(~removed_idx);

    % ---------- info ----------
    info = struct('meanR',meanR, 'maxR',maxR, 'z',z, ...
                  'muZ',muZ, 'sdZ',sdZ, 'zThr',zThr, ...
                  'outlierTrack',outlierTrack, 'rAccept',rAccept, ...
                  'centers',centers, 'correctedIdx',{correctedIdx}, ...
                  'removedEndPts',{removedEndPts}, ...
                  'removed_idx',removed_idx, 'params',opts);
end
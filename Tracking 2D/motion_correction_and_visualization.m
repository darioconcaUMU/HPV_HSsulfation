function [corrected_tracks, avg_motion, idx_lowmob, idx_full_corr] = ...
    motion_correction_and_visualization(tracks, mobility_struct, mobility_thresh, varargin)
% Tracks: cell array, each [N x 2]
% mobility_struct: struct array with field .maxDisp
% mobility_thresh: threshold for reference (low-mobility) tracks
%
% Optional name-value pairs for plotting:
%   'plot_drift'      (default: true)
%   'plot_corrected'  (default: true)
%   'plot_16corr'     (default: true)
%   'plot_sidebyside' (default: true)

% ------- Parse options -------
p = inputParser;
addParameter(p, 'plot_drift', true);
addParameter(p, 'plot_corrected', true);
addParameter(p, 'plot_16corr', true);
addParameter(p, 'plot_sidebyside', true);
parse(p, varargin{:});
opts = p.Results;

% ------- 1. Find reference tracks (full-length & low-mobility) -------
track_lengths = cellfun(@(trk) size(trk,1), tracks);
max_length = max(track_lengths);
idx_full = find(track_lengths == max_length);
idx_lowmob = idx_full([mobility_struct(idx_full).maxDisp] <= mobility_thresh);

if isempty(idx_lowmob)
    error('No tracks meet both criteria. Adjust threshold.');
end

tracks_full = tracks(idx_lowmob);

% ------- 2. Center reference tracks and compute average (drift) motion -------
num_good = numel(tracks_full);
centered_tracks = cell(num_good,1);
for i = 1:num_good
    trk = tracks_full{i}(:,2:3);
    trk_centered = trk - mean(trk,1);
    centered_tracks{i} = trk_centered;
end

X = zeros(max_length, 2, num_good);
for i = 1:num_good
    X(:,:,i) = centered_tracks{i};
end
avg_motion = squeeze(mean(X,3)); % [frames x 2]

% ------- 3. Plot average (drift) trajectory, color-coded by time -------
if opts.plot_drift
    figure('Name','Average Reference (Drift) Trajectory','Color','w');
    hold on
    cmap = jet(max_length);
    for i = 1:max_length-1
        plot(avg_motion(i:i+1,1), avg_motion(i:i+1,2), '-', 'Color', cmap(i,:), 'LineWidth', 2);
    end
    plot(0,0,'ro','MarkerFaceColor','r');
    xlabel('X'); ylabel('Y');
    title('Average (Drift) Trajectory (Color = Time)');
    axis equal; grid on;
    colormap(cmap);
    cb = colorbar; cb.Ticks = linspace(0,1,5); cb.TickLabels = round(linspace(1,max_length,5));
    ylabel(cb,'Frame');
    hold off
end

% ------- 4. Drift-correct all tracks -------
corrected_tracks = cell(size(tracks));
for i = 1:numel(tracks)
    trk = tracks{i};
    frame_idx = trk(:,1)+1;          % frame numbers (assume starts at 0, or as given)
    correction = avg_motion(frame_idx, :); % lookup correction for each frame
    corrected_tracks{i} = [frame_idx, trk(:,2:3) - correction];
end

% ------- 5. Plot average drift-corrected trajectory of all full-length tracks -------
idx_full_corr = find(track_lengths == max_length);
tracks_full_corr = corrected_tracks(idx_full_corr);
num_full_corr = numel(tracks_full_corr);

if opts.plot_corrected
    X_corr = zeros(max_length,2,num_full_corr);
    for i = 1:num_full_corr
    temp_track = tracks_full_corr{i}(:,2:3);
        X_corr(:,:,i) = temp_track-mean(temp_track,1);
    end
    mean_corr = squeeze(mean(X_corr,3));

    figure('Name','Average Drift-Corrected Trajectory','Color','w');
    hold on
    for i = 1:max_length-1
        plot(mean_corr(i:i+1,1), mean_corr(i:i+1,2), '-', 'Color', cmap(i,:), 'LineWidth', 2);
    end
    plot(0,0,'ro','MarkerFaceColor','r');
    xlabel('X'); ylabel('Y');
    title('Average Drift-Corrected Trajectory (Color = Time)');
    axis equal; grid on;
    colormap(cmap); cb = colorbar;
    cb.Ticks = linspace(0,1,5); cb.TickLabels = round(linspace(1,max_length,5));
    ylabel(cb,'Frame');
    hold off
end

% ------- 6. Plot first 16 corrected full-length tracks in subplots -------
n_show = min(16, num_full_corr);
if opts.plot_16corr
    figure('Name','First 16 Drift-Corrected Tracks','Color','w');
    tiledlayout(4,4,'Padding','compact','TileSpacing','compact');
    for k = 1:n_show
        nexttile;
        trk = tracks_full_corr{k}(:,2:3);
        N = size(trk,1);
        cmap = jet(N);
        for i = 1:N-1
            plot(trk(i:i+1,1), trk(i:i+1,2), '-', 'Color', cmap(i,:), 'LineWidth', 1.2);
            hold on;
        end
        plot(trk(1,1),trk(1,2),'ro','MarkerFaceColor','r');
        axis equal; grid on;
        set(gca, 'XTick', [], 'YTick', []);
        title(['Track ' num2str(k)]);
    end
    % Shared colorbar
    colormap(jet(max_length));
    cb = colorbar('Position',[0.92 0.11 0.02 0.815]);
    cb.Ticks = linspace(0,1,5);
    cb.TickLabels = round(linspace(1,max_length,5));
    ylabel(cb, 'Frame');
end

% ------- 7. Plot not corrected first 16 full-length tracks -------
n_show = min(16, num_full_corr);
if opts.plot_16corr
    figure('Name','First 16 Raw Tracks','Color','w');
    tiledlayout(4,4,'Padding','compact','TileSpacing','compact');
    for k = 1:n_show
        nexttile;
        trk = tracks{k}(:,2:3);
        N = size(trk,1);
        cmap = jet(N);
        for i = 1:N-1
            plot(trk(i:i+1,1), trk(i:i+1,2), '-', 'Color', cmap(i,:), 'LineWidth', 1.2);
            hold on;
        end
        plot(trk(1,1),trk(1,2),'ro','MarkerFaceColor','r');
        axis equal; grid on;
        set(gca, 'XTick', [], 'YTick', []);
        title(['Track ' num2str(k)]);
    end
    % Shared colorbar
    colormap(jet(max_length));
    cb = colorbar('Position',[0.92 0.11 0.02 0.815]);
    cb.Ticks = linspace(0,1,5);
    cb.TickLabels = round(linspace(1,max_length,5));
    ylabel(cb, 'Frame');
end

end

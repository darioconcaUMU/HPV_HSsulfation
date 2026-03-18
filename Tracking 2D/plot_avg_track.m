function plot_avg_track(track_cell)

% tracks: cell array, each cell is an [N x 2] (x,y) matrix

tracks = track_cell;

% 1. Find tracks with max length
track_lengths = cellfun(@(trk) size(trk,1), tracks);
max_length = max(track_lengths);
idx_full = find(track_lengths == max_length);

% 2. Keep only those tracks
tracks_full = tracks(idx_full);

% 3. Center each track by its mean position
centered_tracks = cell(size(tracks_full));
for i = 1:numel(tracks_full)
    trk = tracks_full{i};
    trk_centered = trk - mean(trk,1); % subtract mean (x and y)
    centered_tracks{i} = trk_centered;
end

% 4. Stack to get average position at each frame
num_full = numel(centered_tracks);
X = zeros(max_length, 2, num_full);
for i = 1:num_full
    X(:,:,i) = centered_tracks{i};
end
mean_positions = squeeze(mean(X,3)); % [frames x 2]

figure; hold on

% Plot all individual centered tracks as faint lines
for i = 1:num_full
    plot(centered_tracks{i}(:,1), centered_tracks{i}(:,2), 'Color', [0.7 0.7 0.7 0.3]);
end

% Plot mean trajectory with color by time
cmap = jet(max_length); % or parula, hot, etc.
for i = 1:max_length-1
    plot(mean_positions(i:i+1,1), mean_positions(i:i+1,2), '-', ...
         'Color', cmap(i,:), 'LineWidth', 2);
end

plot(0,0,'ro','MarkerFaceColor','r'); % origin
xlabel('Mean X');
ylabel('Mean Y');
title('Average Centered Position vs Frame (Color = Time)');
axis equal; grid on;

% Add colorbar to indicate time/frame
colormap(cmap);
cb = colorbar;
cb.Ticks = linspace(0,1,5);
cb.TickLabels = round(linspace(1,max_length,5));
ylabel(cb, 'Frame');

hold off;

plot_first16_colored_tracks(track_cell)

end

function plot_first16_colored_tracks(track_cell)
% Plots first 16 tracks, color-coded by time, after mean-centering

tracks = track_cell;
track_lengths = cellfun(@(trk) size(trk,1), tracks);
max_length = max(track_lengths);
idx_full = find(track_lengths == max_length);

tracks_full = tracks(idx_full);

n_show = min(16, numel(tracks_full));

figure('Name','First 16 Centered Tracks (color = time)', 'Color','w');
tiledlayout(4,4, 'Padding','compact', 'TileSpacing','compact');

for k = 1:n_show
    nexttile;
    trk = tracks_full{k};
    trk_centered = trk - mean(trk,1);
    N = size(trk_centered,1);
    cmap = jet(N);
    for i = 1:N-1
        plot(trk_centered(i:i+1,1), trk_centered(i:i+1,2), '-', ...
            'Color', cmap(i,:), 'LineWidth', 1);
        hold on;
    end
    plot(0,0,'ro','MarkerFaceColor','r'); % origin
    axis equal; grid on;
    set(gca, 'XTick', [], 'YTick', []);
    title(['Track ' num2str(k)]);
end

% Add a single colorbar for time at the right of all subplots
cb = colorbar('Position',[0.92 0.11 0.02 0.815]);
colormap(jet(max_length));
cb.Ticks = linspace(0,1,5);
cb.TickLabels = round(linspace(1,max_length,5));
ylabel(cb, 'Frame');

end
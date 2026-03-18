% Define the number of binders to simulate
binders = round([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 17, 21, ...
                 26, 33, 41, 51, 64, 80, 100]);

% Loop over all binder numbers
for i = 1:length(binders)

    % Repeat each condition 50 times
    for j = 1:50

        % Run Brownian tether simulation for the current number of binders
        % Ignore the first output and store the position array
        [~, posArray{i,j}] = brownianTether3DMulti(binders(i));
    end

    % Display current binder number as progress indicator
    binders(i)
end

% Loop again to analyze the simulated trajectories
for i = 1:length(binders)
    for j = 1:50

        % Compute the mean particle position over the trajectory
        mean_pos = mean(posArray{i,j}, 1);

        % Replicate the mean position so it matches the trajectory length
        mean_pos = repmat(mean_pos, length(posArray{i,j}), 1);

        % Center the trajectory by subtracting its mean position
        pos_avg{i,j} = posArray{i,j} - mean_pos;

        % Compute radial distance in the XY plane from the centered trajectory
        distanceXY(i,j,:) = sqrt((pos_avg{i,j}(:,1)).^2 + (pos_avg{i,j}(:,2)).^2);
    end

    % Display current binder number as progress indicator
    binders(i)

    % Collect all XY distances for this binder number across repeats
    temp = distanceXY(i,:,:);
    dist_all(i,:) = temp(:);
end
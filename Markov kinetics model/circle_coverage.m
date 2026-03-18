clc; clear; %close all;
simulate = 1;

%% Parameters (arrays: one case per element)
%lambdatot = [38.1, 1.6]*1e3*0.035/100;  % points per µm^2
%lambdatot = 38.1*1e3*[80.86, 8.89, 0.56, 0.035]/100;  % points per µm^2
lambdatot = [1.1, 1.3, 1.5, 1.6]*1e3;  % points per µm^2
%lambdatot = [1.1, 1.3, 1.5, 1.6]*1e3*4.0672;  % points per µm^2
R         = sqrt([1726, 1726, 1726, 1655]/pi)/1000; % µm;
% R         = sqrt([1.818655106633252e+03, 1.813055846046759e+03, ...
%                   1.768420097624555e+03, 5.793722030157644e+02]/pi)/1000; % µm
N         = 2;         % Minimum number of points required inside circle
L         = 10;         % Domain size (µm)
N_max     = 20;        % Maximum N considered for distance maps
res       = 300;       % Grid resolution for simulation plots

d_th_tot  = R;                     % Distance threshold == circle radius
A_tot     = pi .* R.^2;            % Circle areas
Kvals     = 0:N_max;               % For tabulating PDFs

%% Poisson probabilities per case
P_eq_k   = zeros(numel(lambdatot), numel(Kvals)); % P(K = k)
P_ge_N   = zeros(numel(lambdatot), 1);            % P(K >= N)
lambdaAs = lambdatot(:) .* A_tot(:);              % mean count in circle

for ii = 1:numel(lambdatot)
    mu = lambdaAs(ii);
    % PDF table for k=0..N_max (for display if desired)
    P_eq_k(ii, :) = poisspdf(Kvals, mu);
    % Fraction of area where circle covers >= N points:
    % (in a homogeneous Poisson field, this equals probability at a random location)
    P_ge_N(ii) = 1 - poisscdf(N-1, mu);
end

%% Display summary
fprintf('=== Poisson results per case (N >= %d) ===\n', N);
for ii = 1:numel(lambdatot)
    fprintf('Case %d: lambda = %.1f /um^2, R = %.4f um  ->  P>=N = %.4f (%.2f%%)\n', ...
        ii, lambdatot(ii), R(ii), P_ge_N(ii), 100*P_ge_N(ii));
end

%% (Optional) simulation for one selected case
if simulate
    for idx_sim = 1:length(lambdatot)      
        lambda  = lambdatot(idx_sim);
        d_th    = d_th_tot(idx_sim);

        % Generate Poisson-distributed points in LxL box
        num_points = poissrnd(lambda * L^2);
        x = L * rand(num_points, 1);
        y = L * rand(num_points, 1);

        % Grid
        [Xg, Yg] = meshgrid(linspace(0, L, res), linspace(0, L, res));

        % Distances to nearest N_max neighbors for each grid point
        D_all = zeros(res, res, N_max);
        for lin = 1:numel(Xg)
            [r, c] = ind2sub(size(Xg), lin);
            dx = x - Xg(lin);
            dy = y - Yg(lin);
            d  = sqrt(dx.^2 + dy.^2);
            d  = sort(d, 'ascend');
            upto = min(N_max, numel(d));
            D_all(r, c, 1:upto) = d(1:upto);
            if upto < N_max
                D_all(r, c, upto+1:end) = inf; % if fewer points exist
            end
        end

        % Thresholded area for N_selected
        N_selected = N;
        thresholded_area = D_all(:,:,N_selected) <= d_th; % true where at least N points within d_th

        % Cumulative map: how many neighbors within radius d_th (0..N_max)
        cumulative_map = zeros(size(Xg));
        for n = 1:N_max
            cumulative_map = cumulative_map + (D_all(:,:,n) <= d_th);
        end
        max_level = max([max(cumulative_map(:)), 1]);

        % Area percentages for each level (0..max_level)
        total_area = numel(Xg);
        area_percentages = zeros(1, max_level+1);
        for n = 0:max_level
            area_percentages(n+1) = sum(cumulative_map(:) == n)/total_area*100;
        end

        fprintf('\nArea percentages (simulation, case %d, radius = %.0f nm):\n', idx_sim, d_th*1000);
        for n = 0:max_level
            fprintf('  #points within radius = %2d : %6.2f %%\n', n, area_percentages(n+1));
        end

        % Plot 1: single N-th distance thresholded region
        figure; hold on; axis equal tight;
        colormap(parula(2));
        contourf(Xg, Yg, thresholded_area, 1, 'LineColor', 'none');
        scatter(x, y, 4, 'r', 'filled');
        title(sprintf('\\ge %d points within radius %.0f nm', N_selected, d_th*1000));
        xlabel('X (µm)'); ylabel('Y (µm)');

        % Plot 2: cumulative map with discrete colors
        figure; hold on; axis equal tight;
        colormap(parula(max_level+1));
        contourf(Xg, Yg, cumulative_map, 'LineColor', 'none');
        clim([0 max_level]);
        cb = colorbar('Ticks', 0:max_level, 'TickLabels', string(0:max_level));
        cb.Label.String = '# points within radius';
        scatter(x, y, 4, 'r', 'filled');
        title(sprintf('# points within radius %.0f nm (max %d)', d_th*1000, N_max));
        xlabel('X (µm)'); ylabel('Y (µm)');
    end
end




% %%
% %% Generate Poisson-Distributed Points
% num_points = poissrnd(lambda * L^2); % Number of points follows Poisson distribution
% x = L * rand(num_points, 1); % Random x-coordinates
% y = L * rand(num_points, 1); % Random y-coordinates
%
% %% Plot the Poisson Distribution
% figure; hold on; axis equal;
% scatter(x, y, 20, 'b', 'filled'); % Scatter plot of points
% xlabel('X'); ylabel('Y');
% title(['Poisson-Distributed Points (\lambda = ' num2str(lambda) ' per unit area)']);
% grid on;
% %
% d_th = R;   % Distance threshold
%
% %% Generate Poisson-Distributed Points
% num_points = poissrnd(lambda * L^2); % Poisson-distributed number of points
% x = L * rand(num_points, 1); % Random x-coordinates
% y = L * rand(num_points, 1); % Random y-coordinates
%
% %% Define a grid for computing distances
% res = 100;  % Resolution of the grid
% [Xg, Yg] = meshgrid(linspace(0, L, res), linspace(0, L, res));
% D = zeros(size(Xg));  % Store distance to the N-th closest point
%
% %% Compute distance to the N-th closest point for each grid location
% for i = 1:numel(Xg)
%     % Compute distances from grid point to all scattered points
%     distances = sqrt((x - Xg(i)).^2 + (y - Yg(i)).^2);
%     % Sort distances and get the N-th smallest one
%     sorted_distances = sort(distances);
%     D(i) = sorted_distances(N);
% end
%
% %% Identify the area where D > d_th
% for i = 1: length(d_th)
% thresholded_area = D > d_th(i);
%
% %% Plot results
% figure; hold on; axis equal;
% contourf(Xg, Yg, thresholded_area, 1, 'LineColor', 'none', 'FaceAlpha', 0.5); % Highlight area
% colorbar;
% title(['Regions where distance to ' num2str(N) '-th closest point > ' num2str(d_th)]);
% scatter(x, y, 2, 'b', 'filled'); % Scatter plot of points
% xlabel('X'); ylabel('Y');
% end
%% Analyze Simulation Results
% Assume that you have already computed your simulation results, stored in a
% cell array "track" with dimensions: [length(k_on) x length(binder) x 10]
% Each cell track{i,j,k} is a structure with the following fields:
%   - positions: an N x 2 matrix (nm) of particle positions over time
%   - bonds: a vector (length N) with the number of bonds at each time step
%
% Also, the simulation parameters are in the structure "time" with field time.dt.
%
% This script computes:
% 1. Total path length (sum of consecutive distances)
% 2. Maximum distance (largest pairwise distance on the convex hull)
% 3. Diffusion coefficient (from MSD vs. time)
% 4. Average number of bonds

% Preallocate a results cell array with the same size as "track"
results = cell(size(track));

% Loop over each simulation
for i = 1:size(track,1)
    for j = 1:size(track,2)
        for k = 1:size(track,3)
            % Extract trajectory and bond history for this simulation
            traj = track{i,j,k}.positions;  % [N x 2] matrix of positions
            bonds = track{i,j,k}.bonds;       % vector of bond numbers

            % Number of recorded positions
            N = size(traj,1);

            % 1. Total Path Length: sum of Euclidean distances between consecutive positions
            dists = sqrt(sum(diff(traj,1,1).^2,2));
            total_length = sum(dists);

            % 2. Maximum Distance: longest distance among points on the convex hull
            try
                % Compute convex hull indices using convhull
                [hull_idx, hull_area] = convhull(traj(:,1), traj(:,2));
                hull_points = traj(hull_idx, :);
                % Compute all pairwise distances between hull points
                pd = pdist(hull_points);
                max_distance = max(pd);
            catch
                if max([length(unique(traj(:,1))),length(unique(traj(:,2)))])>1
                    clear temp_traj
                    if length(unique(traj(:,1)))>length(unique(traj(:,2)))
                        [temp_traj(:,1), temp] = unique(traj(:,1));
                        temp_traj(:,2) = traj(temp,2);
                    else
                        [temp_traj(:,2), temp] = unique(traj(:,2));
                        temp_traj(:,1) = traj(temp,1);
                    end
                    pd = pdist([temp_traj(:,1) temp_traj(:,2)]);
                    max_distance = max(pd);
                    hull_area = 0;
                else
                    max_distance = 0;
                    hull_area = 0;
                end
            end

            % % 3. Diffusion Coefficient using MSD
            % % Compute the time-averaged MSD for different lag times
            % max_tau = floor(N/3);
            % if max_tau<20
            %     D = NaN;
            %     gof.adjrsquare = NaN;
            % else
            %     msd = zeros(max_tau,1);
            %     for tau = 1:max_tau
            %         % Displacements for the given lag time tau
            %         delta = traj(1+tau:end,:) - traj(1:end-tau,:);
            %         msd(tau) = mean(sum(delta.^2,2));
            %     end
            %     % Fit a line to the first few points (e.g., first 10 or fewer if not available)
            %     t_fit = ( (1:max_tau)' * time.dt );
            % 
            %     % Assume your data is loaded as vectors t (time) and MSD (mean squared displacement)
            %     % Define your fitting function
            %     ft = fittype('A*(1-exp(-b*t))+C','independent','t','coefficients',{'A','b','C'});
            % 
            %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            %     opts.Display = 'Off';
            %     opts.Lower = [0 0 0];
            %     opts.Robust = 'Bisquare';
            %     opts.StartPoint = [max(msd) mean(t_fit)/5 0];
            % 
            %     % Fit the MSD data with curve fitting toolbox
            %     try
            %         [fit_result, gof] = fit( t_fit, msd, ft, opts );
            %     catch
            %         try
            %             opts.StartPoint(2) = opts.StartPoint(2)/10;
            %             [fit_result, gof] = fit( t_fit, msd, ft, opts );
            %         catch
            %             gof.adjrsquare = NaN;
            %             fit_result.A = NaN;
            %             fit_result.b = NaN;
            %         end
            %     end
            % 
            %     % Extract parameters
            %     A = fit_result.A;
            %     tau = 1/fit_result.b;
            % 
            %     % Estimate diffusion coefficient (D) from parameters (for confined diffusion, D can be derived as below)
            %     D = A/(4*tau);
            % 
            %     % Define threshold for good fit (adjust as needed, here 0.9)
            %     rsquare_thresh = 0.3;
            % 
            %     % Check goodness of fit, reject poor fits
            %     if gof.adjrsquare < rsquare_thresh
            %         disp('Fit rejected due to low adjusted R-square. Setting D to NaN');
            %         D = NaN;
            %     end
            % 
            %     % fit_points = min(10, max_tau);
            %     % t_fit = ( (1:fit_points)' * time.dt );
            %     % p = polyfit(t_fit, msd(1:fit_points), 1);
            %     % slope = p(1);
            %     % % In 2D, MSD = 4*D*t so the diffusion coefficient D is:
            %     % D = slope / 4
            % 
            % end

            % 4. Average Number of Bonds over the simulation
            avg_bonds = mean(bonds);

            % Save metrics in a structure
            res.Nsteps = N;
            res.length = total_length;
            res.max_distance = max_distance;
            res.hull_area = hull_area;
            res.D = 0;%D;
            res.D_qual = 0;%gof.adjrsquare;
            res.avg_bonds = avg_bonds;

            % Store in the results cell array
            results{i,j,k} = res;
        end
    end
end

%% Assume you already have the cell array "results" from your simulation
% results{i,j,k} is a structure with fields:
%   .length, .max_distance, .D, .avg_bonds
%
% Let:
%   - num_kon = length(k_on) (i dimension)
%   - num_binder = length(binder) (j dimension)
%   - num_runs = 10 (k dimension)

num_kon = size(results,1);
num_binder = size(results,2);
num_runs = size(results,3);

% Preallocate matrices for averaged metrics
avg_length      = zeros(num_kon, num_binder);
avg_time = zeros(num_kon, num_binder);
avg_max_dist    = zeros(num_kon, num_binder);
avg_D           = zeros(num_kon, num_binder);
avg_bonds       = zeros(num_kon, num_binder);
avg_speed       = zeros(num_kon, num_binder);
avg_radius      = zeros(num_kon, num_binder);
avg_D_ext       = zeros(num_kon, num_binder);

% Loop over k_on (i) and binder (j) indices, and average over runs (k)
for i = 1:num_kon
    for j = 1:num_binder
        temp_length = zeros(num_runs,1);
        temp_time   = zeros(num_runs,1);
        temp_max    = zeros(num_runs,1);
        temp_D      = zeros(num_runs,1);
        temp_bonds  = zeros(num_runs,1);
        temp_speed  = zeros(num_runs,1);
        temp_radius = zeros(num_runs,1);
        temp_D_ext  = zeros(num_runs,1);
        for k = 1:num_runs
            temp_length(k) = results{i,j,k}.length;
            temp_time(k) = results{i,j,k}.Nsteps;
            temp_max(k)    = results{i,j,k}.max_distance;
            temp_D(k)      = results{i,j,k}.D;
            temp_bonds(k)  = results{i,j,k}.avg_bonds;
            temp_speed(k) = temp_length(k)/temp_time(k);
            temp_radius(k)  = sqrt(results{i,j,k}.hull_area)/2;
            temp_D_ext(k) = results{i,j,k}.hull_area/temp_time(k);
        end
        avg_length(i,j) = mean(temp_length);
        avg_time(i,j) = mean(temp_time);
        avg_max_dist(i,j) = mean(temp_max);
        avg_D(i,j) = mean(temp_D);
        avg_bonds(i,j) = mean(temp_bonds);
        avg_speed(i,j) = mean(temp_speed);
        avg_radius(i,j) = mean(temp_radius);
        avg_D_ext(i,j) = mean(temp_D_ext);
    end
end

%% Create imagesc plots for each averaged metric

[X, Y] = meshgrid(binder, kinetics.k_on);

figure;
scaleX = 'log'; %'lin' or 'log'
scaleY = 'log'; %'lin' or 'log'

subplot(3,2,1);
pcolor(X, Y, avg_length);
shading flat;
colorbar;
title('Average Path Length');
xlabel('Binder concentration');
ylabel('k_{on}');
set(gca, 'YDir', 'normal');  % Ensure the y-axis is in the normal direction
set(gca, 'XScale', scaleX);
set(gca, 'YScale', scaleY);

subplot(3,2,2);
pcolor(X, Y, avg_max_dist);
shading flat;
colorbar;
title('Average Maximum Distance');
xlabel('Binder concentration');
ylabel('k_{on}');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', scaleX);
set(gca, 'YScale', scaleY);

% subplot(2,3,3);
% pcolor(X, Y, avg_D);
% shading flat;
% colorbar;
% title('Average Diffusion Coefficient');
% xlabel('Binder concentration');
% ylabel('k_{on}');
% set(gca, 'YDir', 'normal');
% set(gca, 'XScale', scaleX);
% set(gca, 'YScale', scaleY);

subplot(3,2,3);
pcolor(X, Y, avg_D_ext);
shading flat;
colorbar;
title('Average Diffusion Coefficient');
xlabel('Binder concentration');
ylabel('k_{on}');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', scaleX);
set(gca, 'YScale', scaleY);

subplot(3,2,4);
pcolor(X, Y, avg_bonds);
shading flat;
colorbar;
title('Average Number of Bonds');
xlabel('Binder concentration');
ylabel('k_{on}');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', scaleX);
set(gca, 'YScale', scaleY);

subplot(3,2,5);
pcolor(X, Y, avg_time);
shading flat;
colorbar;
title('Average Time track');
xlabel('Binder concentration');
ylabel('k_{on}');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', scaleX);
set(gca, 'YScale', scaleY);

subplot(3,2,6);
pcolor(X, Y, avg_speed);
shading flat;
colorbar;
title('Average Speed track');
xlabel('Binder concentration');
ylabel('k_{on}');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', scaleX);
set(gca, 'YScale', scaleY);


% Optionally, if you want a combined figure with tight layout:
% set(gcf, 'Position', [100 100 1200 800]);

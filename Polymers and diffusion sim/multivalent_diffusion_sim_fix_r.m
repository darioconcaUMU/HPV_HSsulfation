%% Particle Binding Simulation in Matlab
% This script simulates a particle that moves by forming and breaking bonds 
% to randomly distributed binders on a surface.
% The particle's position is always the center-of-mass of its attached bonds.
% The simulation stops when the particle has no bonds.

function [track, binders_loc] = multivalent_diffusion_sim_fix_r(binders, kinetics, time)

% all dimensions in nm
% %% Parameters
% L = 1000;                    % Domain size: [0,L] x [0,L]
% binder_concentration = 8e-2;%1.6e-3;  % Binder density (number per unit area)
% N_binders = round(binder_concentration * L^2);  % Total number of binders
% r_int = 20;               % Interaction radius (only binders within this distance can bind)
% k_on = 3e3;                  % Attachment rate (1/s)
% k_off = 3e2;               % Detachment rate (1/s)
% dt = 0.01;                 % Time step (s)
% max_time = 100;            % Maximum simulation time (s)

% Dimension fraction of virion interaction area and 1/100 of koff
L = 20;                    % Domain size: [0,L] x [0,L]
binder_concentration = binders.concentration;  % Binder density (number per unit area)
N_binders = round(binder_concentration * L^2);  % Total number of binders
r_int = sqrt(1/pi);               % Interaction radius (only binders within this distance can bind)
k_on = kinetics.k_on;                  % Attachment rate (1/s)
k_off = kinetics.k_off;               % Detachment rate (1/s)
dt = time.dt;                 % Time step (s)
max_time = time.max;            % Maximum simulation time (s)

% r_int = (r_int*25 + 11.519*(1:200).^-0.3203)/25;

%% Generate random binder positions
binders_loc = [L*rand(N_binders,1), L*rand(N_binders,1)];

%% Initialize particle: start at a random binder
start_idx = randi(N_binders);
particle_bonds = start_idx;        % Stores indices of binders attached to the particle
particle_pos = binders_loc(start_idx, :);% Initial particle position set to the binder's position

%% Data storage for analysis
trajectory = particle_pos;         % Store particle positions over time
bondHistory = length(particle_bonds);% Number of bonds at each time point
times = 0;                         % Time vector
time = 0;

%% Simulation Loop
while ~isempty(particle_bonds) && time < max_time
    time = time + dt;
    new_bonds = [];

    % ----- Attachment Step -----
    % Check free binders that are within the interaction radius.
    free_binders = setdiff(1:N_binders, particle_bonds);
    if ~isempty(free_binders)
        free_positions = binders_loc(free_binders, :);
        dists = sqrt(sum((free_positions - particle_pos).^2, 2));
        potential = free_binders(dists < r_int);
        % % For each potential binder, try to attach with probability (k_on * dt)
        % new_bonds = potential(rand(1,length(potential)) < (1 - exp(-k_on * dt)));

        temp = rand(length(potential), 1);
        new_bonds = potential(temp<(1 - exp(-k_on * dt)));

    end

    % ----- Detachment Step -----
    % Each attached bond has a chance to detach.

    detached = [];

    temp = rand(length(particle_bonds), 1);
        detached = particle_bonds(temp<1 - exp(-k_off * dt));

    % Remove bonds that have detached (they become free again)
    particle_bonds = setdiff(particle_bonds, detached);
    particle_bonds = [particle_bonds; new_bonds'];
    
    % ----- Update Particle Position -----
    if ~isempty(particle_bonds)
        pos = binders_loc(particle_bonds, :);
        % New position is the geometrical (center-of-mass) of attached bonds
        particle_pos = mean(pos, 1);
        trajectory = [trajectory; particle_pos];
        times = [times; time];
        bondHistory = [bondHistory; length(particle_bonds)];
    else
        % End simulation if no bonds remain
        break;
    end
end

track.positions = trajectory;
track.bonds = bondHistory;

%% Plotting Results

% % Plot binder positions and the particle trajectory
% figure;
% plot(binders_loc(:,1), binders_loc(:,2), 'k.', 'MarkerSize', 10); hold on;
% plot(trajectory(:,1), trajectory(:,2), 'r-', 'LineWidth', 2);
% xlabel('X Position');
% ylabel('Y Position');
% title('Particle Trajectory via Binding/Unbinding');
% legend('Binders','Trajectory');
% 
% % Plot the number of bonds as a function of time
% figure;
% plot(times, bondHistory, 'b-', 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Number of Bonds');
% title('Number of Bonds Over Time');

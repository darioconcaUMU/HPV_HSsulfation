clear

% Unbinding rate constant
kinetics.k_off = 1;

% Time settings for the simulation
time.dt = 0.1;       % time step
time.max = 500;      % total simulation time

% Previous linear parameter sweeps (kept here for reference)
% kinetics.k_on = 0.2:0.2:5;
% binder = [0.2:0.2:2, 3:10, 12:2:30];

% Define logarithmically spaced k_on values from 0.1 to 10000
temp = log([0.1 10000]);
kinetics.k_on = exp(linspace(temp(1), temp(2), 35));

% Define logarithmically spaced binder concentrations from 0.5 to 100
temp = log([0.5 100]);
binder = exp(linspace(temp(1), temp(2), 25));

% Loop over all k_on values
for i = 1:length(kinetics.k_on)

    % Loop over all binder concentrations
    for j = 1:length(binder)

        % Repeat each parameter combination 25 times
        for k = 1:25

            % Create a temporary kinetics structure for this run
            temp_kinetics = kinetics;
            temp_kinetics.k_on = kinetics.k_on(i);

            % Set binder concentration for this run
            binders.concentration = binder(j);

            % Run simulation and store outputs
            [track{i,j,k}, binders_loc{i,j,k}] = multivalent_diffusion_sim(binders, temp_kinetics, time);
        end
    end

    % Display progress in command window
    i
end

% Save all simulation results to a MAT file
save('Diffution fix radius 5000 steps.mat')
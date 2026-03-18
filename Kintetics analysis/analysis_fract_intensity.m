function [fraction, parameters] = analysis_fract_intensity(parameters, intensity_fraction, track_props, spot_record, metadata)

% Store intensity fraction threshold in parameters
parameters.intensityFraction = intensity_fraction;

% Filter intensity fractions greater than the previous threshold
parameters.intensityFraction = parameters.intensityFraction(parameters.intensityFraction >= parameters.prev_high);

% Define intensity sections based on provided thresholds
sections = [parameters.prev_low * parameters.prev_high, parameters.intensityFraction, ceil(max(spot_record(:,5)) + 1)];

% Loop through each intensity section to perform detailed analysis
for k = 2:length(sections)

    % Extract initial prominence values from track properties
    start_int = [track_props(:).initial_prominence];
    
    % Create a label describing the current intensity section
    label = [' - Intensities: ', num2str(sections(k-1)), '-', num2str(sections(k))];

    % Select tracks whose initial prominence falls within the current section
    selection = start_int > sections(k-1) & start_int < sections(k);

    % Calculate particle frame fraction, surface coverage, and its slope
    [fraction{k}.particle_frame_fraction{k}, fraction{k}.frame_props, ...
        fraction{k}.surf_coverage, fraction{k}.surf_cov_slope] = ...
        surface_coverage_stuck(track_props(selection), metadata, label);

    % Perform detailed analysis on track duration, cumulative duration, and arrival times
    [fraction{k}.duration, fraction{k}.cum_duration, fraction{k}.arrival, fraction{k}.cum_arrival, fraction{k}.gof_analysis, ...
        fraction{k}.fit_dissociation, fraction{k}.fit_association] = track_analysis(track_props(selection), metadata, label);

    % Optional: Perform double exponential fitting on dissociation data if parameter is enabled
    if isfield(parameters, 'double_exp')
        if parameters.double_exp == 1
            [fraction{k}.fit_dissociation_2exp, fraction{k}.gof_2edx] = doubleExpFit(...
                fraction{k}.cum_duration, fraction{k}.duration, metadata, fraction{k}.fit_dissociation);
        end
    end

end

fraction(1) = [];  % Removes first element because it is always empty

% % recreated graphs for all particles
% 
% % Create a label describing the current intensity section
% label = 'All particles';
% 
% % Select tracks whose initial prominence falls within the current section
% selection = start_int > 0 & start_int < sections(end);
% 
% % Calculate particle frame fraction, surface coverage, and its slope
% [fraction{k}.particle_frame_fraction{k}, fraction{k}.frame_props, ...
%     fraction{k}.surf_coverage, fraction{k}.surf_cov_slope] = ...
%     surface_coverage_stuck(track_props(selection), metadata, label);
% 
% % Perform detailed analysis on track duration, cumulative duration, and arrival times
% [fraction{k}.duration, fraction{k}.cum_duration, fraction{k}.arrival, fraction{k}.cum_arrival, fraction{k}.gof_analysis, ...
%     fraction{k}.fit_dissociation, fraction{k}.fit_association] = track_analysis(track_props(selection), metadata, label);
% 
% % Optional: Perform double exponential fitting on dissociation data if parameter is enabled
% if isfield(parameters, 'double_exp')
%     if parameters.double_exp == 1
%         [fraction{k}.fit_dissociation_2exp, fraction{k}.gof_2edx] = doubleExpFit(...
%             fraction{k}.cum_duration, fraction{k}.duration, metadata, fraction{k}.fit_dissociation);
%     end
% end
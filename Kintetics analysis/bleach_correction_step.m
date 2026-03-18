function [tracks, track_props, bleach_tracks, bleach_props, step_amplitude] = ...
    bleach_correction_step(tracks, track_props, params)

% This function determines which spots have disapeared because the particle
% detached and which have bleached instead. Compared to the previous
% versions of "bleach_correction", the sorting mechanism has been
% completely changed. Sorting is done by evaluating the size of the
% intensity drop at the end of the track to the noise level. If the spot is
% significantly higher than the noise, the particle is considered detached,
% otherwise it is considered bleached.

if isempty(track_props)
    warning('No tracks detected, so no bleaching analysis was perfored');
    track_props(1).bleached = [];
    bleach_tracks = tracks;
    bleach_props = track_props;
    step_amplitude = [];
    return
end

% Makes sure that at least one detachment event is observed. Otherwise the
% analysis is not performed

valid = false;

% decides if to perform a complete step fit of the intenisty drop or simply
% evaluate the intensity before and after the particle disapearence. The
% first option is more accurate but significantly slower. In most
% application the second ("false") should be used.
complete_fit = false;

if complete_fit
    
    % In the case complete_fit is true, a fit of the trace intensity around
    % the drop is performed.
    
    % Cycles through all the detected tracks
    for i=1:length(tracks)
        
        % Extracts the track of interest
        current = tracks{i};
        
        % Calculated the mean intensity after the particle disappearance
        mean_int(i) = mean(current(end-params.buffer+1:end,5));
        
        % The sixth column of tracks is 0 when the particle is present and
        % then counts up to 10 frames after detachment. The track is
        % considered in the analysis if the count up reaches the value in
        % params.buffer.
        
        if max(tracks{i}(:,6))==params.buffer
            
            % Selects a segment of the trace around the particle
            % disappearance which is 4 times as long as the intensity
            % recorded after detachment.
            % If the track is shorter than that, the whole track is
            % selected.
            temp = min([length(current) params.buffer*4])-1;
            
            segment = current(end-temp+1:end,5);
            
            % Defines starting parameters for the following fit. (a0:
            % offset, b0: ampliture step, c0: x_position).
            a0 = mean(segment(1:end-params.buffer));
            b0 = a0-mean(segment(end-params.buffer+1:end));
            mu0 = temp - params.buffer;
            
            % Define fitting function (continous step) and options.
            ft = fittype('a - b*normcdf(x,mu,sig)','indep','x');
            Fopts = fitoptions(ft);
            Fopts.Lower = [0 0 mu0-3 0];
            Fopts.Upper = [max(segment) max(segment) mu0+3 3];
            Fopts.StartPoint = [a0 b0 mu0 1];
            
            % Performes the fit
            [mdl gof] = fit((1:length(segment))',segment,ft,Fopts);
            
            % Saves the sdjusted R squared for subsequent track sorting
            adjrsq(i) = gof.adjrsquare;
            %         plot(mdl, (1:length(segment))',segment)
            
            % saves mean intensity and standard deviation of the pixel
            % value after detachment. 
            
            bkg_int(i) = mdl.a-mdl.b; % Does this line really make sense??
            bkg_std(i) = std(segment(end-params.buffer+1:end));
            step_amplitude(i) = mdl.b;
            valid = true;
        else
            % If no proper detachment is detected NaN is used.
            bkg_int(i) = NaN;
            bkg_std(i) = NaN;
            step_amplitude(i) = NaN;
            adjrsq(i) = NaN;
        end
    end
    
else
    
    % To only find the difference before and after the detachment event
    
    for i=1:length(tracks)
        
        % Estracts the part of the track before detachent
        track_left = tracks{i}(tracks{i}(:,6)==0, 5);
        % Finds the minimum between the length of the track and the buffer
        temp = min([length(track_left) params.buffer])-1;
        % Consider a "buffer" amount of points before the detachment (final segment).
        % If the track is too short, it considers the whole track.
        segment_left = track_left(end-temp:end);
        % Calculates the average intensity of the whole track
        mean_int(i) = mean(track_left);
        % Calculates the average intensity of the final segment
        end_int(i) = mean(segment_left);
        % If the particle disappeared and a segment as long as "buffer" was
        % collected after desappearance, then the difference between the
        % intensity pre- and post-detachment is calculated
        if max(tracks{i}(:,6))==params.buffer
            track_right = tracks{i}(tracks{i}(:,6)~=0, 5);
            segment_right = track_right(1:params.buffer);
            bkg_int(i) = mean(segment_right);
            bkg_std(i) = std(segment_right);
            step_amplitude(i) = end_int(i) - bkg_int(i);
            valid = true;
        else
            % Else NaN is used.
            bkg_int(i) = NaN;
            bkg_std(i) = NaN;
            step_amplitude(i) = NaN;
        end
    end
    
end

% If no detachment is detected, the analysis is not performed
if ~valid
    warning('No detachment event observed, so no bleaching analysis was performed');
    for i=1:size(track_props,2)
        track_props(i).bleached = 0;
    end
    bleach_tracks = [];%tracks;
    bleach_props = [];%track_props;
    step_amplitude = [];
    return
end
    

% If the minumum amplitude of the intensity drop is not set, 3 times the
% average (over all tracks) of the standard deviation of the background 
% (intensity after detachment is used.

if isnan(params.minStep)
    min_step = 3*nanmean(bkg_std);
else
    min_step = params.minStep;
end

% figure
% plot(mean_int,slope,'.')
% figure
% plot(mean_int,end_int,'.')
% figure
% plot(1:i,slope,'.')

% Finds bleached particles as the ones in which the intensity drop is below
% the threshold.

bleach_ind = step_amplitude<min_step;

% Creats variables containing the tracks and the track info for all the
% particles considered to be bleaching.
if any(bleach_ind)
    bleach_tracks = tracks(bleach_ind);
    bleach_props = track_props(bleach_ind);
else
    disp('No bleaching events detected');
    bleach_tracks = {};
    bleach_props = [];
end

% Adds a column to track_props array indicating if the track has bleached.
% This will used to ignored bleached tracks in the kinetics analysis.
for i=1:size(track_props,2)
    track_props(i).bleached = bleach_ind(i);
end

% Plots a graph off all particles that desappeared during the video in a
% space with average intensity of the particle at the end of the track on
% the y-axis and the intensity drop on the x-axis. A vertical balck line
% indicates the threhold used to setect bleaching and the tracks
% disregarded as bleached are shown in red.
figure
semilogy(step_amplitude(~bleach_ind), mean_int(~bleach_ind),'.');
hold on
semilogy(step_amplitude(bleach_ind), mean_int(bleach_ind),'.r','MarkerSize',12);

x_min = min([step_amplitude min_step/2])/1.1;
x_max = max(step_amplitude)*1.1;
lim_x = [x_min, x_max];

x_guide = linspace(x_min,x_max,50);
line_guide = x_guide+nanmean(bkg_int);

plot(x_guide, line_guide, '--k');

mean_int_detech = mean_int(~isnan(step_amplitude));

lim_y = prctile(mean_int, [2 98]);
lim_y = [lim_y(1)/1.5 lim_y(2)*1.2];

semilogy([min_step, min_step], lim_y,'-k');

%semilogy(lim_x, [min_int min_int],'-k');
xlim(lim_x);
ylim(lim_y);

xlabel('Intenisty drop after detachment');
ylabel('Intensity end track');
title('Tracks bleaching');

drawnow



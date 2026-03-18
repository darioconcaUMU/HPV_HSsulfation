% This script should be used to perform the kinetic analysis of a single
% TIRF video

clear
close all

% Do you want to process all files in the selected directory?
batch_processing = true;
% Do you want to divide image in several sub-areas to speed up the analysis
% This might affect the analysis.
% Insert a single number n to divide the image in nxn subimages, insert a
% vector [n, m] to divide as nxm. If no subdivision is wanted, enter 1.
parallel_areas = 4;

% smart_param
smart_param = false;

% Double expoential analysis
double_exp = 1;

%% Parameters initialisation
%% Image analysis and particle tracking.
% Parameters that should be adjusted by the user depending on the particles
% and conditions used.

% Directory the script will open when run
defaultDirectory = 'C:\Users\daco0028\Documents\Experiments\Data\2022\Microscopy\03032022\Registered';

% (real) time between consecutive frames [unit: s].
param2.dt = 1;

% The following parameters are probably fine for multiple applications.
% however, they might need a tweek in special occasions.

% Starting frame
param2.istart = 1;

% Min track length
param2.minTrackLength = 3;

% Maximum distance allowed for linking particles in two successive frames
% in pixels
param2.maxJumpDistance = 3;

% Minimum area of a detected spot (in pixels)
param2.minSpotSize = 3;

% Number of frames recorded after detachment
param2.buffer = 10;

%% Blinking correction

% Maximum alloed distance between the centroids of the tracks
param2.maxDistanceLink = 2;

% Maximum time inteval allowed between the end of the first trace and the
% start of the second one
param2.maxDelay = 5;%round(150)/param.dt;

% Currently NOT in use. Compares the intensity of the two traces to
% undestand is they are created by the same particle
param2.intensityLink = 0.5;

%% Bleaching correction

% Maximum normalised slope required for a track to be excluded
param2.bleachSlope = -0.02;

% Maximum normalised intensity at the end of the track required for a track
% to be excluded
param2.bleachIntensity = 10;

param2.minStep = NaN;

%% Calls functions and loops through files

old_dir_2 = cd;

if batch_processing
    
    % Display a dialog box from where to the image-file to be analyzed is chosen.
    % Default directory is User-dependent.
    
    selPath = uigetdir(defaultDirectory, 'Select folder');
    if selPath == 0
        warning("No folder was selected. Analisis is cancelled");
        return
    end
    selPath = strcat(selPath, '\');
    cd(selPath);
    
    files = dir('*SpotSize*.mat');
    fileNameMAT = {files(:).name};
    
    if isempty(fileNameMAT)
        warning("No file was selected. Analisis is cancelled");
        return
    end
    
else
    [file,selPath] = ...
       uigetfile('*SpotSize*','Select the TIFF image-file', defaultDirectory);
    if file == 0
       warning("No file was selected. Analisis is cancelled");
       return
    end
    fileNameMAT{1} = file;
end

for i=1:length(fileNameMAT)
    
    cd(selPath);
    
    close all
    clear spot_record metadata tracks track_props duration cum_duration...
       arrival cum_arrival particle_frame frame_props surf_coverage
    
    currentfileMAT = fileNameMAT{i};
    
    load(currentfileMAT);
    
    selpath = selPath;
    param = param2;
    
    metadata.deltaT = param.dt;
    
    cd(old_dir_2);
    
    [tracks, track_props] = separate_tracks(spot_record);
    
    % Remove too short tracks
    [tracks, track_props] = filter_duration(tracks, track_props, param.minTrackLength);
    
    % Calculates surface coverage
    [particle_frame, frame_props, surf_coverage] = surface_coverage_stuck(track_props, metadata);
    
    % Blinking correction
    [tracks, track_props, links] = blink_correction_4(tracks, track_props, ...
        [param.maxDistanceLink param.maxDelay param.intensityLink]);
    
    % tracks doubling correction
    [tracks, track_props] = remove_duplicates_2(tracks, track_props, param.maxDistanceLink*2);
    
    % Bleaching correction
    [tracks, track_props, bleach_tracks, bleach_props, step_size] = ...
        bleach_correction_step(tracks, track_props, param);
    
    % Analyze remaining tracks
     [duration,cum_duration,arrival,cum_arrival,gof_analysis,...
         fit_dissociation,fit_association] = track_analysis(track_props,metadata);
     
    % Analysis dissociaiton with double exponential as well
    if double_exp == 1
        if length(cum_duration)>5
            [fit_dissociation_2exp, gof_2edx] = doubleExpFit(cum_duration, duration, metadata, fit_dissociation);
        else
            disp('Insufficient data points for double exponential fit');
        end
    end
    
    cd(selpath);
    
    clear selPath param2
    
    %save(strcat(matFileName, '.mat'));
    save(currentfileMAT);
    
    selPath = selpath;
    param2 = param;
    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        %FigName = get(FigHandle, 'Name');
        savefig(FigHandle, strcat(currentfile(1:end-4), '_', num2str(iFig), '.fig'));
        saveas(FigHandle, strcat(currentfile(1:end-4), '_', num2str(iFig), '.png'));
    end
    
    cd(old_dir)
end

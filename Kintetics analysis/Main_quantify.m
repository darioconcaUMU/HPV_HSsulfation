% This script should be used to perform the kinetic analysis of a single
% TIRF video

clear
close all

% Do you want to process all files in the selected directory?
 batch_processing = true;
% batch_processing = true;

% Do you want to divide image in several sub-areas to speed up the analysis
% This might affect the analysis.
% Insert a single number n to divide the image in nxn subimages, insert a
% vector [n, m] to divide as nxm. If no subdivision is wanted, enter 1.

 parallel_areas = 1;

%% Parameters initialisation
%% Image analysis and particle tracking.
% Parameters that should be adjusted by the user depending on the particles
% and conditions used.

% Directory the script will open when run
defaultDirectory = 'C:\Users\k_sei\Documents\MATLAB\TIRF analysis -new version\Input';

% (real) time between consecutive frames [unit: s].
 param.dt = 0.2;

% High threshold used to detect new particles

 param.threshold_high = 30;


% Low threshold used to detect particle detachment

 param.threshold_low = 20;

% The following parameters are probably fine for multiple applications.
% however, they might need a tweek in special occasions.

% Starting frame
param.istart = 1;

% Min track length
param.minTrackLength = 10;

% Maximum distance allowed for linking particles in two successive frames
% in pixels
param.maxJumpDistance = 6;

% Minimum area of a detected spot (in pixels)
param.minSpotSize = 3;

%% Blinking correction

% Maximum allowed distance between the centroids of the tracks
param.maxDistanceLink = 5;

% Maximum time inteval allowed between the end of the first trace and the
% start of the second one
param.maxDelay = 5;

% Currently NOT in use. Compares the intensity of the two traces to
% undestand is they are created by the same particle
param.intensityLink = 0.5;

%% Calls functions and loops through files

old_dir = cd;

if batch_processing
    
    % Display a dialog box from where to the image-file to be analyzed is chosen.
    % Default directory is User-dependent.
    
    selpath = uigetdir(defaultDirectory, 'Select folder');
    if selpath == 0
        warning("No folder was selected. Analisis is cancelled");
        return
    end
    selpath = strcat(selpath, '/');
    cd(selpath);
    
    files = dir('*.tif');
    fileName = {files(:).name};
    
    if isempty(fileName)
        warning("No file was selected. Analisis is cancelled");
        return
    end
    
else
    [file,selpath] = ...
        uigetfile('*.tif','Select the TIFF image-file', defaultDirectory);
    if file == 0
        warning("No file was selected. Analisis is cancelled");
        return
    end
    fileName{1} = file;
end

cd(old_dir);

for i=1:length(fileName)
    
    close all
    clear spot_record metadata tracks track_props duration cum_duration...
        arrival cum_arrival particle_frame frame_props surf_coverage
    
    currentfile = fileName{i};
    
    % Image analysis
    [spot_record, metadata, matFileName] = ...
        tirf_analyzer_2v4_quantification(param, currentfile, selpath, parallel_areas);
    cd(old_dir);
    %metadata.deltaT = param.dt;
    [tracks, track_props] = separate_tracks(spot_record);
    
%     % Blinking correction
%     [tracks, track_props, links] = blink_correction_2(tracks, track_props, ...
%         [param.maxDistanceLink param.maxDelay param.intensityLink]);
    
    % Remove too long tracks
    [tracks, track_props] = filter_duration_short(tracks, track_props, param.minTrackLength);
    
    for j=1:length(track_props)
        track_props(j).bleached = 0;
    end
    
    % Analyze remaining tracks
    [duration,cum_duration,arrival,cum_arrival] = track_analysis(track_props,metadata);
    
    cd(selpath);
    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        %FigName = get(FigHandle, 'Name');
        savefig(FigHandle, strcat(currentfile(1:end-4), '_', num2str(iFig), '.fig'));
        saveas(FigHandle, strcat(currentfile(1:end-4), '_', num2str(iFig), '.png'));
    end
    
    close all
    
    save(strcat(matFileName, '.mat'));
    
    cd(old_dir);
end
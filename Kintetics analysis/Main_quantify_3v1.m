% This script should be used to qunatify the particle concentration from
% the number of particles "bouncing" on the surface, thus entering and
% exiting the TIRF volume.
% The script works similarly to Main_v3_0, but only considers short lived
% tracks instead of long ones. No surface coverage or blinking is performed
% as only short lived spots are considered.
% The script produces a graph for attachment and one for dissociation. The
% former has to be interpreted as cumulative number of tracks detected and
% it is linearly proportional to the particle concentraiton in solution, as
% the probability for particles to move into the TIRF volume is dependent
% on the concentration. The latter should be used as a control that the
% analysis parameters are correct. The length of the detected tracks should
% decrease exponentially and be negligible by 3-4 frames as the probability
% that a particle stays in the TIRF volume while diffusing for extended
% period of time is very small. This also depends on the acquisition frame
% rate (suggested 5-10 fps).

clear
close all

% Do you want to process all files in the selected directory?
% batch_processing = false;
batch_processing = true;

% Do you want to divide image in several sub-areas to speed up the analysis
% This might affect the analysis.
% Insert a single number n to divide the image in nxn subimages, insert a
% vector [n, m] to divide as nxm. If no subdivision is wanted, enter 1.

parallel_areas = 3;

% smart_param
smart_param = false;

%% Parameters initialisation
%% Image analysis and particle tracking.
% Parameters that should be adjusted by the user depending on the particles
% and conditions used.

% Directory the script will open when run
defaultDirectory = 'D:\Bouncing HPV16488 on POPCPEG';

% (real) time between consecutive frames [unit: s].
param.dt = 1;

% High threshold used to detect new particles
prominence_high = 50;% 30 40 50 70 100 150 200 300 400 500 700];
adjust_iter = 1; % Automatically adjusts the threshold if too low

% Low threshold used to detect particle detachment
param.prominence_low = 0.4;

% The following parameters are probably fine for multiple applications.
% however, they might need a tweek in special occasions.

% Starting frame
param.istart = 1;
param.iend = [];

% Min track length
param.minTrackLength = 10;

% Maximum distance allowed for linking particles in two successive frames
% in pixels
param.maxJumpDistance = 6;

% Minimum area of a detected spot (in pixels)
param.minSpotSize = 3;

% Number of frames recorded after detachment
param.buffer = 10;

%% Blinking correction

% Maximum allowed distance between the centroids of the tracks
param.maxDistanceLink = 2;

% Maximum time inteval allowed between the end of the first trace and the
% start of the second one
param.maxDelay = 2;

% Currently NOT in use. Compares the intensity of the two traces to
% undestand is they are created by the same particle
param.intensityLink = 0.5;

%% Bleaching correction

% Maximum normalised slope required for a track to be excluded
param.bleachSlope = -0.02;

% Maximum normalised intensity at the end of the track required for a track
% to be excluded
param.bleachIntensity = 10;

param.minStep = NaN;

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

    files = [dir('*.tif'); dir('*.tiff'); dir('*.nd2')];
    fileName = {files(:).name};

    if isempty(fileName)
        warning("No file was selected. Analisis is cancelled");
        return
    end

else
    [file,selpath] = uigetfile( ...
        {'*.tif;*.tiff;*.nd2', 'Image files (*.tif, *.tiff, *.nd2)'; ...
        '*.tif;*.tiff',      'TIFF files (*.tif, *.tiff)'; ...
        '*.nd2',             'ND2 files (*.nd2)'}, ...
        'Select image file', defaultDirectory);
    if file == 0
        warning("No file was selected. Analisis is cancelled");
        return
    end
    fileName{1} = file;
end

cd(old_dir);

for i=1:length(fileName)

    currentfile = fileName{i};

    if strcmp(file_in(end-2:2),'nd2')
        reader = bfGetReader(file_in);
    end

    for j = prominence_high

        close all
        clear spot_record metadata tracks track_props duration cum_duration...
            arrival cum_arrival particle_frame frame_props surf_coverage

        param.prominence_high = j;

        

        if smart_param == true

            [~,~,ext] = fileparts(currentfile);
            if strcmp(ext,'.tif') || strcmp(ext,'.tiff')

            file_in=strcat(selpath,currentfile);
            I_int = imread(file_in,1);

            h=fspecial('gaussian',[3 3],1);
            I_int = imfilter(double(I_int),h,'same','conv');

            prom_level = zeros(100,1);
            num_peaks = prom_level;
            wid = prom_level;
            pro = prom_level;
            for k = 1:100
                prom_level(k) = 5*param.prominence_high/100*k;
                [~,~,wdt,pmn] = find_peaks_2D(I_int, prom_level(k));
                num_peaks(k) = length(wdt);
                wid(k) = mean(wdt);
                pro(k) = mean(pmn);
            end
            figure
            subplot(3,1,1)
            semilogy(prom_level, num_peaks, '-o');
            subplot(3,1,2)
            plot(prom_level, wid, '-o');
            subplot(3,1,3)
            plot(prom_level, pro, '-o');

            param.prominence_high = input('Please type the prominence required for new particles (Higher threshold): ');
            temp_prom = input('Please type the prominence required for an already detected particle (Lower threshold): ');
            param.prominence_low = temp_prom/param.prominence_high;
            close

            else
                disp('Smart parameters only available for TIFF files');
            end

        end

        % Temporary patch

        file_in=strcat(selpath,currentfile);
        if strcmp(file_in(end-2:end),'nd2')
            I_int = bfGetPlane(reader, 1);
        else
            I_int = imread(file_in,1);
        end
        I_int = imgaussfilt(double(I_int),1,'Padding','symmetric');
       
        if adjust_iter
            max_iterations = 30;
        else
            max_iterations = 1;
        end

        flag_high = 1;
        flag_low = 1;
        iteraction = 1;
    
        while iteraction <= max_iterations

            [~,~,wdt,pmn] = find_peaks_2D(I_int, j);
            num_peaks = length(wdt);
            if num_peaks>150
                param.prominence_high(i) = j+10;
                j = param.prominence_high(i);
            else
                flag_low = 0;
            end
            if max(I_int(:))-mean(I_int(:)) < j/3
                param.prominence_high(i) = round(j/20)*10;
                j = param.prominence_high(i);
            else
                flag_high = 0;
            end
            if ~flag_low && ~flag_high
                iteraction = 100;
            else
                iteraction = iteraction+1;
            end
        end

        if flag_low
            display(['Too low threshold, skip ' currentfile '_prom_' num2str(j)]);
            continue
        elseif flag_high
            display(['Too high threshold, skip ' currentfile '_prom_' num2str(j)]);
            continue
        else
            display(['File: ' currentfile '. Threshold used: ' num2str(j)]);
        end
        
        % Image analysis
        [spot_record, metadata, matFileName] = ...
            tirf_analyzer_3v1(param, currentfile, selpath, parallel_areas);
        cd(old_dir);

        if ~isempty(spot_record)

            %metadata.deltaT = param.dt;
            [tracks, track_props] = separate_tracks(spot_record);

            % Blinking correction
            [tracks, track_props, links] = blink_correction_4(tracks, track_props, ...
                [param.maxDistanceLink param.maxDelay param.intensityLink]);

            % tracks doubling correction
            [tracks, track_props] = remove_duplicates_2(tracks, track_props, param.maxDistanceLink*2);

            % Remove too long tracks
            [tracks, track_props] = filter_duration_short(tracks, track_props, param.minTrackLength);

            for k=1:length(track_props)
                track_props(k).bleached = 0;
            end

            % Analyze remaining tracks
            [duration,cum_duration,arrival,cum_arrival] = track_analysis(track_props,metadata);

            cd(selpath);

            save(strcat(matFileName, '.mat'));

            FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
            for iFig = 1:length(FigList)
                FigHandle = FigList(iFig);
                %FigName = get(FigHandle, 'Name');
                savefig(FigHandle, strcat(currentfile(1:end-4), '_prom_', num2str(j), '_' , num2str(iFig), '.fig'));
                saveas(FigHandle, strcat(currentfile(1:end-4), '_prom_', num2str(j), '_', num2str(iFig), '.png'));
            end
            cd(old_dir)
        end
    end
end
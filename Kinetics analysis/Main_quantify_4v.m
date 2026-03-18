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
prominence_high = 20;

% Low threshold used to detect particle detachment
param.prominence_low = 0.4;

% The following parameters are probably fine for multiple applications.
% however, they might need a tweek in special occasions.

% Starting frame
param.istart = 1;
param.iend = 100;

% Min track length
param.minTrackLength = 10;

% Maximum distance allowed for linking particles in two successive frames
% in pixels
param.maxJumpDistance = 3;

% Minimum area of a detected spot (in pixels)
param.minSpotSize = 3;

% Number of frames recorded after detachment
param.buffer = 0;

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
    if strcmp(file_in(end-2:2),'nd2')
        reader = bfGetReader(file_in);
    end

    for j = prominence_high

        close all
        clear spot_record metadata tracks track_props duration cum_duration...
            arrival cum_arrival particle_frame frame_props surf_coverage

        param.prominence_high = j;

        currentfile = fileName{i};

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

        [~,~,wdt,pmn] = find_peaks_2D(I_int, j);
        num_peaks = length(wdt);

        if num_peaks>150
            display(['Too low threshold, skip ' currentfile '_prom_' num2str(j)]);
            continue
        end
        if max(I_int(:))-mean(I_int(:)) < j/3
            display(['Too high threshold, skip ' currentfile '_prom_' num2str(j)]);
            continue
        end

        % Image analysis
        [spot_record, metadata, matFileName] = ...
            tirf_analyzer_3v1(param, currentfile, selpath, parallel_areas);
        cd(old_dir);

        spot_record = spot_record(spot_record(:,6) == 0, :);

        if ~isempty(spot_record)

            %metadata.deltaT = param.dt;
            [tracks, track_props] = separate_tracks(spot_record);

            % Blinking correction
            [tracks, track_props, links] = blink_correction_4(tracks, track_props, ...
                [param.maxDistanceLink param.maxDelay param.intensityLink]);

            % tracks doubling correction
            [tracks, track_props] = remove_duplicates_2(tracks, track_props, param.maxDistanceLink*2);

            %% ---- Long-track classification + per-spot table + plots (NEW) ----

            % Track lengths (in frames / detections)
            trackLengths = zeros(length(tracks),1);
            for kk = 1:length(tracks)
                trackLengths(kk) = size(tracks{kk},1);
            end

            % Define "long" tracks (surface-bound etc.)
            % NOTE: param.minTrackLength is used in this script as the long-track cutoff
            isLongTrack = trackLengths > param.minTrackLength;

            % Build a spot-level table:
            % Each row = one detection (one spot in one frame)
            % Uses the track identity AFTER blink-correction / dedup (kk index)
            allFirst  = [];
            allFrame  = [];
            allX      = [];
            allY      = [];
            allI      = [];
            allDrop   = [];
            allTrack  = [];
            allProm   = [];
            allIsLong = [];

            for kk = 1:length(tracks)
                tr = tracks{kk};

                % Defensive: handle older spot_record formats (7 cols) vs newer (8 cols)
                if size(tr,2) >= 8
                    prom = tr(:,8);
                else
                    prom = NaN(size(tr,1),1);
                end

                allFirst  = [allFirst; tr(:,1)]; %#ok<AGROW>
                allFrame  = [allFrame; tr(:,2)]; %#ok<AGROW>
                allX      = [allX;     tr(:,3)]; %#ok<AGROW>
                allY      = [allY;     tr(:,4)]; %#ok<AGROW>
                allI      = [allI;     tr(:,5)]; %#ok<AGROW>
                allDrop   = [allDrop;  tr(:,6)]; %#ok<AGROW>
                allTrack  = [allTrack; repmat(kk, size(tr,1), 1)]; %#ok<AGROW>
                allProm   = [allProm;  prom]; %#ok<AGROW>
                allIsLong = [allIsLong; repmat(isLongTrack(kk), size(tr,1), 1)]; %#ok<AGROW>
            end

            spotTable = table( ...
                allFrame, allFirst, allX, allY, allI, allProm, allTrack, logical(allIsLong), ...
                'VariableNames', {'frame','firstFrame','x','y','intensity','prominence','trackId','isLongTrack'} );

            % Counts per frame: total vs excluding long tracks
            minF = min(spotTable.frame);
            maxF = max(spotTable.frame);
            frames = (minF:maxF)';

            nTotal = accumarray(spotTable.frame - minF + 1, 1, [numel(frames),1], @sum, 0);
            idxShort = ~spotTable.isLongTrack;
            nExclLong = accumarray(spotTable.frame(idxShort) - minF + 1, 1, [numel(frames),1], @sum, 0);

            countsTable = table(frames, nTotal, nExclLong, ...
                'VariableNames', {'frame','nSpotsTotal','nSpotsExcludingLongTracks'});

            %% --- Plot 1: Spots per frame before/after excluding long tracks ---
            fig1 = figure('Name','SpotsPerFrame_BeforeAfterLongTrackExclusion');
            plot(countsTable.frame, countsTable.nSpotsTotal, '-'); hold on;
            plot(countsTable.frame, countsTable.nSpotsExcludingLongTracks, '-');
            xlabel('Frame'); ylabel('# spots / frame');
            legend({'Total','Excluding long tracks'}, 'Location','best');
            title(sprintf('Spots per frame (prom=%g): before vs after long-track exclusion', param.prominence_high));
            grid on;

            %% --- Plot 2: Track-length distribution ---
            fig2 = figure('Name','TrackLengthDistribution');
            histogram(trackLengths);
            xlabel('Track length (frames)');
            ylabel('Count');
            title(sprintf('Track-length distribution (prom=%g)', param.prominence_high));
            grid on;

            %% ---- Save outputs (PNG + FIG + MAT) ----
            cd(selpath);

            baseName = sprintf('%s_prom_%g', currentfile(1:end-4), param.prominence_high);

            % Save plots
            saveas(fig1, [baseName '_spotsPerFrame.png']);
            savefig(fig1, [baseName '_spotsPerFrame.fig']);

            saveas(fig2, [baseName '_trackLengthHist.png']);
            savefig(fig2, [baseName '_trackLengthHist.fig']);

            % Save results (single .mat containing tables + key outputs)
            out = struct();
            out.currentfile = currentfile;
            out.param = param;
            out.metadata = metadata;
            out.matFileNameFromAnalyzer = matFileName;

            out.spot_record = spot_record;     % raw analyzer output
            out.tracks = tracks;               % after separate/blink/dedup (as in your pipeline)
            out.track_props = track_props;     % from your pipeline
            out.trackLengths = trackLengths;
            out.isLongTrack = isLongTrack;

            out.spotTable = spotTable;         % per-spot table (frame/prominence/isLong etc.)
            out.countsTable = countsTable;     % per-frame counts before/after long exclusion

            save([baseName '_MainQuantifyOut.mat'], 'out', '-v7.3');

            cd(old_dir);

            % (Optional) if you still want your existing track_analysis outputs,
            % you can compute them here WITHOUT removing long tracks, e.g.:
            % [duration,cum_duration,arrival,cum_arrival] = track_analysis(track_props,metadata);
            % and then save to out as well.

            cd(old_dir)
        end
    end
end
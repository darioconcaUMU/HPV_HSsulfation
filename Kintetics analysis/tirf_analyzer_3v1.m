%% Introduction and file initialization
% Program for kinetic analysis of TIRF data. Written by Olov Wahlsten
% (2015-11-24), and extensively rewritten in 2021 in Bally group.
% Inspiration is taken from Peter J??nsson's software. (2007-10-25)
% Mokhtar Mapar found and fixed a bug relating to spot movements. (2016-06-21)
% The program is capable of:
% - Detecting vesicles bound to a surface.
% - Connect particles in differnt frames and track position and intensity.
%   Although some movement of the vesicles on the surface can be dealt
%   with, the software is not ment for sigle particle tracking.
%
% This script deals with the analysis of the .tif file containing the TIRF
% images. The following function should be called to clean the data and
% perform the kinetic analysis.

function [spot_record, metadata, matFileName] = tirf_analyzer_3v1(param, fileName, pathName, num_sec)

tic

if nargin == 3
    num_sec = 1;
end

script_version = 'tirf_analyser_3v1';

cd(pathName);

dataFileName = fileName(1:end-4);% for naming generated files/images
file_in=strcat(pathName,fileName);

%Booleans for troubleshooting
frameByFramepng = false;    %If true, Saves the current frame (with spots labelled) as a .png file under the folder "processed".
frameByFrameMATs = false;       %If true, the script will save the Workspace (as a .MAT) for EVERY frame.

% Dynamic prominence. If true it addapts the required prominence to the average
% intenisity of the image. If we assume equilibrium, this should reflect
% overall bleaching

dynamicThreshold = false;

%% Variable initialization
% -------------------------------------------------------------------------
% Parameters specified by user:
% (dt, tmin, istart, iend, th_high, th_medium, th_low, spot_rad, rem_size)
% -------------------------------------------------------------------------

dt = param.dt; % dt = (real) time between consecutive frames [unit: s].

% This only works for visualisation purposes. All tracks, even shorter, are
% saved and are then eliminated in a separate function: "filter_duration.m"
tmin = param.minTrackLength; %tmin = minimum bound time for a spot to be considered as a vesicle [unit: frames].

istart = param.istart; %istart = start frame [unit: frame].

iend = length(imfinfo(file_in));

if isfield(param, 'iend') 
    if iend>=param.iend
        iend = param.iend;
    elseif ~isnan(iend) && ~isempty(iend)
        warning('Fewer frames detected than the in the input. Video processed till the last frame');
    end
end

N_frames=iend-istart+1; % Number of frames to analyze.
t=(0:N_frames-1); % Time points for all frames, in frames.

% New spots are considered only if the produce a pick with prominence
% exciding this value

prom_high = param.prominence_high;

% The particle is considered to be detached if the prominence goes below
% this percentage of prominence_high;
prom_low = prom_high*param.prominence_low;

% The maximum distance between spots in two consecutive frames to still be considered as
% the same spot. This needs to be lowered for high density films.
maximumJumpDistance = param.maxJumpDistance;

% This should not be needed might have a look at the width of the peak...
% % Disregards spots that are composed by less than minimumSpotSize.
% % minimumSpotSize=2 seems to not function properly. Choose at least minimumSpotSize=3. I recommend 4 (KT)
minimumSpotSize = param.minSpotSize;

% If delayedDrop = true, the script does not immediately drops particles if
% they fall below the low threshold. It instead waits for a consecutive
% delay_offset frames. This helps in maitaining paticles with strong
% fluctuation in the intensity.

delayedDrop = true;
delay_offset = 1;

% Records the intensity several frames after dropping a particle to make
% sure a detachment and not bleaching occurred

delay_buffer = param.buffer;

% Divide the image in smaller section to speed up the analysis, which can
% be quite slow for high-coverage. Input here the number of sections you
% want the image to be divided into in the x and y direction (Max 9 section
% total)

if length(num_sec) == 1
    sub_sec_x = num_sec;
    sub_sec_y = num_sec;
else
    sub_sec_x = num_sec(1);
    sub_sec_y = num_sec(2);
end

num_seg = sub_sec_x*sub_sec_y;

% Set parameter for image analysis
SE = strel('square',2);
%h=fspecial('gaussian',[3 3],1);

% -------------------------------------------------------------------------
%% Matrix initialization and comments

% The information about all the particles detected is contained in the
% structure "spot_record". Each row contains the information regarding a
% single spot in a single frame. The columns are organised in the following
% way:
% 1) Time the particle has first been detected in frames;
% 2) Current frame;
% 3) X position of the spot centroid in the current frame;
% 4) Y position of the spot centroid in the current frame;
% 5) Intensity at the centroid pixel in the current frame;
% 6) Number of consecutive frames the particle centroid has been below the
%    lower threshold previous to the currently analysed. This is used to
%    allow the particle to fall below the threshold and still be considered
%    in the analysis, in case of high noise. Only used if "delayedDropDoes =
%    true".
% 7) Track index. Indicates at which track the spot is attributed to.

% "spot_past", "spot_present" and "spot_dropped are auxiliary matrix with the same
% structure of "spot_record" and are used to combine spots in consecutive
% frames to obtain the tracks.

% spot_last=[];
%
% spot_present=[];
%
% spot_dropped = [];
%
spot_record = [];

%% Different thing

threshold_log = zeros(N_frames,2);

if strcmp(file_in(end-2),'nd2')
    reader = bfGetReader(file_in);
end

for i=1:(iend-istart+1)
    
    % Reads the frame
    
    frame_num = i+istart-1;
    if strcmp(file_in(end-2),'nd2')
        I_raw = bfGetPlane(reader, frame_num);
    else
        I_raw = imread(file_in,frame_num);
    end
    
    if dynamicThreshold
        if i==1
            start_int = sum(I_raw(:));
        else
            prom_high = prom_high*sum(I_raw(:))/start_int;
            prom_low = prom_high*param.prominence_low;
        end
    end
    
    threshold_log(i,:) = [prom_high prom_low];
    
    % Applies the Gaussian filter to the image
    I = imgaussfilt(double(I_raw),1,'Padding','symmetric');
    
    % Divides the image into subsections.
    % Subsections might not be needed but for now I keep them in
    
    if i==1
        
        [xmax, ymax] = size(I);
        ymax_sec = floor(ymax / sub_sec_y);
        ymax = sub_sec_y * ymax_sec;
        
        xmax_sec = floor(xmax / sub_sec_x);
        xmax = sub_sec_x * xmax_sec;
        
        for j=1:sub_sec_y
            for k=1:sub_sec_x
                y_offset(k,j) = (j-1)*ymax_sec+1;
                x_offset(k,j) = (k-1)*xmax_sec+1;
            end
        end
        x_offset = x_offset(:);
        y_offset = y_offset(:);
        
        spot_last_seg = cell(1,length(x_offset));
        
    end
    
    clear I_sec
    
    for j=1:sub_sec_y
        for k=1:sub_sec_x
                %I_sec{k,j} = I((j-1)*ymax_sec+1:j*ymax_sec, (k-1)*xmax_sec+1:k*xmax_sec);
                I_sec{k,j} = I((k-1)*xmax_sec+1:k*xmax_sec, (j-1)*ymax_sec+1:j*ymax_sec);
        end
    end
    
    % Makes it a 1D cell array
    I_sec = I_sec(:);
    
    % I_int does not represent the original image anymore, it is just the
    % int version of I, so after applying the Gaussian filter
    I_int = int16(I);
    
    % Cycles through the segments
    for k = 1:num_seg
        
        spot_pres = [];
        detached = [];
        
        % Pulls some variables out of cell arrays to make the code easier to
        % read. It might improve performance but not in a noticeble way
        
        spot_last = spot_last_seg{k};
        
        % Finds the spots that have a counter equal or above delay
        % offset. This should not be considered in the analysis. Their
        % counter is simply increased to record the background value.
        
        if ~isempty(spot_last)
            
            detached_pos = spot_last(:,6)>=delay_offset;
            detached = spot_last(detached_pos,:);
            spot_last(detached_pos,:) = [];
            
            % Remove particle in which the counter (column 6) has reached the
            % delay buffer. Enough values have been acquired for the
            % background.
            
            detached(detached(:,6)>=delay_buffer,:) = [];
            
        end
    
        [row,column,width,prominence,intensity_pk] = find_peaks_2D(I_sec{k}, prom_low);
        
        select_high = prominence>prom_high;
        centroids = [column(select_high),row(select_high)];
        prominence_spots = prominence(select_high);
        intensity_spots = intensity_pk(select_high);
        
        
        if ~isempty(centroids)

            full_data = [centroids prominence_spots intensity_spots];
           
            % Removes centroids that would cause a out of bound error in the
            % analysis
            
            invalid_centroid = unique([find(full_data(:,1)<=0.5 | full_data(:,1)>=(xmax_sec-0.5)); ...
                find(full_data(:,2)<=0.5 | full_data(:,2)>=(ymax_sec-0.5))]);
            
            full_data(invalid_centroid,:) = [];
            
            if ~isempty(spot_last)
                [spot_pres, full_data, spot_last] = match_spots(full_data, spot_last, maximumJumpDistance, t(i));
            end
            
            % Calculates the number of spots in spot_present
            count = size(spot_pres,1) +1;
            
            % Creates a new track for all the spots that have not been asigned
            % to an already existing one.
           
            % This is used to create a unique identifier for each spot
            % FRAME, SECTION, sequential number for the spots.
            counter = 0;
            
            % Particles are not ordered sequencially as it would violate parfor
            % rules, instead a unique indentifier is used. This should not
            % impact further analysis.
            
            % prominence_spots = prominence_spots(prominence>prom_high);
            
            for j = 1:size(full_data,1)
                
                spot_pres(count,:) = [t(i), t(i), full_data(j,1), full_data(j,2),...
                    full_data(j,4), 0, str2num([num2str(i) num2str(k,'%03.f')...
                    num2str(counter,'%05.f')]), full_data(j,3)];
                counter = counter+1;
                count = count+1;
            end
        end 
            
        %% ADDING DIM SPOTS TO EXISTING TRACKS
        
        % Checks all the spots of the required size that have intensity higher
        % than threshold low.
        % This time only particles already present in the previous frame will
        % be accepted.
        
        % Since this part of the software can only append spots to an
        % existing traces, if no spots are present already, this part
        % can be skipped
        
        if ~isempty(spot_last)
            
            % this time we only look at the peaks that have a
            % prominence lower than the prominence threshold
            centroids = [column(~select_high), row(~select_high)];
            prominence_spots = prominence(~select_high);
            intensity_spots = intensity_pk(~select_high);

            if ~isempty(centroids)

                full_data = [centroids prominence_spots intensity_spots];
                
                % Removes centroids that would cause a out of bound error in the
                % analysis
                invalid_centroid = unique([find(full_data(:,1)<=0.5 | full_data(:,1)>=(xmax_sec-0.5)); ...
                    find(full_data(:,2)<=0.5 | full_data(:,2)>=(ymax_sec-0.5))]);
                full_data(invalid_centroid,:) = [];
                
                [spot_temp,~,spot_last] = match_spots(full_data, spot_last, maximumJumpDistance, t(i));
                
                % Appends the dim particles to spot_present
                spot_pres = [spot_pres; spot_temp];
            end
            
            % Appends the spots which have surpassed delay offset, which
            % had been eliminated from the analysis in the beginning of the
            % loop
            
            spot_last = [spot_last; detached];
            
            % At this point "spot_last" should be composed by all the
            % particles that have not been linked to a spot in the current
            % image
          
        else  
            spot_last = detached;
        end
        
        if ~isempty(spot_last)

            % Int will here contain the current pixel values at each current centroid position (rounded).
            rr = round(spot_last(:,4));   % y (row)
            cc = round(spot_last(:,3));   % x (col)

            % clamp to bounds to avoid out-of-range
            rr = max(1, min(size(I_sec{k},1), rr));
            cc = max(1, min(size(I_sec{k},2), cc));

            Int = I_sec{k}(sub2ind(size(I_sec{k}), rr, cc));


            remaining_spots = size(spot_last,1);

            spot_last(:,6) = spot_last(:,6)+1;
            spot_last(:,2) = t(i)*ones(remaining_spots,1);
            spot_last(:,5) = Int;
            
            dropped = spot_last(:,6) == delay_offset;
        
            spot_pres = [spot_pres; spot_last];
            
            % Deletes spots that are below threshold_low from the spots
            % recorded in the previous frame.
            spot_dropped_seg{k} = spot_last(dropped,:);
            
        else          
            spot_dropped_seg{k} = [];
        end
        % Resets spot_present and spot_last
        spot_last_seg{k} = spot_pres;
        spot_pres = [];
        
    end  
    
    % Collect all particle position in a single variable, adjusting the
    % position
    
    spot_present = [];
    spot_dropped = [];
    
    for k = 1:num_seg
        temp_present = spot_last_seg{k};
        if ~isempty(temp_present)
            temp_present(:,3) = temp_present(:,3) + y_offset(k);
            temp_present(:,4) = temp_present(:,4) + x_offset(k);
            spot_present = [spot_present; temp_present];
        end
        
        temp_dropped = spot_dropped_seg{k};
        if ~isempty(temp_dropped)
            temp_dropped(:,3) = temp_dropped(:,3) + y_offset(k);
            temp_dropped(:,4) = temp_dropped(:,4) + x_offset(k);
            spot_dropped = [spot_dropped; temp_dropped];
        end
    end
    
    % Appends spot present to the full list of spots detected in the
    % whole video
    
    spot_record = [spot_record; spot_present];  
    
    %% Commands to display the frame and the analysis.
    
    % Displays which frame is being analysed on the command window
    if mod(i,50)==0
        disp(['Analysed frame: ',num2str(i),' of ',num2str(N_frames)]);
    end
    
    
    % Shows the current image with the identified spots with a cross over
    % their respective centroid position.
    
    figure(1); imshow(I_int,mean(I_int(:))+[-prom_low, prom_high])
    title(sprintf(['Frame number: ',num2str(istart+i-1)]))
    hold on
    
    % In case spots were detected in the current frame. It plots them as
    % crosses. Red if firmly attached (observed for more than tmin frames)
    % and yellow otherwise
    
    if ~isempty(spot_present)
        
        spot_present = spot_present(spot_present(:,6)<delay_offset,:);
        % t_bound is a column vector stating the time present spots have
        % been bound to the surface so far.
        t_bound_present=spot_present(:,2)-spot_present(:,1);
        
        % If the time bound (t_bound) is shorter than (or equal to) tmin
        % it is considered as loosely attached and is put in the
        % spots_loose_present matrix.
        
        spots_loose_present=spot_present(t_bound_present<=t(tmin),3:4);
        
        % If the time bound (t_bound) is longer than tmin,
        % it is considered as firmly attached and is put in the spots_firm_present
        % matrix.
        
        spots_firm_present=spot_present(t_bound_present>t(tmin),3:4);
        
        % Detected spots that have been attached to the surface shorter than tmin
        % (loosely attached) are marked with yellow crosses.
        % Detected spots that have been attached to the surface longer than tmin
        % (firmly attached) are marked with red crosses.
        plot(spots_loose_present(:,1),spots_loose_present(:,2),'y+')
        plot(spots_firm_present(:,1),spots_firm_present(:,2),'r+')
    end
    
    if ~isempty(spot_dropped) % This is true if spot_past contains info.
        
        % t_bound_past is a column vector stating the time the spot that
        % now was dropped has staied in the video
        
        t_bound_past=spot_dropped(:,2)-spot_dropped(:,1);
        
        % Find the loosely attached spots that have detached.
        % The time above th_medium has to be smaller or equal to tmin
        spots_off_loose=spot_dropped(t_bound_past<=t(tmin),3:4);
        
        % Find the firmly attached spots that have detached.
        % The time above th_medium has to be greater than tmin AND the
        % intensity of the spot before detaching was higher than
        % (threshold_high-threshold_low)/4+threshold_low
        
        bleach_threshold = mean(I_int(:))+(prom_high-prom_low)/4+prom_low;
        spots_off_firm=spot_dropped(t_bound_past>t(tmin) & spot_dropped(:,5) > bleach_threshold,3:4);
        
        % Give an indication of possible bleaching tracks. This should be
        % considered only an indication for the user. This is NOT used in
        % the final analysis in any way.
        % The time above th_medium has to be greater than tmin AND the
        % intensity of the spot before detaching was lower than:
        % (threshold_high-threshold_low)/4+threshold_low
        
        spots_off_bleach=spot_dropped(t_bound_past>t(tmin) & spot_dropped(:,5) <= bleach_threshold,3:4);
        
        % Loosely attached spots that just detached are marked with yellow rings.
        plot(spots_off_loose(:,1),spots_off_loose(:,2),'yo','LineWidth',2)
        % Firmly attached spots that just detached are marked with red rings.
        plot(spots_off_firm(:,1),spots_off_firm(:,2),'ro','LineWidth',2)
        % Firmly attached spots that bleached are marked with green rings.
        plot(spots_off_bleach(:,1),spots_off_bleach(:,2),'go','LineWidth',2)
    end
    hold off
    drawnow
    
    %%%	Saves the current frame (with spots labelled) as a .png file under the folder "processed".
    if frameByFramepng == true
        file_str=strcat('./processed/',dataFileName,' frame',num2str(i),'.png');
        print(file_str,'-dpng')
    end
    
    spot_present = [];
    spot_dropped = [];
    
end

close;

%% Saving .mat files

% Saves useful variables in a mat file.
% Gives the MAT-file a name including the sensitivity and the spotsize (KT)

Sens_high = num2str(prom_high);
Sens_low = num2str(prom_low);
Sensitivity = strcat(Sens_high,", ",Sens_low);
SpotSize = num2str(minimumSpotSize);
matFileName = strcat(dataFileName,' - Sens ',Sensitivity,' - SpotSize ',SpotSize);

metadata.scriptVersion = script_version;
metadata.analysed = file_in;
metadata.deltaT = dt;
metadata.frameRate = 1/dt;
metadata.startFrame = istart;
metadata.totalFrames = iend-istart+1;
metadata.totalFrameVideo = iend;
metadata.saveImageFrame = frameByFramepng;
metadata.saveMATFrame = frameByFrameMATs;
metadata.thresholdHigh = prom_high;
metadata.thresholdLow = prom_low;
metadata.thresholdFrame = threshold_log;
metadata.spotSize = minimumSpotSize;
metadata.maxJump = maximumJumpDistance;
metadata.connectivity = [];
metadata.delayDrop = delay_offset;

save(strcat(matFileName,'.mat'),'metadata','spot_record')

toc

end   


%% Auxiliary function for spots addition to tracks

function [spot_present, full_data, spot_last] = match_spots(full_data, spot_last, max_dist,t)

%%%%%% GENERATE MATRIX WITH DISTANCE OF ALL POINTS TO ALL POINTS
%%%%%% FROM FRAME BEFORE. THEN PICK THE MINIMUM, CONNECT THE TWO
%%%%%% AND DISCARD THEM IN THE MATRIX. REPEAT TILL MIN DIST>THRESHOLD
%centroids = full_data(:,1:2);  %cnetroids are the first two columns of
%full_data

if ~isempty(spot_last)
    
    % Preallocates memory to improve performance
    spot_present = NaN(size(full_data,1), 8);
    count = 1;
    
    % Calculated distance between the all new spots and all the spots
    % in the previous frame.
    % columns are the new spots, rows the old
    [pos_new_x, pos_last_x] = meshgrid(full_data(:,1),spot_last(:,3));
    [pos_new_y, pos_last_y] = meshgrid(full_data(:,2),spot_last(:,4));
    temp_index = spot_last(:,7);
    indeces = NaN(length(temp_index),1);
    dist_matrix = hypot(pos_new_x-pos_last_x, pos_new_y-pos_last_y);
    while true
        
        % For now in case of two spot being exactly at the same
        % distance the software will use the first one encountered.
        % I don't think this is a major issue but we might want to
        % consider the intensity as well in the future..
        
        if isempty(dist_matrix)
            break
        end
        
        % Find all the absolute minimum distance
        if size(dist_matrix,1)>1
            [temp_min, min_pos] = min(dist_matrix);
            [min_value, min_new] = min(temp_min);
            min_last = min_pos(min_new);
        else
            min_pos = 1;
            [min_value, min_new] = min(dist_matrix);
            min_last = 1;
        end
        
        % If the minimum distance is larger than the distance
        % threshold, then no more spots need to be connected and we
        % escape the loop
        if min_value > max_dist
            break
        end
        % intensity_centroid = I(round(centroids(min_new,2)),round(centroids(min_new,1)));
        pos_index = spot_last(:,7) == temp_index(min_last);
        spot_present(count,:) = [spot_last(pos_index,1), t, full_data(min_new,1), full_data(min_new,2),...
            full_data(min_new,4), 0, temp_index(min_last),  full_data(min_new,3)];
        indeces(count) = temp_index(min_last);
        count = count+1;
        
        % Remove spots from distance matrix so it will not be
        % double counted
        
        dist_matrix(min_last,:) = [];
        dist_matrix(:,min_new) = [];
        full_data(min_new,:) = [];
        temp_index(min_last) = [];
    end
    
    % Remove extra rows that were not used
    spot_present(all(isnan(spot_present), 2), :) = [];
    spot_last(ismember(spot_last(:,7),indeces),:) = [];
else
    % If no spots were present in the previous frame, an empty matrix
    % is returned
    spot_present = [];
end

end

function [row,column,width,prominence,intensity] = find_peaks_2D(data, min_prom)

data_size = size(data);
data2 = data';
% creates a 1D array of the matrix in both directions
[pks1, locs1, w1, p1] = findpeaks(double(data(:)), 'MinPeakHeight',min_prom, 'MinPeakProminence', min_prom);  % Do we want to keep the peak height?
[pks2, locs2, w2, p2] = findpeaks(double(data2(:)), 'MinPeakHeight',min_prom, 'MinPeakProminence', min_prom);
% since the second matrix is rotated it converts the position to be
% compatible to the peaks position found in the first matrix
[c2, r2] = ind2sub(fliplr(data_size), locs2);
locs2 = sub2ind(data_size, r2, c2);
% Only allows peaks that are found in both x and y direction
[ind, pos1, pos2] = intersect(locs1, locs2);
width = min([w1(pos1) w2(pos2)],[],2);
prominence = min([p1(pos1) p2(pos2)],[],2);
%prominence = mean([p1(pos1) p2(pos2)]);
% Finds x and y positon of the peaks
[row, column] = ind2sub(data_size, ind);

% --- peak intensity at the same coordinates ---
intensity = double(data(sub2ind(data_size, row, column)));

end


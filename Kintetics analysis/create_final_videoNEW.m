function create_final_videoNEW(track_props, tracks, spot_record, params, file_name, directory)

% This function allows to create an AVI video of the analysis of each
% frame.
% Input:
% - track_props, tracks, spot_record: variables which are created during
% the analysis of the TIFF file.
% - file_name: name of the TIFF file that was analysed
% - directory: directory where the TIFF file is saved.
% The AVI file created by the script is saved in the same folder as the
% script with the name "[file_name]_video.avi". The video shows the
% original video to the left and the video with overlayed signs to show
% which particles are considered in the analysis to the right. The simbols
% used are:
% - cross: detected particle
% - star: indicates that the detected particle will non be present in the
% next frame.
% - circle: last recorded position of a particle that detached on this
% frame
% Color:
% - Red: particle considered as attached in the analysis
% - Yellow: particle ignored during the analysis as the track was too short
% or a duplicate of another already considered particle.
% - Green: particle ignored in the detachment analysis because it bleached
% instad of detaching.

frames = min([track_props(:).beginning]):max([track_props(:).ending]);

track_index = [track_props(:).index];
track_beginning = [track_props(:).beginning];
track_ending = [track_props(:).ending];
track_bleached = [track_props(:).bleached];

% create a table connecting spot indexes with track indexes
index_comp = [];
for i=1:length(track_props)
    for j=1:length(track_props(i).spot_index)
        index_comp = [index_comp; [track_props(i).spot_index(j) track_index(i) track_bleached(i)]];
    end
end

% creates a single table with all spots after filtering
spot_filter = [];
spot_end = [];
for i=1:length(tracks)
    tracks{i} = tracks{i}(tracks{i}(:,6) == 0,:);
    % tracks{i}(:,tracks{i}(:,6)~=0) = [];
    spot_filter = [spot_filter; tracks{i}(1:end-1,:)];
    spot_end = [spot_end; tracks{i}(end,:)];
end

% Create a cell array with all the spots that have been filtered out
% (belonging to too short tracks)

spot_real = spot_record(spot_record(:,6) == 0,:);
index_excluded = setxor(spot_real(:,7), spot_filter(:,7));
spot_excluded = spot_real(ismember(spot_real(:,7),index_excluded),:);

% create two cell arrays with all spots in each frame, both filtered and
% excluded
for i=1:length(frames)
    temp_frame = spot_filter(spot_filter(:,2)==frames(i),:);
    bleached_frame = index_comp(ismember(index_comp(:,1),temp_frame(:,7)),3);
    spot_frame{i}.good = temp_frame(~bleached_frame,:);
    spot_frame{i}.bleach = temp_frame(logical(bleached_frame),:);
    
    spot_frame{i}.exc = spot_excluded(spot_excluded(:,2)==frames(i),:);
    
    temp_end = spot_end(spot_end(:,2)==frames(i),:);
    bleached_end = index_comp(ismember(index_comp(:,1),temp_end(:,7)),3);
    spot_frame{i}.ending = temp_end(~bleached_end,:);
    spot_frame{i}.end_bleach = temp_end(logical(bleached_end),:);
    
    % Find only new arrivals
    new_arrivals_temp = [track_props.beginning] == i-1;
    new_arrival_index = find(new_arrivals_temp);
    new_arrival_frame{i} = [[track_props(new_arrival_index).pos_mean_x];... 
            [track_props(new_arrival_index).pos_mean_y]]';
    
end
spot_frame = spot_frame';

% Create cell array with missing spots for each frame. The position is the
% average position of the track.

% for i=1:length(frames)
%     cond_time = track_beginning<=frames(i) & track_ending>=frames(i);
%     temp = track_index(cond_time);
%     all_particle_frame = index_comp(ismember(index_comp(:,2), temp),1);
%     cond_miss = ~ismember(all_particle_frame,spot_frame{i}.good(:,7));
%     missing_spot_index = all_particle_frame(cond_miss);
%     missing_track_index = index_comp(ismember(index_comp(:,1), missing_spot_index),2);
%     temp_track_props = track_props(ismember([track_props(:).index], missing_track_index));
%     temp_miss = [[temp_track_props(:).pos_mean_x]' [temp_track_props(:).pos_mean_y]'];
%     spot_frame{i}.miss_bleach = temp_miss([temp_track_props(:).bleached]==1,:);
%     spot_frame{i}.miss = temp_miss([temp_track_props(:).bleached]==0,:);
% end

% Diplays the points on the image and save it as an avi file


v = VideoWriter(strcat(file_name, '_video.avi'),'Uncompressed AVI');
v.FrameRate = 5;
open(v)
f = figure(1);
f.Position = [50 50 1400 650];

% Calculate range image
or_dir = cd;
cd(directory)
I_int=imread(file_name,i);
mod_int = mode(I_int(:));
interval_int = [mod_int/2 mod_int+2*params.prominence_high];

for i=1:100%length(frames)
    
    or_dir = cd;
    cd(directory)
    I_int=imread(file_name,i);
    subplot(1,2,1)
    %imshow(I_int, [prctile(I_int(:),30) prctile(I_int(:),95)*5]);
    imshow(I_int, interval_int);
    title(sprintf(['Frame number: ',num2str(frames(i))]))
    
    hold on
    try
        plot(new_arrival_frame{i}(:,1),new_arrival_frame{i}(:,2),'+r');
    end
    if i<length(frames)-1
        try
            plot(new_arrival_frame{i+1}(:,1),new_arrival_frame{i+1}(:,2),'or');
        end
    end
    
    subplot(1,2,2)
%     int_track = [track_props(:).avg_intensity];
%     imshow(I_int,'DisplayRange', [min(int_track)/10 prctile(int_track,90)*5]);
%     log_I = log(I_int - min(I_int(:))+1);
%     imshow(log_I, [0 max(log_I(:))]);
    %imshow(I_int, [prctile(I_int(:),30) prctile(I_int(:),95)*5]);
    imshow(I_int, interval_int);
    title('Detected particles')
    hold on
    try
    plot(spot_frame{i}.bleach(:,3),spot_frame{i}.bleach(:,4),'+g');
    end
    try
    plot(spot_frame{i}.exc(:,3),spot_frame{i}.exc(:,4),'+y');
    end
    try
    plot(spot_frame{i}.good(:,3),spot_frame{i}.good(:,4),'+r');
    end
%     try
%     plot(spot_frame{i}.miss(:,1),spot_frame{i}.miss(:,2),'.r');
%     end
%     try
%     plot(spot_frame{i}.miss_bleach(:,1),spot_frame{i}.miss_bleach(:,2),'.g');
%     end
    try
    plot(spot_frame{i}.ending(:,3),spot_frame{i}.ending(:,4),'*r');
    end
    try
    plot(spot_frame{i}.end_bleach(:,3),spot_frame{i}.end_bleach(:,4),'*g');
    end
    try
    plot(spot_frame{i-1}.ending(:,3),spot_frame{i-1}.ending(:,4),'or');
    end
    try
    plot(spot_frame{i-1}.end_bleach(:,3),spot_frame{i-1}.end_bleach(:,4),'og');
    end
    hold off
    drawnow
    
    cd(or_dir);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);


    
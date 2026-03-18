function props = calculate_properties(track_raw, tr_index, bleached)

% This function takes as input "track", a cell array in which each cell is
% a complete track, i.e. it contains the information of the particle in all
% frames it was detected. See separate tracks for more info.
% The output is "props", a structure containing the most importan general
% information for each track. The fields are:
% 1. "index": An integer unequivocally indentifying the track.
% 2-3. "pos_mean_x" and "pos_mean_y": The mean position of the particle
%   along the whole track in x and y.
% 4. "pos_std": Standard deviation of the displacement of the particle from
%   the first frame.
% 5. "max_disp": The maximum displacement of the particle from the position
%   in the first frame.
% 6-7. "start" and "end": The frame in which the particle was first and last
%   detected
% 8. "duration": Lenght of the track in FRAMES.
% 9. "avg.intensity": Average fluorescence intensity at the particle
%   centroid.

track = track_raw(track_raw(:,6)==0,:);

props.index = tr_index;
props.spot_index = track(1,7);
props.pos_mean_x = mean(track(:,3));
props.pos_mean_y = mean(track(:,4));
displacement = sqrt((track(1,3)-track(:,3)).^2 + ...
    (track(1,4)-track(:,4)).^2);
props.pos_std = std(displacement);
props.max_displ = max(displacement);
props.beginning = track(1,1);
props.ending = track(end,2);
props.duration = props.ending - props.beginning + 1;
props.avg_intensity = mean(track(:,5));
temp = min(size(track,1),3);                                    %% It is calcualting all tracks even if too short. Change when you have time
props.initial_intensity = mean(track(1:temp,5));
props.avg_prominence = mean(track(:,8));
props.initial_prominence = mean(track(1:temp,8));
if nargin>2
   props.bleaching = bleached;
end

function [tracks, track_props] = separate_tracks(spot_record)

track_index = unique(spot_record(:,7));
number_tracks = length(track_index);
tracks = cell(number_tracks, 1);

for i=1:number_tracks
    tracks{i} = spot_record(spot_record(:,7) == track_index(i),:);
    track_props(i) = calculate_properties(tracks{i},i);
end

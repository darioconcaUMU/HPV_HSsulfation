function tracks = extract_track_position(tracks_MSS)

for i=1:size(tracks_MSS,1)
    temp1 = tracks_MSS(i,1:8:end);
    nan_pos = isnan(temp1);
    temp1 = temp1(~nan_pos);
    temp2 = tracks_MSS(i,2:8:end);
    temp2 = temp2(~nan_pos);
    tracks{i} = [temp1' temp2'];
end
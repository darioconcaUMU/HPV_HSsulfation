function summary_segments = compile_subtracks(transDiffResults)

n_tracks = length(transDiffResults);

cumulative_table = [];

for i=1:n_tracks
    data_subtrack = transDiffResults(i).segmentClass.momentScalingSpectrum;
    n_subtracks = size(data_subtrack,1);
    track_ID = i*ones(n_subtracks,1);
    data_subtrack = [data_subtrack track_ID];
    cumulative_table = [cumulative_table; data_subtrack];
end

summary_segments.direct = cumulative_table(cumulative_table(:,3)==3,:);
summary_segments.free = cumulative_table(cumulative_table(:,3)==2,:);
summary_segments.confined = cumulative_table(cumulative_table(:,3)==1,:);
summary_segments.immobile = cumulative_table(cumulative_table(:,3)==0,:);

summary_segments.mean_std_count(:,1:3) = summerise_data(summary_segments.direct);
summary_segments.mean_std_count(:,4:6) = summerise_data(summary_segments.free);
summary_segments.mean_std_count(:,7:9) = summerise_data(summary_segments.confined);
summary_segments.mean_std_count(:,10:12) = summerise_data(summary_segments.immobile);

end

function summary = summerise_data(data)

num_seg = size(data,1);
summary(1,:) = [mean(data(:,2) - data(:,1)), std(data(:,2) - data(:,1)), num_seg]; %duration
summary(2,:) = [mean(data(:,4)) std(data(:,4)), num_seg]; %slope
summary(3,:) = [mean(data(:,19)) std(data(:,19)), num_seg]; % Diffusion coefficient
summary(4,:) = [mean(data(:,20)) std(data(:,20)), num_seg]; % Confinement radius

end

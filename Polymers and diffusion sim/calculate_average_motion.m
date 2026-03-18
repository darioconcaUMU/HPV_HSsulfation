function average_motion = calculate_average_motion(analysis_result)

for i = 1:size(analysis_result,1)
    for j = 1:size(analysis_result,2)
        temp = analysis_result{i,j};
        for k = 1:length(temp)
            temp_motionType = temp(k).segmentClass.momentScalingSpectrum(:,3);
            temp_length = temp(k).segmentClass.momentScalingSpectrum(:,2)-temp(k).segmentClass.momentScalingSpectrum(:,1);
            average_motion(i,j,k) = sum(temp_motionType(~isnan(temp_motionType)).*temp_length(~isnan(temp_motionType)))/sum(temp_length(~isnan(temp_motionType)));
        end
    end
end
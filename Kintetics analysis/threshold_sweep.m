function [num_part, int_range] = threshold_sweep(file_name, interval, resolution)

% Test sevaral threshold values to show how this impacts the number of
% detected particles
% Input:
% - file_name: name of the tiff file to analyse
% - interval: interval of threshold values to be tested
% - resolution: steps in the threshold intensity used for the interval
% sweep
% Output:
% - num_particles: array with the number of particles detected for each
% threshold level tested
% - threslevels: array containing the threshold levels tested
% The script also generates a semi-log plot of num_particles vs threslevels

select_frame = 50;

I_int = imread(file_name,select_frame);

% % In case frame averaging is needed
% I_int = 0;
% for i=1:5
%     temp = imread(file_name,select_frame+i);
%     I_int = I_int + temp;
% end
% I_int = I_int/5;

I = imgaussfilt(double(I_int),1,'Padding','symmetric');

int_range = interval(1):resolution:interval(2);

for i=1:length(int_range)
    temp = find_peaks_2D(I, int_range(i));
    num_part(i) = length(temp);
end

semilogy(int_range, num_part, '.');
xlabel('Threshold value');
ylabel('Number of particles detected');

end


function [row,column,width,prominence] = find_peaks_2D(data, min_prom)

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
% Finds x and y positon of the peaks
[row, column] = ind2sub(data_size, ind);

end
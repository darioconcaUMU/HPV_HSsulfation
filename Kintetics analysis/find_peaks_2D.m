function [row,column,width,prominence] = find_peaks_2D(data, min_prom)

data_size = size(data);
data2 = data';
% creates a 1D array of the matrix in both directions
[pks1, locs1, w1, p1] = findpeaks(double(data(:)), 'MinPeakHeight',min_prom, 'MinPeakProminence',min_prom);  % Do we want to keep the peak height?
[pks2, locs2, w2, p2] = findpeaks(double(data2(:)), 'MinPeakHeight',min_prom, 'MinPeakProminence',min_prom);
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

% % Plots to show that they match
% imagesc(data)
% hold on
% plot(column, row, 'or');

end
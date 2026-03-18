[file,selpath] = uigetfile('*.tif','Select the TIFF image-file');
file_in=strcat(selpath,file);
I_int=imread(file_in,1);

h=fspecial('gaussian',[3 3],1);
SE = strel('square',2);
I = imfilter(double(I_int),h,'same','conv');
min_perc = prctile(I(:), 99);
max_perc = prctile(I(:), 99.995);
k = 1;

for i = linspace(min_perc, max_perc, 100)
    spots_logical = I>i;
    spots_logical = bwareaopen(spots_logical,3,8);
    Lh = imclose(spots_logical,SE);
    s = regionprops(Lh,I,'Area');
    num_part(k) = length(s);
    area(k) = mean([s.Area]);
    thresh(k) = i;
    k = k+1;    
end

figure
subplot(2,1,1)
semilogy(thresh, num_part,'.')
hold on
semilogy(thresh, movmean(num_part,5))
subplot(2,1,2)
semilogy(thresh, area,'.')
hold on
semilogy(thresh, movmean(area,5))

a = 1;


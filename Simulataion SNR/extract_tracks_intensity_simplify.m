function [tracks, intensities] = extract_tracks_intensity_simplify(xmlFile, tiffFile, saveMAT)
% extract_tracks_intensity  Load TrackMate XML + TIFF movie and return
%                           per-track intensity time series at particle
%                           centres (bilinear interpolated).
%
%   [tracks, intensities] = extract_tracks_intensity(xmlFile, tiffFile)
%   returns:
%       tracks      ? 1?N cell array, each cell M?3  [t x y]
%       intensities ? 1?N cell array, each cell M?1  (same order as tracks)
%
%   [...] = extract_tracks_intensity(..., saveMAT) saves a MAT file with
%   variables tracks and intensities when saveMAT is true (default false).
%
%   Requires Image Processing Toolbox for imread of multi?page TIFF.
%
% -----------------------------------------------------------------------
if nargin < 3, saveMAT = false; end
assert(isfile(xmlFile),  'XML file not found.');
assert(isfile(tiffFile), 'TIFF file not found.');

%% Step 1: obtain tracks (reuse previous function if available)
if exist('extract_tracks_hull', 'file') == 2
    tracks = extract_tracks_hull(xmlFile, false);
else
    error('Function extract_tracks_hull.m must be on path.');
end
N = numel(tracks);

%% Step 2: read TIFF metadata to know dimensions
info = imfinfo(tiffFile);
numFrames = numel(info);
imgW = info(1).Width; imgH = info(1).Height;

%% Step 3: load entire stack into memory (uint16/8) ? may be large but fast
fprintf('Reading %d-frame TIFF...\n', numFrames);
stack = zeros(imgH, imgW, numFrames, class(imread(tiffFile, 1)) );
for k = 1:numFrames
    stack(:,:,k) = imread(tiffFile, k);
end

%% Step 4: sample intensity at each detection
intensities = cell(1, N);
for i = 1:N
    T = tracks{i};
    n = size(T,1);
    intensities{i} = zeros(n,1, class(stack));
    for j = 1:n
        
        t  = T(j,1)+1;          % frame index (TrackMate uses 0?based t)
        x  = min(round(T(j,2)/0.183),imgW);            % sub?pixel position
        y  = min(round(T(j,3)/0.183),imgH);
        % x  = min(round(T(j,2)),imgW);            % sub?pixel position
        % y  = min(round(T(j,3)),imgH);
        % Bilinear interpolation (double precision)
        slice = imgaussfilt(stack(:,:,t),1);
        x_low = max(1,x-3);
        x_high = min(x+3,imgW);
        y_low = max(1,y-3);
        y_high = min(y+3,imgH);
        area = slice(y_low:y_high,x_low:x_high);
        intensities{i}(j) = max(area(:));
        % frameDbl = double(stack(:,:,t));
        % intensities{i}(j) = interp2(Xgrid, Ygrid, frameDbl, x, y, 'linear', NaN);
    end
end

%% Optional save
if saveMAT
    [p,f,~] = fileparts(xmlFile);
    matName = fullfile(p, [f '_intensity.mat']);
    save(matName, 'tracks', 'intensities');
    fprintf('Saved intensities to %s\n', matName);
end
end

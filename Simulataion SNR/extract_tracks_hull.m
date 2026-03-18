function [tracks, hullArea, hullMaxDist] = extract_tracks_hull(xmlFile, saveMAT)
% extract_tracks_hull  Extracts particle tracks from a TrackMate XML file
%                      and computes the convex‑hull area and maximal
%                      inter‑point distance for every track.
%
%   [tracks, hullArea, hullMaxDist] = extract_tracks_hull(xmlFile)
%   reads *xmlFile* (TrackMate/track export), returns:
%       tracks      – 1×N cell array, each cell is an M×3 matrix [t x y]
%       hullArea    – 1×N vector, area of 2‑D convex hull (units²)
%       hullMaxDist – 1×N vector, maximum pairwise distance within hull
%
%   [...] = extract_tracks_hull(xmlFile, saveMAT) additionally saves the
%   results to a MAT file when saveMAT is true (default = false). The MAT
%   filename mirrors the XML name, e.g.  "...Tracks.mat".
%
%   The function assumes <particle> elements contain <detection> children
%   with numeric attributes t, x, y. z is ignored.
%
%   Example:
%       [T, A, D] = extract_tracks_hull('simulation_video_SNR3p61_Tracks.xml', true);
%       disp(A(1:5));   % first five hull areas
%
%   Requires: Statistics and Machine Learning Toolbox (for pdist) or uses
%   a fallback to a double loop if pdist is unavailable.
%
% -----------------------------------------------------------------------
% Author: ChatGPT             Date: 29‑Apr‑2025
% -----------------------------------------------------------------------

if nargin < 2, saveMAT = false; end
assert(isfile(xmlFile), 'File not found: %s', xmlFile);

%% Parse XML
try
    doc = xmlread(xmlFile);
catch ME
    error('Failed to read XML: %s', ME.message);
end

particles = doc.getElementsByTagName('particle');
numTracks = particles.getLength;

tracks      = cell(1, numTracks);
hullArea    = zeros(1, numTracks);
hullMaxDist = zeros(1, numTracks);

for i = 0:numTracks-1   % Java NodeList is 0‑based
    pNode = particles.item(i);
    detections = pNode.getElementsByTagName('detection');
    nSpots  = detections.getLength;
    T = zeros(nSpots, 3);  % [t x y]

    for j = 0:nSpots-1
        dNode = detections.item(j);
        T(j+1,1) = str2double(dNode.getAttribute('t'));
        T(j+1,2) = str2double(dNode.getAttribute('x'));
        T(j+1,3) = str2double(dNode.getAttribute('y'));
    end

    tracks{i+1} = T;

    % %% Convex‑hull metrics (skip if < 3 detections)
    % if nSpots >= 3
    %     k = convhull(T(:,2), T(:,3));
    %     hullArea(i+1) = polyarea(T(k,2), T(k,3));
    %     % Maximum distance between hull vertices (farthest pair)
    %     V = [T(k,2) T(k,3)];
    %     if exist('pdist', 'file') == 2
    %         hullMaxDist(i+1) = max(pdist(V));
    %     else  % fallback O(N^2)
    %         m = size(V,1);
    %         dmax = 0;
    %         for a = 1:m-1
    %             for b = a+1:m
    %                 d = hypot(V(a,1)-V(b,1), V(a,2)-V(b,2));
    %                 if d > dmax, dmax = d; end
    %             end
    %         end
    %         hullMaxDist(i+1) = dmax;
    %     end
    % else
    %     hullArea(i+1) = 0;
    %     hullMaxDist(i+1) = 0;
    % end
end

% %% Optionally save
% if saveMAT
%     [p,f,~] = fileparts(xmlFile);
%     matName = fullfile(p, [f '.mat']);
%     save(matName, 'tracks', 'hullArea', 'hullMaxDist');
%     fprintf('Saved results to %s\n', matName);
% end

end % function

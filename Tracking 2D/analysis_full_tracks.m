function analysis = analysis_full_tracks(tracks)

for i=1:length(tracks)
    temp_track = tracks{i}(~isnan(tracks{i}));
    analysis(i) = analyseTrack(tracks{i}, 10, 0.3);
end

end


function stats = analyseTrack(track, dt, maxLagFrac)
% analyseTrackConfined  MSD analysis with free- and confined-diffusion fits
%
%   stats = analyseTrackConfined(track)
%   stats = analyseTrackConfined(track, dt)
%   stats = analyseTrackConfined(track, dt, maxLagFrac)
%
%   INPUT
%       track        N×2 matrix of x- and y-coordinates
%       dt           (optional) sampling interval, default = 1
%       maxLagFrac   (optional) longest lag to include in MSD (fraction
%                    of track length), default = 0.25
%
%   OUTPUT (structure *stats*)
%       .duration     number of points track
%       .maxDisp      maximum pairwise displacement
%       .avgDisp      mean displacement
%       .hullArea     convex-hull area
%       .tau          lag-time vector used for MSD calculation
%       .msd          MSD values
%
%       ─ Free-diffusion fit ─
%       .D_free       diffusion coefficient (µm² s⁻¹)
%       .Rsq_free     R² of the linear fit
%
%       ─ Confined-diffusion fit ─
%       .D_conf       diffusion coefficient (µm² s⁻¹)
%       .L_conf       confinement length scale L  (µm)  (MSD plateau = L²)
%       .Rsq_conf     R² of the confined fit
%
%   The confined model assumes 2-D Brownian motion in a circular domain
%   with radius L / √2, giving MSD plateau L².

% ----------------------------------------- Input handling
if nargin < 2 || isempty(dt),         dt = 1;      end
if nargin < 3 || isempty(maxLagFrac), maxLagFrac = 0.25; end

if size(track,2) > 3
    track = track';
end

if size(track,2) == 3
    track = track(:,2:3);
end

assert(size(track,2) == 2,'track must be N×2 matrix');

N = size(track,1);

% ----------------------------------------- Basic statistics
stats.duration = N;
stats.maxDisp  = max(pdist(track));
stats.avgDisp = mean(pdist(track));
% ── Convex-hull area (handle degenerate cases)
if size(track,1) < 3                       % fewer than 3 points → no hull
    stats.hullArea = 0;
else
    if min([length(unique(track(:,1))),length(unique(track(:,2)))]) < 3                        % all points colinear
        stats.hullArea = 0;
    else
        k = convhull(track(:,1), track(:,2));  % indices of hull vertices
        stats.hullArea = polyarea(track(k,1), track(k,2));
    end
end

% ----------------------------------------- MSD calculation
maxLag = floor(maxLagFrac*N);
stats.tau = (1:maxLag).' * dt;     % column vector
stats.msd = zeros(maxLag,1);

for n = 1:maxLag
    disp_n       = track(1+n:end,:) - track(1:end-n,:);
    stats.msd(n) = mean(sum(disp_n.^2,2));
end

% =========================================================
% A) FREE-DIFFUSION FIT  (MSD = 4 D τ)
% =========================================================
p_lin          = stats.tau \ stats.msd;   % slope by least squares
stats.D_free   = p_lin / 4;               % 2-D: slope = 4 D
msdFit_free    = p_lin * stats.tau;
stats.Rsq_free = 1 - sum((stats.msd-msdFit_free).^2) / ...
                    sum((stats.msd-mean(stats.msd)).^2);

% =========================================================
% B) CONFINED-DIFFUSION FIT  (non-linear)
%     MSD(τ) = L² [1 - exp(-4 D τ / L²)]
% =========================================================
% Initial guesses:  D from free fit,  L² = plateau estimate (mean of
% last 20 % points)
idxPlateau     = round(maxLag*0.8):maxLag;
L2_init        = mean(stats.msd(idxPlateau));
par0           = [max(L2_init,eps)    max(stats.D_free,eps)]; % [L²  D]

% Objective to minimise
modelFun = @(p,t) p(1)*(1-exp(-4*p(2)*t/p(1)));          % MSD model
costFun  = @(p)  modelFun(p,stats.tau) - stats.msd;

% Use fminsearch (no toolboxes required)
opts      = optimset('Display','off');
p_est     = fminsearch(@(p) sum(costFun(p).^2), par0, opts);

L2        = p_est(1);
stats.L_conf = sqrt(L2);
stats.D_conf = p_est(2);

msdFit_conf  = modelFun(p_est,stats.tau);
stats.Rsq_conf = 1 - sum((stats.msd-msdFit_conf).^2) / ...
                     sum((stats.msd-mean(stats.msd)).^2);
end

